use std::collections::{HashMap, HashSet};
use std::mem;

use crate::{
    intersections::{SegmentId, SegmentIntersection},
    polygon::{calculate_region_area, Contour, Polygon},
    vector2::{Vector2f, Vector2fBits},
};

pub struct ProximityMerger {
    epsilon: f32,
    encountered_points: Vec<Vector2f>,
}

impl ProximityMerger {
    pub fn new(epsilon: f32) -> ProximityMerger {
        ProximityMerger {
            epsilon,
            encountered_points: Vec::new(),
        }
    }

    pub fn map(&mut self, point: Vector2f) -> Vector2f {
        for &encountered_point in self.encountered_points.iter() {
            if (point - encountered_point).length_squared() <= self.epsilon * self.epsilon {
                return encountered_point;
            }
        }

        self.encountered_points.push(point);
        return point;
    }
}

#[derive(Debug, Clone)]
pub struct Node {
    position: Vector2f,
    // node is marked dirty whenever a new edge is inserted
    dirty: bool,
    // edges, ordered CCW when not dirty
    edges: Vec<usize>,
    // look up index of edge by edge value (=neighbor index)
    // only valid when not dirty
    edge_to_edge_idx: HashMap<usize, usize>,
    // when tracing regions this stores whether the edge was traced along
    // note: this is per direction as the connected node will point back
    traced_edges: Vec<bool>,
}

impl Node {
    fn new(position: Vector2f) -> Node {
        Node {
            position,
            dirty: false,
            edges: Vec::new(),
            edge_to_edge_idx: HashMap::new(),
            traced_edges: Vec::new(),
        }
    }

    fn add_edge(&mut self, edge: usize) -> bool {
        if !self.edges.contains(&edge) {
            self.dirty = true;
            self.edges.push(edge);
            return true;
        }
        return false;
    }

    fn next_edge_ccw(&self, current_idx: usize) -> usize {
        debug_assert!(!self.dirty);
        if current_idx >= self.edges.len() {
            panic!("Edge index out of bounds")
        }

        return (current_idx + self.edges.len() + 1) % self.edges.len();
    }

    pub fn edges(&self) -> &Vec<usize> {
        &self.edges
    }

    pub fn position(&self) -> Vector2f {
        self.position
    }
}

#[derive(Debug, Clone)]
pub struct Graph {
    #[cfg(debug_assertions)]
    island: bool,
    nodes: Vec<Node>,
    position_lookup: HashMap<Vector2fBits, usize>,
}

impl Graph {
    pub fn new() -> Self {
        Self {
            #[cfg(debug_assertions)]
            island: false,
            nodes: Vec::new(),
            position_lookup: HashMap::new(),
        }
    }

    pub fn insert_segment(&mut self, from: Vector2f, to: Vector2f) {
        if from == to {
            return;
        }

        let from_idx = self.get_or_insert_node(from);
        let to_idx = self.get_or_insert_node(to);

        self.nodes[from_idx].add_edge(to_idx);
        self.nodes[to_idx].add_edge(from_idx);
    }

    fn get_or_insert_node(&mut self, position: Vector2f) -> usize {
        let position_bits = Vector2f::to_bits(position);
        if let Some(index) = self.position_lookup.get(&position_bits) {
            return *index;
        }

        debug_assert!(!self.island);
        let index = self.nodes.len();
        self.nodes.push(Node::new(position));
        self.position_lookup.insert(position_bits, index);

        return index;
    }

    pub fn split_by_islands(mut self) -> Vec<Graph> {
        let mut node_moved_out = vec![false; self.nodes.len()];
        let mut subgraphs = Vec::new();

        for node_idx in 0..node_moved_out.len() {
            if node_moved_out[node_idx] {
                continue;
            }

            let island = self.find_island(node_idx);

            // move out nodes
            let mut graph_nodes = Vec::with_capacity(island.len());
            for &i in &island {
                node_moved_out[i] = true;

                // set position so that we can use it to look up the index later
                let placeholder_node = Node::new(self.nodes[i].position);

                // use mem::replace so that we can keep the allocated edge Vec
                graph_nodes.push(mem::replace(&mut self.nodes[i], placeholder_node));
            }

            // create position lookup
            let mut graph_position_lookup: HashMap<Vector2fBits, usize> = HashMap::new();
            for (i, node) in graph_nodes.iter().enumerate() {
                graph_position_lookup.insert(Vector2f::to_bits(node.position), i);
            }

            // fix node indices
            for node in &mut graph_nodes {
                for neighbor_idx in &mut node.edges {
                    *neighbor_idx = *graph_position_lookup
                        .get(&Vector2f::to_bits(self.nodes[*neighbor_idx].position))
                        .unwrap();
                }
            }

            subgraphs.push(Graph {
                #[cfg(debug_assertions)]
                island: true,
                nodes: graph_nodes,
                position_lookup: graph_position_lookup,
            });
        }

        debug_assert!(self.nodes.len() == subgraphs.iter().map(|g| g.nodes.len()).sum());

        #[cfg(debug_assertions)]
        for subgraph in &subgraphs {
            subgraph.validate();
        }

        return subgraphs;
    }

    fn validate(&self) {
        for (i, node) in self.nodes.iter().enumerate() {
            debug_assert!(!node.edges.is_empty(), "node has no edges");
            debug_assert!(
                self.position_lookup.get(&Vector2f::to_bits(node.position)) == Some(&i),
                "node has wrong lookup"
            );
            debug_assert!(
                node.edges.len() == node.edges.iter().collect::<HashSet<_>>().len(),
                "node has duplicate edges"
            );

            for &edge in &node.edges {
                let neighbor = &self.nodes[edge];
                debug_assert!(neighbor.edges.contains(&i), "neighbor has no back-edge");
            }
        }
    }

    fn find_island(&self, start_node_idx: usize) -> HashSet<usize> {
        let mut island_node_indices = HashSet::new();

        let mut stack = Vec::new();
        island_node_indices.insert(start_node_idx);
        stack.extend(&self.nodes[start_node_idx].edges);

        while !stack.is_empty() {
            if let Some(&node_idx) = stack.last() {
                if !island_node_indices.contains(node_idx) {
                    island_node_indices.insert(*node_idx);
                    stack.extend(&self.nodes[*node_idx].edges)
                } else {
                    stack.pop();
                }
            }
        }

        return island_node_indices;
    }

    /// Traces the regions (faces) of the planar graph. The first traced region is the outline.
    /// The outline always runs CCW, the interior regions run CW.
    /// The result of this function is only valid if there are no disjoint islands.
    pub fn trace_regions(&mut self) -> Island {
        // clean up nodes: sort edges in CCW order and create edge reverse lookup
        let nodes_ptr = self.nodes.as_ptr();
        for node in &mut self.nodes {
            if !node.dirty {
                continue;
            }

            node.edges.sort_by(|n1, n2| {
                let n1 = unsafe { &*(nodes_ptr.add(*n1)) };
                let n2 = unsafe { &*(nodes_ptr.add(*n2)) };

                let a1 = (n1.position - node.position).angle();
                let a2 = (n2.position - node.position).angle();

                return f32::total_cmp(&a1, &a2);
            });

            node.edge_to_edge_idx.clear();
            for (i, &edge) in node.edges.iter().enumerate() {
                node.edge_to_edge_idx.insert(edge, i);
            }

            node.dirty = false;
        }

        // angle 0 points in +X direction, with edges going CCW
        // -> start tracing with node of lowest y-value to ensure the first region is the outline
        let node_idx_min_y = (0..self.nodes.len())
            .min_by(|&i, &j| f32::total_cmp(&self.nodes[i].position.y, &self.nodes[j].position.y))
            .unwrap();

        // reset traced edges
        for node in &mut self.nodes {
            node.traced_edges.clear();
            node.traced_edges.resize(node.edges.len(), false);
        }

        // since we start with node with min y, the first region is the outline
        let outline = self.trace_region(node_idx_min_y, 0);

        let mut interior = Vec::new();
        for node_idx in 0..self.nodes.len() {
            for edge_idx in 0..self.nodes[node_idx].traced_edges.len() {
                if self.nodes[node_idx].traced_edges[edge_idx] {
                    continue;
                }

                let region = self.trace_region(node_idx, edge_idx);
                interior.push(region);
            }
        }

        let island = Island { outline, interior };
        // all edges were traced
        debug_assert!(self.nodes.iter().all(|n| n.traced_edges.iter().all(|&t| t)));
        debug_assert!(island_contour_windings_are_valid(&island));
        debug_assert!(region_edge_count_matches_graph(&self, &island));
        return island;
    }

    fn trace_region(&mut self, start_node_idx: usize, start_outgoing_edge_idx: usize) -> Contour {
        let mut region = Contour::new();

        let mut current_node_idx = start_node_idx;
        let mut current_outgoing_edge_idx = start_outgoing_edge_idx;

        loop {
            // add point to contour
            let current_node = &mut self.nodes[current_node_idx];
            region.push(current_node.position);

            // mark the edge as traced
            if current_node.traced_edges[current_outgoing_edge_idx] {
                panic!("Edge traced twice");
            }
            current_node.traced_edges[current_outgoing_edge_idx] = true;

            // traverse the edge
            let previous_node_idx = current_node_idx;
            current_node_idx = current_node.edges[current_outgoing_edge_idx];

            // find incoming edge on connecting node
            let current_node = &self.nodes[current_node_idx];
            let incoming_edge_idx = *current_node
                .edge_to_edge_idx
                .get(&previous_node_idx)
                .unwrap();

            // new outgoing edge is next to incomine one, rotating CCW, i.e.,
            // we need to take the next edge in order of edge index
            current_outgoing_edge_idx = current_node.next_edge_ccw(incoming_edge_idx);

            // if we reached original node and edge, we are done
            if current_node_idx == start_node_idx
                && current_outgoing_edge_idx == start_outgoing_edge_idx
            {
                return region;
            }
        }
    }

    pub fn nodes(&self) -> &Vec<Node> {
        &self.nodes
    }
}

#[derive(Clone)]
pub struct Island {
    // runs CCW
    pub outline: Contour,
    // runs CW
    pub interior: Vec<Contour>,
}

fn island_contour_windings_are_valid(island: &Island) -> bool {
    if calculate_region_area(&island.outline) <= 0.0 {
        return false;
    }

    for interior in &island.interior {
        if calculate_region_area(interior) >= 0.0 {
            return false;
        }
    }

    return true;
}

fn region_edge_count_matches_graph(graph: &Graph, island: &Island) -> bool {
    // counts each edge twice because it's bidirectional
    let graph_edge_count: usize = graph.nodes.iter().map(|n| n.edges.len()).sum();
    // counts each edge twice because each edge is traced in both directions
    let region_edge_count =
        island.outline.len() + island.interior.iter().map(|r| r.len()).sum::<usize>();
    println!("graph: {}, island: {}", graph_edge_count, region_edge_count);
    return graph_edge_count == region_edge_count;
}

fn strip_intersection_type(
    intersections: &HashMap<SegmentId, Vec<SegmentIntersection>>,
) -> HashMap<SegmentId, Vec<Vector2f>> {
    let mut stripped = HashMap::new();

    for (segment_id, intersections) in intersections.iter() {
        let mut stripped_intersections = Vec::new();
        for intersection in intersections.iter() {
            match intersection {
                SegmentIntersection::Point(p) => stripped_intersections.push(*p),
                SegmentIntersection::Segment(p1, p2) => {
                    stripped_intersections.push(*p1);
                    stripped_intersections.push(*p2);
                }
                _ => panic!("Unexpected intersection type"),
            }
        }
        stripped.insert(*segment_id, stripped_intersections);
    }

    return stripped;
}

fn merge_nearby_points(
    proximity_merger: &mut ProximityMerger,
    intersections: &HashMap<SegmentId, Vec<Vector2f>>,
) -> HashMap<SegmentId, Vec<Vector2f>> {
    let mut merged = intersections.clone();

    for (_, points) in merged.iter_mut() {
        for point in points.iter_mut() {
            *point = proximity_merger.map(*point);
        }
    }

    return merged;
}

pub fn build_graph(
    polygon: &Polygon,
    intersections: &HashMap<SegmentId, Vec<SegmentIntersection>>,
    epsilon: f32,
) -> (Graph, ProximityMerger) {
    let stripped_intersections = strip_intersection_type(intersections);

    let mut proximity_merger = ProximityMerger::new(epsilon);
    let mut merged_intersections =
        merge_nearby_points(&mut proximity_merger, &stripped_intersections);

    let mut graph = Graph::new();

    for (contour_idx, contour) in polygon.contours().iter().enumerate() {
        for (segment_idx, start_point) in contour.iter().enumerate() {
            let end_point = contour[(segment_idx + 1) % contour.len()];

            let seg_id = SegmentId {
                contour: contour_idx,
                segment: segment_idx,
            };
            // calling unwrap because every line must have an intersection at the start and endpoint
            let intersections_opt = merged_intersections.get_mut(&seg_id);
            if intersections_opt.is_none() {
                panic!(
                    "Segment with ID ({}, {}) has no intersections",
                    seg_id.contour, seg_id.segment
                );
            }
            let intersections = intersections_opt.unwrap();

            // sort intersection points along the segment going from start_point to end_point
            let segment_vec = end_point - *start_point;
            intersections.sort_by(|p1, p2| {
                let v1 = segment_vec.dot(*p1 - *start_point);
                let v2 = segment_vec.dot(*p2 - *start_point);
                return f32::total_cmp(&v1, &v2);
            });

            // each pair in order is a sub-segment
            for i in 0..intersections.len() - 1 {
                graph.insert_segment(intersections[i], intersections[i + 1]);
            }
        }
    }

    return (graph, proximity_merger);
}
