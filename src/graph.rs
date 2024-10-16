use std::collections::{HashMap, HashSet};
use std::mem;

use crate::{
    polygon::{calculate_region_area, Contour},
    vector2::{Vector2f, Vector2fBits},
};

#[derive(Debug, Clone)]
pub struct Node {
    position: Vector2f,
    // edges, arbirary order, when tracing regions it is reodered to go CCW
    edges: Vec<usize>,
    // look up index of edge by edge value (=neighbor index)
    // is only initialized before tracing regions
    edge_to_edge_idx: HashMap<usize, usize>,
    // when tracing regions this stores whether the edge was traced along
    // note: this is per direction as the connected node will point back
    traced_edges: Vec<bool>,
}

impl Node {
    fn new(position: Vector2f) -> Node {
        Node {
            position,
            edges: Vec::new(),
            edge_to_edge_idx: HashMap::new(),
            traced_edges: Vec::new(),
        }
    }

    fn add_edge(&mut self, edge: usize) -> bool {
        if !self.edges.contains(&edge) {
            self.edges.push(edge);
            return true;
        }
        return false;
    }

    fn remove_edge(&mut self, edge: usize) -> bool {
        if let Some(idx) = self.edges.iter().position(|&e| e == edge) {
            self.edges.remove(idx);
            return true;
        }
        return false;
    }

    fn next_edge_ccw(&self, current_idx: usize) -> usize {
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
    nodes: Vec<Node>,
    position_lookup: HashMap<Vector2fBits, usize>,
}

impl Graph {
    pub fn new() -> Self {
        Self {
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

        let added_from = self.nodes[from_idx].add_edge(to_idx);
        let added_to = self.nodes[to_idx].add_edge(from_idx);

        debug_assert!(added_from == added_to);
    }

    fn get_or_insert_node(&mut self, position: Vector2f) -> usize {
        let position_bits = Vector2f::to_bits(position);
        if let Some(index) = self.position_lookup.get(&position_bits) {
            return *index;
        }

        let index = self.nodes.len();
        self.nodes.push(Node::new(position));
        self.position_lookup.insert(position_bits, index);

        return index;
    }

    pub fn remove_segment(&mut self, from: Vector2f, to: Vector2f) -> bool {
        if from == to {
            return false;
        }

        let from_idx = self.get_or_insert_node(from);
        let to_idx = self.get_or_insert_node(to);

        let removed_from = self.nodes[from_idx].remove_edge(to_idx);
        let removed_to = self.nodes[to_idx].remove_edge(from_idx);

        debug_assert!(removed_from == removed_to);

        return removed_from;
    }

    /// Trims off nodes with only one edge, and nodes that will have
    /// one edge after trimming of those nodes.
    pub fn trim_pendant_paths(&mut self) -> bool {
        let mut removed = false;

        for i in 0..self.nodes.len() {
            if self.trim_pendant_path(i) {
                removed = true;
            }
        }

        return removed;
    }

    fn trim_pendant_path(&mut self, mut node_idx: usize) -> bool {
        let mut removed = false;

        while self.nodes[node_idx].edges.len() == 1 {
            let neighbor_idx = self.nodes[node_idx].edges[0];
            self.remove_segment(
                self.nodes[node_idx].position,
                self.nodes[neighbor_idx].position,
            );
            node_idx = neighbor_idx;

            removed = true;
        }

        return removed;
    }

    pub fn split_by_islands(mut self) -> Vec<Graph> {
        let mut node_moved_out = vec![false; self.nodes.len()];
        let mut subgraphs = Vec::new();

        let num_empty_nodes = self.nodes.iter().filter(|n| n.edges.is_empty()).count();

        for node_idx in 0..node_moved_out.len() {
            // skip nodes that had all their edges deleted
            if self.nodes[node_idx].edges.is_empty() {
                continue;
            }

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
                nodes: graph_nodes,
                position_lookup: graph_position_lookup,
            });
        }

        // subgraphs should have as many nodes as original graph minus nodes that had all their edges removed
        debug_assert!(
            self.nodes.len() - num_empty_nodes == subgraphs.iter().map(|g| g.nodes.len()).sum()
        );

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
    pub fn trace_regions(&mut self) -> (Contour, Vec<Contour>) {
        // clean up nodes: sort edges in CCW order and create edge reverse lookup
        let nodes_ptr = self.nodes.as_ptr();
        for node in &mut self.nodes {
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
        }

        // angle 0 points in +X direction, with edges going CCW
        // -> start tracing with node of lowest y-value to ensure the first region is the outline
        let node_idx_min_y = (0..self.nodes.len())
            // filter out nodes that had all their edges deleted
            .filter(|&i| !&self.nodes[i].edges.is_empty())
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

        // all edges were traced
        debug_assert!(self.nodes.iter().all(|n| n.traced_edges.iter().all(|&t| t)));
        debug_assert!(contour_windings_are_valid(&outline, &interior));
        debug_assert!(region_edge_count_matches_graph(&self, &outline, &interior));
        return (outline, interior);
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

fn contour_windings_are_valid(outline: &Contour, interior: &Vec<Contour>) -> bool {
    if calculate_region_area(outline) <= 0.0 {
        return false;
    }

    for interior in interior {
        if calculate_region_area(interior) >= 0.0 {
            return false;
        }
    }

    return true;
}

fn region_edge_count_matches_graph(
    graph: &Graph,
    outline: &Contour,
    interior: &Vec<Contour>,
) -> bool {
    // counts each edge twice because it's bidirectional
    let graph_edge_count: usize = graph.nodes.iter().map(|n| n.edges.len()).sum();
    // counts each edge twice because each edge is traced in both directions
    let region_edge_count = outline.len() + interior.iter().map(|r| r.len()).sum::<usize>();
    return graph_edge_count == region_edge_count;
}

pub fn build_graph(subdivided_contours: &[Contour]) -> Graph {
    let mut graph = Graph::new();

    for contour in subdivided_contours {
        for current in 0..contour.len() {
            let next = (current + 1) % contour.len();
            graph.insert_segment(contour[current], contour[next]);
        }
    }

    // trim pendant paths as they will pose a problem with region tracing otherwise
    graph.trim_pendant_paths();

    return graph;
}
