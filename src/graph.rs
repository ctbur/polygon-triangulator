use std::collections::HashMap;

use crate::{
    intersections::{SegmentId, SegmentIntersection},
    polygon::{Countour, Polygon},
    vector2::Vector2f,
};

#[derive(Debug, Clone)]
pub struct Node {
    pub position: Vector2f,
    pub edges: Vec<usize>,
}

#[derive(Debug, Clone)]
pub struct Graph {
    nodes: Vec<Node>,
    position_lookup: HashMap<(i64, i64), usize>,
    epsilon: f32,
}

impl Graph {
    pub fn new(epsilon: f32) -> Self {
        Self {
            nodes: Vec::new(),
            position_lookup: HashMap::new(),
            epsilon,
        }
    }

    fn insert_segment(&mut self, start: Vector2f, end: Vector2f) {
        if start == end {
            return;
        }

        let start_index = self.get_or_insert_node(start);
        let end_index = self.get_or_insert_node(end);

        if !self.nodes[start_index].edges.contains(&start_index) {
            self.nodes[start_index].edges.push(end_index);
        }
        if !self.nodes[end_index].edges.contains(&end_index) {
            self.nodes[end_index].edges.push(start_index);
        }
    }

    fn get_or_insert_node(&mut self, position: Vector2f) -> usize {
        let position_i64 = self.position_to_i64(position);
        if let Some(index) = self.position_lookup.get(&position_i64) {
            return *index;
        }

        let index = self.nodes.len();
        self.nodes.push(Node {
            position,
            edges: Vec::new(),
        });
        self.position_lookup.insert(position_i64, index);

        return index;
    }

    fn position_to_i64(&self, position: Vector2f) -> (i64, i64) {
        (
            ((position.x as f64) / (self.epsilon as f64)).floor() as i64,
            ((position.y as f64) / (self.epsilon as f64)).floor() as i64,
        )
    }

    pub fn nodes(&self) -> &Vec<Node> {
        &self.nodes
    }

    pub fn into_nodes(self) -> Vec<Node> {
        self.nodes
    }
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
    intersections: &HashMap<SegmentId, Vec<Vector2f>>,
    epsilon: f32,
) -> HashMap<SegmentId, Vec<Vector2f>> {
    let mut merged = intersections.clone();
    let mut encountered_points = Vec::new();

    for (_, points) in merged.iter_mut() {
        for point in points.iter_mut() {
            // search encountered_points for a point that is closer than epsilon
            let mut found_neighbor_point: Option<&Vector2f> = None;
            for encountered_point in encountered_points.iter() {
                if (*point - *encountered_point).length_squared() <= epsilon * epsilon {
                    found_neighbor_point = Some(encountered_point);
                    break;
                }
            }

            if let Some(neighbor_point) = found_neighbor_point {
                // replace point with neighbor_point
                *point = *neighbor_point;
            } else {
                // add point to encountered_points
                encountered_points.push(*point);
            }
        }
    }

    return merged;
}

pub fn build_graph(
    polygon: &Polygon,
    intersections: &HashMap<SegmentId, Vec<SegmentIntersection>>,
    epsilon: f32,
) -> Graph {
    let stripped_intersections = strip_intersection_type(intersections);
    let mut merged_intersections = merge_nearby_points(&stripped_intersections, epsilon);

    let mut graph = Graph::new(epsilon);

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

    return graph;
}

/// Traces the regions (faces) of the planar graph. The first region is the outline.
fn trace_regions(graph: Graph) -> Vec<Countour> {
    let mut nodes = graph.into_nodes();

    // sort edges in CCW order
    let nodes_ptr = nodes.as_ptr();
    for node in &mut nodes {
        node.edges.sort_by(|n1, n2| {
            let n1: &Node = unsafe { &*(nodes_ptr.add(*n1)) };
            let n2: &Node = unsafe { &*(nodes_ptr.add(*n2)) };

            let a1 = (n1.position - node.position).angle();
            let a2 = (n2.position - node.position).angle();

            return f32::total_cmp(&a1, &a2);
        });
    }

    let mut nodes_map = HashMap::new();
    for (node_idx, node) in nodes.into_iter().enumerate() {
        nodes_map.insert(node_idx, node);
    }

    let mut regions = Vec::new();
    loop {
        // angle 0 points in +X, with edges going CCW
        // -> find point with lowest y value to start tracing with outline region
        let mut node_idx_lowest_y = None;
        for (node_idx, node) in &nodes_map {
            match node_idx_lowest_y {
                None => node_idx_lowest_y = Some(*node_idx),
                Some(idx) => {
                    if node.position.y < nodes_map.get(&idx).unwrap().position.y {
                        node_idx_lowest_y = Some(idx);
                    }
                }
            }
        }

        if node_idx_lowest_y.is_none() {
            // no nodes are left
            break;
        }

        let mut region = Countour::new();
        let start_node_idx = node_idx_lowest_y.unwrap();
        let mut current_node_idx = start_node_idx;
        let mut outgoing_edge = 0;
        loop {
            // add point to contour
            let current_node = nodes_map.get_mut(&current_node_idx).unwrap();
            region.push(current_node.position);

            // remove the edge, and the node if there are no edges left
            current_node.edges.remove(outgoing_edge);
            if current_node.edges.is_empty() {
                nodes_map.remove(&current_node_idx);
            }

            // traverse the edge
            let previous_node_idx = current_node_idx;
            current_node_idx = outgoing_edge;

            // if we reached original node, we are done
            if current_node_idx == start_node_idx {
                break;
            }

            // find incoming edge on connecting node
            let next_node = nodes_map.get(&outgoing_edge).unwrap();
            let mut incoming_edge_idx_opt = None;
            for (edge_idx, edge) in next_node.edges.iter().enumerate() {
                // if the edge points to the current node, that's where we come from
                if *edge == previous_node_idx {
                    incoming_edge_idx_opt = Some(edge_idx);
                    break;
                }
            }
            let incoming_edge_idx = incoming_edge_idx_opt.unwrap();

            // new outgoing edge is next to incomine one, rotating CCW
            let current_node = nodes_map.get(&current_node_idx).unwrap();
            outgoing_edge = current_node.edges[(incoming_edge_idx + 1) % current_node.edges.len()];
        }

        regions.push(region);
    }

    return regions;
}
