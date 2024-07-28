use std::collections::HashMap;

use crate::{
    intersections::{SegmentId, SegmentIntersection},
    polygon::Polygon,
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
