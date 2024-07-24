use std::collections::HashMap;

use crate::{
    intersections::{SegmentId, SegmentIntersection},
    polygon::Polygon,
    vector2::Vector2f,
};

pub struct Node {
    position: Vector2f,
    edges: Vec<usize>,
}

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

    pub fn insert_segment(&mut self, start: Vector2f, end: Vector2f) {
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
        if let Some(index) = self.position_lookup.get(&position) {
            return *index;
        }

        let index = self.nodes.len();
        self.nodes.push(Node {
            position,
            edges: Vec::new(),
        });
        self.position_lookup.insert(position, index);

        return index;
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
    let merged_intersections = merge_nearby_points(&stripped_intersections, epsilon);
    let mut graph = Graph::new(epsilon);
    return graph;
}
