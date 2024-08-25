use std::collections::HashMap;

use crate::{polygon::Contour, vector2::Vector2f};

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

    pub fn insert_segment(&mut self, start: Vector2f, end: Vector2f) {
        if start == end {
            return;
        }

        let start_index = self.get_or_insert_node(start);
        let end_index = self.get_or_insert_node(end);

        if !self.nodes[start_index].edges.contains(&end_index) {
            self.nodes[start_index].edges.push(end_index);
        }
        if !self.nodes[end_index].edges.contains(&start_index) {
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

pub fn build_graph(subdivided_contours: &Vec<Contour>, epsilon: f32) -> Graph {
    let mut graph = Graph::new(epsilon);
    for contour in subdivided_contours {
        for i in 0..contour.len() {
            let (from, to) = (contour[i], contour[(i + 1) % contour.len()]);
            graph.insert_segment(from, to);
        }
    }

    return graph;
}
