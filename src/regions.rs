use std::collections::{BTreeMap, HashMap};

use crate::graph::Graph;
use crate::polygon::Contour;
use crate::vector2::Vector2f;

struct Node {
    position: Vector2f,
    // whether the edge was traced along
    // note: this is per direction as the connected node will point back
    traced_edges: Vec<bool>,
    // edges are ordered CCW
    edges: Vec<usize>,
    // look up index of edge by edge value (=neighbor index)
    edge_to_edge_idx: HashMap<usize, usize>,
}

impl Node {
    fn first_untraced_index(&self) -> Option<usize> {
        for (idx, traced) in self.traced_edges.iter().enumerate() {
            if !traced {
                return Some(idx);
            }
        }
        return None;
    }

    fn next_edge_ccw(&self, current_idx: usize) -> usize {
        if current_idx >= self.edges.len() {
            panic!("Edge index out of bounds")
        }

        return (current_idx + self.edges.len() + 1) % self.edges.len();
    }
}

/// Traces the regions (faces) of the planar graph. The first region is the outline.
pub fn trace_regions(graph: Graph) -> Vec<Contour> {
    let mut graph_nodes = graph.into_nodes();

    // sort edges in CCW order
    let nodes_ptr = graph_nodes.as_ptr();
    for node in &mut graph_nodes {
        node.edges.sort_by(|n1, n2| {
            let n1 = unsafe { &*(nodes_ptr.add(*n1)) };
            let n2 = unsafe { &*(nodes_ptr.add(*n2)) };

            let a1 = (n1.position - node.position).angle();
            let a2 = (n2.position - node.position).angle();

            return f32::total_cmp(&a1, &a2);
        });
    }

    let mut nodes = Vec::new();
    for (idx, node) in graph_nodes.into_iter().enumerate() {
        let mut edge_to_edge_idx = HashMap::new();

        for (edge_idx, neighbor) in node.edges.iter().enumerate() {
            edge_to_edge_idx.insert(*neighbor, edge_idx);
        }

        nodes.push(Node {
            position: node.position,
            traced_edges: vec![false; node.edges.len()],
            edges: node.edges,
            edge_to_edge_idx,
        });
    }

    let mut regions = Vec::new();
    loop {
        // angle 0 points in +X, with edges going CCW
        // -> find point with lowest y value to start tracing with outline region
        let mut node_idx_lowest_y = None;
        for (node_idx, node) in nodes.iter().enumerate() {
            if node.first_untraced_index().is_none() {
                continue;
            }

            match node_idx_lowest_y {
                None => node_idx_lowest_y = Some(node_idx),
                Some(idx) => {
                    if node.position.y < nodes[idx].position.y {
                        node_idx_lowest_y = Some(idx);
                    }
                }
            }
        }

        if node_idx_lowest_y.is_none() {
            // no nodes are left
            break;
        }

        let mut region = Contour::new();
        let start_node_idx = node_idx_lowest_y.unwrap();
        let start_outgoing_edge_idx = nodes[start_node_idx].first_untraced_index().unwrap();

        let mut current_node_idx = start_node_idx;
        let mut current_outgoing_edge_idx = start_outgoing_edge_idx;
        loop {
            // add point to contour
            let current_node = &mut nodes[current_node_idx];
            region.push(current_node.position);

            // mark the edge as traced
            if current_node.traced_edges[current_outgoing_edge_idx] {
                //panic!("Edge traced twice");
            }
            current_node.traced_edges[current_outgoing_edge_idx] = true;

            // traverse the edge
            let previous_node_idx = current_node_idx;
            current_node_idx = current_node.edges[current_outgoing_edge_idx];

            // if we reached original node, we are done
            if current_node_idx == start_node_idx {
                break;
            }

            // find incoming edge on connecting node
            let current_node = &nodes[current_node_idx];
            let incoming_edge_idx = *current_node
                .edge_to_edge_idx
                .get(&previous_node_idx)
                .unwrap();

            // new outgoing edge is next to incomine one, rotating CCW, i.e.,
            // we need to take the next edge in order of edge index
            current_outgoing_edge_idx = current_node.next_edge_ccw(incoming_edge_idx);
        }

        regions.push(region);
    }

    return regions;
}
