use std::collections::{BTreeMap, HashMap};

use crate::graph::Graph;
use crate::polygon::Contour;
use crate::vector2::Vector2f;

struct Node {
    position: Vector2f,
    // edge_idx is the index of the edge in a list where the edges are ordered CCW
    neighbor_to_edge_idx: HashMap<usize, usize>,
    edges_by_idx: BTreeMap<usize, usize>,
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

    let mut nodes = HashMap::new();
    for (idx, node) in graph_nodes.into_iter().enumerate() {
        let mut neighbor_to_edge_idx = HashMap::new();
        let mut edges_by_idx = BTreeMap::new();

        for (edge_idx, neighbor) in node.edges.iter().enumerate() {
            neighbor_to_edge_idx.insert(*neighbor, edge_idx);
            edges_by_idx.insert(edge_idx, *neighbor);
        }

        nodes.insert(
            idx,
            Node {
                position: node.position,
                edges_by_idx,
                neighbor_to_edge_idx,
            },
        );
    }

    let mut regions = Vec::new();
    loop {
        // angle 0 points in +X, with edges going CCW
        // -> find point with lowest y value to start tracing with outline region
        let mut node_idx_lowest_y = None;
        for (node_idx, node) in &nodes {
            match node_idx_lowest_y {
                None => node_idx_lowest_y = Some(*node_idx),
                Some(idx) => {
                    if node.position.y < nodes.get(&idx).unwrap().position.y {
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
        let mut current_node_idx = start_node_idx;
        let mut outgoing_edge_idx = *nodes
            .get_mut(&current_node_idx)
            .unwrap()
            .edges_by_idx
            .first_entry()
            .unwrap()
            .key();
        loop {
            // add point to contour
            let mut current_node = nodes.get_mut(&current_node_idx).unwrap();
            region.push(current_node.position);

            // remove the edge, and the node if there are no edges left
            let outgoing_edge = current_node
                .edges_by_idx
                .remove(&outgoing_edge_idx)
                .unwrap();
            if current_node.edges_by_idx.is_empty() {
                nodes.remove(&current_node_idx);
            }

            // traverse the edge
            let previous_node_idx = current_node_idx;
            current_node_idx = outgoing_edge;

            // if we reached original node, we are done
            if current_node_idx == start_node_idx {
                break;
            }

            // find incoming edge on connecting node
            current_node = nodes.get_mut(&current_node_idx).unwrap();
            let incoming_edge_idx = current_node
                .neighbor_to_edge_idx
                .get(&previous_node_idx)
                .unwrap();

            // new outgoing edge is next to incomine one, rotating CCW, i.e.,
            // we need to take the next edge in order of edge index
            let outgoing_edge_idx_opt = current_node
                .edges_by_idx
                .range(incoming_edge_idx + 1..)
                .next();
            outgoing_edge_idx = match outgoing_edge_idx_opt {
                // there is an edge with larger index than the incoming one
                Some((idx, _)) => *idx,
                // we need to wrap around
                None => *current_node.edges_by_idx.first_entry().unwrap().key(),
            };
        }

        regions.push(region);
    }

    return regions;
}
