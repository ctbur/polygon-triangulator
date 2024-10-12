use std::collections::{HashMap, HashSet};

use crate::graph::Graph;
use crate::polygon::{calculate_region_area, Contour};
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

#[derive(Clone)]
pub struct Island {
    // runs CCW
    pub outline: Contour,
    // runs CW
    pub interior: Vec<Contour>,
}

/// Traces the regions (faces) of the planar graph. The first region is the outline.
/// The outline always runs CCW, the interior regions run CW.
pub fn trace_regions(graph: Graph) -> Vec<Island> {
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

    // convert graph::Node into Node
    let mut nodes = Vec::new();
    for node in graph_nodes.into_iter() {
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

    // angle 0 points in +X direction, with edges going CCW
    // -> trace islands, in order of y-value to ensure that the first region is the outline
    let mut node_indices_by_y = Vec::from_iter(0..nodes.len());
    node_indices_by_y.sort_by(|&i, &j| f32::total_cmp(&nodes[i].position.y, &nodes[j].position.y));

    let mut islands = Vec::new();
    for i in node_indices_by_y {
        if nodes[i].first_untraced_index().is_none() {
            continue;
        }

        let island = trace_island(&mut nodes, i);
        debug_assert!(island_contour_windings_are_valid(&island));
        islands.push(island);
    }

    return islands;
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

fn trace_island(nodes: &mut Vec<Node>, node_idx_lowest_y: usize) -> Island {
    let outline = trace_region(nodes, node_idx_lowest_y);

    let mut island_node_indices = find_island(nodes, node_idx_lowest_y);
    let mut interior = Vec::new();
    while let Some(&node_idx) = island_node_indices.iter().next() {
        if nodes[node_idx].first_untraced_index().is_none() {
            island_node_indices.remove(&node_idx);
            continue;
        }

        let region = trace_region(nodes, node_idx);
        interior.push(region);
    }

    return Island { outline, interior };
}

fn find_island(nodes: &mut Vec<Node>, start_node_idx: usize) -> HashSet<usize> {
    let mut island_node_indices = HashSet::new();

    let mut stack = Vec::new();
    island_node_indices.insert(start_node_idx);
    stack.extend(&nodes[start_node_idx].edges);

    while !stack.is_empty() {
        if let Some(&node_idx) = stack.last() {
            if !island_node_indices.contains(node_idx) {
                island_node_indices.insert(*node_idx);
                stack.extend(&nodes[*node_idx].edges)
            } else {
                stack.pop();
            }
        }
    }

    return island_node_indices;
}

fn trace_region(nodes: &mut Vec<Node>, start_node_idx: usize) -> Contour {
    let mut region = Contour::new();
    let start_outgoing_edge_idx = nodes[start_node_idx].first_untraced_index().unwrap();

    let mut current_node_idx = start_node_idx;
    let mut current_outgoing_edge_idx = start_outgoing_edge_idx;

    loop {
        // add point to contour
        let current_node = &mut nodes[current_node_idx];
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
        let current_node = &nodes[current_node_idx];
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
