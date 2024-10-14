use crate::{
    graph::Island,
    polygon::{self, Contour},
};

pub struct IslandNode {
    pub island: usize,
    pub interior: Vec<InteriorNode>,
}

pub struct InteriorNode {
    pub children: Vec<IslandNode>,
}

/// Creates a hierarchy of regions that describes which regions contain which other ones.
/// `regions` is a sequence of CCW regions each followed by one or more CW regions.
/// The regions are provided in increasing min-y-coordinate.
pub fn build_region_hierarchy(islands: &[Island]) -> Vec<IslandNode> {
    let mut island_nodes = Vec::with_capacity(islands.len());
    for (island_idx, island) in islands.iter().enumerate() {
        let mut interior_nodes = Vec::with_capacity(island.interior.len());
        for _ in 0..island.interior.len() {
            interior_nodes.push(InteriorNode {
                children: Vec::new(),
            })
        }
        island_nodes.push(IslandNode {
            island: island_idx,
            interior: interior_nodes,
        })
    }

    build_hierarchy_from_island_nodes(islands, &mut island_nodes);
    return island_nodes;
}

fn build_hierarchy_from_island_nodes(islands: &[Island], island_nodes: &mut Vec<IslandNode>) {
    // check all island nodes with each other and move inside if contained
    let mut container_idx = 0;
    while container_idx < island_nodes.len() {
        let mut containee_idx = 0;
        while containee_idx < island_nodes.len() {
            if container_idx == containee_idx {
                containee_idx += 1;
                continue;
            }

            if let Some(interior_node_idx) = find_interior_node_containing_island_node(
                islands,
                &island_nodes[container_idx],
                &island_nodes[containee_idx],
            ) {
                // move containee into container node
                let containee = island_nodes.remove(containee_idx);
                island_nodes[container_idx].interior[interior_node_idx]
                    .children
                    .push(containee);

                // update iterator
                if container_idx > containee_idx {
                    container_idx -= 1;
                }
            }

            containee_idx += 1;
        }

        container_idx += 1;
    }

    // repeat build hierarchy for all children for multiple nesting
    for island_node in island_nodes {
        for interior_node in &mut island_node.interior {
            build_hierarchy_from_island_nodes(islands, &mut interior_node.children);
        }
    }
}

fn find_interior_node_containing_island_node(
    islands: &[Island],
    container: &IslandNode,
    containee: &IslandNode,
) -> Option<usize> {
    if !are_regions_nested(
        &islands[container.island].outline,
        &islands[containee.island].outline,
    ) {
        return None;
    }

    for idx in 0..container.interior.len() {
        let region = &islands[container.island].interior[idx];
        let outline = &islands[containee.island].outline;
        if are_regions_nested(region, outline) {
            return Some(idx);
        }
    }

    panic!("Node is contained by outline but not by interior region");
}

fn are_regions_nested(container: &Contour, containee: &Contour) -> bool {
    return polygon::calculate_winding_number(container, containee[0]) != 0;
}

#[cfg(test)]
mod tests {
    use crate::{regions, showcase};

    use super::*;

    #[test]
    fn todo() {}
}
