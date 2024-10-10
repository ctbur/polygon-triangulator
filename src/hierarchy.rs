use crate::{
    polygon::{self, Contour},
    regions::Island,
};

pub struct IslandNode {
    // runs in CCW order
    pub outline: Contour,
    pub interior: Vec<InteriorNode>,
}

pub struct InteriorNode {
    // runs in CW order
    pub region: Contour,
    pub children: Vec<IslandNode>,
}

/// Creates a hierarchy of regions that describes which regions contain which other ones.
/// `regions` is a sequence of CCW regions each followed by one or more CW regions.
/// The regions are provided in increasing min-y-coordinate.
pub fn build_region_hierarchy(islands: Vec<Island>) -> Vec<IslandNode> {
    let mut island_nodes = Vec::with_capacity(islands.len());
    for island in islands {
        let mut interior_nodes = Vec::with_capacity(island.interior.len());
        for region in island.interior {
            interior_nodes.push(InteriorNode {
                region,
                children: Vec::new(),
            })
        }
        island_nodes.push(IslandNode {
            outline: island.outline,
            interior: interior_nodes,
        })
    }

    build_hierarchy_from_island_nodes(&mut island_nodes);
    return island_nodes;
}

fn build_hierarchy_from_island_nodes(island_nodes: &mut Vec<IslandNode>) {
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
            build_hierarchy_from_island_nodes(&mut interior_node.children);
        }
    }
}

fn find_interior_node_containing_island_node(
    container: &IslandNode,
    containee: &IslandNode,
) -> Option<usize> {
    if !are_regions_nested(&container.outline, &containee.outline) {
        return None;
    }

    for (idx, interior_node) in container.interior.iter().enumerate() {
        if are_regions_nested(&interior_node.region, &containee.outline) {
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
