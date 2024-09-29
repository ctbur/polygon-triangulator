use std::mem;

use crate::polygon::{self, Contour};

pub struct OuterNode {
    // runs in CCW order
    pub outline: Contour,
    pub interior: Vec<InnerNode>,
}

pub struct InnerNode {
    // runs in CW order
    pub region: Contour,
    pub children: Vec<OuterNode>,
}

/// Creates a hierarchy of regions that describes which regions contain which other ones.
/// `regions` is a sequence of CCW regions each followed by one or more CW regions.
/// The regions are provided in increasing min-y-coordinate.
pub fn create_region_hierarchy(regions: Vec<Contour>) -> OuterNode {
    let mut outer_nodes = Vec::new();

    for region in regions {
        let is_ccw = polygon::calculate_region_area(&region) > 0.0;
        if is_ccw {
            outer_nodes.push(OuterNode {
                outline: region,
                interior: Vec::new(),
            })
        } else {
            let outer_node = outer_nodes.last_mut().expect("First region is not CCW");
            outer_node.interior.push(InnerNode {
                region,
                children: Vec::new(),
            })
        }
    }

    // check for all outlines if they are contained in any interior region

    panic!()
}

fn ligma(mut outer_nodes: Vec<OuterNode>) -> Vec<OuterNode> {
    let mut current = outer_nodes.remove(0);
    let mut outer_nodes_outside = Vec::new();

    // iterate over outer_nodes once in order to make sure the order is preserved
    for outer_node in outer_nodes {
        let mut new_parent_idx_opt = None;
        for (idx, inner_node) in current.interior.iter_mut().enumerate() {
            if are_regions_nested(&inner_node.region, &outer_node.outline) {
                new_parent_idx_opt = Some(idx);
            }
        }

        if let Some(new_parent_idx) = new_parent_idx_opt {
            current.interior[new_parent_idx].children.push(outer_node);
        } else {
            outer_nodes_outside.push(outer_node);
        }
    }

    for inner_node in &mut current.interior {
        let children = mem::replace(&mut inner_node.children, Vec::new());
        inner_node.children = ligma(children);
    }

    let mut result = vec![current];
    if !outer_nodes_outside.is_empty() {
        result.append(&mut ligma(outer_nodes_outside));
    }
    return result;
}

fn are_regions_nested(container: &Contour, containee: &Contour) -> bool {
    panic!()
}

#[cfg(test)]
mod tests {
    use crate::{regions, showcase};

    use super::*;

    #[test]
    fn todo() {}
}
