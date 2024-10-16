use crate::{
    debugp,
    graph::Graph,
    polygon::{self, Contour},
    vector2::{self, Vector2f},
};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::ops;

const DBG_WN: &str = "winding numbers";

pub enum WindingRule {
    Odd,
    NonZero,
    Positive,
    Negative,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FillMode {
    /// If a region describes an area inside the polygon
    Filled,
    /// If a region describes an area outside the polygon
    Hole,
}

impl ops::Not for FillMode {
    type Output = Self;

    fn not(self) -> Self::Output {
        match self {
            FillMode::Hole => FillMode::Filled,
            FillMode::Filled => FillMode::Hole,
        }
    }
}

impl WindingRule {
    fn get_fill_mode(&self, winding_number: i32) -> FillMode {
        let inside = match self {
            Self::Odd => winding_number % 2 == 1,
            Self::NonZero => winding_number != 0,
            Self::Positive => winding_number > 0,
            Self::Negative => winding_number < 0,
        };

        return if inside {
            FillMode::Filled
        } else {
            FillMode::Hole
        };
    }
}

struct RegionEntry {
    upper_from: Vector2f,
    upper_to: Vector2f,
    island: usize,
    region: usize,
}

#[derive(Debug, Eq, PartialEq, Clone, Copy, Hash)]
struct PointId {
    contour: usize,
    point: usize,
}

struct Edge {
    contour: usize,
    from: usize,
    to: usize,
}

struct SweepState<'a> {
    contours: &'a [Contour],
    edges: Vec<Edge>,
}

impl<'a> SweepState<'a> {
    fn start_edge_pair(&mut self, from: PointId) {
        self.edges.push(Edge {
            contour: from.contour,
            from: from.point,
            to: (from.point + 1) % self.contours[from.contour].len(),
        });
        self.edges.push(Edge {
            contour: from.contour,
            from: from.point,
            to: (from.point + self.contours[from.contour].len() - 1)
                % self.contours[from.contour].len(),
        });
        debugp!(
            DBG_WN,
            "Contour {}: insert edges {}->{}, {}->{}",
            from.contour,
            from.point,
            (from.point + 1) % self.contours[from.contour].len(),
            from.point,
            (from.point + self.contours[from.contour].len() - 1)
                % self.contours[from.contour].len()
        );
    }

    fn end_edge_pair(&mut self, to: PointId) {
        let mut indices_to_delete = Vec::new();
        for (idx, edge) in self.edges.iter().enumerate() {
            if edge.contour == to.contour && edge.to == to.point {
                indices_to_delete.push(idx);
            }
        }

        for &idx in &indices_to_delete {
            debugp!(
                DBG_WN,
                "Contour {}: delete edge {}->{}",
                self.edges[idx].contour,
                self.edges[idx].from,
                self.edges[idx].to,
            );
        }

        debug_assert!(indices_to_delete.len() == 2);

        indices_to_delete.reverse();
        for i in indices_to_delete {
            self.edges.remove(i);
        }
    }

    fn advance_chain(&mut self, to: PointId) {
        let mut edge_idx_opt = None;
        for (idx, edge) in self.edges.iter().enumerate() {
            if edge.contour == to.contour && edge.to == to.point {
                edge_idx_opt = Some(idx);
                break;
            }
        }

        let edge_idx = edge_idx_opt.unwrap();
        let edge = &mut self.edges[edge_idx];
        let (old_from, old_to) = (edge.from, edge.to);

        let contour_len = self.contours[edge.contour].len();
        let point_prev = (to.point + contour_len - 1) % contour_len;
        let point_next = (to.point + 1) % contour_len;

        // "to" of new edge becomes next point in contour that "from" of old edge was not
        edge.to = if edge.from == point_prev {
            point_next
        } else if edge.from == point_next {
            point_prev
        } else {
            panic!(
                "Edge {}->{} doesn't use consecutive points in contour",
                old_from, old_to
            )
        };

        // "from" of new edge becomes "to" of old edge
        edge.from = to.point;

        debugp!(
            DBG_WN,
            "Contour {}: replace edge {}->{} with {}->{}",
            edge.contour,
            old_from,
            old_to,
            edge.from,
            edge.to
        );
    }

    fn calculate_winding_number(&self, upper_from: Vector2f, upper_to: Vector2f) -> i32 {
        let mut winding_number = 0;

        for edge in &self.edges {
            let edge_from = self.contours[edge.contour][edge.from];
            let edge_to = self.contours[edge.contour][edge.to];
            if self.is_edge_equal_or_higher(edge_from, edge_to, upper_from, upper_to) {
                // indices in contour direction
                let (contour_from, contour_to) = if edge.from < edge.to {
                    (edge.from, edge.to)
                } else {
                    (edge.to, edge.from)
                };
                let is_cw = vector2::comp_points_x_dir(
                    self.contours[edge.contour][contour_from],
                    self.contours[edge.contour][contour_to],
                )
                .is_lt();
                winding_number += if is_cw { -1 } else { 1 };
            }
        }

        return winding_number;
    }

    fn is_edge_equal_or_higher(
        &self,
        edge_from: Vector2f,
        edge_to: Vector2f,
        base_from: Vector2f,
        base_to: Vector2f,
    ) -> bool {
        if edge_from == base_from {
            return (base_to - base_from).cross(edge_to - edge_from) >= 0.0;
        }

        // TODO: do we need to handle the case where the edge is vertical?
        let v = edge_to - edge_from;
        let y_dist = edge_from.y + v.y * (base_from.x - edge_from.x) / v.x - base_from.y;
        return y_dist > 0.0;
    }
}

fn calculate_regions_inside(
    islands: &mut [Island],
    contours: &[Contour],
    winding_rule: WindingRule,
) {
    debugp!(DBG_WN, "# Start calculate_regions_inside");

    // find leftmost point for each region and put in lookup map
    let mut region_by_min_x = HashMap::new();
    for (island_idx, island) in islands.iter().enumerate() {
        for (region_idx, region) in island.interior.iter().enumerate() {
            let min_x_idx = (0..region.len())
                .min_by(|&i, &j| vector2::comp_points_x_dir(region[i], region[j]))
                .unwrap();

            let region_entry = RegionEntry {
                upper_from: region[min_x_idx],
                upper_to: region[(min_x_idx + 1) % region.len()],
                island: island_idx,
                region: region_idx,
            };
            region_by_min_x.insert(Vector2f::to_bits(region[min_x_idx]), region_entry);
        }
    }

    // create segment IDs
    let mut point_ids = Vec::with_capacity(contours.iter().map(|c| c.len()).sum());
    for (contour_idx, contour) in contours.iter().enumerate() {
        for point_idx in 0..contour.len() {
            point_ids.push(PointId {
                contour: contour_idx,
                point: point_idx,
            });
        }
    }
    // sort by x-direction
    point_ids.sort_by(|a, b| {
        vector2::comp_points_x_dir(contours[a.contour][a.point], contours[b.contour][b.point])
    });

    // run plane sweep to calculate winding numbers
    let mut sweep_state = SweepState {
        contours: &contours,
        edges: Vec::new(),
    };
    let mut prev_point_opt = None;
    for point_id in point_ids {
        let contour = &contours[point_id.contour];

        let prev = contour[(point_id.point + contour.len() - 1) % contour.len()];
        let point = contour[point_id.point];
        let next = contour[(point_id.point + 1) % contour.len()];

        if let Some(prev_point) = prev_point_opt {
            // if we just handled all of the events at a given point
            if prev_point != point {
                // and if there is a region start there
                if let Some(ref region_entry) = region_by_min_x.get(&Vector2f::to_bits(prev_point))
                {
                    // then we calculate the winding number
                    let winding_number = sweep_state
                        .calculate_winding_number(region_entry.upper_from, region_entry.upper_to);
                    islands[region_entry.island].interior_fill_mode[region_entry.region] =
                        Some(winding_rule.get_fill_mode(winding_number));
                    debugp!(
                        DBG_WN,
                        "Region ({}, {}) has winding number {} and is inside: {:?}",
                        region_entry.island,
                        region_entry.region,
                        winding_number,
                        winding_rule.get_fill_mode(winding_number)
                    );
                }
            }
        }
        prev_point_opt = Some(point);

        // update sweep state
        match categorize_point(prev, point, next) {
            PointType::Start => sweep_state.start_edge_pair(point_id),
            PointType::End => sweep_state.end_edge_pair(point_id),
            PointType::Chain => sweep_state.advance_chain(point_id),
        };
    }

    for island in &mut *islands {
        println!("Fill mode: {:?}", island.interior_fill_mode);
    }
    debug_assert!(islands
        .iter()
        .all(|i| i.interior_fill_mode.iter().all(|o| o.is_some())));
    debugp!(DBG_WN, "# End calculate_regions_inside");
}

#[derive(Debug, PartialEq, Eq)]
enum PointType {
    Start,
    End,
    Chain,
}

fn categorize_point(prev: Vector2f, point: Vector2f, next: Vector2f) -> PointType {
    let prev_comp = vector2::comp_points_x_dir(point, prev);
    let next_comp = vector2::comp_points_x_dir(point, next);

    return match (prev_comp, next_comp) {
        (Ordering::Less, Ordering::Less) => PointType::Start,
        (Ordering::Greater, Ordering::Greater) => PointType::End,
        (Ordering::Greater, Ordering::Less) => PointType::Chain,
        (Ordering::Less, Ordering::Greater) => PointType::Chain,
        (_, _) => panic!("Region contains repeated points"),
    };
}

struct Island {
    graph: Graph,
    // runs CCW
    outline: Contour,
    // runs CW
    interior: Vec<Contour>,
    interior_fill_mode: Vec<Option<FillMode>>,
}

struct IslandNode {
    island: Island,
    interior: Vec<InteriorNode>,
}

struct InteriorNode {
    children: Vec<IslandNode>,
}

/// Creates a hierarchy of regions that describes which regions contain which other ones.
/// `regions` is a sequence of CCW regions each followed by one or more CW regions.
/// The regions are provided in increasing min-y-coordinate.
fn build_region_hierarchy(islands: Vec<Island>) -> Vec<IslandNode> {
    let mut island_nodes = Vec::with_capacity(islands.len());
    for island in islands {
        let mut interior_nodes = Vec::with_capacity(island.interior.len());
        for _ in 0..island.interior.len() {
            interior_nodes.push(InteriorNode {
                children: Vec::new(),
            })
        }
        island_nodes.push(IslandNode {
            island,
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
    if !are_regions_nested(&container.island.outline, &containee.island.outline) {
        return None;
    }

    for idx in 0..container.interior.len() {
        let region = &container.island.interior[idx];
        let outline = &containee.island.outline;
        if are_regions_nested(region, outline) {
            return Some(idx);
        }
    }

    panic!("Node is contained by outline but not by interior region");
}

fn are_regions_nested(container: &Contour, containee: &Contour) -> bool {
    return polygon::calculate_winding_number(container, containee[0]) != 0;
}

fn decompose_hierarchy(island_node: IslandNode) {
    // for each region that is filled
    // find all interior regions that are not visible
    // get the outline of those
}

pub fn decompose(graph: Graph, subdivided_contours: &[Contour], winding_rule: WindingRule) {
    let island_graphs = graph.clone().split_by_islands();
    let mut islands: Vec<_> = island_graphs
        .into_iter()
        .map(|mut graph| {
            let (outline, interior) = graph.trace_regions();
            Island {
                graph,
                outline,
                interior_fill_mode: vec![None; interior.len()],
                interior,
            }
        })
        .collect();

    calculate_regions_inside(&mut islands, &subdivided_contours, winding_rule);

    let hierarchy = build_region_hierarchy(islands);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wn_triangle() {
        let contour = vec![
            Vector2f::new(0.0, 0.0),
            Vector2f::new(2.0, 0.0),
            Vector2f::new(1.0, -1.0),
        ];

        let mut outline = contour.clone();
        outline.reverse();

        let island = Island {
            outline,
            interior: vec![contour.clone()],
        };

        calculate_regions_inside(&[island], &[contour], WindingRule::Odd);
    }
}
