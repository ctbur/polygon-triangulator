use crate::debugp;
use crate::{
    graph::Island,
    graph::ProximityMerger,
    polygon::Contour,
    vector2::{self, Vector2f},
};
use std::cmp::Ordering;
use std::collections::HashMap;

const DBG_WN: &str = "winding numbers";

pub enum WindingRule {
    Odd,
    NonZero,
    Positive,
    Negative,
}

impl WindingRule {
    fn is_inside(&self, winding_number: i32) -> bool {
        match self {
            Self::Odd => winding_number % 2 == 1,
            Self::NonZero => winding_number != 0,
            Self::Positive => winding_number > 0,
            Self::Negative => winding_number < 0,
        }
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

pub struct RegionId {
    island: usize,
    region: usize,
}

pub fn calculate_regions_inside(
    islands: &[Island],
    contours: &[Contour],
    // TODO: move this elsewhere
    proximity_merger: &mut ProximityMerger,
    winding_rule: WindingRule,
) -> Vec<RegionId> {
    debugp!(DBG_WN, "# Start calculate_regions_inside");

    // align contour with region points
    let mut aligned_contours = Vec::with_capacity(contours.len());
    for contour in contours {
        let mut aligned_contour = Vec::new();

        for &point in contour {
            let aligned_point = proximity_merger.map(point);

            // don't add duplicate points
            if let Some(&prev) = aligned_contour.last() {
                if aligned_point == prev {
                    continue;
                }
            }

            aligned_contour.push(aligned_point);
        }

        // ensure start and end are not duplicates either
        while aligned_contour[0] == *aligned_contour.last().unwrap() {
            aligned_contour.pop();
        }

        aligned_contours.push(aligned_contour);
    }

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
    let mut point_ids = Vec::with_capacity(aligned_contours.iter().map(|c| c.len()).sum());
    for (contour_idx, contour) in aligned_contours.iter().enumerate() {
        for point_idx in 0..contour.len() {
            point_ids.push(PointId {
                contour: contour_idx,
                point: point_idx,
            });
        }
    }
    // sort by x-direction
    point_ids.sort_by(|a, b| {
        vector2::comp_points_x_dir(
            aligned_contours[a.contour][a.point],
            aligned_contours[b.contour][b.point],
        )
    });

    // run plane sweep to calculate winding numbers
    let mut regions_inside = Vec::new();
    let mut sweep_state = SweepState {
        contours: &aligned_contours,
        edges: Vec::new(),
    };
    let mut prev_point_opt = None;
    for point_id in point_ids {
        let contour = &aligned_contours[point_id.contour];

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
                    let inside = winding_rule.is_inside(winding_number);
                    if inside {
                        regions_inside.push(RegionId {
                            island: region_entry.island,
                            region: region_entry.region,
                        });
                    }
                    debugp!(
                        DBG_WN,
                        "Region ({}, {}) has winding number {} and is inside: {}",
                        region_entry.island,
                        region_entry.region,
                        winding_number,
                        inside
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

    debugp!(DBG_WN, "# End calculate_regions_inside");
    return regions_inside;
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

fn find_center_of_ear(region: Contour) -> Vector2f {
    for i in 0..region.len() {
        let (a, b, c) = (i, (i + 1) % region.len(), (i + 2) % region.len());

        let mut ear = true;
        for (j, &p) in region.iter().enumerate() {
            if j == a || j == b || j == c {
                continue;
            }

            if is_point_in_triangle(p, [region[a], region[b], region[c]]) {
                ear = false;
                break;
            }
        }

        if ear {
            return (region[a] + region[b] + region[c]) / 3.0;
        }
    }

    // mathematically there should always be an ear
    panic!()
}

fn is_point_in_triangle(p: Vector2f, t: [Vector2f; 3]) -> bool {
    let e0 = (t[1] - t[0]).cross(p);
    let e1 = (t[2] - t[1]).cross(p);
    let e2 = (t[0] - t[2]).cross(p);

    return e0 == e1 && e1 == e2;
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

        let mut proximity_merger = ProximityMerger::new(0.001);

        calculate_regions_inside(
            &[island],
            &[contour],
            &mut proximity_merger,
            WindingRule::Odd,
        );
    }
}
