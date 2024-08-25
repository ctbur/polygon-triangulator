use std::{
    cmp::Ordering,
    collections::{HashMap, HashSet},
};

use crate::{
    polygon::{self, Contour},
    vector2::Vector2f,
};

fn comp_points_x_dir(a: Vector2f, b: Vector2f) -> Ordering {
    return f32::total_cmp(&a.x, &b.x).then_with(|| f32::total_cmp(&a.y, &b.y));
}

#[derive(Clone, Copy, Debug)]
struct Point {
    contour: usize,
    idx: usize,
}

#[derive(Debug)]
struct Edge {
    contour: usize,
    from_idx: usize,
    to_idx: usize,
}

struct SweepLine<'a> {
    subdivided_contours: &'a Vec<Contour>,
    edges: Vec<(Edge, i32)>,
}

impl<'a> SweepLine<'a> {
    fn start_region(&mut self, lower: Edge, upper: Edge, winding_number: i32) {
        debug_assert_eq!(lower.from_idx, upper.from_idx);

        let from_position = self.subdivided_contours[lower.contour][lower.from_idx];
        let (edge_idx, prev_winding_number) =
            if let Some(edge_idx) = self.find_edge_idx_above(from_position) {
                (edge_idx, self.edges[edge_idx].1)
            } else {
                (self.edges.len(), 0)
            };

        println!("Insert edges: {:?}, {:?}", lower, upper);

        self.edges.insert(edge_idx, (lower, winding_number));
        self.edges
            .insert(edge_idx + 1, (upper, prev_winding_number));
    }

    fn end_region(&mut self, to_point: Point) {
        let edge_idx = self.find_edge_by_to_point(to_point);

        // remove lower
        let (lower, _) = self.edges.remove(edge_idx);

        // remove upper
        let (upper, _) = self.edges.remove(edge_idx);
        debug_assert_eq!(upper.contour, to_point.contour);
        debug_assert_eq!(upper.to_idx, to_point.idx);

        println!("Remove edges: {:?}, {:?}", lower, upper);
    }

    fn update_edge(&mut self, to_point: Point, new_to_point_idx: usize) {
        let edge_idx = self.find_edge_by_to_point(to_point);

        // old `to` point becomes new `from` point
        self.edges[edge_idx].0.from_idx = to_point.idx;
        self.edges[edge_idx].0.to_idx = new_to_point_idx;
    }

    fn query_winding_number(&self, position: Vector2f) -> i32 {
        if let Some(edge_idx) = self.find_edge_idx_above(position) {
            return self.edges[edge_idx].1;
        } else {
            return 0;
        }
    }

    fn find_edge_by_to_point(&self, to_point: Point) -> usize {
        for (edge_idx, (edge, _)) in self.edges.iter().enumerate() {
            if edge.contour == to_point.contour && edge.to_idx == to_point.idx {
                return edge_idx;
            }
        }

        panic!("Edges: {:?}, to_point: {:?}", self.edges, to_point);
    }

    fn find_edge_idx_above(&self, position: Vector2f) -> Option<usize> {
        // go from highest edge to lowest, return on first positive y-dist
        for (idx, (edge, _)) in self.edges.iter().rev().enumerate() {
            let edge_from = self.subdivided_contours[edge.contour][edge.from_idx];
            let edge_to = self.subdivided_contours[edge.contour][edge.to_idx];

            // TODO: do we need to handle the case where the edge is vertical?
            let v = edge_to - edge_from;
            let y_dist = edge_from.y + v.y * (position.x - edge_from.x) / v.x - position.y;

            if y_dist > 0.0 {
                return Some(idx);
            }
        }

        return None;
    }
}

pub fn calculate_winding_numbers(
    subdivided_contours: &Vec<Contour>,
    regions: &Vec<Contour>,
) -> Vec<i32> {
    // TODO: find clear distinction between outline and interior regions
    let regions: Vec<_> = regions
        .iter()
        .filter(|r| polygon::calculate_region_area(r) < 0.0)
        .map(|r| r.clone())
        .collect();
    // TODO: remove
    for region in &regions {
        let mut points = HashSet::new();
        for &p in region {
            points.insert(Vector2f::to_bits(p));
        }
        assert_eq!(
            region.len(),
            points.len(),
            "Region really has duplicate points!"
        );
    }

    let total_points = subdivided_contours.iter().map(|c| c.len()).sum();
    let mut contour_points_order = vec![Point { contour: 0, idx: 0 }; total_points];

    let mut idx = 0;
    for (contour_idx, contour) in subdivided_contours.iter().enumerate() {
        for (point_idx, _) in contour.iter().enumerate() {
            contour_points_order[idx].contour = contour_idx;
            contour_points_order[idx].idx = point_idx;
            idx += 1;
        }
    }

    contour_points_order.sort_by(|id1, id2| {
        comp_points_x_dir(
            subdivided_contours[id1.contour][id1.idx],
            subdivided_contours[id2.contour][id2.idx],
        )
    });

    let mut sweep_line = SweepLine {
        subdivided_contours,
        edges: Vec::new(),
    };

    let mut region_lookup = HashMap::new();
    for region in &regions {
        for i in 0..region.len() {
            let (from, to) = (region[i], region[(i + 1) % region.len()]);
            // each region is uniquely identifiable by two consecutive points
            region_lookup.insert((Vector2f::to_bits(from), Vector2f::to_bits(to)), i);
        }
    }

    let mut region_winding_numbers = vec![None; regions.len()];

    for p in contour_points_order {
        let contour = &subdivided_contours[p.contour];
        let point = contour[p.idx];

        let prev_idx = (p.idx + contour.len() - 1) % contour.len();
        let prev = contour[prev_idx];

        let next_idx = (p.idx + 1) % contour.len();
        let next = contour[next_idx];

        println!(
            "\nPoint type: {:?}, p: {:?}, edges/regions: {:?}",
            categorize_point(prev, point, next),
            p,
            sweep_line.edges,
        );
        for c in contour {
            print!("({},{}),", c.x, c.y);
        }
        println!();
        match categorize_point(prev, point, next) {
            PointType::Start { downward } => {
                let winding_number = sweep_line.query_winding_number(point);
                let edge1 = Edge {
                    contour: p.contour,
                    from_idx: p.idx,
                    to_idx: prev_idx,
                };
                let edge2 = Edge {
                    contour: p.contour,
                    from_idx: p.idx,
                    to_idx: next_idx,
                };

                let (lower, upper, winding_number) = if downward {
                    (edge2, edge1, winding_number + 1)
                } else {
                    (edge1, edge2, winding_number - 1)
                };
                sweep_line.start_region(lower, upper, winding_number);

                // region runs in CW order
                if let Some(&region_idx) =
                    region_lookup.get(&(Vector2f::to_bits(point), Vector2f::to_bits(prev)))
                {
                    if let Some(wn) = region_winding_numbers[region_idx] {
                        debug_assert_eq!(wn, winding_number);
                    }
                    region_winding_numbers[region_idx] = Some(winding_number);
                }
            }
            PointType::End {} => {
                sweep_line.end_region(p);
            }
            PointType::Chain { upper } => {
                let new_to_point_idx = if upper { next_idx } else { prev_idx };
                sweep_line.update_edge(p, new_to_point_idx);
            }
        }
    }

    println!("Result: {:?}", region_winding_numbers);
    let region_winding_numbers_result: Vec<_> = region_winding_numbers
        .iter()
        .map(|wn| wn.unwrap())
        .collect();
    return region_winding_numbers_result;
}

#[derive(Debug)]
enum PointType {
    Start { downward: bool },
    End {},
    Chain { upper: bool },
}

fn categorize_point(prev: Vector2f, point: Vector2f, next: Vector2f) -> PointType {
    let prev_comp = comp_points_x_dir(point, prev);
    let next_comp = comp_points_x_dir(point, next);

    match (prev_comp, next_comp) {
        (Ordering::Less, Ordering::Less) => {
            if (prev - point).cross(next - point) > 0.0 {
                PointType::Start { downward: false }
            } else {
                PointType::Start { downward: true }
            }
        }
        (Ordering::Greater, Ordering::Greater) => PointType::End {},
        // assuming clockwise order
        (Ordering::Less, Ordering::Greater) => PointType::Chain { upper: false },
        (Ordering::Greater, Ordering::Less) => PointType::Chain { upper: true },
        (_, _) => panic!("Contour contains consecutive equal points"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interior_diamond() {
        let mut diamond = vec![
            Vector2f::new(1.0, 1.0),
            Vector2f::new(2.0, 2.0),
            Vector2f::new(3.0, 1.0),
            Vector2f::new(2.0, 0.0),
        ];

        let wns = calculate_winding_numbers(&vec![diamond.clone()], &vec![diamond.clone()]);
        assert_eq!(wns, vec![-1]);

        diamond.reverse();
        let wns = calculate_winding_numbers(
            &vec![diamond.clone(), diamond.clone(), diamond.clone()],
            &vec![diamond.clone()],
        );
        assert_eq!(wns, vec![3]);
    }
}
