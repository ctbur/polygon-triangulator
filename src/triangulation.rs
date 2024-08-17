use core::f32;
use std::cmp::Ordering;

use crate::polygon::Contour;
use crate::vector2::Vector2f;

pub type Triangle = [Vector2f; 3];

fn comp_points_x_dir(a: Vector2f, b: Vector2f) -> Ordering {
    return f32::total_cmp(&a.x, &b.x).then_with(|| f32::total_cmp(&a.y, &b.y));
}

/// Uses the shoelace formulate to calculate the signed polygon area.
/// The sign is positive if the region runs CCW.
fn calculate_region_area(region: &Contour) -> f32 {
    if region.len() < 3 {
        return 0.0;
    }

    let mut area = 0.0;

    for i in 0..region.len() - 1 {
        let p_i = region[i];
        let p_ni = region[i + 1];
        area += (p_i.y + p_ni.y) * (p_i.x - p_ni.x);
    }

    let p_i = region.last().unwrap();
    let p_ni = region[0];
    area += (p_i.y + p_ni.y) * (p_i.x - p_ni.x);

    return area / 2.0;
}

struct Edge {
    from: usize,
    to: usize,
    helper: usize,
}

struct SweepLineState<'a> {
    region: &'a Contour,
    edges: Vec<Edge>,
}

impl<'a> SweepLineState<'a> {
    fn insert_edge_pair(&mut self, lower: Edge, upper: Edge) {
        println!(
            "Insert edge pair: ({}->{}), ({}->{})",
            lower.from, lower.to, upper.from, upper.to
        );
        // find edge immediately above upper.from
        // insert below
        self.edges.push(lower);
        self.edges.push(upper);
    }

    fn delete_edge_pair(&mut self, to_point_idx: usize) -> (Edge, Edge) {
        let mut upper_idx = None;
        for (idx, edge) in self.edges.iter().enumerate() {
            if edge.to == to_point_idx {
                upper_idx = Some(idx);
            }
        }
        let upper = self.edges.remove(upper_idx.unwrap());

        let mut lower_idx = None;
        for (idx, edge) in self.edges.iter().enumerate() {
            if edge.to == to_point_idx {
                lower_idx = Some(idx);
            }
        }
        let lower = self.edges.remove(lower_idx.unwrap());

        return (upper, lower);
    }

    fn replace_edge(&mut self, to_point_idx: usize, new_edge: Edge) -> Edge {
        println!(
            "Replace edge to {} with ({}->{})",
            to_point_idx, new_edge.from, new_edge.to
        );
        let mut edge_idx = None;
        for (idx, edge) in self.edges.iter().enumerate() {
            if edge.to == to_point_idx {
                edge_idx = Some(idx);
            }
        }
        return std::mem::replace(&mut self.edges[edge_idx.unwrap()], new_edge);
    }

    fn find_edge_above(&mut self, point_idx: usize) -> &mut Edge {
        println!("Find edge above {}", point_idx);
        // find edge immediately above point_idx
        let mut immediate_upper_edge = None;
        let mut lowest_y_dist = f32::INFINITY;

        if self.edges.is_empty() {
            panic!();
        }

        for (idx, edge) in self.edges.iter_mut().enumerate() {
            println!("idx: {} ({}->{})", idx, edge.from, edge.to);
            // skip the lower edges
            if idx % 2 == 0 {
                continue;
            }

            // find distance from point to the end
            let edge_from = self.region[edge.from];
            let edge_to = self.region[edge.to];
            let point = self.region[point_idx];
            println!("{}, {}->{}", point, edge_from, edge_to);

            // TODO: do we need to handle the case where the edge is vertical?
            let v = edge_to - edge_from;
            let y_dist = edge_from.y + v.y * (point.x - edge_from.x) / v.x - point.y;

            println!("{} > 0.0 && {} < {}", y_dist, y_dist, lowest_y_dist);
            if y_dist > 0.0 && y_dist < lowest_y_dist {
                lowest_y_dist = y_dist;
                immediate_upper_edge = Some(edge);
            }
        }

        return immediate_upper_edge.unwrap();
    }
}

fn prev_point_idx(region: &Contour, point_idx: usize, ccw: bool) -> usize {
    return if ccw {
        (point_idx + region.len() - 1) % region.len()
    } else {
        (point_idx + 1) % region.len()
    };
}

fn next_point_idx(region: &Contour, point_idx: usize, ccw: bool) -> usize {
    return if ccw {
        (point_idx + 1) % region.len()
    } else {
        (point_idx + region.len() - 1) % region.len()
    };
}

/// Triangulates a region, i.e., a simple polygon (concave, no intersections)
/// Uses the plane-sweep algorithm described in https://www.cs.umd.edu/class/fall2021/cmsc754/Lects/lect05-triangulate.pdf
pub fn triangulate_region(region: &Contour) -> Vec<Triangle> {
    if region.len() < 3 {
        return Vec::new();
    }

    // check if the region is running CCW
    let ccw = calculate_region_area(region) > 0.0;

    // contains indices of region which are going to be sorted
    let mut points_order = Vec::from_iter(0..region.len());
    points_order.sort_by(|a, b| comp_points_x_dir(region[*a], region[*b]));

    for point_idx in &points_order {
        let point = region[*point_idx];
        println!("point: {}, idx: {}", point, point_idx);
    }

    // sweep across plane from x-direction
    let mut sweep_line = SweepLineState {
        region,
        edges: Vec::new(),
    };
    let mut i = 0;
    for point_idx in points_order {
        i += 1;

        let prev_idx = prev_point_idx(region, point_idx, ccw);
        let next_idx = next_point_idx(region, point_idx, ccw);
        println!(
            "\nindex: ccw: {} - prev: {}, point: {}, next: {}",
            ccw, prev_idx, point_idx, next_idx
        );

        let prev = region[prev_idx];
        let point = region[point_idx];
        let next = region[next_idx];
        println!(
            "coords: ccw: {} - prev: {}, point: {}, next: {}",
            ccw, prev, point, next
        );

        println!(
            "Looking at point {}, type {:?}, iteration {}",
            point_idx,
            categorize_point(prev, point, next),
            i,
        );

        match categorize_point(prev, point, next) {
            PointType::Start => {
                // insert both edges
                let upper = Edge {
                    from: point_idx,
                    to: prev_idx,
                    helper: point_idx,
                };
                let lower = Edge {
                    from: point_idx,
                    to: next_idx,
                    helper: point_idx, // TODO: check
                };
                sweep_line.insert_edge_pair(lower, upper);
            }
            PointType::End => {
                // delete upper ending at point_idx sweep line state
                // delete lower ending at point_idx from sweep line state
                let (upper, _) = sweep_line.delete_edge_pair(point_idx);
                fix_up(region, ccw, point_idx, &upper);
            }
            PointType::Split => {
                // top_edge = find edge immediately above v
                // diagonal connect point_idx to top_edge.helper
                let top_edge = sweep_line.find_edge_above(point_idx);
                split_polygon_across(point_idx, top_edge.helper);

                // insert both edges
                let upper = Edge {
                    from: point_idx,
                    to: next_idx,
                    helper: point_idx,
                };
                let lower = Edge {
                    from: point_idx,
                    to: prev_idx,
                    helper: point_idx,
                };
                sweep_line.insert_edge_pair(lower, upper);
            }
            PointType::Merge => {
                // _, lower = delete edges ending at point_idx
                let (_, lower) = sweep_line.delete_edge_pair(point_idx);

                // top_edge = find edge immediately above v
                let top_edge = sweep_line.find_edge_above(point_idx);

                fix_up(region, ccw, point_idx, top_edge);
                fix_up(region, ccw, point_idx, &lower);
                // top_edge.helper = point_idx - persist in state
                top_edge.helper = point_idx;
            }
            PointType::UpperChain => {
                // edge = replace edge ending at point_idx with Edge point_idx -> prev_idx, helper: point_idx
                let edge = Edge {
                    from: point_idx,
                    to: prev_idx,
                    helper: point_idx,
                };
                fix_up(region, ccw, point_idx, &edge);
                sweep_line.replace_edge(point_idx, edge);
            }
            PointType::LowerChain => {
                // edge = replace edge ending at point_idx with Edge point_idx -> next_idx, helper: point_idx
                let edge = Edge {
                    from: point_idx,
                    to: next_idx,
                    helper: point_idx,
                };
                sweep_line.replace_edge(point_idx, edge);
            }
        }
    }

    return Vec::new();
}

fn fix_up(region: &Contour, ccw: bool, point_idx: usize, edge: &Edge) {
    let prev = region[prev_point_idx(region, edge.helper, ccw)];
    let next = region[next_point_idx(region, edge.helper, ccw)];

    if categorize_point(prev, region[edge.helper], next) == PointType::Merge {
        split_polygon_across(edge.helper, point_idx);
    }
}

fn split_polygon_across(a: usize, b: usize) {
    println!("Splitting polygon with edge from {} to {}", a, b);
}

#[derive(Debug, PartialEq, Eq)]
enum PointType {
    Start,
    End,
    Split,
    Merge,
    UpperChain,
    LowerChain,
}

fn categorize_point(prev: Vector2f, point: Vector2f, next: Vector2f) -> PointType {
    let prev_comp = comp_points_x_dir(point, prev);
    let next_comp = comp_points_x_dir(point, next);

    return match (prev_comp, next_comp) {
        (Ordering::Less, Ordering::Less) => {
            if is_point_concave(prev, point, next) {
                PointType::Split
            } else {
                PointType::Start
            }
        }
        (Ordering::Greater, Ordering::Greater) => {
            if is_point_concave(prev, point, next) {
                PointType::Merge
            } else {
                PointType::End
            }
        }
        (Ordering::Greater, Ordering::Less) => PointType::LowerChain,
        (Ordering::Less, Ordering::Greater) => PointType::UpperChain,
        (_, _) => panic!("Region contains duplicate points"),
    };
}

fn is_point_concave(prev: Vector2f, point: Vector2f, next: Vector2f) -> bool {
    // we assume that the region is running CCW
    // thus the interior angle is on the left side of the above point chain
    return (prev - point).cross(next - point) > 0.0;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculate_region_area() {
        let mut region = vec![
            Vector2f::new(1.0, 1.0),
            Vector2f::new(2.0, 2.0),
            Vector2f::new(1.0, 3.0),
            Vector2f::new(0.0, 2.0),
        ];

        let area = calculate_region_area(&region);
        assert_eq!(area, 2.0);

        region.reverse();

        let rev_area = calculate_region_area(&region);
        assert_eq!(rev_area, -2.0);
    }

    #[test]
    fn test_categorize_point() {
        let top = Vector2f::new(5.0, 11.0);
        let left = Vector2f::new(4.0, 10.0);
        let bottom = Vector2f::new(5.0, 9.0);
        let right = Vector2f::new(6.0, 10.0);

        assert_eq!(categorize_point(top, left, bottom), PointType::Start);
        assert_eq!(categorize_point(bottom, right, top), PointType::End);
        assert_eq!(categorize_point(bottom, left, top), PointType::Split);
        assert_eq!(categorize_point(top, right, bottom), PointType::Merge);
        assert_eq!(categorize_point(right, top, left), PointType::UpperChain);
        assert_eq!(categorize_point(left, bottom, right), PointType::LowerChain);
    }

    #[test]
    fn test_diamond() {
        let diamond = vec![
            Vector2f::new(1.0, 1.0),
            Vector2f::new(2.0, 2.0),
            Vector2f::new(1.0, 3.0),
            Vector2f::new(0.0, 2.0),
        ];
        triangulate_region(&diamond);
    }

    #[test]
    fn test_rect() {
        let rect = vec![
            Vector2f::new(1.0, 1.0),
            Vector2f::new(1.0, 2.0),
            Vector2f::new(2.0, 2.0),
            Vector2f::new(2.0, 1.0),
        ];
        triangulate_region(&rect);
    }

    #[test]
    fn test_arrow() {
        let arrow = vec![
            Vector2f::new(0.0, 20.0),
            Vector2f::new(5.0, 15.0),
            Vector2f::new(0.0, 10.0),
            Vector2f::new(20.0, 15.0),
        ];
        triangulate_region(&arrow);
    }

    #[test]
    fn test_hourglass() {
        let hourglass = vec![
            Vector2f::new(0.0, 0.0),
            Vector2f::new(10.0, 0.0),
            Vector2f::new(8.0, 5.0),
            Vector2f::new(10.0, 10.0),
            Vector2f::new(0.0, 10.0),
            Vector2f::new(2.0, 5.0),
        ];
        triangulate_region(&hourglass);
    }
}
