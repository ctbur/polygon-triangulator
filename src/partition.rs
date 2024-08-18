use core::f32;
use std::cmp::Ordering;

use crate::graph::Graph;
use crate::polygon::{self, Contour};
use crate::vector2::Vector2f;

fn comp_points_x_dir(a: Vector2f, b: Vector2f) -> Ordering {
    return f32::total_cmp(&a.x, &b.x).then_with(|| f32::total_cmp(&a.y, &b.y));
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

    fn delete_edge_pair(
        &mut self,
        to_point_idx: usize,
        is_merge_point: bool,
        ccw: bool,
    ) -> (Edge, Edge) {
        println!("Delete edge pair to {}", to_point_idx);

        let mut first_idx = None;
        for (idx, edge) in self.edges.iter().enumerate() {
            if edge.to == to_point_idx {
                first_idx = Some(idx);
            }
        }
        let first = self.edges.remove(first_idx.unwrap());

        let mut second_idx = None;
        for (idx, edge) in self.edges.iter().enumerate() {
            if edge.to == to_point_idx {
                second_idx = Some(idx);
            }
        }
        let second = self.edges.remove(second_idx.unwrap());

        let (predecessor, successor) = if first.from == (to_point_idx + 1) % self.region.len() {
            (second, first)
        } else {
            (first, second)
        };

        // ccw: merge -> successor is lower
        // ccw: end -> predecessor is upper
        let tmp = if ccw && is_merge_point {
            (successor, predecessor)
        } else {
            (predecessor, successor)
        };

        println!(
            "Lower: ({}->{}), Upper. ({}->{})",
            tmp.0.from, tmp.0.to, tmp.1.from, tmp.1.to
        );

        return tmp;
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
        // find edge immediately above point_idx
        let mut immediate_upper_edge = None;
        let mut lowest_y_dist = f32::INFINITY;

        if self.edges.is_empty() {
            panic!();
        }

        for edge in self.edges.iter_mut() {
            // find distance from point to the end
            let edge_from = self.region[edge.from];
            let edge_to = self.region[edge.to];
            let point = self.region[point_idx];

            // TODO: do we need to handle the case where the edge is vertical?
            let v = edge_to - edge_from;
            let y_dist = edge_from.y + v.y * (point.x - edge_from.x) / v.x - point.y;

            if y_dist > 0.0 && y_dist < lowest_y_dist {
                lowest_y_dist = y_dist;
                immediate_upper_edge = Some(edge);
            }
        }
        println!(
            "Find edge above {}: ({}->{})",
            point_idx,
            immediate_upper_edge.as_ref().unwrap().from,
            immediate_upper_edge.as_ref().unwrap().to
        );
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

/// Partitions a region, i.e., a simple polygon (concave, no intersections)
/// Uses the plane-sweep algorithm described in https://www.cs.umd.edu/class/fall2021/cmsc754/Lects/lect05-triangulate.pdf
pub fn partition_region(graph: &mut Graph, region: &Contour) {
    if region.len() < 3 {
        return;
    }

    // check if the region is running CCW
    let ccw = polygon::calculate_region_area(region) > 0.0;

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

    for &point_idx in &points_order {
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
            "Looking at point {}, type {:?}",
            point_idx,
            categorize_point(prev, point, next),
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
                let (_, upper) = sweep_line.delete_edge_pair(point_idx, false, ccw);
                fix_up(graph, region, ccw, point_idx, &upper);
            }
            PointType::Split => {
                // top_edge = find edge immediately above v
                // diagonal connect point_idx to top_edge.helper
                let top_edge = sweep_line.find_edge_above(point_idx);
                println!(
                    "Splitting polygon with edge from {} to {}",
                    point_idx, top_edge.helper
                );
                graph.insert_segment(region[point_idx], region[top_edge.helper]);

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
                let (lower, _) = sweep_line.delete_edge_pair(point_idx, true, ccw);

                // top_edge = find edge immediately above v
                let top_edge = sweep_line.find_edge_above(point_idx);

                fix_up(graph, region, ccw, point_idx, top_edge);
                fix_up(graph, region, ccw, point_idx, &lower);
                // top_edge.helper = point_idx - persist in state
                println!(
                    "Set edge ({}->{}) helper to {}",
                    top_edge.from, top_edge.to, point_idx
                );
                top_edge.helper = point_idx;
            }
            PointType::UpperChain => {
                let edge = Edge {
                    from: point_idx,
                    to: prev_idx,
                    helper: point_idx,
                };
                let prev_edge = sweep_line.replace_edge(point_idx, edge);
                fix_up(graph, region, ccw, point_idx, &prev_edge);
            }
            PointType::LowerChain => {
                let top_edge = sweep_line.find_edge_above(point_idx);
                fix_up(graph, region, ccw, point_idx, &top_edge);
                let edge = Edge {
                    from: point_idx,
                    to: next_idx,
                    helper: point_idx,
                };
                sweep_line.replace_edge(point_idx, edge);
            }
        }
    }
}

fn fix_up(graph: &mut Graph, region: &Contour, ccw: bool, point_idx: usize, edge: &Edge) {
    let prev = region[prev_point_idx(region, edge.helper, ccw)];
    let next = region[next_point_idx(region, edge.helper, ccw)];

    println!(
        "Fix up: P={}, E=({}->{}), H={}",
        point_idx, edge.from, edge.to, edge.helper
    );
    if categorize_point(prev, region[edge.helper], next) == PointType::Merge {
        println!(
            "Splitting polygon with edge from {} to {}",
            edge.helper, point_idx
        );
        graph.insert_segment(region[edge.helper], region[point_idx]);
    }
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
            if is_point_reflex(prev, point, next) {
                PointType::Split
            } else {
                PointType::Start
            }
        }
        (Ordering::Greater, Ordering::Greater) => {
            if is_point_reflex(prev, point, next) {
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

fn is_point_reflex(prev: Vector2f, point: Vector2f, next: Vector2f) -> bool {
    // we assume that the region is running CCW
    // thus the interior angle is on the left side of the above point chain
    return (prev - point).cross(next - point) > 0.0;
}

#[cfg(test)]
mod tests {
    use crate::{regions, showcase};

    use super::*;

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

    fn get_graph(region: &Contour) -> Graph {
        let mut graph = Graph::new(0.01);
        for i in 0..region.len() {
            graph.insert_segment(region[i], region[(i + 1) % region.len()]);
        }
        return graph;
    }

    fn validate_partitioning(graph: Graph) {
        let monotone_polygons = regions::trace_regions(graph);

        // check if polygons really are x-monotone
        for region in monotone_polygons.iter().skip(1) {
            assert!(polygon::is_region_x_monotone(region))
        }
    }

    #[test]
    fn test_partition_diamond() {
        let diamond = vec![
            Vector2f::new(1.0, 1.0),
            Vector2f::new(2.0, 2.0),
            Vector2f::new(1.0, 3.0),
            Vector2f::new(0.0, 2.0),
        ];
        let mut graph = get_graph(&diamond);
        partition_region(&mut graph, &diamond);
        validate_partitioning(graph);
    }

    #[test]
    fn test_partition_rect() {
        let rect = vec![
            Vector2f::new(1.0, 1.0),
            Vector2f::new(1.0, 2.0),
            Vector2f::new(2.0, 2.0),
            Vector2f::new(2.0, 1.0),
        ];
        let mut graph = get_graph(&rect);
        partition_region(&mut graph, &rect);
        validate_partitioning(graph);
    }

    #[test]
    fn test_partition_arrow() {
        let arrow = vec![
            Vector2f::new(0.0, 20.0),
            Vector2f::new(5.0, 15.0),
            Vector2f::new(0.0, 10.0),
            Vector2f::new(20.0, 15.0),
        ];
        let mut graph = get_graph(&arrow);
        partition_region(&mut graph, &arrow);
        validate_partitioning(graph);
    }

    #[test]
    fn test_partition_hourglass() {
        let hourglass = vec![
            Vector2f::new(0.0, 0.0),
            Vector2f::new(10.0, 0.0),
            Vector2f::new(8.0, 5.0),
            Vector2f::new(10.0, 10.0),
            Vector2f::new(0.0, 10.0),
            Vector2f::new(2.0, 5.0),
        ];
        let mut graph = get_graph(&hourglass);
        partition_region(&mut graph, &hourglass);
        validate_partitioning(graph);
    }

    #[test]
    fn test_partition_contour() {
        let contour = vec![
            Vector2f::new(400.0, 450.0),
            Vector2f::new(400.0, 700.0),
            Vector2f::new(700.0, 700.0),
            Vector2f::new(700.0, 400.0),
            Vector2f::new(450.0, 400.0),
            Vector2f::new(450.0, 450.0),
        ];
        let mut graph = get_graph(&contour);
        partition_region(&mut graph, &contour);
        validate_partitioning(graph);
    }

    #[test]
    fn test_partition_star() {
        let star = showcase::get_star().contours()[0].clone();
        let mut graph = get_graph(&star);
        partition_region(&mut graph, &star);
        validate_partitioning(graph);
    }
}
