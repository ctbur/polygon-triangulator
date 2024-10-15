use std::collections::HashMap;
use std::mem;

use crate::polygon::{Contour, Polygon};
use crate::vector2::Vector2f;

#[derive(Debug, Eq, PartialEq, Clone, Copy, Hash)]
pub struct SegmentId {
    pub contour: usize,
    pub segment: usize,
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum SegmentIntersection {
    None,
    Point(Vector2f),
    Segment(Vector2f, Vector2f),
}

pub fn find_intersections(
    polygon: &Polygon,
    epsilon: f32,
) -> HashMap<SegmentId, Vec<SegmentIntersection>> {
    let mut intersections = HashMap::new();

    let contours = polygon.contours();
    for i in 0..contours.len() {
        for j in i..contours.len() {
            let contour = &contours[i];
            let other_contour = &contours[j];

            for k in 0..contour.len() {
                let p0 = contour[k];
                let p1 = contour[(k + 1) % contour.len()];

                for l in 0..other_contour.len() {
                    if i == j && k == l {
                        continue;
                    }

                    let p2 = other_contour[l];
                    let p3 = other_contour[(l + 1) % other_contour.len()];

                    match intersect_segments(p0, p1, p2, p3, epsilon) {
                        SegmentIntersection::None => {}
                        intersection => {
                            let seg_id1 = SegmentId {
                                contour: i,
                                segment: k,
                            };
                            intersections
                                .entry(seg_id1)
                                .or_insert_with(Vec::new)
                                .push(intersection);

                            let seg_id2 = SegmentId {
                                contour: j,
                                segment: l,
                            };
                            intersections
                                .entry(seg_id2)
                                .or_insert_with(Vec::new)
                                .push(intersection);
                        }
                    }
                }
            }
        }
    }

    intersections
}

fn intersect_segments(
    mut p0: Vector2f,
    mut p1: Vector2f,
    mut p2: Vector2f,
    mut p3: Vector2f,
    epsilon: f32,
) -> SegmentIntersection {
    // function uses parametric equations to find the intersection point
    let div = (p0.x - p1.x) * (p2.y - p3.y) - (p0.y - p1.y) * (p2.x - p3.x);

    if div.abs() < epsilon {
        // lines are parallel or degenerate
        let is_point_0 = points_are_equal(p0, p1, epsilon);
        let is_point_1 = points_are_equal(p2, p3, epsilon);

        if is_point_0 && is_point_1 {
            // degenerate case 1: both lines are points
            if (p0.x - p2.x).abs() < epsilon && (p0.y - p2.y).abs() < epsilon {
                return SegmentIntersection::Point(p0);
            } else {
                return SegmentIntersection::None;
            }
        } else if is_point_0 {
            // degenerate case 2: one of the lines is a point
            if intersect_segment_and_point(p2, p3, p0, epsilon) {
                return SegmentIntersection::Point(p0);
            } else {
                return SegmentIntersection::None;
            }
        } else if is_point_1 {
            // degenerate case 2: one of the lines is a point
            if intersect_segment_and_point(p0, p1, p2, epsilon) {
                return SegmentIntersection::Point(p2);
            } else {
                return SegmentIntersection::None;
            }
        }

        // parallel case
        // pick shorter line segment and check both endpoints

        // p0 and p1 should be the shorter line segment
        if (p0 - p1).length_squared() > (p2 - p3).length_squared() {
            // swap p0 and p1 with p2 and p3
            mem::swap(&mut p0, &mut p2);
            mem::swap(&mut p1, &mut p3);
        }

        // both segments should point in the same direction
        if (p1 - p0).dot(p3 - p2) < 0.0 {
            mem::swap(&mut p2, &mut p3);
        }

        let p0_on_segment = intersect_segment_and_point(p2, p3, p0, epsilon);
        let p1_on_segment = intersect_segment_and_point(p2, p3, p1, epsilon);

        if p0_on_segment && p1_on_segment {
            return SegmentIntersection::Segment(p0, p1);
        } else if p0_on_segment {
            if points_are_equal(p0, p3, epsilon) {
                return SegmentIntersection::Point(p0);
            }
            return SegmentIntersection::Segment(p0, p3);
        } else if p1_on_segment {
            if points_are_equal(p2, p1, epsilon) {
                return SegmentIntersection::Point(p1);
            }
            return SegmentIntersection::Segment(p2, p1);
        }

        return SegmentIntersection::None;
    } else {
        // general case
        let t1 = ((p0.x - p2.x) * (p2.y - p3.y) - (p0.y - p2.y) * (p2.x - p3.x)) / div;
        let t2 = ((p0.y - p1.y) * (p0.x - p2.x) - (p0.x - p1.x) * (p0.y - p2.y)) / div;

        if t1 >= 0.0 && t1 <= 1.0 && t2 >= 0.0 && t2 <= 1.0 {
            // calculate the intersection point using parametric equation of the first line
            return SegmentIntersection::Point(Vector2f::new(
                p0.x + t1 * (p1.x - p0.x),
                p0.y + t1 * (p1.y - p0.y),
            ));
        } else {
            return SegmentIntersection::None;
        }
    }
}

fn points_are_equal(p0: Vector2f, p1: Vector2f, epsilon: f32) -> bool {
    (p0 - p1).length_squared() < epsilon * epsilon
}

fn intersect_segment_and_point(p0: Vector2f, p1: Vector2f, point: Vector2f, epsilon: f32) -> bool {
    // check if the points are all aligned
    let cross = (p1 - p0).cross(point - p0);
    if cross.abs() > epsilon {
        return false;
    }

    // check if the point is inside the segment
    let dot = (p1 - p0).dot(point - p0);
    if dot < 0.0 {
        return false;
    }

    if dot > (p1 - p0).length_squared() {
        return false;
    }

    return true;
}

/// Struct used to merge points that are closer than epsilon into a single point.
/// This is used that mitigate errors from floating-point calculations.
struct ProximityMerger {
    epsilon: f32,
    encountered_points: Vec<Vector2f>,
}

impl ProximityMerger {
    pub fn new(epsilon: f32) -> ProximityMerger {
        ProximityMerger {
            epsilon,
            encountered_points: Vec::new(),
        }
    }

    fn map(&mut self, point: Vector2f) -> Vector2f {
        for &encountered_point in self.encountered_points.iter() {
            if (point - encountered_point).length_squared() <= self.epsilon * self.epsilon {
                return encountered_point;
            }
        }

        self.encountered_points.push(point);
        return point;
    }
}

fn strip_intersection_type(
    intersections: &HashMap<SegmentId, Vec<SegmentIntersection>>,
) -> HashMap<SegmentId, Vec<Vector2f>> {
    let mut stripped = HashMap::new();

    for (segment_id, intersections) in intersections.iter() {
        let mut stripped_intersections = Vec::new();
        for intersection in intersections.iter() {
            match intersection {
                SegmentIntersection::Point(p) => stripped_intersections.push(*p),
                SegmentIntersection::Segment(p1, p2) => {
                    stripped_intersections.push(*p1);
                    stripped_intersections.push(*p2);
                }
                _ => panic!("Unexpected intersection type"),
            }
        }
        stripped.insert(*segment_id, stripped_intersections);
    }

    return stripped;
}

fn merge_nearby_points(
    proximity_merger: &mut ProximityMerger,
    intersections: &HashMap<SegmentId, Vec<Vector2f>>,
) -> HashMap<SegmentId, Vec<Vector2f>> {
    let mut merged = intersections.clone();

    for (_, points) in merged.iter_mut() {
        for point in points.iter_mut() {
            *point = proximity_merger.map(*point);
        }
    }

    return merged;
}

pub fn subdivide_contours_at_intersections(polygon: &Polygon, epsilon: f32) -> Vec<Contour> {
    let intersections = find_intersections(&polygon, epsilon);
    let stripped_intersections = strip_intersection_type(&intersections);

    let mut proximity_merger = ProximityMerger::new(epsilon);
    let mut merged_intersections =
        merge_nearby_points(&mut proximity_merger, &stripped_intersections);

    let mut subdivided_contours = Vec::with_capacity(polygon.contours().len());

    for (contour_idx, contour) in polygon.contours().iter().enumerate() {
        let mut subdivided_contour = Contour::new();

        for (segment_idx, start_point) in contour.iter().enumerate() {
            let end_point = contour[(segment_idx + 1) % contour.len()];

            let seg_id = SegmentId {
                contour: contour_idx,
                segment: segment_idx,
            };
            // calling unwrap because every line must have an intersection at the start and endpoint
            let intersections_opt = merged_intersections.get_mut(&seg_id);
            if intersections_opt.is_none() {
                panic!(
                    "Segment with ID ({}, {}) has no intersections",
                    seg_id.contour, seg_id.segment
                );
            }
            let intersections = intersections_opt.unwrap();

            // sort intersection points along the segment going from start_point to end_point
            let segment_vec = end_point - *start_point;
            intersections.sort_by(|p1, p2| {
                let v1 = segment_vec.dot(*p1 - *start_point);
                let v2 = segment_vec.dot(*p2 - *start_point);
                return f32::total_cmp(&v1, &v2);
            });

            // insert all points, but skip duplicates
            for p in intersections {
                if subdivided_contour.last() != Some(p) {
                    subdivided_contour.push(*p);
                }
            }
        }

        // ensure start and end are not duplicates either
        while subdivided_contour[0] == *subdivided_contour.last().unwrap() {
            subdivided_contour.pop();
        }

        subdivided_contours.push(subdivided_contour);
    }

    return subdivided_contours;
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPSILON: f32 = 0.0001;

    #[test]
    fn test_intersect_segments_general_case_no_intersection() {
        let p0 = Vector2f::new(0.0, 0.0);
        let p1 = Vector2f::new(2.0, 2.0);
        let p2 = Vector2f::new(3.0, 3.0);
        let p3 = Vector2f::new(4.0, 5.0);
        assert_eq!(
            intersect_segments(p0, p1, p2, p3, EPSILON),
            SegmentIntersection::None
        );
    }

    #[test]
    fn test_intersect_segments_general_case_with_intersection() {
        let p0 = Vector2f::new(0.0, 0.0);
        let p1 = Vector2f::new(4.0, 4.0);
        let p2 = Vector2f::new(2.0, 0.0);
        let p3 = Vector2f::new(0.0, 2.0);
        assert_eq!(
            intersect_segments(p0, p1, p2, p3, EPSILON),
            SegmentIntersection::Point(Vector2f::new(1.0, 1.0)),
        );
    }

    #[test]
    fn test_intersect_segments_parallel_no_overlap() {
        let p0 = Vector2f::new(0.0, 0.0);
        let p1 = Vector2f::new(2.0, 2.0);
        let p2 = Vector2f::new(0.0, 1.0);
        let p3 = Vector2f::new(2.0, 3.0);
        assert_eq!(
            intersect_segments(p0, p1, p2, p3, EPSILON),
            SegmentIntersection::None
        );
    }

    #[test]
    fn test_intersect_segments_parallel_with_overlap() {
        let p0 = Vector2f::new(0.0, 0.0);
        let p1 = Vector2f::new(2.0, 2.0);
        let p2 = Vector2f::new(1.0, 1.0);
        let p3 = Vector2f::new(3.0, 3.0);
        assert_eq!(
            intersect_segments(p0, p1, p2, p3, EPSILON),
            SegmentIntersection::Segment(p2, p1),
        );
    }

    #[test]
    fn test_intersect_segments_parallel_with_overlap_point() {
        let p0 = Vector2f::new(0.0, 0.0);
        let p1 = Vector2f::new(2.0, 2.0);
        let p2 = Vector2f::new(2.0, 2.0);
        let p3 = Vector2f::new(3.0, 3.0);
        assert_eq!(
            intersect_segments(p0, p1, p2, p3, EPSILON),
            SegmentIntersection::Point(Vector2f::new(2.0, 2.0))
        );
    }

    #[test]
    fn test_intersect_segments_degenerate_case1_overlap() {
        let p0 = Vector2f::new(0.0, 0.0);
        let p1 = Vector2f::new(0.0, 0.0);
        let p2 = Vector2f::new(0.0, 0.0);
        let p3 = Vector2f::new(0.0, 0.0);
        assert_eq!(
            intersect_segments(p0, p1, p2, p3, EPSILON),
            SegmentIntersection::Point(Vector2f::new(0.0, 0.0))
        );
    }

    #[test]
    fn test_intersect_segments_degenerate_case1_no_overlap() {
        let p0 = Vector2f::new(0.0, 0.0);
        let p1 = Vector2f::new(0.0, 0.0);
        let p2 = Vector2f::new(2.0, 0.0);
        let p3 = Vector2f::new(2.0, 0.0);
        assert_eq!(
            intersect_segments(p0, p1, p2, p3, EPSILON),
            SegmentIntersection::None
        );
    }

    #[test]
    fn test_intersect_segments_degenerate_case2() {
        let p0 = Vector2f::new(0.0, 0.0);
        let p1 = Vector2f::new(1.0, 1.0);
        let p2 = Vector2f::new(0.0, 0.0);
        let p3 = Vector2f::new(0.0, 0.0);
        assert_eq!(
            intersect_segments(p0, p1, p2, p3, EPSILON),
            SegmentIntersection::Point(Vector2f::new(0.0, 0.0))
        );
    }

    #[test]
    fn test_intersect_segments_degenerate_case3() {
        let p0 = Vector2f::new(1.0, 1.0);
        let p1 = Vector2f::new(1.0, -10.0);
        let p2 = Vector2f::new(1.0, 1.0);
        let p3 = Vector2f::new(1.0, -10.0);
        assert_eq!(
            intersect_segments(p0, p1, p2, p3, EPSILON),
            SegmentIntersection::Segment(Vector2f::new(1.0, 1.0), Vector2f::new(1.0, -10.0))
        );
    }
}
