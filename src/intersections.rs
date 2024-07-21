use std::mem;

use crate::polygon::Polygon;
use crate::vector2::Vector2f;

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum SegmentIntersection {
    None,
    Point(Vector2f),
    Segment(Vector2f, Vector2f),
}

pub fn find_intersections(polygon: &Polygon, epsilon: f32) -> Vec<SegmentIntersection> {
    let mut intersections = Vec::new();

    let contours = polygon.contours();
    for i in 0..contours.len() {
        for j in i + 1..contours.len() {
            let contour = &contours[i];
            let other_contour = &contours[j];

            for k in 0..contour.len() {
                let p0 = contour[k];
                let p1 = contour[(k + 1) % contour.len()];

                for l in 0..other_contour.len() {
                    let p2 = other_contour[l];
                    let p3 = other_contour[(l + 1) % other_contour.len()];

                    match intersect_segments(p0, p1, p2, p3, epsilon) {
                        SegmentIntersection::None => {}
                        intersection => intersections.push(intersection),
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
