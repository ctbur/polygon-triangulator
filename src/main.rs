mod vector2;

use std::mem::{self, swap};

use self::vector2::Vector2f;
use raylib::prelude::*;

fn main() {
    let (mut rl, thread) = raylib::init()
        .size(2500, 1500)
        .title("Hello, World")
        .build();

    while !rl.window_should_close() {
        let mut d = rl.begin_drawing(&thread);

        d.clear_background(Color::WHITE);
        d.draw_text("Hello, world!", 12, 12, 100, Color::BLACK);

        let mut polygon = Polygon::new();

        polygon.move_to(100.0, 100.0);
        polygon.line_to(500.0, 100.0);
        polygon.line_to(500.0, 500.0);
        polygon.line_to(100.0, 800.0);

        // draw a triangle inside the previous polygon
        polygon.move_to(200.0, 200.0);
        polygon.line_to(400.0, 200.0);
        polygon.line_to(300.0, 400.0);

        // draw a rectangle with rounded corners
        let mut rect = Polygon::new();
        rect.move_to(100.0, 100.0);
        rect.line_to(500.0, 100.0);
        rect.quad_bezier_curve_to(600.0, 100.0, 600.0, 200.0, 50);
        rect.line_to(600.0, 500.0);
        rect.quad_bezier_curve_to(600.0, 600.0, 500.0, 600.0, 50);
        rect.line_to(100.0, 600.0);
        rect.quad_bezier_curve_to(0.0, 600.0, 0.0, 500.0, 50);
        rect.line_to(0.0, 200.0);
        rect.quad_bezier_curve_to(0.0, 100.0, 100.0, 100.0, 50);

        draw_polygon(&mut d, &polygon, Vector2f::new(0.0, 0.0));
        draw_polygon(&mut d, &rect, Vector2f::new(700.0, 0.0));

        // discretize polygons and draw them a bit further down
        let scale_factor = 10.0;
        let discrete_polygon = discretize(&polygon, scale_factor);
        let discrete_rect = discretize(&rect, scale_factor);

        draw_discrete_polygon(
            &mut d,
            &discrete_polygon,
            Vector2::new(0.0, 700.0),
            scale_factor,
        );
        draw_discrete_polygon(
            &mut d,
            &discrete_rect,
            Vector2::new(700.0, 700.0),
            scale_factor,
        );
    }
}

fn draw_polygon(d: &mut RaylibDrawHandle, polygon: &Polygon, offset: Vector2f) {
    for contour in &polygon.contours {
        for i in 0..contour.len() {
            let start = contour[i];
            let end = contour[(i + 1) % contour.len()];
            d.draw_line_ex(start + offset, end + offset, 1.0, Color::BLACK);
        }
    }
}

fn draw_discrete_polygon(
    d: &mut RaylibDrawHandle,
    polygon: &DiscretePolygon,
    offset: Vector2,
    scale_factor: f32,
) {
    for contour in &polygon.contours {
        for i in 0..contour.len() {
            let (x0, y0) = contour[i];
            let (x1, y1) = contour[(i + 1) % contour.len()];
            d.draw_line_ex(
                Vector2::new(x0 as f32, y0 as f32) / scale_factor + offset,
                Vector2::new(x1 as f32, y1 as f32) / scale_factor + offset,
                1.0,
                Color::BLACK,
            );
        }
    }
}

struct Polygon {
    contours: Vec<Vec<Vector2f>>,
}

impl Polygon {
    fn new() -> Self {
        Self {
            contours: Vec::new(),
        }
    }

    fn move_to(&mut self, x: f32, y: f32) {
        // if there is a previous contour that has less then 3 points, remove it
        if let Some(contour) = self.contours.last() {
            if contour.len() < 3 {
                self.contours.pop();
            }
        }

        self.contours.push(vec![Vector2f::new(x, y)]);
    }

    fn line_to(&mut self, x: f32, y: f32) {
        if let Some(contour) = self.contours.last_mut() {
            assert!(contour.len() >= 1);

            contour.push(Vector2f::new(x, y));
        } else {
            panic!("line_to called without a previous move_to");
        }
    }

    fn quad_bezier_curve_to(&mut self, x1: f32, y1: f32, x2: f32, y2: f32, steps: u32) {
        if let Some(contour) = self.contours.last_mut() {
            assert!(contour.len() >= 1);

            let p0 = *contour.last().unwrap();
            let (x0, y0) = (p0.x, p0.y);

            for i in 1..=steps {
                let t = i as f32 / steps as f32;
                let x = (1.0 - t).powi(2) * x0 + 2.0 * (1.0 - t) * t * x1 + t.powi(2) * x2;
                let y = (1.0 - t).powi(2) * y0 + 2.0 * (1.0 - t) * t * y1 + t.powi(2) * y2;
                contour.push(Vector2f::new(x, y));
            }
        } else {
            panic!("bezier_curve_to called without a previous move_to");
        }
    }
}

struct DiscretePolygon {
    contours: Vec<Vec<(i64, i64)>>,
}

fn discretize(polygon: &Polygon, scale_factor: f32) -> DiscretePolygon {
    let mut contours = Vec::new();

    for contour in &polygon.contours {
        let mut discrete_contour = Vec::new();

        for p in contour {
            let xd = (p.x * scale_factor) as i64;
            let yd = (p.y * scale_factor) as i64;
            discrete_contour.push((xd, yd));
        }

        contours.push(discrete_contour);
    }

    DiscretePolygon { contours }
}

#[derive(Debug, PartialEq, Clone, Copy)]
enum SegmentIntersection {
    None,
    Point(Vector2f),
    Segment(Vector2f, Vector2f),
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
