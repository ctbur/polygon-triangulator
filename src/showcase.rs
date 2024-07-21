use raylib::prelude::*;

use crate::intersections::{self, SegmentIntersection};
use crate::{vector2::Vector2f, Polygon};

pub struct Showcase {
    raylib_handle: RaylibHandle,
    raylib_thread: RaylibThread,
    polygon: Polygon,
    intersections: Vec<SegmentIntersection>,
}

impl Showcase {
    pub fn new(
        raylib_handle: RaylibHandle,
        raylib_thread: RaylibThread,
        polygon: Polygon,
    ) -> Showcase {
        let intersections = intersections::find_intersections(&polygon, 0.001);
        Showcase {
            raylib_handle,
            raylib_thread,
            polygon,
            intersections,
        }
    }

    pub fn process(&mut self) -> bool {
        if self.raylib_handle.window_should_close() {
            return false;
        }

        let mut d = self.raylib_handle.begin_drawing(&self.raylib_thread);

        d.clear_background(Color::WHITE);
        draw_polygon(&mut d, &self.polygon, Vector2f::new(0.0, 0.0));
        draw_intersections(&mut d, &self.intersections, Vector2f::new(0.0, 0.0));

        return true;
    }
}

fn draw_polygon(d: &mut RaylibDrawHandle, polygon: &Polygon, offset: Vector2f) {
    for contour in polygon.contours() {
        for i in 0..contour.len() {
            let start = contour[i];
            let end = contour[(i + 1) % contour.len()];
            d.draw_line_ex(start + offset, end + offset, 4.0, Color::BLACK);
        }
    }
}

fn draw_intersections(
    d: &mut RaylibDrawHandle,
    intersections: &[SegmentIntersection],
    offset: Vector2f,
) {
    let scale = 4.0;
    for intersection in intersections {
        match intersection {
            &SegmentIntersection::Point(p) => {
                d.draw_ring(
                    p + offset,
                    4.0 * scale,
                    5.0 * scale,
                    0.0,
                    360.0,
                    20,
                    Color::GREEN.alpha(0.5),
                );
            }
            &SegmentIntersection::Segment(p0, p1) => {
                d.draw_ring(
                    p0 + offset,
                    4.0 * scale,
                    5.0 * scale,
                    0.0,
                    360.0,
                    20,
                    Color::BLUE.alpha(0.5),
                );
                d.draw_ring(
                    p1 + offset,
                    4.0 * scale,
                    5.0 * scale,
                    0.0,
                    360.0,
                    20,
                    Color::BLUE.alpha(0.5),
                );
                d.draw_line_ex(p0 + offset, p1 + offset, 1.0 * scale, Color::BLUE);
            }
            _ => {
                panic!();
            }
        }
    }
}

/*fn draw_discrete_polygon(
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
}*/
