use raylib::prelude::*;

use crate::{find_intersections, vector2::Vector2f, Polygon, SegmentIntersection};

use super::Showcase;

pub struct IntersectionsShowcase {
    polygon: Polygon,
    intersections: Vec<SegmentIntersection>,
}

impl IntersectionsShowcase {
    pub fn new(polygon: Polygon) -> IntersectionsShowcase {
        let intersections = find_intersections(&polygon, 0.001);
        IntersectionsShowcase {
            polygon,
            intersections,
        }
    }
}

impl Showcase for IntersectionsShowcase {
    fn render(&mut self, d: &mut RaylibDrawHandle) {
        d.clear_background(Color::WHITE);

        draw_polygon(d, &self.polygon, Vector2f::new(0.0, 0.0));
        draw_intersections(d, &self.intersections, Vector2f::new(0.0, 0.0));
    }
}

fn draw_polygon(d: &mut RaylibDrawHandle, polygon: &Polygon, offset: Vector2f) {
    for contour in &polygon.contours {
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
                panic!("reeear");
            }
        }
    }
}
