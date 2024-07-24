use std::collections::HashMap;

use raylib::prelude::*;

use crate::graph::{self, Graph};
use crate::intersections::{self, SegmentId, SegmentIntersection};
use crate::{vector2::Vector2f, Polygon};

pub struct Showcase {
    raylib_handle: RaylibHandle,
    raylib_thread: RaylibThread,
    polygons: Vec<Polygon>,
    intersections: Vec<HashMap<SegmentId, Vec<SegmentIntersection>>>,
    graphs: Vec<Graph>,
    selected_polygon: usize,
    camera: Camera2D,
}

impl Showcase {
    pub fn build(
        raylib_handle: RaylibHandle,
        raylib_thread: RaylibThread,
        polygons: Vec<Polygon>,
    ) -> Showcase {
        // find intersections for all polygons
        let mut intersections = Vec::new();
        for polygon in &polygons {
            intersections.push(intersections::find_intersections(polygon, 0.01));
        }

        // build graphs for all polygons
        let mut graphs = Vec::new();
        for (polygon, intersections) in polygons.iter().zip(intersections.iter()) {
            graphs.push(graph::build_graph(&polygon, &intersections, 0.01));
        }

        let camera = Camera2D {
            offset: Vector2::new(0.0, 0.0),
            target: Vector2::new(0.0, 0.0),
            rotation: 0.0,
            zoom: 1.0,
        };
        Showcase {
            raylib_handle,
            raylib_thread,
            polygons,
            intersections,
            graphs,
            selected_polygon: 0,
            camera,
        }
    }

    pub fn process(&mut self) -> bool {
        if self.raylib_handle.window_should_close() {
            return false;
        }

        self.handle_input();
        self.render();
        return true;
    }

    fn handle_input(&mut self) {
        // select next/previous polygon
        if self.raylib_handle.is_key_pressed(KeyboardKey::KEY_PAGE_UP) {
            self.selected_polygon = (self.selected_polygon + 1) % self.polygons.len();
        }

        if self
            .raylib_handle
            .is_key_pressed(KeyboardKey::KEY_PAGE_DOWN)
        {
            self.selected_polygon =
                (self.selected_polygon + self.polygons.len() - 1) % self.polygons.len();
        }

        // select polygons with number keys
        for i in 0..self.polygons.len() {
            let key_by_number = match i {
                0 => KeyboardKey::KEY_ONE,
                1 => KeyboardKey::KEY_TWO,
                2 => KeyboardKey::KEY_THREE,
                3 => KeyboardKey::KEY_FOUR,
                4 => KeyboardKey::KEY_FIVE,
                5 => KeyboardKey::KEY_SIX,
                6 => KeyboardKey::KEY_SEVEN,
                7 => KeyboardKey::KEY_EIGHT,
                8 => KeyboardKey::KEY_NINE,
                9 => KeyboardKey::KEY_ZERO,
                _ => break,
            };

            if self.raylib_handle.is_key_pressed(key_by_number) {
                self.selected_polygon = i;
            }
        }

        // move camera with arrow keys
        let mut offset = Vector2::new(0.0, 0.0);
        if self.raylib_handle.is_key_down(KeyboardKey::KEY_LEFT) {
            offset.x += 10.0;
        }
        if self.raylib_handle.is_key_down(KeyboardKey::KEY_RIGHT) {
            offset.x -= 10.0;
        }
        if self.raylib_handle.is_key_down(KeyboardKey::KEY_UP) {
            offset.y += 10.0;
        }
        if self.raylib_handle.is_key_down(KeyboardKey::KEY_DOWN) {
            offset.y -= 10.0;
        }

        // zoom in/out with mouse wheel
        let zoom = self.raylib_handle.get_mouse_wheel_move();
        self.camera.zoom += zoom * 0.1;

        // reset camera with space
        if self.raylib_handle.is_key_pressed(KeyboardKey::KEY_SPACE) {
            self.camera.offset = Vector2::new(0.0, 0.0);
            self.camera.zoom = 1.0;
        }

        self.camera.offset += offset * self.raylib_handle.get_frame_time() * 60.0;
    }

    fn render(&mut self) {
        let mut d = self.raylib_handle.begin_drawing(&self.raylib_thread);
        d.clear_background(Color::WHITE);
        d.draw_fps(10, 10);

        let mut c = d.begin_mode2D(self.camera);
        draw_polygon(
            &mut c,
            &self.polygons[self.selected_polygon],
            Vector2f::new(0.0, 0.0),
        );
        draw_intersections(
            &mut c,
            &self.intersections[self.selected_polygon],
            Vector2f::new(0.0, 0.0),
        );
    }
}

fn draw_polygon<'a, T>(d: &mut RaylibMode2D<'a, T>, polygon: &Polygon, offset: Vector2f) {
    let colors = [
        Color::RED,
        Color::BLUE,
        Color::GREEN,
        Color::ORANGE,
        Color::PURPLE,
        Color::PINK,
    ];

    for (i, contour) in polygon.contours().iter().enumerate() {
        for j in 0..contour.len() {
            let start = contour[j];
            let end = contour[(j + 1) % contour.len()];
            d.draw_line_ex(start + offset, end + offset, 4.0, colors[i % colors.len()]);
        }
    }
}

fn draw_intersections<'a, T>(
    d: &mut RaylibMode2D<'a, T>,
    intersections: &HashMap<SegmentId, Vec<SegmentIntersection>>,
    offset: Vector2f,
) {
    let scale = 4.0;
    for (_, intersections) in intersections {
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
                        Color::BLACK.alpha(0.3),
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
                        Color::BLACK.alpha(0.3),
                    );
                    d.draw_ring(
                        p1 + offset,
                        4.0 * scale,
                        5.0 * scale,
                        0.0,
                        360.0,
                        20,
                        Color::BLACK.alpha(0.3),
                    );
                    d.draw_line_ex(
                        p0 + offset,
                        p1 + offset,
                        1.0 * scale,
                        Color::BLACK.alpha(0.6),
                    );
                }
                _ => {
                    panic!();
                }
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

pub fn get_self_overlapping_contour() -> Polygon {
    let mut polygon = Polygon::new();

    polygon.move_to(100.0, 100.0);
    polygon.line_to(500.0, 100.0);
    polygon.line_to(100.0, 500.0);
    polygon.line_to(500.0, 500.0);

    polygon
}

pub fn get_arbitrary_polygon() -> Polygon {
    let mut polygon = Polygon::new();

    polygon.move_to(100.0, 100.0);
    polygon.line_to(500.0, 100.0);
    polygon.line_to(500.0, 500.0);
    polygon.line_to(100.0, 800.0);

    // draw a triangle inside the previous polygon
    polygon.move_to(200.0, 200.0);
    polygon.line_to(400.0, 200.0);
    polygon.line_to(300.0, 1000.0);

    // draw another triangle
    polygon.move_to(500.0, 200.0);
    polygon.line_to(500.0, 300.0);
    polygon.line_to(700.0, 250.0);

    // draw another triangle
    polygon.move_to(100.0, 600.0);
    polygon.line_to(100.0, 1000.0);
    polygon.line_to(700.0, 800.0);

    polygon
}

pub fn get_rounded_rect() -> Polygon {
    let mut polygon = Polygon::new();

    let width = 800.0;
    let height = 800.0;
    let radius = 100.0;

    polygon.move_to(radius, 0.0);
    polygon.quad_bezier_curve_to(0.0, 0.0, 0.0, radius, 10);
    polygon.line_to(0.0, height - radius);
    polygon.quad_bezier_curve_to(0.0, height, radius, height, 10);
    polygon.line_to(width - radius, height);
    polygon.quad_bezier_curve_to(width, height, width, height - radius, 10);
    polygon.line_to(width, radius);
    polygon.quad_bezier_curve_to(width, 0.0, width - radius, 0.0, 10);
    polygon.line_to(radius, 0.0);

    polygon
}

pub fn get_star() -> Polygon {
    let mut polygon = Polygon::new();

    let points = 5;
    let radius = 100.0;
    let center = Vector2f::new(400.0, 400.0);

    for i in 0..points * 2 {
        let angle = 2.0 * std::f32::consts::PI / (points * 2) as f32 * i as f32;
        let length = if i % 2 == 0 { radius } else { radius / 2.0 };
        let x = center.x + angle.cos() * length;
        let y = center.y + angle.sin() * length;

        if i == 0 {
            polygon.move_to(x, y);
        } else {
            polygon.line_to(x, y);
        }
    }

    polygon
}

pub fn get_grid() -> Polygon {
    let mut polygon = Polygon::new();

    let cell_size = 50.0;
    let width = 800.0;
    let height = 800.0;

    let num_rows = (height / cell_size) as usize;
    let num_cols = (width / cell_size) as usize;
    polygon.move_to(0.0, 0.0);
    for row in 0..num_rows {
        for col in 0..num_cols {
            let x = col as f32 * cell_size;
            let y = row as f32 * cell_size;
            polygon.line_to(x, y);
            polygon.line_to(x + cell_size, y);
            polygon.line_to(x + cell_size, y + cell_size);
            polygon.line_to(x, y + cell_size);
        }
    }

    polygon
}
