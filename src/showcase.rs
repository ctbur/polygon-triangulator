use raylib::prelude::*;

use crate::intersections::{self, SegmentIntersection};
use crate::{vector2::Vector2f, Polygon};

pub struct Showcase {
    raylib_handle: RaylibHandle,
    raylib_thread: RaylibThread,
    polygons: Vec<Polygon>,
    intersections: Vec<Vec<SegmentIntersection>>,
    selected_polygon: usize,
    camera: Camera2D,
}

impl Showcase {
    pub fn build(
        raylib_handle: RaylibHandle,
        raylib_thread: RaylibThread,
        polygons: Vec<Polygon>,
    ) -> Showcase {
        let mut intersections = Vec::new();
        for polygon in &polygons {
            intersections.push(intersections::find_intersections(polygon, 0.01));
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

        self.camera.offset += offset;
    }

    fn render(&mut self) {
        let mut d = self.raylib_handle.begin_drawing(&self.raylib_thread);
        d.clear_background(Color::WHITE);

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
    for contour in polygon.contours() {
        for i in 0..contour.len() {
            let start = contour[i];
            let end = contour[(i + 1) % contour.len()];
            d.draw_line_ex(start + offset, end + offset, 4.0, Color::BLACK);
        }
    }
}

fn draw_intersections<'a, T>(
    d: &mut RaylibMode2D<'a, T>,
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
    polygon.quad_bezier_curve_to(radius, 0.0, 0.0, radius, 10);
    polygon.line_to(0.0, height - radius);
    polygon.quad_bezier_curve_to(0.0, height - radius, radius, height, 10);
    polygon.line_to(width - radius, height);
    polygon.quad_bezier_curve_to(width - radius, height, width, height - radius, 10);
    polygon.line_to(width, radius);
    polygon.quad_bezier_curve_to(width, radius, width - radius, 0.0, 10);
    polygon.line_to(radius, 0.0);

    polygon
}

pub fn get_star() -> Polygon {
    let mut polygon = Polygon::new();

    let points = 5;
    let radius = 100.0;
    let center = Vector2f::new(400.0, 400.0);

    polygon.move_to(center.x, center.y - radius);
    for i in 0..points * 2 {
        let angle = 2.0 * std::f32::consts::PI / (points * 2) as f32 * i as f32;
        let length = if i % 2 == 0 { radius } else { radius / 2.0 };
        let x = center.x + angle.cos() * length;
        let y = center.y + angle.sin() * length;
        polygon.line_to(x, y);
    }

    polygon
}

pub fn get_grid() -> Polygon {
    let mut polygon = Polygon::new();

    let width = 800.0;
    let height = 800.0;
    let cell_size = 50.0;

    for x in (0..(width as i32)).step_by(cell_size as usize) {
        polygon.move_to(x as f32, 0.0);
        polygon.line_to(x as f32, height);
    }

    for y in (0..(height as i32)).step_by(cell_size as usize) {
        polygon.move_to(0.0, y as f32);
        polygon.line_to(width, y as f32);
    }

    polygon
}
