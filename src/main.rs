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

        draw_polygon(&mut d, &polygon, Vector2::new(0.0, 0.0));
        draw_polygon(&mut d, &rect, Vector2::new(700.0, 0.0));

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

fn draw_polygon(d: &mut RaylibDrawHandle, polygon: &Polygon, offset: Vector2) {
    for contour in &polygon.contours {
        for i in 0..contour.len() {
            let (x0, y0) = contour[i];
            let (x1, y1) = contour[(i + 1) % contour.len()];
            d.draw_line_ex(
                Vector2::new(x0, y0) + offset,
                Vector2::new(x1, y1) + offset,
                1.0,
                Color::BLACK,
            );
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
    contours: Vec<Vec<(f32, f32)>>,
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

        self.contours.push(vec![(x, y)]);
    }

    fn line_to(&mut self, x: f32, y: f32) {
        if let Some(contour) = self.contours.last_mut() {
            assert!(contour.len() >= 1);

            contour.push((x, y));
        } else {
            panic!("line_to called without a previous move_to");
        }
    }

    fn quad_bezier_curve_to(&mut self, x1: f32, y1: f32, x2: f32, y2: f32, steps: u32) {
        if let Some(contour) = self.contours.last_mut() {
            assert!(contour.len() >= 1);

            let (x0, y0) = *contour.last().unwrap();

            for i in 1..=steps {
                let t = i as f32 / steps as f32;
                let x = (1.0 - t).powi(2) * x0 + 2.0 * (1.0 - t) * t * x1 + t.powi(2) * x2;
                let y = (1.0 - t).powi(2) * y0 + 2.0 * (1.0 - t) * t * y1 + t.powi(2) * y2;
                contour.push((x, y));
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

        for (x, y) in contour {
            let xd = (x * scale_factor) as i64;
            let yd = (y * scale_factor) as i64;
            discrete_contour.push((xd, yd));
        }

        contours.push(discrete_contour);
    }

    DiscretePolygon { contours }
}
