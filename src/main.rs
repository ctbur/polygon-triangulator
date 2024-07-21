mod intersections;
mod polygon;
mod showcase;
mod vector2;

use polygon::Polygon;

fn main() {
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

    let (mut rl, thread) = raylib::init()
        .size(2500, 1500)
        .title("Hello, World")
        .build();

    let mut showcase = showcase::Showcase::new(rl, thread, polygon);

    while showcase.process() {}
}
