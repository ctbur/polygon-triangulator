mod intersections;
mod polygon;
mod showcase;
mod vector2;

use polygon::Polygon;

fn main() {
    let (rl, thread) = raylib::init()
        .size(2500, 1500)
        .title("Hello, World")
        .build();

    let polygons = vec![
        showcase::get_arbitrary_polygon(),
        showcase::get_star(),
        showcase::get_rounded_rect(),
        showcase::get_grid(),
    ];
    let mut showcase = showcase::Showcase::build(rl, thread, polygons);

    while showcase.process() {}
}
