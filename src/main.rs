mod graph;
mod intersections;
mod polygon;
mod regions;
mod showcase;
mod vector2;

fn main() {
    let (rl, thread) = raylib::init()
        .size(2500, 1500)
        .title("Hello, World")
        .build();

    let polygons = vec![
        showcase::get_self_overlapping_contour(),
        showcase::get_arbitrary_polygon(),
        showcase::get_star(),
        showcase::get_rounded_rect(),
        showcase::get_overlapping_rects(),
        showcase::get_grid(),
    ];
    let mut showcase = showcase::Showcase::build(rl, thread, polygons);

    while showcase.process() {}
}
