mod graph;
mod hierarchy;
mod intersections;
mod partition;
mod polygon;
mod regions;
mod showcase;
mod triangulation;
mod vector2;

fn main() {
    let (mut rl, thread) = raylib::init()
        .size(2500, 1500)
        .title("Hello, World")
        .build();
    rl.set_target_fps(60);

    let polygons = vec![
        showcase::get_self_overlapping_contour(),
        showcase::get_arbitrary_polygon(),
        showcase::get_star(),
        showcase::get_rounded_rect(),
        showcase::get_overlapping_rects(),
        showcase::get_grid(),
        showcase::get_nested_boxes(),
    ];
    let mut showcase = showcase::Showcase::build(rl, thread, polygons);

    while showcase.process() {}
}
