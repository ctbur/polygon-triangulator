use ::core::f32;
use std::collections::{HashMap, HashSet};
use std::time::{SystemTime, UNIX_EPOCH};

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use raylib::prelude::*;

use crate::decomposition::WindingRule;
use crate::graph::{self, Graph};
use crate::intersections::{self, SegmentId, SegmentIntersection};
use crate::polygon::{Contour, Polygon};
use crate::triangulation::{self, Triangle};
use crate::vector2::Vector2f;
use crate::{decomposition, partition};

#[derive(PartialEq, Eq)]
enum DrawingMode {
    Contours,
    Graph,
    Regions,
    MonotoneRegions,
    Triangulation,
}

const COLOR_SEQUENCE: [Color; 6] = [
    Color::RED,
    Color::BLUE,
    Color::GREEN,
    Color::ORANGE,
    Color::PURPLE,
    Color::PINK,
];

struct TriangulationStages {
    polygon: Polygon,
    intersections: HashMap<SegmentId, Vec<SegmentIntersection>>,
    graph: Graph,
    islands: Vec<(Contour, Vec<Contour>)>,
    monotone_regions: Vec<(Contour, Vec<Contour>)>,
    triangulations: Vec<(Contour, Vec<Triangle>)>,
}

pub struct Showcase {
    raylib_handle: RaylibHandle,
    raylib_thread: RaylibThread,
    stages: Vec<TriangulationStages>,
    selected_polygon: usize,
    drawing_mode: DrawingMode,
    camera: Camera2D,
}

impl Showcase {
    pub fn build(
        raylib_handle: RaylibHandle,
        raylib_thread: RaylibThread,
        polygons: Vec<Polygon>,
    ) -> Showcase {
        let epsilon = 0.01;
        let mut stages = Vec::new();
        for polygon in polygons {
            let intersections = intersections::find_intersections(&polygon, epsilon);
            let subdivided_contours =
                intersections::subdivide_contours_at_intersections(&polygon, epsilon);
            let graph = graph::build_graph(&subdivided_contours);

            let mut island_graphs = graph.clone().split_by_islands();
            let islands: Vec<_> = island_graphs
                .iter_mut()
                .map(|g| g.trace_regions())
                .collect();
            decomposition::decompose(graph.clone(), &subdivided_contours, WindingRule::Odd);

            // islands where all regions are monotone
            let mut monotone_regions = Vec::new();
            for (i, (_, ref interior)) in islands.iter().enumerate() {
                for region in interior {
                    partition::partition_region(&mut island_graphs[i], region);
                }
                monotone_regions.push(island_graphs[i].trace_regions());
            }

            let mut triangulations = Vec::new();
            for (_, ref interior) in &islands {
                for region in interior {
                    let mut triangles = Vec::new();
                    triangulation::triangulate_monotone_region(&mut triangles, region);
                    triangulations.push((region.clone(), triangles));
                }
            }

            let stage = TriangulationStages {
                polygon,
                intersections,
                graph,
                islands,
                monotone_regions,
                triangulations,
            };
            stages.push(stage);
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
            stages,
            selected_polygon: 0,
            drawing_mode: DrawingMode::Contours,
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
            self.selected_polygon = (self.selected_polygon + 1) % self.stages.len();
        }

        if self
            .raylib_handle
            .is_key_pressed(KeyboardKey::KEY_PAGE_DOWN)
        {
            self.selected_polygon =
                (self.selected_polygon + self.stages.len() - 1) % self.stages.len();
        }

        // select polygons with number keys
        for i in 0..self.stages.len() {
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

        // select drawing mode
        if self.raylib_handle.is_key_pressed(KeyboardKey::KEY_C) {
            self.drawing_mode = DrawingMode::Contours;
        } else if self.raylib_handle.is_key_pressed(KeyboardKey::KEY_G) {
            self.drawing_mode = DrawingMode::Graph;
        } else if self.raylib_handle.is_key_pressed(KeyboardKey::KEY_R) {
            self.drawing_mode = DrawingMode::Regions;
        } else if self.raylib_handle.is_key_pressed(KeyboardKey::KEY_M) {
            self.drawing_mode = DrawingMode::MonotoneRegions;
        } else if self.raylib_handle.is_key_pressed(KeyboardKey::KEY_T) {
            self.drawing_mode = DrawingMode::Triangulation;
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

        match self.drawing_mode {
            DrawingMode::Contours => {
                draw_contours(&mut c, &self.stages[self.selected_polygon].polygon)
            }
            DrawingMode::Graph => draw_graph(&mut c, &self.stages[self.selected_polygon].graph),
            DrawingMode::Regions => {
                draw_regions(&mut c, &self.stages[self.selected_polygon].islands)
            }
            DrawingMode::MonotoneRegions => {
                draw_regions(&mut c, &self.stages[self.selected_polygon].monotone_regions)
            }
            DrawingMode::Triangulation => {
                draw_triangulations(&mut c, &self.stages[self.selected_polygon].triangulations);
            }
        }

        if self.drawing_mode != DrawingMode::Regions {
            draw_intersections(
                &mut c,
                &self.stages[self.selected_polygon].intersections,
                Vector2f::new(0.0, 0.0),
            );
        }
    }
}

fn draw_contours<'a, T>(d: &mut RaylibMode2D<'a, T>, polygon: &Polygon) {
    for (i, contour) in polygon.contours().iter().enumerate() {
        draw_contour(
            d,
            contour,
            Vector2f::ZERO,
            4.0,
            COLOR_SEQUENCE[i % COLOR_SEQUENCE.len()],
        );
    }
}

fn draw_contour<'a, T>(
    d: &mut RaylibMode2D<'a, T>,
    contour: &Contour,
    offset: Vector2f,
    thickness: f32,
    color: Color,
) {
    for j in 0..contour.len() {
        let start = contour[j] + offset;
        let end = contour[(j + 1) % contour.len()] + offset;
        d.draw_line_ex(start, end, thickness, color);
    }
}

fn draw_graph<'a, T>(d: &mut RaylibMode2D<'a, T>, graph: &Graph) {
    let mut rng = StdRng::seed_from_u64(1234567823);
    let mut visited_nodes = HashSet::new();
    let mut node_stack = Vec::new();

    let nodes = graph.nodes();
    let mut color_counter = 0;
    for i in 0..nodes.len() {
        if visited_nodes.contains(&i) {
            continue;
        }

        node_stack.push(i);
        while let Some(node_idx) = node_stack.pop() {
            let node = &nodes[node_idx];

            for other_node_idx in node.edges() {
                if visited_nodes.contains(other_node_idx) {
                    continue;
                }

                let other_node = &nodes[*other_node_idx];
                let magnitude = rng.gen_range(0.0..10.0);
                let angle = rng.gen_range(0.0..2.0 * f32::consts::PI)
                    + ((SystemTime::now()
                        .duration_since(UNIX_EPOCH)
                        .unwrap()
                        .as_millis()
                        % 2000) as f32)
                        / 2000.0
                        * 2.0
                        * f32::consts::PI;
                let offset = Vector2f::new(angle.sin(), angle.cos()) * magnitude;
                d.draw_line_ex(
                    node.position() + offset,
                    other_node.position() + offset,
                    4.0,
                    COLOR_SEQUENCE[color_counter % COLOR_SEQUENCE.len()],
                );
                color_counter += 1;

                node_stack.push(*other_node_idx);
            }

            visited_nodes.insert(node_idx);
        }
    }
}

fn draw_regions<'a, T>(d: &mut RaylibMode2D<'a, T>, regions: &Vec<(Contour, Vec<Contour>)>) {
    let mut rng = StdRng::seed_from_u64(1234567890);
    for (ref outline, ref interior) in regions {
        draw_contour(d, &outline, Vector2f::ZERO, 8.0, Color::BLACK.alpha(0.3));
        for (i, region) in interior.iter().enumerate() {
            let color = COLOR_SEQUENCE[i % COLOR_SEQUENCE.len()].alpha(0.5);

            let magnitude = rng.gen_range(0.0..10.0);
            let angle = rng.gen_range(0.0..2.0 * f32::consts::PI)
                + ((SystemTime::now()
                    .duration_since(UNIX_EPOCH)
                    .unwrap()
                    .as_millis()
                    % 2000) as f32)
                    / 2000.0
                    * 2.0
                    * f32::consts::PI;
            let offset = Vector2f::new(angle.sin(), angle.cos()) * magnitude;
            draw_contour(d, region, offset, 4.0, color);
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

fn draw_triangulations<'a, T>(
    d: &mut RaylibMode2D<'a, T>,
    triangluations: &Vec<(Contour, Vec<Triangle>)>,
) {
    println!("draw");
    let mut color_idx = 0;
    for (c, triangles) in triangluations {
        for t in triangles {
            d.draw_triangle(
                c[t[0]],
                c[t[1]],
                c[t[2]],
                COLOR_SEQUENCE[color_idx % COLOR_SEQUENCE.len()],
            );
        }
        color_idx += 1;
    }

    for (c, triangles) in triangluations {
        for t in triangles {
            d.draw_line_ex(c[t[0]], c[t[1]], 3.0, Color::BLACK);
            d.draw_line_ex(c[t[1]], c[t[2]], 3.0, Color::BLACK);
            d.draw_line_ex(c[t[2]], c[t[0]], 3.0, Color::BLACK);
        }
        color_idx += 1;
    }
}

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

pub fn get_overlapping_rects() -> Polygon {
    let mut polygon = Polygon::new();

    polygon.move_to(50.0, 50.0);
    polygon.line_to(800.0, 50.0);
    polygon.line_to(800.0, 800.0);
    polygon.line_to(50.0, 800.0);

    polygon.move_to(150.0, 150.0);
    polygon.line_to(150.0, 450.0);
    polygon.line_to(450.0, 450.0);
    polygon.line_to(450.0, 150.0);

    polygon.move_to(400.0, 400.0);
    polygon.line_to(400.0, 700.0);
    polygon.line_to(700.0, 700.0);
    polygon.line_to(700.0, 400.0);

    return polygon;
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

fn add_box(polygon: &mut Polygon, v1: Vector2f, v2: Vector2f) {
    polygon.move_to(v1.x, v1.y);
    polygon.line_to(v2.x, v1.y);
    polygon.line_to(v2.x, v2.y);
    polygon.line_to(v1.x, v2.y);
}

pub fn get_nested_boxes() -> Polygon {
    let mut polygon = Polygon::new();

    add_box(
        &mut polygon,
        Vector2f::new(50.0, 50.0),
        Vector2f::new(850.0, 850.0),
    );
    add_box(
        &mut polygon,
        Vector2f::new(100.0, 100.0),
        Vector2f::new(500.0, 500.0),
    );
    add_box(
        &mut polygon,
        Vector2f::new(600.0, 150.0),
        Vector2f::new(700.0, 300.0),
    );
    add_box(
        &mut polygon,
        Vector2f::new(200.0, 200.0),
        Vector2f::new(300.0, 300.0),
    );
    add_box(
        &mut polygon,
        Vector2f::new(150.0, 600.0),
        Vector2f::new(700.0, 750.0),
    );
    add_box(
        &mut polygon,
        Vector2f::new(200.0, 550.0),
        Vector2f::new(650.0, 700.0),
    );

    return polygon;
}
