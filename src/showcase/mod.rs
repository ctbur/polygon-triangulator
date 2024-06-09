use raylib::drawing::RaylibDrawHandle;

mod intersections;

pub use intersections::*;

pub trait Showcase {
    fn render(&mut self, r: &mut RaylibDrawHandle);
}
