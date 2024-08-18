use core::f32;

use crate::polygon::{self, Contour};

pub type Triangle = [usize; 3];

struct MonotoneRegionTriangulator<'a> {
    triangles: &'a mut Vec<Triangle>,
    region: &'a Contour,
    // current reflex chain, first entry is start_idx of the currently untriangulated region
    reflex_chain: Vec<usize>,
    // if the reflex chain is on the upper or lower chain
    is_reflex_chain_upper: bool,
}

impl<'a> MonotoneRegionTriangulator<'a> {
    fn start(
        triangles: &'a mut Vec<Triangle>,
        region: &'a Contour,
        start_idx: usize,
    ) -> MonotoneRegionTriangulator<'a> {
        println!("Start: {}", start_idx);
        MonotoneRegionTriangulator {
            triangles,
            region,
            reflex_chain: vec![start_idx],
            is_reflex_chain_upper: false,
        }
    }

    fn process_point(&mut self, point_idx: usize, is_upper_point: bool) {
        if self.reflex_chain.len() >= 2 {
            if is_upper_point != self.is_reflex_chain_upper {
                let new_start_idx = *self.reflex_chain.last().unwrap();

                self.process_all(point_idx);

                self.reflex_chain.push(new_start_idx);
            } else {
                self.process_visible(point_idx);
            }
        }

        self.reflex_chain.push(point_idx);
        self.is_reflex_chain_upper = is_upper_point;
    }

    fn end(&mut self, point_idx: usize) {
        println!("End: {} - reflex chain: {:?}", point_idx, self.reflex_chain);
        self.process_all(point_idx);
    }

    fn process_all(&mut self, point_idx: usize) {
        println!("Case 1: drain");
        // Case 1: process the reflex chain fully
        let mut prev_ref_chain_idx = self.reflex_chain.pop().unwrap();
        while let Some(ref_chain_idx) = self.reflex_chain.pop() {
            self.push_triangle([point_idx, prev_ref_chain_idx, ref_chain_idx]);
            prev_ref_chain_idx = ref_chain_idx;
        }
    }

    fn process_visible(&mut self, point_idx: usize) {
        println!("Case 2: Reflex or not");
        // Case 2: process the reflex chain until reflex_chain[-2] is not visible from point
        // i.e., point -> reflex_chain[-1] -> reflex_chain[-2] is reflex
        while self.reflex_chain.len() >= 2 && !self.is_top_point_reflex(point_idx) {
            let triangle = [
                point_idx,
                self.reflex_chain.pop().unwrap(),
                *self.reflex_chain.last().unwrap(),
            ];
            self.push_triangle(triangle);
        }
    }

    fn is_top_point_reflex(&self, point_idx: usize) -> bool {
        let point = self.region[point_idx];
        let prev = self.region[self.reflex_chain[self.reflex_chain.len() - 1]];
        let prev_prev = self.region[self.reflex_chain[self.reflex_chain.len() - 2]];

        return (prev_prev - prev).cross(point - prev) < 0.0;
    }

    // Adds a triangle, in CCW order
    fn push_triangle(&mut self, mut t: Triangle) {
        let (a, b, c) = (self.region[t[0]], self.region[t[1]], self.region[t[2]]);
        if (b - a).cross(c - a) > 0.0 {
            t.reverse();
        }
        self.triangles.push(t);
    }
}

/// Triangulates an x-monotone polygon running in CW order with no duplicate points
/// Uses the plane-sweep algorithm described in https://www.cs.umd.edu/class/fall2021/cmsc754/Lects/lect05-triangulate.pdf
pub fn triangulate_monotone_region(triangles: &mut Vec<Triangle>, region: &Contour) {
    if region.len() < 3 {
        return;
    }

    if polygon::calculate_region_area(region) > 0.0 {
        //panic!("Region needs to run in CW order");
        println!("Skipping region not running in CW order");
        return;
    }

    // find start point, i.e., leftmost point
    let mut start_idx = 0;
    let mut min_x = f32::INFINITY;
    for (i, point) in region.iter().enumerate() {
        if min_x > point.x {
            min_x = point.x;
            start_idx = i;
        }
    }

    let mut triangulator = MonotoneRegionTriangulator::start(triangles, region, start_idx);

    let mut upper_idx = (start_idx + 1) % region.len();
    let mut lower_idx = (start_idx + region.len() - 1) % region.len();

    // sweep across points in x-direction
    while upper_idx != lower_idx {
        // if the point to be processed is on the upper or lower chain
        let is_upper_point = region[upper_idx].x < region[lower_idx].x;

        let point_idx;
        if is_upper_point {
            point_idx = upper_idx;
            upper_idx = (upper_idx + 1) % region.len();
        } else {
            point_idx = lower_idx;
            lower_idx = (lower_idx + region.len() - 1) % region.len();
        }
        println!(
            "Next event: {}, {} - reflex chain: {:?}, upper_idx: {}, lower_idx: {}",
            point_idx,
            if is_upper_point { "upper" } else { "lower" },
            triangulator.reflex_chain,
            lower_idx,
            upper_idx
        );

        triangulator.process_point(point_idx, is_upper_point);
    }

    triangulator.end(upper_idx);
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use crate::{polygon::Contour, vector2::Vector2f};

    use super::*;

    fn validate_triangulation(region: &Contour, triangles: &Vec<Triangle>) {
        let mut edge_count = HashMap::new();
        fn insert_edge(edge_count: &mut HashMap<[usize; 2], usize>, from: usize, to: usize) {
            let mut edge = [from, to];
            edge.sort();

            let entry = edge_count.entry(edge).or_default();
            *entry += 1;
        }

        for triangle in triangles {
            insert_edge(&mut edge_count, triangle[0], triangle[1]);
            insert_edge(&mut edge_count, triangle[1], triangle[2]);
            insert_edge(&mut edge_count, triangle[2], triangle[0]);
        }

        for i in 0..region.len() {
            let mut edge = [i, (i + 1) % region.len()];
            edge.sort();
            assert_eq!(
                edge_count[&edge], 1,
                "Contour edge appears {} times instead of once",
                edge_count[&edge]
            );
            edge_count.remove(&edge);
        }

        for (_, count) in &edge_count {
            assert_eq!(
                *count, 2,
                "Internal edge appears {} times instead of twice",
                count
            );
        }
    }

    #[test]
    fn test_triangulate_triangle_upper() {
        let upper = vec![
            Vector2f::new(0.0, 0.0),
            Vector2f::new(1.0, 1.0),
            Vector2f::new(2.0, 0.0),
        ];
        let mut triangles = Vec::new();
        triangulate_monotone_region(&mut triangles, &upper);
        validate_triangulation(&upper, &triangles);
    }

    #[test]
    fn test_triangulate_triangle_lower() {
        let lower = vec![
            Vector2f::new(0.0, 0.0),
            Vector2f::new(2.0, 0.0),
            Vector2f::new(1.0, -1.0),
        ];
        let mut triangles = Vec::new();
        triangulate_monotone_region(&mut triangles, &lower);
        validate_triangulation(&lower, &triangles);
    }

    #[test]
    fn test_triangulate_diamond() {
        let diamond = vec![
            Vector2f::new(0.0, 2.0),
            Vector2f::new(1.0, 3.0),
            Vector2f::new(2.0, 2.0),
            Vector2f::new(1.0, 1.0),
        ];
        let mut triangles = Vec::new();
        triangulate_monotone_region(&mut triangles, &diamond);
        validate_triangulation(&diamond, &triangles);
    }

    #[test]
    fn test_triangulate_region() {
        let region = vec![
            // start
            Vector2f::new(0.0, 0.0),
            // upper chain
            Vector2f::new(1.0, 1.0),
            Vector2f::new(2.0, 3.0),
            Vector2f::new(3.0, 6.0),
            Vector2f::new(4.0, 5.0),
            Vector2f::new(5.0, 1.0),
            Vector2f::new(6.0, 3.0),
            // end
            Vector2f::new(7.0, 0.0),
            // lower chain
            Vector2f::new(5.0, 0.0),
        ];
        let mut triangles = Vec::new();
        triangulate_monotone_region(&mut triangles, &region);
        validate_triangulation(&region, &triangles);
    }
}
