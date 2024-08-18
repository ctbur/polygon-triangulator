use crate::vector2::Vector2f;

pub type Contour = Vec<Vector2f>;

pub struct Polygon {
    contours: Vec<Contour>,
}

impl Polygon {
    pub fn new() -> Self {
        Self {
            contours: Vec::new(),
        }
    }

    pub fn move_to(&mut self, x: f32, y: f32) {
        // if there is a previous contour that has less then 3 points, remove it
        if let Some(contour) = self.contours.last() {
            if contour.len() < 3 {
                self.contours.pop();
            }
        }

        self.contours.push(vec![Vector2f::new(x, y)]);
    }

    pub fn line_to(&mut self, x: f32, y: f32) {
        if let Some(contour) = self.contours.last_mut() {
            assert!(contour.len() >= 1);

            contour.push(Vector2f::new(x, y));
        } else {
            panic!("line_to called without a previous move_to");
        }
    }

    pub fn quad_bezier_curve_to(&mut self, x1: f32, y1: f32, x2: f32, y2: f32, steps: u32) {
        if let Some(contour) = self.contours.last_mut() {
            assert!(contour.len() >= 1);

            let p0 = *contour.last().unwrap();
            let (x0, y0) = (p0.x, p0.y);

            for i in 1..=steps {
                let t = i as f32 / steps as f32;
                let x = (1.0 - t).powi(2) * x0 + 2.0 * (1.0 - t) * t * x1 + t.powi(2) * x2;
                let y = (1.0 - t).powi(2) * y0 + 2.0 * (1.0 - t) * t * y1 + t.powi(2) * y2;
                contour.push(Vector2f::new(x, y));
            }
        } else {
            panic!("bezier_curve_to called without a previous move_to");
        }
    }

    pub fn contours(&self) -> &Vec<Contour> {
        &self.contours
    }
}

/// Uses the shoelace formulate to calculate the signed polygon area.
/// The sign is positive if the region runs CCW.
pub fn calculate_region_area(region: &Contour) -> f32 {
    if region.len() < 3 {
        return 0.0;
    }

    let mut area = 0.0;

    for i in 0..region.len() - 1 {
        let p_i = region[i];
        let p_ni = region[i + 1];
        area += (p_i.y + p_ni.y) * (p_i.x - p_ni.x);
    }

    let p_i = region.last().unwrap();
    let p_ni = region[0];
    area += (p_i.y + p_ni.y) * (p_i.x - p_ni.x);

    return area / 2.0;
}

pub fn is_region_x_monotone(region: &Contour) -> bool {
    if region.len() <= 3 {
        return true;
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

    // find first chain end
    let mut current = start_idx;
    while region[current].x <= region[(current + 1) % region.len()].x {
        current = (current + 1) % region.len();

        if current == start_idx {
            // all region points have the same x coordinate
            return false;
        }
    }

    // check if second chain reaches end
    let end = current;
    current = start_idx;
    while current != end {
        let next = (current + region.len() - 1) % region.len();

        if region[current].x > region[next].x {
            // x position in chain decreases before reaching `end`
            return false;
        }

        current = next;
    }

    return true;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculate_region_area() {
        let mut region = vec![
            Vector2f::new(1.0, 1.0),
            Vector2f::new(2.0, 2.0),
            Vector2f::new(1.0, 3.0),
            Vector2f::new(0.0, 2.0),
        ];

        let area = calculate_region_area(&region);
        assert_eq!(area, 2.0);

        region.reverse();

        let rev_area = calculate_region_area(&region);
        assert_eq!(rev_area, -2.0);
    }
}
