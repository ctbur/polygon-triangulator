use crate::vector2::Vector2f;

pub type Countour = Vec<Vector2f>;

pub struct Polygon {
    contours: Vec<Countour>,
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

    pub fn contours(&self) -> &Vec<Countour> {
        &self.contours
    }
}

struct DiscretePolygon {
    contours: Vec<Vec<(i64, i64)>>,
}

fn discretize(polygon: &Polygon, scale_factor: f32) -> DiscretePolygon {
    let mut contours = Vec::new();

    for contour in &polygon.contours {
        let mut discrete_contour = Vec::new();

        for p in contour {
            let xd = (p.x * scale_factor) as i64;
            let yd = (p.y * scale_factor) as i64;
            discrete_contour.push((xd, yd));
        }

        contours.push(discrete_contour);
    }

    DiscretePolygon { contours }
}
