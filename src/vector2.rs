use std::cmp::Ordering;

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vector2f {
    pub x: f32,
    pub y: f32,
}

pub type Vector2fBits = (u32, u32);

impl Vector2f {
    pub const ZERO: Vector2f = Vector2f { x: 0.0, y: 0.0 };

    pub fn new(x: f32, y: f32) -> Vector2f {
        Vector2f { x, y }
    }

    pub fn to_bits(v: Vector2f) -> Vector2fBits {
        (f32::to_bits(v.x), f32::to_bits(v.y))
    }

    pub fn from_bits(v: Vector2fBits) -> Vector2f {
        Vector2f::new(f32::from_bits(v.0), f32::from_bits(v.1))
    }

    pub fn dot(&self, other: Vector2f) -> f32 {
        self.x * other.x + self.y * other.y
    }

    pub fn cross(&self, other: Vector2f) -> f32 {
        self.x * other.y - self.y * other.x
    }

    pub fn length(&self) -> f32 {
        (self.x * self.x + self.y * self.y).sqrt()
    }

    pub fn angle(&self) -> f32 {
        self.y.atan2(self.x)
    }

    pub fn length_squared(&self) -> f32 {
        self.x * self.x + self.y * self.y
    }

    pub fn normalized(&self) -> Vector2f {
        let length = self.length();
        Vector2f {
            x: self.x / length,
            y: self.y / length,
        }
    }
}

impl std::ops::Add for Vector2f {
    type Output = Vector2f;

    fn add(self, other: Vector2f) -> Vector2f {
        Vector2f {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl std::ops::Sub for Vector2f {
    type Output = Vector2f;

    fn sub(self, other: Vector2f) -> Vector2f {
        Vector2f {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl std::ops::Mul<f32> for Vector2f {
    type Output = Vector2f;

    fn mul(self, scalar: f32) -> Vector2f {
        Vector2f {
            x: self.x * scalar,
            y: self.y * scalar,
        }
    }
}

impl std::ops::Div<f32> for Vector2f {
    type Output = Vector2f;

    fn div(self, scalar: f32) -> Vector2f {
        Vector2f {
            x: self.x / scalar,
            y: self.y / scalar,
        }
    }
}

impl std::ops::Neg for Vector2f {
    type Output = Vector2f;

    fn neg(self) -> Vector2f {
        Vector2f {
            x: -self.x,
            y: -self.y,
        }
    }
}

impl std::ops::AddAssign for Vector2f {
    fn add_assign(&mut self, other: Vector2f) {
        *self = Vector2f {
            x: self.x + other.x,
            y: self.y + other.y,
        };
    }
}

impl std::ops::SubAssign for Vector2f {
    fn sub_assign(&mut self, other: Vector2f) {
        *self = Vector2f {
            x: self.x - other.x,
            y: self.y - other.y,
        };
    }
}

impl std::ops::MulAssign<f32> for Vector2f {
    fn mul_assign(&mut self, scalar: f32) {
        *self = Vector2f {
            x: self.x * scalar,
            y: self.y * scalar,
        };
    }
}

impl std::ops::DivAssign<f32> for Vector2f {
    fn div_assign(&mut self, scalar: f32) {
        *self = Vector2f {
            x: self.x / scalar,
            y: self.y / scalar,
        };
    }
}

impl std::fmt::Display for Vector2f {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "({}, {})", self.x, self.y)
    }
}

impl Into<raylib::core::math::Vector2> for Vector2f {
    fn into(self) -> raylib::core::math::Vector2 {
        raylib::core::math::Vector2::new(self.x, self.y)
    }
}

impl Into<raylib::ffi::Vector2> for Vector2f {
    fn into(self) -> raylib::ffi::Vector2 {
        raylib::ffi::Vector2 {
            x: self.x,
            y: self.y,
        }
    }
}

pub fn comp_points_x_dir(a: Vector2f, b: Vector2f) -> Ordering {
    return f32::total_cmp(&a.x, &b.x).then_with(|| f32::total_cmp(&a.y, &b.y));
}
