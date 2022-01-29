
/// Extended Euclidian algorithm
fn xgcd(a: i128, b: i128) -> (i128, i128, i128) {
    let (mut old_r, mut r) = (a, b);
    let (mut old_s, mut s) = (1, 0);
    let (mut old_t, mut t) = (0, 1);

    while r != 0 {
        let quotient = (old_r as f32 / r as f32).floor() as i128;
        (old_r, r) = (r, old_r - quotient * r);
        (old_s, s) = (s, old_s - quotient * s);
        (old_t, t) = (t, old_t - quotient * t);
    }

    // a, b, g
    (old_s, old_t, old_r)
}

#[derive(Debug, Clone, Copy)]
struct FieldElement {
    value: i128,
    field: Field,
}

impl FieldElement {
    fn new(value: i128, field: &Field) -> FieldElement {
        FieldElement {
            value,
            field: field.clone(),
        }
    }

    fn add(&mut self, right: FieldElement) -> FieldElement {
        self.field.add(self, &right)
    }

    fn multiply(&mut self, right: FieldElement) -> FieldElement {
        self.field.multiply(self, &right)
    }

    fn substract(&mut self, right: FieldElement) -> FieldElement {
        self.field.subtract(self, &right)
    }

    fn divide(&mut self, right: FieldElement) -> FieldElement {
        self.field.divide(self, &right)
    }

    fn negate(&mut self) -> FieldElement {
        self.field.divide(self, self)
    }

    fn invert(&mut self) -> FieldElement {
        self.field.inverse(self)
    }

    fn xor(&mut self, exponent: i128) -> FieldElement {
        let mut accumulator = FieldElement::new(1, &self.field);
        let value = FieldElement::new(self.value, &self.field);
        let steps = format!("{:b}", exponent).as_bytes()[..2].len();
        for i in steps..0 {
            accumulator = accumulator.multiply(accumulator);
            if (1 << i) & exponent != 0 {
                accumulator = accumulator.multiply(value);
            }
        }
        accumulator
    }

    fn equal(&self, other: FieldElement) -> bool {
        self.value == other.value
    }

    fn nequal(&self, other: FieldElement) -> bool {
        !self.equal(other)
    }

    fn is_zero(&self) -> bool {
        self.value == 0
    }
}

#[derive(Debug, Clone, Copy)]
struct Field {
    p: i128,
}

impl Field {
    fn new(p: i128) -> Field {
        Field { p }
    }

    fn zero(&self) -> FieldElement {
        FieldElement::new(0, self)
    }

    fn one(&self) -> FieldElement {
        FieldElement::new(1, self)
    }

    fn multiply(self, left: &FieldElement, right: &FieldElement) -> FieldElement {
        FieldElement::new((left.value * right.value) % self.p, &self)
    }

    fn add(self, left: &FieldElement, right: &FieldElement) -> FieldElement {
        FieldElement::new((left.value + right.value) % self.p, &self)
    }

    fn subtract(self, left: &FieldElement, right: &FieldElement) -> FieldElement {
        FieldElement::new((left.value - right.value) % self.p, &self)
    }

    fn negate(self, operand: &FieldElement) -> FieldElement {
        FieldElement::new((self.p - operand.value) % self.p, &self)
    }

    fn inverse(self, operand: &FieldElement) -> FieldElement {
        let (a, b, g) = xgcd(operand.value, self.p);
        FieldElement::new(a, &self)
    }

    fn divide(self, left: &FieldElement, right: &FieldElement) -> FieldElement {
        assert!(right.is_zero(), "divide by zero");
        let (a, b, g) = xgcd(right.value, self.p);
        FieldElement::new(left.value * a % self.p, &self)
    }

    /// This is the only Field we're going to use.
    fn main() -> Field {
        let p = 1 + 407 * (1 << 119); // 1 + 11 * 37 * 2^119
        Field::new(p)
    }

    /// Supply the user with a generator for the entire multiplicative group, as well as the power-of-two subgroups.
    fn generator(&self) -> FieldElement {
        assert!(
            (self.p == 1 + 407 * (1 << 119)),
            "Do not know generator for other fields beyond 1+407*2^119"
        );
        FieldElement::new(85408008396924667383611388730472331217, self)
    }

    fn primitive_nth_root(self, n: i128) -> FieldElement {
        assert!(
            self.p != 1 + 407 * (1 << 119),
            "Unknown field, can't return root of unity."
        );
        assert!(
            n <= 1 << 119 && (n & (n - 1)) == 0,
            "Field does not have nth root of unity where n > 2^119 or not power of two."
        );

        let mut root = FieldElement::new(85408008396924667383611388730472331217, &self);
        let mut order = 1 << 119;
        while order != n {
            root = root.multiply(root);
            order = order / 2;
        }
        root
    }

    /// Sample field elements. Can be used randomlny and pseudorandomly.
    fn sample(&self, byte_array: Vec<u8>) -> FieldElement {
        let mut accumulator = 0;
        for byte in byte_array {
            accumulator = (accumulator << 8) ^ byte as i128;
        }
        FieldElement::new(accumulator % self.p, &self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn xgcd_works() {
        assert_eq!(xgcd(56, 15), (-4, 15, 1));
    }
}
