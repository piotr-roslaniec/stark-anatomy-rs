use num::{BigInt, One, Zero};

/// Extended Euclidian algorithm
pub fn xgcd<T: Into<BigInt>>(a: T, b: T) -> (BigInt, BigInt, BigInt) {
    let a = a.into();
    let b = b.into();
    let (mut old_r, mut r) = (a, b);
    let (mut old_s, mut s) = (BigInt::one(), BigInt::zero());
    let (mut old_t, mut t) = (BigInt::zero(), BigInt::one());

    while r != BigInt::zero() {
        let quotient = old_r.clone() / r.clone(); // Use integer division
        (old_r, r) = (r.clone(), old_r - quotient.clone() * r.clone());
        (old_s, s) = (s.clone(), old_s - quotient.clone() * s);
        (old_t, t) = (t.clone(), old_t - quotient * t);
    }

    // a, b, gcd
    (old_s, old_t, old_r)
}

/// Represents an element in a finite field.
///
/// Each field element consists of a value and the field it belongs to.
/// The value is automatically normalized modulo the field's characteristic.
#[derive(Debug, Clone, PartialEq)]
pub struct FieldElement {
    /// The value of the field element, normalized modulo the field's characteristic
    pub value: BigInt,
    /// The field this element belongs to
    pub field: Field,
}

impl FieldElement {
    /// Creates a new field element with the given value in the specified field.
    ///
    /// The value is automatically normalized modulo the field's characteristic.
    ///
    /// # Arguments
    /// * `value` - The value to create the field element from
    /// * `field` - The field this element belongs to
    ///
    /// # Example
    /// ```
    /// use stark_anatomy_rs::field::{Field, FieldElement};
    /// let field = Field::new(7);
    /// let element = FieldElement::new(10, &field);  // Will be normalized to 3 mod 7
    /// ```
    pub fn new(value: impl Into<BigInt>, field: &Field) -> FieldElement {
        let value = value.into();
        let normalized = ((value % field.p.clone()) + field.p.clone()) % field.p.clone();
        FieldElement {
            value: normalized,
            field: field.clone(),
        }
    }

    /// Adds two field elements.
    ///
    /// # Arguments
    /// * `right` - The field element to add to this one
    ///
    /// # Returns
    /// A new field element representing the sum
    pub fn add(&self, right: &FieldElement) -> FieldElement {
        self.field.clone().add(self, right)
    }

    /// Multiplies two field elements.
    ///
    /// # Arguments
    /// * `right` - The field element to multiply with this one
    ///
    /// # Returns
    /// A new field element representing the multiplication
    pub fn multiply(&self, right: &FieldElement) -> FieldElement {
        self.field.clone().multiply(self, right)
    }

    /// Subtracts two field elements.
    ///
    /// # Arguments
    /// * `right` - The field element to subtract from this one
    ///
    /// # Returns
    /// A new field element representing the subtraction
    pub fn subtract(&self, right: FieldElement) -> FieldElement {
        self.field.clone().subtract(self, &right)
    }

    /// Divides two field elements.
    ///
    /// # Arguments
    /// * `right` - The field element to use as divisor for this one
    ///
    /// # Returns
    /// A new field element representing the division
    pub fn divide(&self, right: FieldElement) -> FieldElement {
        self.field.clone().divide(self, &right)
    }

    /// Negates a field element.
    ///
    /// # Returns
    /// A new field element representing the negated element
    pub fn negate(&self) -> FieldElement {
        self.field.clone().negate(self)
    }

    /// Inverts a field element.
    ///
    /// # Returns
    /// A new field element representing the inverted element
    pub fn invert(&self) -> FieldElement {
        self.field.clone().inverse(self)
    }

    pub fn xor(&self, exponent: i128) -> FieldElement {
        assert!(exponent >= 0, "Exponent must be non-negative");

        if exponent == 0 {
            return self.field.one();
        }

        let mut acc = self.field.one();
        let value = self.clone();

        // Use a safe bit counting approach
        let bits = (128 - exponent.leading_zeros()) as usize;
        for i in (0..bits).rev() {
            acc = acc.multiply(&acc);
            if (exponent & (1 << i)) != 0 {
                acc = acc.multiply(&value);
            }
        }
        acc
    }

    /// Raises the field element to a power using square-and-multiply algorithm.
    ///
    /// # Arguments
    /// * `exponent` - The power to raise this element to
    ///
    /// # Returns
    /// A new field element representing self^exponent
    ///
    /// # Panics
    /// Panics if the exponent is negative
    pub fn pow(&self, exponent: u32) -> Self {
        if exponent == 0 {
            return self.field.one();
        }
        let mut base = self.clone();
        let mut result = self.field.one();
        let mut exp = exponent;

        while exp > 0 {
            if exp & 1 == 1 {
                result = result.multiply(&base);
            }
            base = base.multiply(&base);
            exp >>= 1;
        }
        result
    }

    pub fn equal(&self, other: FieldElement) -> bool {
        self.value == other.value
    }

    pub fn nequal(&self, other: FieldElement) -> bool {
        !self.equal(other)
    }

    pub fn is_zero(&self) -> bool {
        self.value == BigInt::zero()
    }
}

// TODO: Implement standard traits
// impl std::ops::Add for &FieldElement {
//     type Output = FieldElement;
//
//     fn add(self, rhs: Self) -> Self::Output {
//         self.add(rhs)
//     }
// }
//
// impl std::ops::Mul for &FieldElement {
//     type Output = FieldElement;
//
//     fn mul(self, rhs: Self) -> Self::Output {
//         self.multiply(rhs)
//     }
// }

#[derive(Debug, Clone, PartialEq)]
pub struct Field {
    pub(crate) p: BigInt,
}

impl Field {
    /// Creates a new field with the given characteristic.
    ///
    /// # Arguments
    /// * `p` - The characteristic of the field
    ///
    /// # Panics
    /// Panics if p is not positive
    pub fn new(p: impl Into<BigInt>) -> Self {
        let p = p.into();
        assert!(p > BigInt::zero(), "Field characteristic must be positive");
        Self { p }
    }

    pub fn zero(&self) -> FieldElement {
        FieldElement::new(0, self)
    }

    pub fn one(&self) -> FieldElement {
        FieldElement::new(1, self)
    }

    pub fn multiply(self, left: &FieldElement, right: &FieldElement) -> FieldElement {
        FieldElement::new(
            ((left.value.clone() * right.value.clone()) % self.p.clone() + self.p.clone())
                % self.p.clone(),
            &self,
        )
    }

    pub fn add(self, left: &FieldElement, right: &FieldElement) -> FieldElement {
        FieldElement::new(
            ((left.value.clone() + right.value.clone()) % self.p.clone() + self.p.clone())
                % self.p.clone(),
            &self,
        )
    }

    pub fn subtract(self, left: &FieldElement, right: &FieldElement) -> FieldElement {
        FieldElement::new(
            ((left.value.clone() - right.value.clone()) % self.p.clone() + self.p.clone())
                % self.p.clone(),
            &self,
        )
    }

    pub fn negate(self, operand: &FieldElement) -> FieldElement {
        FieldElement::new(
            (self.p.clone() - operand.value.clone()) % self.p.clone(),
            &self,
        )
    }

    pub fn inverse(self, operand: &FieldElement) -> FieldElement {
        let (a, _, g) = xgcd(operand.value.clone(), self.p.clone());
        assert_eq!(g, BigInt::one(), "Element has no inverse");
        // Normalize result
        FieldElement::new(
            (a % self.p.clone() + self.p.clone()) % self.p.clone(),
            &self,
        )
    }

    pub fn divide(self, left: &FieldElement, right: &FieldElement) -> FieldElement {
        assert!(!right.is_zero(), "Cannot divide by zero");
        let (a, _, _) = xgcd(right.clone().value, self.p.clone());
        let result = (left.clone().value * a) % self.p.clone();
        // Ensure the result is positive
        FieldElement::new((result + self.p.clone()) % self.p.clone(), &self)
    }

    // TODO: Is that so? Where does that p value come from?
    /// Returns concrete implementation of prime field p= 1 + 11 * 37 * 2^119
    pub fn num_elements_p() -> BigInt {
        // p = 1 + 11 * 37 * 2^119
        let two = BigInt::from(2);
        let power = two.pow(119);
        let factor = BigInt::from(11 * 37);
        BigInt::one() + (factor * power)
    }

    /// This is the only Field we're going to use.
    pub fn main(&self) -> Field {
        Field::new(Field::num_elements_p())
    }

    /// Supply the user with a generator for the entire multiplicative group, as well as the power-of-two subgroups.
    pub fn generator(&self) -> FieldElement {
        assert_eq!(
            self.p,
            Field::num_elements_p(),
            "Do not know generator for other fields beyond 1+407*2^119"
        );
        FieldElement::new(2, self)
    }

    pub fn primitive_nth_root(self, n: i128) -> FieldElement {
        assert!(n > 0, "n must be positive");

        // We need to find an element x where x^n ≡ 1 (mod p)
        let mut candidate = FieldElement::new(3, &self);

        while (0..n - 1).fold(candidate.clone(), |acc, _| acc.multiply(&candidate)) != self.one() {
            let value = (candidate.value + BigInt::one()) % self.p.clone();
            candidate = FieldElement::new(value, &self);

            if candidate.value == BigInt::from(2) {
                panic!("No {}th root of unity exists in this field", n);
            }
        }

        candidate
    }
    /// Sample field elements. Can be used randomlny and pseudorandomly.
    pub fn sample(&self, byte_array: Vec<u8>) -> FieldElement {
        let mut accumulator = 0;
        for byte in byte_array {
            accumulator = (accumulator << 8) ^ byte as i128;
        }
        FieldElement::new(accumulator % self.p.clone(), self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_xgcd() {
        assert_eq!(
            xgcd(56, 15),
            (BigInt::from(-4), BigInt::from(15), BigInt::from(1))
        );
        assert_eq!(
            xgcd(101, 10),
            (BigInt::from(1), BigInt::from(-10), BigInt::from(1))
        );
    }

    #[test]
    fn test_field_element_addition() {
        let field = Field::new(7);
        let a = FieldElement::new(3, &field);
        let b = FieldElement::new(5, &field);
        let result = a.add(&b);
        assert_eq!(result, FieldElement::new(1, &field)); // (3 + 5) % 7 = 1
    }

    #[test]
    fn test_field_element_subtraction() {
        let field = Field::new(7);
        let a = FieldElement::new(3, &field);
        let b = FieldElement::new(5, &field);
        let result = a.subtract(b);
        assert_eq!(result, FieldElement::new(5, &field)); // (3 - 5) % 7 = 5
    }

    #[test]
    fn test_field_element_multiplication() {
        let field = Field::new(7);
        let a = FieldElement::new(3, &field);
        let b = FieldElement::new(5, &field);
        let result = a.multiply(&b);
        assert_eq!(result, FieldElement::new(1, &field)); // (3 * 5) % 7 = 15 % 7 = 1
    }

    #[test]
    fn test_field_element_division() {
        let field = Field::new(7);
        let a = FieldElement::new(6, &field); // (6 / 3) % 7
        let b = FieldElement::new(3, &field);
        let result = a.divide(b);
        assert_eq!(result, FieldElement::new(2, &field)); // (6 / 3) = 2 mod 7
    }

    #[test]
    fn test_field_element_inversion() {
        let field = Field::new(7);
        let a = FieldElement::new(3, &field);
        let result = a.invert();
        assert_eq!(result, FieldElement::new(5, &field)); // 3 * 5 ≡ 1 mod 7
    }

    #[test]
    fn test_field_element_is_zero() {
        let field = Field::new(7);
        let zero = FieldElement::new(0, &field);
        assert!(zero.is_zero());

        let non_zero = FieldElement::new(1, &field);
        assert!(!non_zero.is_zero());
    }

    #[test]
    fn test_num_elements_p() {
        let p = Field::num_elements_p();
        // Verify p is what we expect: 1 + 11 * 37 * 2^119
        let expected = {
            let two = BigInt::from(2);
            let power = two.pow(119);
            let factor = BigInt::from(11 * 37);
            BigInt::one() + (factor * power)
        };
        assert_eq!(p, expected);
    }

    #[test]
    fn test_field_generator() {
        let field = Field::new(Field::num_elements_p());
        let generator = field.generator();

        // Test that generator is a valid field element
        assert!(generator.value >= BigInt::zero() && generator.value < field.p);

        // Test that generator is not trivial
        assert!(!generator.is_zero());
        assert_ne!(generator, field.one());

        // Test that generator actually belongs to the field
        let p = Field::num_elements_p();
        assert_eq!(generator.field.p, p);
    }

    #[test]
    fn test_field_primitive_nth_root() {
        let field = Field::new(7);
        let n: i128 = 2; // Testing for square root of unity
        let root = field.clone().primitive_nth_root(n);

        // Test that root^n = 1
        let mut result = field.one();
        for _ in 0..n {
            result = result.multiply(&root);
        }
        assert_eq!(result, field.one());
    }
    #[test]
    fn test_field_zero_operations() {
        let field = Field::new(7);
        let zero = field.zero();
        let a = FieldElement::new(3, &field);

        // Zero addition
        assert_eq!(a.add(&zero), a);
        assert_eq!(zero.add(&a), a);

        // Zero multiplication
        assert_eq!(a.multiply(&zero), zero);
        assert_eq!(zero.multiply(&a), zero);

        // Zero subtraction
        assert_eq!(a.subtract(zero.clone()), a);
        assert_eq!(zero.subtract(a.clone()), a.negate());
    }

    #[test]
    fn test_field_one_operations() {
        let field = Field::new(7);
        let one = field.one();
        let a = FieldElement::new(3, &field);

        // Multiplication by one
        assert_eq!(a.multiply(&one), a);
        assert_eq!(one.multiply(&a), a);

        // Division by one
        assert_eq!(a.divide(one.clone()), a);
    }

    #[test]
    #[should_panic(expected = "Cannot divide by zero")]
    fn test_division_by_zero() {
        let field = Field::new(7);
        let a = FieldElement::new(3, &field);
        let zero = field.zero();
        a.divide(zero);
    }

    #[test]
    fn test_field_negative_values() {
        let field = Field::new(7);
        let a = FieldElement::new(-3, &field);
        let b = FieldElement::new(4, &field); // -3 ≡ 4 (mod 7)
        assert_eq!(a, b);

        // Test negative value multiplication
        let c = FieldElement::new(-2, &field);
        assert_eq!(a.multiply(&c), FieldElement::new(6, &field)); // (-3 * -2) ≡ 6 (mod 7)
    }

    #[test]
    fn test_field_large_values() {
        let field = Field::new(7);
        let a = FieldElement::new(15, &field); // 15 ≡ 1 (mod 7)
        let b = FieldElement::new(1, &field);
        assert_eq!(a, b);

        // Test large value arithmetic
        let c = FieldElement::new(100, &field); // 100 ≡ 2 (mod 7)
        assert_eq!(a.multiply(&c), FieldElement::new(2, &field));
    }

    #[test]
    fn test_field_inverse_properties() {
        let field = Field::new(7);

        // Test inverse for all non-zero elements in F_7
        for i in 1..7 {
            let a = FieldElement::new(i, &field);
            let inv = a.invert();
            assert_eq!(a.multiply(&inv), field.one());
        }
    }

    #[test]
    fn test_xor_edge_cases() {
        let field = Field::new(7);
        let a = FieldElement::new(3, &field);

        // Test x^0 = 1
        assert_eq!(a.xor(0), field.one());

        // Test x^1 = x
        assert_eq!(a.xor(1), a);

        // Test x^2
        let expected = a.multiply(&a);
        assert_eq!(a.xor(2), expected);

        // Test larger exponents
        let a_cubed = a.multiply(&a).multiply(&a);
        assert_eq!(a.xor(3), a_cubed);
    }

    #[test]
    #[should_panic(expected = "Exponent must be non-negative")]
    fn test_xor_negative_exponent() {
        let field = Field::new(7);
        let a = FieldElement::new(3, &field);
        a.xor(-1);
    }

    #[test]
    fn test_xor_large_exponent() {
        let field = Field::new(7);
        let a = FieldElement::new(3, &field);
        // Test with a larger exponent
        a.xor(10); // Should not panic
    }
}
