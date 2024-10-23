use crate::field::FieldElement;
use std::iter;

#[derive(Debug, Clone, PartialEq)]
struct Polynomial {
    coeffs: Vec<FieldElement>,
}

impl Polynomial {
    fn new(coeffs: Vec<FieldElement>) -> Polynomial {
        Polynomial { coeffs }
    }

    /// The degree is defined as the index of the last non-zero coefficient.
    /// The zero polynomial has a degree of -1.
    pub fn degree(&self) -> i32 {
        // No coefficients
        if self.coeffs.is_empty() {
            return -1;
        }

        // All coefficients equal to zero
        let zero = self.coeffs[0].field.zero();
        if self.coeffs.iter().all(|c| c == &zero) {
            return -1;
        }

        // Find coefficient with non-zero value
        let mut max_index = 0;
        for (i, coeff) in self.coeffs.iter().enumerate() {
            if !coeff.is_zero() {
                max_index = i;
            }
        }
        max_index as i32
    }

    pub fn negate(&self) -> Polynomial {
        let mut coeffs = Vec::with_capacity(self.coeffs.len());
        for coeff in &self.coeffs {
            coeffs.push(coeff.clone().negate());
        }
        Polynomial::new(coeffs)
    }
    pub fn add(&self, other: &Polynomial) -> Polynomial {
        if self.degree() == -1 {
            return other.clone();
        } else if other.degree() == -1 {
            return self.clone();
        }

        let field = &self.coeffs[0].field;
        let max_len = std::cmp::max(self.coeffs.len(), other.coeffs.len());
        let mut coeffs = iter::repeat(field.zero()).take(max_len).collect::<Vec<_>>();

        for i in 0..self.coeffs.len() {
            coeffs[i] = coeffs[i].add(&self.coeffs[i].clone());
        }
        for i in 0..other.coeffs.len() {
            coeffs[i] = coeffs[i].add(&other.coeffs[i].clone());
        }

        Polynomial::new(coeffs)
    }

    pub fn subtract(&self, other: &Polynomial) -> Polynomial {
        self.add(&other.negate())
    }

    pub fn multiply(&self, other: &Polynomial) -> Polynomial {
        if self.coeffs.is_empty() || other.coeffs.is_empty() {
            return Polynomial::new(vec![]);
        }

        let zero = self.coeffs[0].field.zero();
        let mut buffer: Vec<FieldElement> = iter::repeat(zero)
            .take(self.coeffs.len() + other.coeffs.len() - 1)
            .collect();
        for i in 0..self.coeffs.len() {
            if self.coeffs[i].is_zero() {
                // Optimization for sparse polynomials, skip zero coefficients
                continue;
            }
            for j in 0..other.coeffs.len() {
                buffer[i + j] =
                    buffer[i + j].add(&self.coeffs[i].multiply(&other.coeffs[j].clone()));
            }
        }
        Polynomial::new(buffer)
    }

    pub fn equal(&self, other: &Polynomial) -> bool {
        if self.degree() != other.degree() {
            return false;
        }

        if self.degree() == -1 {
            return true;
        }

        for i in 0..self.coeffs.len() {
            if self.coeffs[i] != other.coeffs[i] {
                return false;
            }
        }
        true
    }

    pub fn not_equal(&self, other: &Polynomial) -> bool {
        !self.equal(other)
    }

    pub fn is_zero(&self) -> bool {
        self.degree() == -1
    }

    pub fn leading_coefficient(&self) -> FieldElement {
        if self.degree() == -1 {
            return self.coeffs[0].field.zero();
        }
        self.coeffs[self.degree() as usize].clone()
    }

    /// Remove trailing zeroes from `Polynomial.coeffs`.
    fn normalize(&self) -> Polynomial {
        if self.coeffs.is_empty() {
            return self.clone();
        }

        let mut coeffs = self.coeffs.clone();
        while coeffs.len() > 1 && coeffs.last().unwrap().is_zero() {
            coeffs.pop();
        }
        Polynomial::new(coeffs)
    }

    pub fn divide(&self, denominator: &Polynomial) -> Option<(Polynomial, Polynomial)> {
        if denominator.degree() == -1 {
            return None;
        }

        if self.degree() < denominator.degree() {
            return Some((Polynomial::new(vec![]), self.clone()));
        }

        let field = &denominator.coeffs[0].field;
        let mut remainder = Polynomial::new(self.coeffs.clone());
        let mut quotient_coeffs =
            vec![field.zero(); (self.degree() - denominator.degree() + 1) as usize];

        for _ in 0..(self.degree() - denominator.degree() + 1) {
            if remainder.degree() < denominator.degree() {
                break;
            }

            let coefficient = remainder
                .leading_coefficient()
                .divide(denominator.leading_coefficient());
            let shift = (remainder.degree() - denominator.degree()) as usize;

            let mut single_term_coeffs = vec![field.zero(); shift];
            single_term_coeffs.push(coefficient.clone());
            let single_term = Polynomial::new(single_term_coeffs);

            let subtractee = single_term.multiply(denominator);

            let coeffs_len = quotient_coeffs.len();
            quotient_coeffs[coeffs_len - 1 - shift] = coefficient;
            remainder = remainder.subtract(&subtractee);
        }

        quotient_coeffs.reverse();
        // Normalize both quotient and remainder by removing trailing zeros
        let quotient = Polynomial::new(quotient_coeffs).normalize();
        let remainder = remainder.normalize();

        Some((quotient, remainder))
    }

    fn true_div(&self, other: &Polynomial) {
        let (_quotient, remainder) = self.divide(other).unwrap();
        // TODO: Return some here?
        assert!(remainder.is_zero());
    }

    fn modulus(&self, other: &Polynomial) -> Option<Polynomial> {
        if let Some((_quotient, remainder)) = self.divide(other) {
            Some(remainder)
        } else {
            None
        }
    }

    fn xor(&self, exponent: u32) -> Polynomial {
        if self.is_zero() {
            return Polynomial::new(vec![]);
        }
        if exponent == 0 {
            return Polynomial::new(vec![self.coeffs[0].field.one()]);
        }

        let mut acc = Polynomial::new(vec![self.coeffs[0].field.one()]);
        for i in format!("{:b}", exponent).len()..0 {
            acc = acc.multiply(&acc);
            if 1 << i != 0 && exponent != 0 {
                acc = acc.multiply(self);
            }
        }
        acc
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::{Field, FieldElement};

    fn poly_from_ints(coeffs: &[i32], field: &Field) -> Polynomial {
        Polynomial::new(
            coeffs
                .iter()
                .map(|&x| FieldElement::new(x, field))
                .collect(),
        )
    }

    #[test]
    fn test_polynomial_degree_edge_cases() {
        let field = Field::new(7);

        // Empty polynomial
        let empty = Polynomial::new(vec![]);
        assert_eq!(empty.degree(), -1);

        // Zero polynomial
        let zero = poly_from_ints(&[0], &field);
        assert_eq!(zero.degree(), -1);

        // Multiple zeros
        let zeros = poly_from_ints(&[0, 0, 0], &field);
        assert_eq!(zeros.degree(), -1);

        // Leading zeros
        let leading_zeros = poly_from_ints(&[1, 0, 0], &field);
        assert_eq!(leading_zeros.degree(), 0);

        // All coefficients the same
        let all_same = poly_from_ints(&[3, 3, 3], &field);
        assert_eq!(all_same.degree(), 2);
    }

    #[test]
    fn test_polynomial_addition_edge_cases() {
        let field = Field::new(7);

        // Adding zero polynomial
        let p1 = poly_from_ints(&[1, 2, 3], &field);
        let zero = Polynomial::new(vec![]);
        assert_eq!(p1.add(&zero), p1);
        assert_eq!(zero.add(&p1), p1);

        // Adding polynomials of different degrees
        let p2 = poly_from_ints(&[1, 2], &field);
        let p3 = poly_from_ints(&[1, 2, 3], &field);
        let result = p2.add(&p3);
        assert_eq!(result, poly_from_ints(&[2, 4, 3], &field));

        // Adding to get zero coefficients
        let p4 = poly_from_ints(&[1, 2, 3], &field);
        let p5 = poly_from_ints(&[6, 5, 4], &field); // -1, -2, -3 in F_7
        assert_eq!(p4.add(&p5).degree(), -1);
    }

    #[test]
    fn test_polynomial_multiplication_edge_cases() {
        let field = Field::new(7);

        // Multiply by zero
        let p1 = poly_from_ints(&[1, 2, 3], &field);
        let zero = Polynomial::new(vec![]);
        assert_eq!(p1.multiply(&zero), zero);
        assert_eq!(zero.multiply(&p1), zero);

        // Multiply by one
        let one = poly_from_ints(&[1], &field);
        assert_eq!(p1.multiply(&one), p1);

        // Multiply by x
        let x = poly_from_ints(&[0, 1], &field);
        let result = p1.multiply(&x);
        assert_eq!(result, poly_from_ints(&[0, 1, 2, 3], &field));
    }

    #[test]
    fn test_polynomial_division_edge_cases() {
        let field = Field::new(7);

        // Division by zero should return None
        let p1 = poly_from_ints(&[1, 2, 3], &field);
        let zero = Polynomial::new(vec![]);
        assert_eq!(p1.divide(&zero), None);

        // Division where degree(numerator) < degree(denominator)
        let p2 = poly_from_ints(&[1, 2], &field);
        let p3 = poly_from_ints(&[1, 2, 3], &field);
        let (quotient, remainder) = p2.divide(&p3).unwrap();
        assert_eq!(quotient, Polynomial::new(vec![]));
        assert_eq!(remainder, p2);

        // Division by monic polynomial (leading coefficient = 1)
        let num = poly_from_ints(&[1, 2, 1], &field); // x^2 + 2x + 1
        let den = poly_from_ints(&[1, 1], &field); // x + 1
        let (quotient, remainder) = num.divide(&den).unwrap();
        assert_eq!(quotient, poly_from_ints(&[1, 1], &field)); // x + 1
        assert_eq!(remainder, poly_from_ints(&[0], &field)); // 0
    }

    #[test]
    fn test_polynomial_composition() {
        let field = Field::new(7);

        // Test (ax + b)(cx + d) = acx^2 + (ad + bc)x + bd
        let p1 = poly_from_ints(&[2, 3], &field); // 3x + 2
        let p2 = poly_from_ints(&[1, 2], &field); // 2x + 1
        let result = p1.multiply(&p2);

        // Should be 6x^2 + 7x + 2 ≡ 6x^2 + 0x + 2 (mod 7)
        assert_eq!(result, poly_from_ints(&[2, 0, 6], &field));
    }

    #[test]
    fn test_polynomial_identities() {
        let field = Field::new(7);
        let p1 = poly_from_ints(&[1, 2, 3], &field);
        let p2 = poly_from_ints(&[4, 5, 6], &field);

        // Test commutativity of addition
        assert_eq!(p1.add(&p2), p2.add(&p1));

        // Test commutativity of multiplication
        assert_eq!(p1.multiply(&p2), p2.multiply(&p1));

        // Test associativity of addition
        let p3 = poly_from_ints(&[2, 3, 4], &field);
        assert_eq!(p1.add(&p2).add(&p3), p1.add(&p2.add(&p3)));

        // Test distributive property
        assert_eq!(
            p1.multiply(&p2.add(&p3)),
            p1.multiply(&p2).add(&p1.multiply(&p3))
        );
    }

    #[test]
    fn test_polynomial_division_with_remainder() {
        let field = Field::new(7);

        // Test (ax^2 + bx + c) = (dx + e)(fx + g) + r
        // Then multiply back and add remainder to verify
        let numerator = poly_from_ints(&[2, 4, 6], &field); // 6x^2 + 4x + 2
        let denominator = poly_from_ints(&[1, 1], &field); // x + 1

        let (quotient, remainder) = numerator.divide(&denominator).unwrap();

        // Verify that numerator = denominator * quotient + remainder
        let product = denominator.multiply(&quotient);
        let result = product.add(&remainder);
        assert_eq!(result, numerator);
    }
}
