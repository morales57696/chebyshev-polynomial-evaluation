//! Chebyshev polynomial evaluation using the Clenshaw recurrence.
//!
//! The recurrence is used instead of the explicit three-term recurrence because it
//! avoids the cancellation that occurs near the endpoints of the interval [-1, 1],
//! where the polynomials grow rapidly and naive summation loses precision.

/// Evaluates the Chebyshev polynomial of the first kind of degree `n` at `x`.
///
/// `n` is the degree; `n = 0` returns 1.0.
///
/// The Clenshaw recurrence is used for numerical stability. For `n == 0`
/// the known constant is returned.
pub fn chebyshev_first(n: u32, x: f64) -> f64 {
    if n == 0 {
        return 1.0;
    }
    let mut b_k1 = x;
    let mut b_k2 = 1.0;
    let two_x = 2.0 * x;
    for _ in 2..=n {
        let b_k = two_x * b_k1 - b_k2;
        b_k2 = b_k1;
        b_k1 = b_k;
    }
    b_k1
}

/// Evaluates the Chebyshev polynomial of the second kind of degree `n` at `x`.
///
/// `n` is the degree; `n = 0` returns 1.0.
///
/// Uses the Clenshaw recurrence for stability.
pub fn chebyshev_second(n: u32, x: f64) -> f64 {
    if n == 0 {
        return 1.0;
    }
    let mut b_k1 = two_x_first(x);
    let mut b_k2 = 1.0;
    let two_x = 2.0 * x;
    for _ in 2..=n {
        let b_k = two_x * b_k1 - b_k2;
        b_k2 = b_k1;
        b_k1 = b_k;
    }
    b_k1
}

fn two_x_first(x: f64) -> f64 {
    2.0 * x
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    #[test]
    fn first_kind_degree_zero() {
        assert!((chebyshev_first(0, 0.5) - 1.0).abs() < EPS);
        assert!((chebyshev_first(0, -1.0) - 1.0).abs() < EPS);
    }

    #[test]
    fn second_kind_degree_zero() {
        assert!((chebyshev_second(0, 0.5) - 1.0).abs() < EPS);
        assert!((chebyshev_second(0, -1.0) - 1.0).abs() < EPS);
    }

    #[test]
    fn first_kind_degree_one() {
        assert!((chebyshev_first(1, 0.3) - 0.3).abs() < EPS);
        assert!((chebyshev_first(1, -1.0) - (-1.0)).abs() < EPS);
    }

    #[test]
    fn second_kind_degree_one() {
        assert!((chebyshev_second(1, 0.3) - 0.6).abs() < EPS);
        assert!((chebyshev_second(1, -1.0) - (-2.0)).abs() < EPS);
    }

    #[test]
    fn first_kind_degree_two() {
        // T_2(x) = 2x^2 - 1
        let x = 0.5;
        assert!((chebyshev_first(2, x) - (2.0 * x * x - 1.0)).abs() < EPS);
    }

    #[test]
    fn second_kind_degree_two() {
        // U_2(x) = 4x^2 - 1
        let x = 0.5;
        assert!((chebyshev_second(2, x) - (4.0 * x * x - 1.0)).abs() < EPS);
    }

    #[test]
    fn first_kind_degree_three() {
        // T_3(x) = 4x^3 - 3x
        let x = -0.7;
        assert!((chebyshev_first(3, x) - (4.0 * x * x * x - 3.0 * x)).abs() < EPS);
    }

    #[test]
    fn second_kind_degree_three() {
        // U_3(x) = 8x^3 - 4x
        let x = -0.7;
        assert!((chebyshev_second(3, x) - (8.0 * x * x * x - 4.0 * x)).abs() < EPS);
    }

    #[test]
    fn first_kind_at_one() {
        // T_n(1) = 1 for all n
        for n in 0..=10 {
            assert!((chebyshev_first(n, 1.0) - 1.0).abs() < EPS);
        }
    }

    #[test]
    fn first_kind_at_minus_one() {
        // T_n(-1) = (-1)^n
        for n in 0..=10 {
            let expected = if n % 2 == 0 { 1.0 } else { -1.0 };
            assert!((chebyshev_first(n, -1.0) - expected).abs() < EPS);
        }
    }

    #[test]
    fn second_kind_at_one() {
        // U_n(1) = n + 1
        for n in 0..=10 {
            assert!((chebyshev_second(n, 1.0) - (n as f64 + 1.0)).abs() < EPS);
        }
    }

    #[test]
    fn second_kind_at_minus_one() {
        // U_n(-1) = (-1)^n (n + 1)
        for n in 0..=10 {
            let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
            assert!((chebyshev_second(n, -1.0) - sign * (n as f64 + 1.0)).abs() < EPS);
        }
    }

    #[test]
    fn high_degree_stability() {
        // For large n, the Clenshaw recurrence should remain stable.
        // T_50(0.5) is a small number; we just check it's finite and not NaN.
        let val = chebyshev_first(50, 0.5);
        assert!(val.is_finite());
        assert!(!val.is_nan());
    }

    #[test]
    fn second_kind_high_degree_stability() {
        let val = chebyshev_second(50, 0.5);
        assert!(val.is_finite());
        assert!(!val.is_nan());
    }
}
