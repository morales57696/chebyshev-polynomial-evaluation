# chebyshev-polynomial-evaluation

Evaluates Chebyshev polynomials of the first and second kind at real points using the Clenshaw recurrence for numerical stability.

```rust
use chebyshev_polynomial_evaluation::{chebyshev_first, chebyshev_second};

fn main() {
    // T_2(0.5) = 2*(0.5)^2 - 1 = -0.5
    assert!((chebyshev_first(2, 0.5) - (-0.5)).abs() < 1e-12);

    // U_3(0.2) = 8*(0.2)^3 - 4*(0.2) = 0.064 - 0.8 = -0.736
    assert!((chebyshev_second(3, 0.2) - (-0.736)).abs() < 1e-12);
}
```

## Why this exists

Evaluating Chebyshev polynomials with the explicit three-term recurrence
`T_{n+1}(x) = 2x T_n(x) - T_{n-1}(x)` is straightforward but loses precision
near the endpoints of `[-1, 1]`. The Clenshaw recurrence is a backward recurrence
that sums the polynomial series in a way that avoids the worst cancellation. This
crate provides a small, dependency-free implementation of that recurrence for both
kinds.

The trade-off is that Clenshaw is slightly more code than the naive recurrence and
returns the value of the polynomial directly rather than its coefficients. For
applications that need coefficients, this crate is not the right tool.

## Edge cases

- Degree `0` is handled explicitly and returns `1.0` for both kinds.
- The functions accept any real `x`, including values outside `[-1, 1]`.
- For very large degrees the result may overflow to infinity, but the recurrence
  itself does not produce NaN for finite inputs.
