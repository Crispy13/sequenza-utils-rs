/// Rounds like Python 3 `round(value, ndigits)` using banker's rounding
/// (ties to nearest even).
///
/// - `ndigits == 0`: rounds to an integer-valued `f64` with ties-to-even.
/// - `ndigits > 0`: rounds to the right of the decimal point.
/// - `ndigits < 0`: rounds to the left of the decimal point.
/// - `NaN` and `±inf` are returned unchanged.
///
/// # Examples
///
/// ```ignore
/// assert_eq!(py_round(2.5, 0), 2.0);
/// assert_eq!(py_round(3.5, 0), 4.0);
/// assert_eq!(py_round(2.675, 2), 2.67);
/// assert_eq!(py_round(1250.0, -2), 1200.0);
/// ```
pub(crate) fn py_round(value: f64, ndigits: i32) -> f64 {
    if !value.is_finite() {
        return value;
    }

    if ndigits == 0 {
        return value.round_ties_even();
    }

    if ndigits > 0 {
        if ndigits > 308 {
            return value;
        }

        let repr = format!("{:.*}", ndigits as usize, value);
        if let Ok(parsed) = repr.parse::<f64>() {
            return parsed;
        }

        let scale = 10_f64.powi(ndigits);
        if !scale.is_finite() || scale == 0.0 {
            return value;
        }

        return (value * scale).round_ties_even() / scale;
    }

    let scale = 10_f64.powi(-ndigits);
    if !scale.is_finite() || scale == 0.0 {
        return value;
    }

    let shifted = value / scale;
    if !shifted.is_finite() {
        return value;
    }

    shifted.round_ties_even() * scale
}