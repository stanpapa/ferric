pub fn e(i: i16, j: i16, t: i16, center: f64, a: &f64, b: &f64) -> f64 {
    let p = a + b;
    let u = a * b / p;

    if t < 0 || t > i + j {
        return 0.0;
    } else if i == 0 && j == 0 && t == 0 {
        return (-u * center * center).exp();
    } else if j == 0 {
        return (1.0 / (2.0 * p)) * e(i - 1, j, t - 1, center, a, b)
            - (u * center / a) * e(i - 1, j, t, center, a, b)
            + f64::from(t + 1) * e(i - 1, j, t + 1, center, a, b);
    } else {
        return (1.0 / (2.0 * p)) * e(i, j - 1, t - 1, center, a, b)
            + (u * center / b) * e(i, j - 1, t, center, a, b)
            + f64::from(t + 1) * e(i, j - 1, t + 1, center, a, b);
    }
}
