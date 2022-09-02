fn boys(n: i16, t: f64) -> f64 {
    let mut f;

    // ----------------------------------------------
    //   Upward recursion for large arguments
    //   ----------------------------------------------
    if t > 30.0 {
        let pi_sqrt = 0.88622692545275801365;
        f = pi_sqrt / t.sqrt();
        for i in 1..=n {
            f *= (f64::from(i) - 0.5) / t;
        }
    }
    // ----------------------------------------------
    //   Downward recursion for small arguments
    //   ----------------------------------------------
    else {
        let pt5_pm = f64::from(n) + 0.5;
        let mut term = 0.5 / pt5_pm;
        let mut sum = term;
        let accuracy = 1e-12;

        let mut i = 1;
        loop {
            term *= t / (pt5_pm + f64::from(i));
            sum += term;
            i += 1;

            if term < accuracy || i > 1000 {
                break;
            }
        }
        f = (-t).exp() * sum;
    }

    f
}

pub fn r(t: i16, u: i16, v: i16, n: i16, p: f64, p_origin: &[f64; 3], distance: &f64) -> f64 {
    let mut value = 0.0;

    if t == 0 && u == 0 && v == 0 {
        let pr2 = p * distance * distance;
        if n == 0 {
            value += boys(n, pr2);
        } else {
            value += (-2.0 * p).powf(n.into()) * boys(n, pr2);
        }
    } else if t == 0 && u == 0 {
        if v > 1 {
            value += f64::from(v - 1) * r(t, u, v - 2, n + 1, p, p_origin, distance);
        }
        value += p_origin[2] * r(t, u, v - 1, n + 1, p, p_origin, distance);
    } else if t == 0 {
        if u > 1 {
            value += f64::from(u - 1) * r(t, u - 2, v, n + 1, p, p_origin, distance);
        }
        value += p_origin[1] * r(t, u - 1, v, n + 1, p, p_origin, distance);
    } else {
        if t > 1 {
            value += f64::from(t - 1) * r(t - 2, u, v, n + 1, p, p_origin, distance);
        }
        value += p_origin[0] * r(t - 1, u, v, n + 1, p, p_origin, distance);
    }

    value
}
