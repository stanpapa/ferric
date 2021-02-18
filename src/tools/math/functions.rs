pub fn factorial(n: i32) -> i32 {
  (1..=n).product() 
}

pub fn factorial2(n: i32) -> i32{
  (1..=n).rev().step_by(2).product()
}

pub fn binomial_coeff(n: i32, k: i32) -> i32 {
  if k <= 0 || k > n {
    panic!("k is larger than n!");
  }
  factorial(n) / (factorial(k) * factorial(n-k))
}

#[cfg(test)]
mod tests {
    use super::*;

  #[test]
  fn test_factorial() {
    assert_eq!(1, factorial(0));
    assert_eq!(1, factorial(1));
    assert_eq!(2, factorial(2));
    assert_eq!(6, factorial(3));
    assert_eq!(24, factorial(4));
    assert_eq!(120, factorial(5));
  }

  #[test]
  fn test_factorial2() {
    assert_eq!(1, factorial2(0));
    assert_eq!(1, factorial2(1));
    assert_eq!(2, factorial2(2));
    assert_eq!(3, factorial2(3));
    assert_eq!(8, factorial2(4));
    assert_eq!(15, factorial2(5));
    assert_eq!(48, factorial2(6));
  }

  #[test]
  fn test_binomial_coeff() {
    assert_eq!(10, binomial_coeff(5,2));
  }
}