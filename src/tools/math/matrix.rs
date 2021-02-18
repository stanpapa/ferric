pub type FMatrix = Matrix<f64>;
pub type IMatrix = Matrix<i64>;

pub struct Matrix<T> {
  rows: usize,
  columns: usize,
  elem: Vec<Vec<T>>,
}

// figure out how to properly initialize to 0 using T
impl<T> Matrix<T> {
  // pub fn new(n: &usize) -> Matrix<T> {
  //   if *n > 0 {
  //     Matrix {
  //       rows: *n,
  //       columns: *n,
  //       elem: vec![vec![0.0; *n]; *n],
  //     }
  //   } else {
  //     panic!("Matrix dimensions must be larger than 0. ABORT.");
  //   }
  // } // new()

  // pub fn init(&mut self) {
  //   for i in 0..self.rows {
  //     for j in 0..self.columns {
  //       self.elem[i][j] = 0.0;
  //     }
  //   }
  // } // init()
  // pub fn::print(&self) {
  //   for i in 0..self.rows {
  //     for j in i..columns {
  //       if 
  //     }
  //   }
  // }
} // impl<T> Matrix<T>

impl FMatrix {
  pub fn new(n: &usize) -> FMatrix {
    if *n > 0 {
      Matrix {
        rows: *n,
        columns: *n,
        elem: vec![vec![0.0; *n]; *n],
      }
    } else {
      panic!("Matrix dimensions must be larger than 0. ABORT.");
    }
  } // new()

  pub fn init(&mut self) {
    for i in 0..self.rows {
      for j in 0..self.columns {
        self.elem[i][j] = 0.0;
      }
    }
  } // init()
} // impl FMatrix


impl IMatrix {
  pub fn new(n: &usize) -> IMatrix {
    if *n > 0 {
      Matrix {
        rows: *n,
        columns: *n,
        elem: vec![vec![0; *n]; *n],
      }
    } else {
      panic!("Matrix dimensions must be larger than 0. ABORT.");
    }
  } // new()

  pub fn init(&mut self) {
    for i in 0..self.rows {
      for j in 0..self.columns {
        self.elem[i][j] = 0;
      }
    }
  } // init()
} // impl IMatrix