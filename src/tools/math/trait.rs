use std::fmt::Display;

pub trait MyTrait: Copy + Clone + Display + Default {}
impl<T: Copy + Clone + Display + Default> MyTrait for T {}
