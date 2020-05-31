extern crate nalgebra;
extern crate num;

use nalgebra::{DMatrix, Dynamic, SymmetricEigen};
use num::complex::Complex64;

#[derive(Debug)]
pub struct System {
    pub ham: DMatrix<Complex64>,
    pub eigensystem: SymmetricEigen<Complex64, Dynamic>,
}

impl System {
    pub fn new(ham: DMatrix<Complex64>) -> System {
        let eigensystem = ham.clone().symmetric_eigen();

        System { ham, eigensystem }
    }
}
