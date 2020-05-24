extern crate nalgebra;
extern crate num;

use nalgebra::{DVector, RowDVector};
use num::complex::Complex64;

#[derive(Clone, Debug)]
pub enum Braket {
    Bra(RowDVector<Complex64>),
    Ket(DVector<Complex64>),
}

impl Braket {
    pub fn len(&self) -> usize {
        use Braket::*;

        match self {
            Bra(ref vec) => vec.len(),
            Ket(ref vec) => vec.len(),
        }
    }
}

impl PartialEq for Braket {
    fn eq(&self, other: &Self) -> bool {
        use Braket::*;

        match std::mem::discriminant(self) == std::mem::discriminant(other) {
            true => match (self, other) {
                (Bra(vec), Bra(other_vec)) => vec == other_vec,
                (Ket(vec), Ket(other_vec)) => vec == other_vec,
                _ => false,
            },
            false => false,
        }
    }
}
