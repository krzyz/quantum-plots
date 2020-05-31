extern crate nalgebra;
extern crate num;

use crate::basis::Basis;
use nalgebra::DMatrix;
use num::complex::Complex64;

pub type Op<'a> = Basis<'a, DMatrix<Complex64>>;

/*
use nalgebra::DMatrix;
use num::complex::Complex64;
use std::ops::Mul;
use crate::state::*;
use crate::braket::*;


impl Mul<State> for Op {
    type Output = State;

    fn mul(self, rhs: State) -> Self::Output {
        use Basis::*;
        use Braket::*;

        if rhs.is_bra() {
            panic!();
        }

        match (self, rhs) {
            (Problem(ref mat), Problem(ref braket)) => match braket {
                Ket(ref vec) => Problem(Ket(mat * vec)),
                Bra(_) => panic!(),
            }
            (Problem(_), Solution(_)) => self.to_solution() * rhs,
            (Solution(ref mat), Problem(ref vec)) => self.to_problem() * vec,
            (Solution(ref mat), Solution(ref braket)) => match braket {
                Ket(ref vec) => Solution(Ket(mat * vec)),
                Bra(_) => panic!(),
            }
        }
    }


}

*/
