extern crate nalgebra;
extern crate num;

use nalgebra::{convert, ComplexField, DVector, RowDVector};
use num::complex::Complex64;
use simba::scalar::SubsetOf;
use crate::braket::*;

#[derive(Clone, Debug)]
pub enum State {
    ProblemBasis(Braket),
    SolutionBasis(Braket),
}

impl State {
    pub fn len(&self) -> usize {
        use State::*;
        match self {
            ProblemBasis(ref vec) => vec.len(),
            SolutionBasis(ref vec) => vec.len(),
        }
    }

    pub fn get_braket(&self) -> &Braket {
        use State::*;
        match self {
            ProblemBasis(ref vec) => &vec,
            SolutionBasis(ref vec) => &vec,
        }
    }

    pub fn problem_ket(vec: DVector<Complex64>) -> State {
        State::ProblemBasis(Braket::Ket(vec.normalize()))
    }

    pub fn solution_ket(vec: DVector<Complex64>) -> State {
        State::SolutionBasis(Braket::Ket(vec.normalize()))
    }

    pub fn problem_bra(vec: RowDVector<Complex64>) -> State {
        State::ProblemBasis(Braket::Bra(vec.normalize()))
    }

    pub fn solution_bra(vec: RowDVector<Complex64>) -> State {
        State::SolutionBasis(Braket::Bra(vec.normalize()))
    }

    pub fn problem_ket_from_slice<N>(data: &[N]) -> State
    where
        N: SubsetOf<Complex64> + ComplexField,
    {
        State::solution_ket(convert(DVector::from_row_slice(data)))
    }

    pub fn solution_ket_from_slice<N>(data: &[N]) -> State
    where
        N: SubsetOf<Complex64> + ComplexField,
    {
        State::solution_ket(convert(DVector::from_row_slice(data)))
    }
}

impl PartialEq for State {
    fn eq(&self, other: &Self) -> bool {
        use State::*;

        match std::mem::discriminant(self) == std::mem::discriminant(other) {
            true => match (self, other) {
                (SolutionBasis(vec), SolutionBasis(other_vec)) => vec == other_vec,
                (ProblemBasis(vec), ProblemBasis(other_vec)) => vec == other_vec,
                _ => false,
            },
            false => false,
        }
    }
}
