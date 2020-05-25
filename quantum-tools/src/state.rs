extern crate nalgebra;
extern crate num;

use crate::braket::*;
use nalgebra::{convert, ComplexField, DVector, RowDVector};
use num::complex::Complex64;
use simba::scalar::SubsetOf;

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

    pub fn get_probabilities(&self) -> Vec<f64> {
        use Braket::*;
        let braket = self.get_braket();

        match braket {
            Bra(ref vec) => vec.iter().map(|x| x.norm_sqr()).collect(),
            Ket(ref vec) => vec.iter().map(|x| x.norm_sqr()).collect(),
        }
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

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;

    #[fixture]
    pub fn solution_psi_2() -> State {
        State::solution_ket_from_slice(&vec![0.2, 0.5])
    }

    #[rstest]
    fn get_probabilities_test(solution_psi_2: State) {
        let probabilities = solution_psi_2.get_probabilities();
        
        assert_eq!(
            probabilities,
            vec![0.13793103448275856, 0.8620689655172411]
        )
    }
}