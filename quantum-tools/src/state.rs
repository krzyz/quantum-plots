extern crate nalgebra;
extern crate num;

use crate::braket::*;
use crate::basis::Basis;
use nalgebra::{convert, ComplexField, DVector, RowDVector};
use num::complex::Complex64;
use simba::scalar::SubsetOf;

pub type State = Basis<Braket>;

impl State {
    pub fn len(&self) -> usize {
        self.get_value().len()
    }

    pub fn get_braket(&self) -> &Braket {
        self.get_value()
    }

    pub fn problem_ket(vec: DVector<Complex64>) -> State {
       Basis::Problem(Braket::Ket(vec.normalize()))
    }

    pub fn solution_ket(vec: DVector<Complex64>) -> State {
        Basis::Solution(Braket::Ket(vec.normalize()))
    }

    pub fn problem_bra(vec: RowDVector<Complex64>) -> State {
        Basis::Problem(Braket::Bra(vec.normalize()))
    }

    pub fn solution_bra(vec: RowDVector<Complex64>) -> State {
        Basis::Solution(Braket::Bra(vec.normalize()))
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