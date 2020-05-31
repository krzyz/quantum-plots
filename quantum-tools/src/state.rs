extern crate nalgebra;
extern crate num;

use crate::basis::Basis;
use crate::braket::Braket;
use crate::context::Context;
use crate::system::System;
use nalgebra::{convert, ComplexField, DMatrix, DVector, RowDVector};
use num::complex::Complex64;
use simba::scalar::SubsetOf;

pub type State<'a> = Basis<'a, Braket>;

impl<'b> State<'b> {
    pub fn len(&self) -> usize {
        self.get_value().len()
    }

    pub fn get_braket(&self) -> &Braket {
        self.get_value()
    }

    pub fn is_bra(&self) -> bool {
        use Braket::*;
        let braket = self.get_braket();

        match braket {
            Bra(_) => true,
            Ket(_) => false,
        }
    }

    pub fn is_ket(&self) -> bool {
        return !self.is_bra();
    }

    pub fn problem_ket<'a>(vec: DVector<Complex64>, system: Option<&'a System>) -> State<'a> {
        Basis::Problem(Context {
            value: Braket::Ket(vec.normalize()),
            system,
        })
    }

    pub fn solution_ket<'a>(vec: DVector<Complex64>, system: Option<&'a System>) -> State<'a> {
        Basis::Solution(Context {
            value: Braket::Ket(vec.normalize()),
            system,
        })
    }

    pub fn problem_bra<'a>(vec: RowDVector<Complex64>, system: Option<&'a System>) -> State<'a> {
        Basis::Problem(Context {
            value: Braket::Bra(vec.normalize()),
            system,
        })
    }

    pub fn solution_bra<'a>(vec: RowDVector<Complex64>, system: Option<&'a System>) -> State<'a> {
        Basis::Solution(Context {
            value: Braket::Bra(vec.normalize()),
            system,
        })
    }

    pub fn problem_ket_from_slice<'a, N>(data: &[N], system: Option<&'a System>) -> State<'a>
    where
        N: SubsetOf<Complex64> + ComplexField,
    {
        State::solution_ket(convert(DVector::from_row_slice(data)), system)
    }

    pub fn solution_ket_from_slice<'a, N>(data: &[N], system: Option<&'a System>) -> State<'a>
    where
        N: SubsetOf<Complex64> + ComplexField,
    {
        State::solution_ket(convert(DVector::from_row_slice(data)), system)
    }

    pub fn get_probabilities(&self) -> Vec<f64> {
        use Braket::*;
        let braket = self.get_braket();

        match braket {
            Bra(ref vec) => vec.iter().map(|x| x.norm_sqr()).collect(),
            Ket(ref vec) => vec.iter().map(|x| x.norm_sqr()).collect(),
        }
    }

    pub fn to_solution(&self) -> State<'b> {
        use Braket::*;

        match self {
            Basis::Solution(_) => (*self).clone(),
            Basis::Problem(ref context) => {
                let system = &self.get_system().unwrap();

                match context.value {
                    Ket(ref vec) => State::solution_ket(
                        &system.eigensystem.eigenvectors.adjoint() * vec,
                        Some(system),
                    ),
                    Bra(ref vec) => {
                        State::solution_bra(vec * &system.eigensystem.eigenvectors, Some(system))
                    }
                }
            }
        }
    }

    pub fn to_problem(&self) -> State<'b> {
        use Braket::*;

        match self {
            Basis::Problem(_) => (*self).clone(),
            Basis::Solution(ref context) => {
                let system = &self.get_system().unwrap();

                match context.value {
                    Ket(ref vec) => {
                        State::problem_ket(&system.eigensystem.eigenvectors * vec, Some(system))
                    }
                    Bra(ref vec) => State::problem_bra(
                        vec * &system.eigensystem.eigenvectors.adjoint(),
                        Some(system),
                    ),
                }
            }
        }
    }

    pub fn evolve(&self, t: f64) -> Result<State<'b>, String> {
        let system = &self
            .get_system()
            .ok_or("Can't evolve a state not tied to any system!")?;

        if self.len() != system.ham.nrows() {
            return Err(format!(
                "Dimensions of ket({}) and ham({}) don't match!",
                self.len(),
                system.ham.nrows()
            ));
        }

        let result = match self {
            Basis::Solution(ref context) => {
                let ham = &system.eigensystem.eigenvalues;
                let diagonal: Vec<Complex64> = ham
                    .iter()
                    .map(|x| (-Complex64::i() * x * t).exp())
                    .collect();

                let dvec = DVector::from_row_slice(&diagonal);
                let matrix = DMatrix::from_diagonal(&dvec);

                match context.value {
                    Braket::Ket(ref vec) => State::solution_ket(matrix * vec, Some(system)),
                    Braket::Bra(ref vec) => State::solution_bra(-vec * matrix, Some(system)),
                }
            }
            Basis::Problem(_) => {
                let solution_psi = &self.to_solution();
                let result_solution_psi = &solution_psi.evolve(t)?;
                result_solution_psi.to_problem()
            }
        };

        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num::Complex;
    use rstest::*;

    #[fixture]
    pub fn solution_psi_1() -> State<'static> {
        State::solution_ket_from_slice(&vec![1.0], None)
    }

    #[fixture]
    pub fn solution_psi_2() -> State<'static> {
        State::solution_ket_from_slice(&vec![0.1, 0.5], None)
    }
    #[fixture]
    pub fn system2() -> System {
        let ham = DMatrix::from_row_slice(
            2,
            2,
            &[
                Complex::new(0., 0.),
                Complex::i(),
                -Complex::i(),
                Complex::new(0., 0.),
            ],
        );

        System::new(ham)
    }

    #[rstest]
    fn evolve_psi_mismatched_lengths_fail(solution_psi_1: State<'static>, system2: System) {
        let psi = solution_psi_1.with_system(&system2);
        let res_ket = psi.evolve(0.1);

        assert_eq!(
            res_ket,
            Err(String::from("Dimensions of ket(1) and ham(2) don't match!"))
        )
    }

    #[rstest]
    fn evolve_solution_psi_valid_success(solution_psi_2: State<'static>, system2: System) {
        let psi = solution_psi_2.with_system(&system2);
        let res_ket = psi.evolve(0.1).unwrap();

        let correct_ket = State::solution_ket_from_slice(
            &vec![
                Complex64::new(0.1951363713407213, -0.01957894383041598),
                Complex64::new(0.9756818567036065, 0.09789471915207991),
            ],
            Some(&system2),
        );

        assert_eq!(res_ket, correct_ket);
    }

    #[rstest]
    fn get_probabilities_test(solution_psi_2: State) {
        let probabilities = solution_psi_2.get_probabilities();

        assert_eq!(probabilities, vec![0.03846153846153846, 0.9615384615384615])
    }
}
