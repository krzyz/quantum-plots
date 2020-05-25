extern crate nalgebra;
extern crate num;

use nalgebra::{ComplexField, DMatrix, DVector, Dynamic, SymmetricEigen};
use num::complex::Complex64;

use crate::*;
use braket::*;
use state::*;

pub struct System {
    ham: Op,
    eigensystem: Option<SymmetricEigen<Complex64, Dynamic>>,
}

impl System {
    pub fn new(ham: Op) -> System {
        System {
            ham,
            eigensystem: None,
        }
    }

    pub fn eigensystem(&mut self) -> &SymmetricEigen<Complex64, Dynamic> {
        match self.eigensystem {
            Some(ref x) => x,
            None => {
                self.eigensystem = Some(self.ham.clone().symmetric_eigen());
                &self.eigensystem.as_ref().unwrap()
            }
        }
    }

    // TODO: test to_solution
    pub fn to_solution(&mut self, psi: &State) -> State {
        use Braket::*;
        use State::*;

        match psi {
            SolutionBasis(_) => (*psi).clone(),
            ProblemBasis(Ket(vec)) => State::solution_ket(&self.eigensystem().eigenvectors.adjoint() * vec),
            ProblemBasis(Bra(vec)) => State::solution_bra(vec * &self.eigensystem().eigenvectors),
        }
    }

    pub fn to_problem(&mut self, psi: &State) -> State {
        use Braket::*;
        use State::*;

        match psi {
            ProblemBasis(_) => (*psi).clone(),
            SolutionBasis(Ket(vec)) => State::problem_ket(&self.eigensystem().eigenvectors * vec),
            SolutionBasis(Bra(vec)) => State::problem_bra(vec * &self.eigensystem().eigenvectors.adjoint()),
        }
    }

    pub fn evolve_psi(&mut self, psi: &State, t: f64) -> Result<State, String> {
        if psi.len() != self.ham.nrows() {
            return Err(format!(
                "Dimensions of ket({}) and ham({}) don't match!",
                psi.len(),
                self.ham.nrows()
            ));
        }

        let result = match psi {
            State::SolutionBasis(ref sys_vec) => {
                let ham = &self.eigensystem().eigenvalues;
                let diagonal: Vec<Complex64> = ham
                    .iter()
                    .map(|x| (-Complex64::i() * x * t).exp())
                    .collect();

                let dvec = DVector::from_row_slice(&diagonal);
                let matrix = DMatrix::from_diagonal(&dvec);

                match sys_vec {
                    Braket::Ket(ref vec) => State::solution_ket(matrix * vec),
                    Braket::Bra(ref vec) => State::solution_bra(-vec * matrix),
                }
            }
            State::ProblemBasis(_) => {
                let solution_psi = &self.to_solution(psi);
                let result_solution_psi = &self.evolve_psi(solution_psi, t).unwrap();
                self.to_problem(result_solution_psi)
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

    #[fixture]
    pub fn solution_psi_1() -> State {
        State::solution_ket_from_slice(&vec![1.0])
    }

    #[fixture]
    pub fn solution_psi_2() -> State {
        State::solution_ket_from_slice(&vec![0.1, 0.5])
    }

    #[rstest]
    fn evolve_psi_mismatched_lengths_fail(solution_psi_1: State, mut system2: System) {
        let res_ket = system2.evolve_psi(&solution_psi_1, 0.1);

        assert_eq!(
            res_ket,
            Err(String::from("Dimensions of ket(1) and ham(2) don't match!"))
        )
    }

    #[rstest]
    fn evolve_solution_psi_valid_success(solution_psi_2: State, mut system2: System) {
        let res_ket = system2.evolve_psi(&solution_psi_2, 0.1).unwrap();

        let correct_ket = State::solution_ket_from_slice(&vec![
            Complex64::new(0.1951363713407213, -0.01957894383041598),
            Complex64::new(0.9756818567036065, 0.09789471915207991),
        ]);

        assert_eq!(res_ket, correct_ket);
    }
}
