extern crate nalgebra;
extern crate num;
extern crate simba;

use nalgebra::{convert, ComplexField, DMatrix, DVector, Dynamic, RowDVector, SymmetricEigen};
use num::complex::Complex64;
use simba::scalar::SubsetOf;

type Op = DMatrix<Complex64>;

trait VectorType {}

impl<N> VectorType for DVector<N> where N: ComplexField {}
impl<N> VectorType for RowDVector<N> where N: ComplexField {}

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
        State::ProblemBasis(Braket::Ket(vec))
    }

    pub fn solution_ket(vec: DVector<Complex64>) -> State {
        State::SolutionBasis(Braket::Ket(vec))
    }

    pub fn problem_bra(vec: RowDVector<Complex64>) -> State {
        State::ProblemBasis(Braket::Bra(vec))
    }

    pub fn solution_bra(vec: RowDVector<Complex64>) -> State {
        State::SolutionBasis(Braket::Bra(vec))
    }

    pub fn problem_ket_from_slice<N>(data: &[N]) -> State
    where
        N: SubsetOf<Complex64> + ComplexField,
    {
        State::solution_ket(convert(DVector::from_row_slice(data).normalize()))
    }

    pub fn solution_ket_from_slice<N>(data: &[N]) -> State
    where
        N: SubsetOf<Complex64> + ComplexField,
    {
        State::solution_ket(convert(DVector::from_row_slice(data).normalize()))
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
            ProblemBasis(Ket(vec)) => {
                State::solution_ket(&self.eigensystem().eigenvectors.adjoint() * vec)
            }
            SolutionBasis(Ket(_)) => (*psi).clone(),
            ProblemBasis(Bra(vec)) => State::solution_bra(vec * &self.eigensystem().eigenvectors),
            SolutionBasis(Bra(_)) => (*psi).clone(),
        }
    }

    pub fn to_problem(&mut self, psi: &State) -> State {
        use Braket::*;
        use State::*;

        match psi {
            ProblemBasis(Ket(_)) => (*psi).clone(),
            SolutionBasis(Ket(vec)) => State::problem_ket(&self.eigensystem().eigenvectors * vec),
            ProblemBasis(Bra(_)) => (*psi).clone(),
            SolutionBasis(Bra(vec)) => {
                State::problem_bra(vec * &self.eigensystem().eigenvectors.adjoint())
            }
        }
    }

    pub fn evolve_psi(&mut self, psi: &State, t: f64) -> Result<State, String> {
        println!("{:#?}", &self.eigensystem());
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

        let correct_ket: State = State::solution_ket_from_slice(&vec![
            Complex64::new(0.1951363713407213, -0.01957894383041598),
            Complex64::new(0.9756818567036065, 0.09789471915207991),
        ]);

        assert_eq!(res_ket, correct_ket);
    }
}
