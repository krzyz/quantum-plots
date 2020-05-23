extern crate nalgebra;
extern crate num;
extern crate simba;

use nalgebra::{convert, ComplexField, DMatrix, DVector, Dynamic, SymmetricEigen};
use num::complex::Complex64;
use simba::scalar::SubsetOf;

type Op = DMatrix<Complex64>;

#[derive(Debug)]
pub enum SystemVector {
    ProblemBasis(DVector<Complex64>),
    SolutionBasis(DVector<Complex64>),
}

impl SystemVector {
    pub fn len(&self) -> usize {
        self.get_vector().len()
    }

    pub fn get_vector(&self) -> &DVector<Complex64> {
        use SystemVector::*;
        match self {
            &ProblemBasis(ref vec) => &vec,
            &SolutionBasis(ref vec) => &vec,
        }
    }
}

impl PartialEq for SystemVector {
    fn eq(&self, other: &Self) -> bool {
        match std::mem::discriminant(self) == std::mem::discriminant(other) {
            true => self.get_vector() == other.get_vector(),
            false => false
        }
    }
}

#[derive(Debug)]
pub enum State {
    Bra(SystemVector),
    Ket(SystemVector),
}

impl State {
    pub fn len(&self) -> usize {
        return self.get_vector().len()
    }

    pub fn get_system_vector(&self) -> &SystemVector {
        use State::*;
        match self {
            &Bra(ref vec) => &vec,
            &Ket(ref vec) => &vec,
        }
    }

    pub fn get_vector(&self) -> &DVector<Complex64> {
        self.get_system_vector().get_vector()
    }

    pub fn problem_ket(vec: DVector<Complex64>) -> State {
        State::Ket(SystemVector::ProblemBasis(vec))
    }

    pub fn solution_ket(vec: DVector<Complex64>) -> State {
        State::Ket(SystemVector::SolutionBasis(vec))
    }

    pub fn problem_ket_from_slice<N> (data: &[N]) -> State
    where N: SubsetOf<Complex64> + ComplexField {
        State::problem_ket(convert(DVector::from_row_slice(data).normalize()))
    }

    pub fn solution_ket_from_slice<N> (data: &[N]) -> State
    where N: SubsetOf<Complex64> + ComplexField {
        State::solution_ket(convert(DVector::from_row_slice(data).normalize()))
    }
}

impl PartialEq for State {
    fn eq(&self, other: &Self) -> bool {
        match std::mem::discriminant(self) == std::mem::discriminant(other) {
            true => self.get_vector() == other.get_vector(),
            false => false
        }
    }
}

pub struct System
{
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

    pub fn evolve_psi(&mut self, psi: State, t: f64) -> Result<State, String> {
        let ham = &self.eigensystem().eigenvalues;
        if psi.len() != ham.len() {
            return Err(format!(
                "Dimensions of ket({}) and ham({}) don't match!",
                psi.len(),
                ham.len()
            ));
        }
        let diagonal: Vec<Complex64> = ham.iter().map(|x| (Complex64::i() * x * t).exp()).collect();
        let dvec = DVector::from_row_slice(&diagonal);
        let matrix = DMatrix::from_diagonal(&dvec);
        let result = match psi {
            State::Ket(sys_vec) => {
                match sys_vec {
                    SystemVector::SolutionBasis(vec) => {
                        State::solution_ket(matrix * vec)
                    },
                    _ => panic!(),
                }
            },
            _ => panic!(),
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
        let res_ket = system2.evolve_psi(solution_psi_1, 0.1);

        assert_eq!(
            res_ket,
            Err(String::from("Dimensions of ket(1) and ham(2) don't match!"))
        )
    }

    #[rstest]
    fn evolve_solution_psi_valid_success(solution_psi_2: State, mut system2: System) {
        let res_ket = system2.evolve_psi(solution_psi_2, 0.1).unwrap();

        let correct_ket: State = State::solution_ket_from_slice(&vec![
            Complex64::new(0.1951363713407213, 0.01957894383041598),
            Complex64::new(0.9756818567036065, -0.09789471915207991),
        ]);

        assert_eq!(res_ket, correct_ket);
    }
}
