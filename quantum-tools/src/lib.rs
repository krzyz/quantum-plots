extern crate nalgebra as na;
extern crate num;

use na::{DMatrix, DVector, Dynamic, SymmetricEigen};
use num::complex::Complex64;

//type Bra = DVector<Complex64>;
type Op = DMatrix<Complex64>;

pub enum SystemVector {
    Problem
}

pub enum State {
    Ket(SystemVector)
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
        match self.eigensystem.take() {
            Some(x) => {
                self.eigensystem = Some(x);
            },
            None => {
                self.eigensystem = Some(self.ham.clone().symmetric_eigen());
            }
        }

        &self.eigensystem.as_ref().unwrap()
    }

    pub fn evolve_psi(&mut self, psi: Ket, t: f64) -> Result<Ket, String> {
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
        let result = matrix * psi;
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
    pub fn psi_1() -> Ket {
        na::convert(DVector::from_row_slice(&vec![1]))
    }

    #[fixture]
    pub fn psi_2() -> Ket {
        na::convert(DVector::from_row_slice(&vec![0.1, 0.5]).normalize())
    }

    #[rstest]
    fn evolve_psi_mismatched_lengths_fail(psi_1: Ket, mut system2: System) {
        let res_ket = system2.evolve_psi(psi_1, 0.1);

        assert_eq!(
            res_ket,
            Err(String::from("Dimensions of ket(1) and ham(2) don't match!"))
        )
    }

    #[rstest]
    fn evolve_psi_valid_success(psi_2: Ket, mut system2: System) {
        let res_ket = system2.evolve_psi(psi_2, 0.1).unwrap();

        let correct_ket: Ket = na::convert(DVector::from_row_slice(&vec![
            Complex64::new(0.1951363713407213, 0.01957894383041598),
            Complex64::new(0.9756818567036065, -0.09789471915207991),
        ]));

        assert_eq!(res_ket, correct_ket);
    }
}
