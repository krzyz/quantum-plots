extern crate nalgebra;
extern crate num;

use nalgebra::DMatrix;
use num::complex::Complex64;

type Op = DMatrix<Complex64>;


pub mod braket;
pub mod state;
pub mod system;
pub mod basis;

pub use braket::*;
pub use state::*;
pub use system::*;
