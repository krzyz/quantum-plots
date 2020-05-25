extern crate nalgebra;
extern crate num;

use crate::Chart;
use crate::DrawResult;
use plotters::prelude::*;
use quantum_tools as qt;
use nalgebra::DMatrix;
use num::complex::Complex;
use std::collections::LinkedList;

extern crate web_sys;

// A macro to provide `println!(..)`-style syntax for `console.log` logging.
#[macro_use]
macro_rules! log {
    ( $( $t:tt )* ) => {
        web_sys::console::log_1(&format!( $( $t )* ).into());
    }
}

pub struct Problem {
    psi_0: qt::State,
    system: qt::System,
    size: usize,
}

impl Problem {
    pub fn new() -> Problem {
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

        let psi_0_problem = qt::State::problem_ket_from_slice(&vec![1.0, 2.0]);
        let size = psi_0_problem.len();
        let mut system = qt::System::new(ham);
        let psi_0 = system.to_solution(&psi_0_problem);

        Problem {
            psi_0,
            system,
            size,
        }
    }
}

pub struct PlotData {
    t: LinkedList<f32>,
    data: Vec<LinkedList<f32>>,
}

impl PlotData {
    pub fn new(problem: &Problem) -> PlotData {
        let t = LinkedList::new();
        let data = vec![LinkedList::new(); problem.size];

        PlotData { t, data }
    }

    pub fn add_time(&mut self, problem: &mut Problem, time: f64) {
        self.t.push_back(time as f32);
        if self.t.len() > 5000 {
            self.t.pop_front();
        }

        let psi_t = problem.system.evolve_psi(&problem.psi_0, time).unwrap();
        let psi_t_problem = problem.system.to_problem(&psi_t);
        let probabilities = psi_t_problem.get_probabilities();

        for (i, prob) in probabilities.iter().enumerate() {
            self.data[i].push_back(*prob as f32);
            if self.data[i].len() > 5000 {
                self.data[i].pop_front();
            }
        }
    }
}

pub fn draw(
    chart_struct: &mut Chart,
) -> DrawResult<impl Fn((i32, i32)) -> Option<(f32, f32)>> {
    let root = &chart_struct.root;
    let font: FontDesc = ("sans-serif", 20.0).into();

    root.fill(&WHITE)?;
    let xmin = 0f32.max(*chart_struct.plot_data.t.front().unwrap_or(&0.0));
    let xmax = 12f32.max(*chart_struct.plot_data.t.back().unwrap_or(&0.0));

    let mut chart = ChartBuilder::on(&root)
        .caption("quantum", font)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_ranged(xmin..xmax, 0f32..1f32)?;

    chart.configure_mesh().x_labels(3).y_labels(3).draw()?;

    for state_probs in chart_struct.plot_data.data.clone().into_iter() {
        chart.draw_series(LineSeries::new(
            chart_struct.plot_data.t.clone().into_iter().zip(state_probs.into_iter()),
            &RED,
        ))?;
    }

    root.present()?;
    return Ok(chart.into_coord_trans());
}
