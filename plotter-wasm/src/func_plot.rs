use crate::DrawResult;
use crate::Chart;
use plotters::prelude::*;
use std::f32::consts;

pub fn draw(chart: &mut Chart, time: i32) -> DrawResult<impl Fn((i32, i32)) -> Option<(f32, f32)>> {
    let t = time as f32 / 1000.0;
    let root = &chart.root;
    let font: FontDesc = ("sans-serif", 20.0).into();

    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption(format!("y=sin(pi * x - {})", t), font)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_ranged(-1f32..1f32, -1.2f32..1.2f32)?;

    chart.configure_mesh().x_labels(3).y_labels(3).draw()?;

    chart.draw_series(LineSeries::new(
        (-50..=50)
            .map(|x| x as f32 / 50.0)
            .map(|x| (x, (consts::PI * x - t).sin())),
        &RED,
    ))?;

    root.present()?;
    return Ok(chart.into_coord_trans());
}
