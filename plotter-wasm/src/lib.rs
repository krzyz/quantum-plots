use func_plot::{PlotData, Problem};
use plotters::coord::Shift;
use plotters::prelude::*;
use wasm_bindgen::prelude::*;

mod func_plot;
mod utils;

#[global_allocator]
static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;

pub type DrawResult<T> = Result<T, Box<dyn std::error::Error>>;

#[wasm_bindgen]
pub struct Chart {
    root: DrawingArea<CanvasBackend, Shift>,
    convert: Box<dyn Fn((i32, i32)) -> Option<(f64, f64)>>,
    plot_data: PlotData,
    problem: Problem,
}

#[wasm_bindgen]
pub struct Point {
    pub x: f64,
    pub y: f64,
}

#[wasm_bindgen]
impl Chart {
    /// Draw provided power function on the canvas element using it's id.
    /// Return `Chart` struct suitable for coordinate conversion.
    pub fn power(&mut self, time: f64) -> Result<(), JsValue> {
        self.plot_data.add_time(&mut self.problem, time);
        let map_coord = func_plot::draw(self).map_err(|err| err.to_string())?;
        self.convert = Box::new(move |coord| map_coord(coord).map(|(x, y)| (x.into(), y.into())));

        Ok(())
    }
    pub fn new(canvas_id: &str) -> Chart {
        utils::set_panic_hook();
        let backend = CanvasBackend::new(canvas_id).expect("cannot find canvas");
        let root = backend.into_drawing_area();
        let problem = Problem::new();
        let plot_data = PlotData::new(&problem);

        Chart {
            root,
            convert: Box::new(move |(x, y)| Some((x.into(), y.into()))),
            plot_data,
            problem,
        }
    }

    /// This function can be used to convert screen coordinates to
    /// chart coordinates.
    pub fn coord(&self, x: i32, y: i32) -> Option<Point> {
        (self.convert)((x, y)).map(|(x, y)| Point { x, y })
    }
}
