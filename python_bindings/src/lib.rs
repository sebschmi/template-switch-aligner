use std::io;

use lib_tsalign::{
    a_star_aligner::{
        alignment_result::AlignmentResult,
        configurable_a_star_align::{Config, a_star_align},
        template_switch_distance::AlignmentType,
    },
    costs::U64Cost,
};
use lib_tsshow::plain_text::show_template_switches;
use pyo3::{exceptions::PyRuntimeError, prelude::*, types::PyDict};
use pythonize::{depythonize, pythonize};

#[pyclass]
struct TSPairwiseAlignment {
    result: AlignmentResult<AlignmentType, U64Cost>,
}

#[pymethods]
impl TSPairwiseAlignment {
    fn viz_template_switches(&self) -> PyResult<()> {
        show_template_switches(io::stdout(), &self.result, &None)
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(())
    }

    fn stats<'a>(&'a self, py: Python<'a>) -> PyResult<Bound<'a, PyAny>> {
        Ok(pythonize(py, self.result.statistics())?)
    }

    fn cigar(&self) -> Option<String> {
        match &self.result {
            AlignmentResult::WithTarget { alignment, .. } => Some(alignment.cigar()),
            AlignmentResult::WithoutTarget { .. } => None,
        }
    }

    fn alignments<'a>(&'a self, py: Python<'a>) -> PyResult<Option<Bound<'a, PyAny>>> {
        match &self.result {
            AlignmentResult::WithTarget { alignment, .. } => {
                let mut container = Vec::new();
                alignment.iter_compact().for_each(|e| container.push(e));
                Ok(Some(pythonize(py, &container)?))
            }
            AlignmentResult::WithoutTarget { .. } => Ok(None),
        }
    }
}

fn py_to_str(o: Bound<'_, PyAny>) -> PyResult<Vec<u8>> {
    let str = o.str()?.to_str()?.as_bytes().to_vec();
    Ok(str)
}

/// Creates a config object by amending the default by values present in the python dictionary
fn create_config(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Config> {
    let Some(kwargs) = kwargs else {
        return Ok(Config::default());
    };

    let config = depythonize::<Config>(kwargs)?;

    Ok(config)
}

/// Align two sequences, accounting for template switches
///
/// The function takes a reference and a query string, and performs a global alignment on both. The output alignment may contain (short-range) template switches.
/// Optionally, settings can be specified. See the (configuration struct)[lib_tsalign::a_star_aligner::configurable_a_star_align::Config] for the available keys and values.
#[pyfunction]
#[pyo3(signature = (reference, query, **kwargs))]
fn align(
    reference: Bound<'_, PyAny>, // Accepting PyAny instead of PyString to allow using e.g. `Bio.Seq` types and alike. String representation will be used.
    query: Bound<'_, PyAny>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<Option<TSPairwiseAlignment>> {
    let reference = py_to_str(reference)?;
    let query = py_to_str(query)?;
    let config = create_config(kwargs)?;

    let r = a_star_align(&reference, &query, &config);

    match r {
        result @ AlignmentResult::WithTarget { .. } => {
            let ts_alignment = TSPairwiseAlignment { result };
            Ok(Some(ts_alignment))
        }
        AlignmentResult::WithoutTarget { .. } => Ok(None),
    }
}

/// Bindings for the `lib_tsalign` library.
#[pymodule]
fn tsalign(m: &Bound<'_, PyModule>) -> PyResult<()> {
    pyo3_log::init();
    m.add_class::<TSPairwiseAlignment>()?;
    m.add_function(wrap_pyfunction!(align, m)?)?;
    Ok(())
}
