use pyo3::prelude::*;

#[pymodule]
fn siminf_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", "0.1.0")?;
    Ok(())
}