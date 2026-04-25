//! siminf-r: Rcpp FFI bindings to siminf-core
//!
//! All functions use #[no_mangle] + extern "C" and receive/return SEXP
//! to be callable directly from R via .Call().

use libc::{c_char, c_double, c_int, c_void};

type SEXP = *mut usize;

const INTSXP: c_int = 13;
const REALSXP: c_int = 14;
const STRSXP: c_int = 16;
const EXTPTRSXP: c_int = 20;

extern "C" {
    fn Rf_mkString(buf: *const c_char) -> SEXP;
    fn Rf_allocVector(t: c_int, len: usize) -> SEXP;
    fn Rf_allocMatrix(t: c_int, nrow: c_int, ncol: c_int) -> SEXP;
    fn Rf_ScalarInteger(v: c_int) -> SEXP;
    fn Rf_xlength(x: SEXP) -> usize;
    fn Rf_NilValue() -> SEXP;
    fn R_MakeExternalPtr(ptr: *mut c_void, tag: SEXP, prot: SEXP) -> SEXP;
    fn R_ExternalPtrAddr(x: SEXP) -> *mut c_void;
    fn R_PreserveObject(x: SEXP);
    fn R_ReleaseObject(x: SEXP);
    fn REAL(x: SEXP) -> *mut c_double;
    fn INTEGER(x: SEXP) -> *mut c_int;
    fn SET_STRING_ELT(x: SEXP, i: usize, v: SEXP);
    fn STRING_ELT(x: SEXP, i: usize) -> SEXP;
    fn R_CHAR(x: SEXP) -> *const c_char;
    fn Rf_mkChar(s: *const c_char) -> SEXP;
    fn Rf_length(x: SEXP) -> usize;
}

fn null_sexp() -> SEXP {
    unsafe { Rf_NilValue() }
}

fn alloc_matrix(t: c_int, nrow: c_int, ncol: c_int) -> SEXP {
    unsafe { Rf_allocMatrix(t, nrow, ncol) }
}

fn mk_string(s: &str) -> SEXP {
    let c_str = std::ffi::CString::new(s).expect("CString::new failed");
    unsafe { Rf_mkString(c_str.as_ptr()) }
}

fn sexp_to_real_slice(sexp: SEXP) -> Option<&'static [c_double]> {
    if sexp.is_null() { return None; }
    unsafe {
        let len = Rf_xlength(sexp);
        if len == 0 { return None; }
        let ptr = REAL(sexp);
        if ptr.is_null() { return None; }
        Some(std::slice::from_raw_parts(ptr, len))
    }
}

fn sexp_to_int_slice(sexp: SEXP) -> Option<&'static [c_int]> {
    if sexp.is_null() { return None; }
    unsafe {
        let len = Rf_xlength(sexp);
        if len == 0 { return None; }
        let ptr = INTEGER(sexp);
        if ptr.is_null() { return None; }
        Some(std::slice::from_raw_parts(ptr, len))
    }
}

fn sexp_to_usize(sexp: SEXP) -> Option<usize> {
    if sexp.is_null() { return None; }
    unsafe {
        let len = Rf_xlength(sexp);
        if len == 0 { return None; }
        let ptr = INTEGER(sexp);
        Some(*ptr as usize)
    }
}

fn real_to_vec(sexp: SEXP) -> Vec<c_double> {
    sexp_to_real_slice(sexp).map(|s| s.to_vec()).unwrap_or_default()
}

fn int_to_vec(sexp: SEXP) -> Vec<c_int> {
    sexp_to_int_slice(sexp).map(|s| s.to_vec()).unwrap_or_default()
}

fn allocate_model_ptr(model: siminf_core::Model) -> SEXP {
    let boxed = Box::new(model);
    let ptr = Box::into_raw(boxed) as *mut c_void;
    let xp = unsafe { R_MakeExternalPtr(ptr, Rf_NilValue(), Rf_NilValue()) };
    unsafe { R_PreserveObject(xp) };
    xp
}

fn allocate_result_ptr(result: siminf_core::TrajectoryResult) -> SEXP {
    let boxed = Box::new(result);
    let ptr = Box::into_raw(boxed) as *mut c_void;
    let xp = unsafe { R_MakeExternalPtr(ptr, Rf_NilValue(), Rf_NilValue()) };
    unsafe { R_PreserveObject(xp) };
    xp
}

fn model_from_sexp(sexp: SEXP) -> Option<&'static mut siminf_core::Model> {
    if sexp.is_null() { return None; }
    unsafe {
        let ptr = R_ExternalPtrAddr(sexp) as *mut siminf_core::Model;
        if ptr.is_null() { None } else { Some(&mut *ptr) }
    }
}

fn result_from_sexp(sexp: SEXP) -> Option<&'static mut siminf_core::TrajectoryResult> {
    if sexp.is_null() { return None; }
    unsafe {
        let ptr = R_ExternalPtrAddr(sexp) as *mut siminf_core::TrajectoryResult;
        if ptr.is_null() { None } else { Some(&mut *ptr) }
    }
}

fn sir_transitions() -> Vec<siminf_core::Transition> {
    use std::sync::Arc;

    let s_to_i = siminf_core::PropensityFn::Custom {
        name: "S_to_I".to_string(),
        eval: Arc::new(|u: &[i32], _: &[f64], _: &[f64], gdata: &[f64], _t: f64| {
            let s = u[0] as f64;
            let i = u[1] as f64;
            let n = (u[0] + u[1] + u[2]) as f64;
            let beta = gdata[0];
            if n > 0.0 && s > 0.0 && i > 0.0 { beta * s * i / n } else { 0.0 }
        }),
    };
    let i_to_r = siminf_core::PropensityFn::Custom {
        name: "I_to_R".to_string(),
        eval: Arc::new(|u: &[i32], _: &[f64], _: &[f64], gdata: &[f64], _t: f64| {
            let i = u[1] as f64;
            let gamma = gdata[1];
            if i > 0.0 { gamma * i } else { 0.0 }
        }),
    };
    vec![
        siminf_core::Transition {
            name: "S_to_I".to_string(),
            from: vec![siminf_core::CompartmentId(0)],
            to: vec![siminf_core::CompartmentId(1)],
            propensity_fn: s_to_i,
        },
        siminf_core::Transition {
            name: "I_to_R".to_string(),
            from: vec![siminf_core::CompartmentId(1)],
            to: vec![siminf_core::CompartmentId(2)],
            propensity_fn: i_to_r,
        },
    ]
}

#[no_mangle]
pub extern "C" fn siminfr_create_sir(
    gdata_sexp: SEXP,
    u0_sexp: SEXP,
    tspan_sexp: SEXP,
    seed_sexp: SEXP,
) -> SEXP {
    let gdata = real_to_vec(gdata_sexp);
    let u0 = int_to_vec(u0_sexp);
    let tspan = real_to_vec(tspan_sexp);
    let seed = sexp_to_usize(seed_sexp).unwrap_or(42) as u64;

    let num_compartments = 3;
    let num_nodes = u0.len() / num_compartments;

    let model = match siminf_core::Model::builder()
        .compartments(&["S", "I", "R"], &[0; 3])
        .transition(sir_transitions()[0].clone())
        .transition(sir_transitions()[1].clone())
        .global_data("beta", gdata[0])
        .global_data("gamma", gdata[1])
        .u0(u0)
        .num_nodes(num_nodes)
        .tspan_points(tspan)
        .seed(seed)
        .build()
    {
        Ok(m) => m,
        Err(_) => return null_sexp(),
    };

    allocate_model_ptr(model)
}

#[no_mangle]
pub extern "C" fn siminfr_run(model_ptr_sexp: SEXP) -> SEXP {
    let model = match model_from_sexp(model_ptr_sexp) {
        Some(m) => m,
        None => return null_sexp(),
    };

    let model_clone = model.clone();
    let mut solver = siminf_core::Solver::new(model_clone);

    match solver.run() {
        Ok(result) => allocate_result_ptr(result),
        Err(_) => null_sexp(),
    }
}

#[no_mangle]
pub extern "C" fn siminfr_trajectory(result_ptr_sexp: SEXP) -> SEXP {
    let result = match result_from_sexp(result_ptr_sexp) {
        Some(p) => p,
        None => return null_sexp(),
    };

    let num_nodes = result.num_nodes;
    let num_compartments = result.num_compartments;
    let num_timepoints = result.tspan.len();
    let nrow = (num_nodes * num_compartments) as c_int;
    let ncol = num_timepoints as c_int;

    let sexp = unsafe { Rf_allocMatrix(INTSXP, nrow, ncol) };

    unsafe {
        let dst = INTEGER(sexp);
        dst.copy_from_nonoverlapping(result.u.as_ptr() as *const c_int, result.u.len());
    }

    sexp
}

#[no_mangle]
pub extern "C" fn siminfr_result_mean(result_ptr_sexp: SEXP) -> SEXP {
    let result = match result_from_sexp(result_ptr_sexp) {
        Some(p) => p,
        None => return null_sexp(),
    };

    let means = result.mean_compartments();
    let num_compartments = result.num_compartments;
    let num_timepoints = result.tspan.len();
    let nrow = num_compartments as c_int;
    let ncol = num_timepoints as c_int;

    let sexp = unsafe { Rf_allocMatrix(REALSXP, nrow, ncol) };

    unsafe {
        let dst = REAL(sexp);
        for (t_idx, row) in means.iter().enumerate() {
            for (c_idx, &val) in row.iter().enumerate() {
                dst.add(t_idx * num_compartments + c_idx).write(val);
            }
        }
    }

    sexp
}

#[no_mangle]
pub extern "C" fn siminfr_compartment_names(result_ptr_sexp: SEXP) -> SEXP {
    let result = match result_from_sexp(result_ptr_sexp) {
        Some(p) => p,
        None => return null_sexp(),
    };

    let n = result.compartment_names.len();
    let sexp = unsafe { Rf_allocVector(STRSXP, n) };

    for (i, name) in result.compartment_names.iter().enumerate() {
        let c_str = std::ffi::CString::new(name.as_str()).expect("CString::new failed");
        unsafe {
            SET_STRING_ELT(sexp, i, Rf_mkChar(c_str.as_ptr()));
        }
    }

    sexp
}

#[no_mangle]
pub extern "C" fn siminfr_tspan(result_ptr_sexp: SEXP) -> SEXP {
    let result = match result_from_sexp(result_ptr_sexp) {
        Some(p) => p,
        None => return null_sexp(),
    };

    let n = result.tspan.len();
    let sexp = unsafe { Rf_allocVector(REALSXP, n) };

    unsafe {
        let dst = REAL(sexp);
        dst.copy_from_nonoverlapping(result.tspan.as_ptr(), result.tspan.len());
    }

    sexp
}

#[no_mangle]
pub extern "C" fn siminfr_num_nodes(result_ptr_sexp: SEXP) -> SEXP {
    let result = match result_from_sexp(result_ptr_sexp) {
        Some(p) => p,
        None => return null_sexp(),
    };
    unsafe { Rf_ScalarInteger(result.num_nodes as c_int) }
}

#[no_mangle]
pub extern "C" fn siminfr_num_compartments(result_ptr_sexp: SEXP) -> SEXP {
    let result = match result_from_sexp(result_ptr_sexp) {
        Some(p) => p,
        None => return null_sexp(),
    };
    unsafe { Rf_ScalarInteger(result.num_compartments as c_int) }
}

#[no_mangle]
pub extern "C" fn siminfr_delete_model(model_ptr_sexp: SEXP) -> SEXP {
    if !model_ptr_sexp.is_null() {
        unsafe {
            let p = R_ExternalPtrAddr(model_ptr_sexp);
            if !p.is_null() {
                drop(Box::from_raw(p as *mut siminf_core::Model));
                R_ReleaseObject(model_ptr_sexp);
            }
        }
    }
    null_sexp()
}

#[no_mangle]
pub extern "C" fn siminfr_delete_result(result_ptr_sexp: SEXP) -> SEXP {
    if !result_ptr_sexp.is_null() {
        unsafe {
            let p = R_ExternalPtrAddr(result_ptr_sexp);
            if !p.is_null() {
                drop(Box::from_raw(p as *mut siminf_core::TrajectoryResult));
                R_ReleaseObject(result_ptr_sexp);
            }
        }
    }
    null_sexp()
}

#[no_mangle]
pub extern "C" fn siminfr_mparse(
    _transitions_sexp: SEXP,
    compartments_sexp: SEXP,
    gdata_sexp: SEXP,
    _ldata_sexp: SEXP,
    u0_sexp: SEXP,
    _v0_sexp: SEXP,
    tspan_sexp: SEXP,
    seed_sexp: SEXP,
    nd_sexp: SEXP,
) -> SEXP {
    let compartments = if compartments_sexp.is_null() {
        return null_sexp();
    } else {
        let len = unsafe { Rf_xlength(compartments_sexp) } as usize;
        (0..len)
            .map(|i| {
                let elt = unsafe { STRING_ELT(compartments_sexp, i) };
                let cptr = unsafe { R_CHAR(elt) };
                let cstr = unsafe { std::ffi::CStr::from_ptr(cptr) };
                cstr.to_string_lossy().into_owned()
            })
            .collect::<Vec<_>>()
    };

    let num_compartments = compartments.len();
    let gdata = real_to_vec(gdata_sexp);
    let u0 = int_to_vec(u0_sexp);
    let tspan = real_to_vec(tspan_sexp);
    let seed = sexp_to_usize(seed_sexp).unwrap_or(42) as u64;
    let nd = sexp_to_usize(nd_sexp).unwrap_or(0);

    if u0.is_empty() || tspan.is_empty() {
        return null_sexp();
    }

    let num_nodes = u0.len() / num_compartments;
    let initial_counts: Vec<i32> = vec![0; num_compartments];

    let transitions = sir_transitions();
    let mut gd = siminf_core::GlobalData::default();
    for (i, &val) in gdata.iter().enumerate() {
        gd.values.insert(format!("p{}", i), val);
    }

    let model = match siminf_core::Model::builder()
        .compartments(
            &compartments.iter().map(|s| s.as_str()).collect::<Vec<_>>(),
            &initial_counts,
        )
        .transition(transitions[0].clone())
        .transition(transitions[1].clone())
        .gdata_struct(gd)
        .u0(u0)
        .num_nodes(num_nodes)
        .tspan_points(tspan)
        .seed(seed)
        .nd(nd)
        .build()
    {
        Ok(m) => m,
        Err(_) => return null_sexp(),
    };

    allocate_model_ptr(model)
}
