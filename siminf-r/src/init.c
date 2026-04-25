#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

static SEXP nil_value_cache = NULL;

SEXP Rf_NilValue(void) {
    if (nil_value_cache == NULL) {
        nil_value_cache = R_NilValue;
    }
    return nil_value_cache;
}

extern SEXP siminfr_create_sir(SEXP gdata, SEXP u0, SEXP tspan, SEXP seed);
extern SEXP siminfr_run(SEXP model_ptr);
extern SEXP siminfr_trajectory(SEXP result_ptr);
extern SEXP siminfr_result_mean(SEXP result_ptr);
extern SEXP siminfr_compartment_names(SEXP result_ptr);
extern SEXP siminfr_tspan(SEXP result_ptr);
extern SEXP siminfr_num_nodes(SEXP result_ptr);
extern SEXP siminfr_num_compartments(SEXP result_ptr);
extern SEXP siminfr_delete_model(SEXP model_ptr);
extern SEXP siminfr_delete_result(SEXP result_ptr);
extern SEXP siminfr_mparse(SEXP transitions, SEXP compartments, SEXP gdata,
                          SEXP ldata, SEXP u0, SEXP v0, SEXP tspan,
                          SEXP seed, SEXP nd);

static const R_CallMethodDef CallEntries[] = {
    {"_siminfr_create_sir",        (DL_FUNC) &siminfr_create_sir,        4},
    {"_siminfr_run",               (DL_FUNC) &siminfr_run,               1},
    {"_siminfr_trajectory",        (DL_FUNC) &siminfr_trajectory,        1},
    {"_siminfr_result_mean",       (DL_FUNC) &siminfr_result_mean,       1},
    {"_siminfr_compartment_names", (DL_FUNC) &siminfr_compartment_names, 1},
    {"_siminfr_tspan",            (DL_FUNC) &siminfr_tspan,             1},
    {"_siminfr_num_nodes",         (DL_FUNC) &siminfr_num_nodes,         1},
    {"_siminfr_num_compartments",  (DL_FUNC) &siminfr_num_compartments,  1},
    {"_siminfr_delete_model",      (DL_FUNC) &siminfr_delete_model,      1},
    {"_siminfr_delete_result",    (DL_FUNC) &siminfr_delete_result,     1},
    {"_siminfr_mparse",           (DL_FUNC) &siminfr_mparse,            9},
    {NULL, NULL, 0}
};

void R_init_siminfr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
}
