setClass(
    "SimInf_model",
    slots = c(
        ptr          = "externalptr",
        gdata        = "numeric",
        tspan        = "numeric",
        u0           = "array",
        v0           = "matrix",
        events       = "SimInf_events",
        replicates   = "integer",
        num_nodes    = "integer",
        num_compartments = "integer",
        compartment_names = "character"
    )
)

setClass(
    "SimInf_events",
    slots = c(
        E          = "dgCMatrix",
        N          = "matrix",
        event      = "integer",
        time       = "integer",
        node       = "integer",
        dest       = "integer",
        n          = "integer",
        proportion = "numeric",
        select     = "integer",
        shift      = "integer"
    )
)

setClass("SIR", contains = c("SimInf_model"))
setClass("SEIR", contains = c("SimInf_model"))
setClass("SIS", contains = c("SimInf_model"))
