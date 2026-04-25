setGeneric("run", function(model, ...)
    standardGeneric("run"))

run <- function(model, solver = c("ssm", "aem")) {
    cat("DEBUG: in run(), solver=", solver, "\n")
    solver <- match.arg(solver)
    cat("DEBUG: after match.arg\n")
    if (!is(model, "SimInf_model"))
        stop("'model' must be a SimInf_model object", call. = FALSE)

    if (is.null(model@ptr))
        stop("Model pointer is NULL", call. = FALSE)

    result_ptr <- .Call(`_siminfr_run`, model@ptr)

    if (is.null(result_ptr))
        stop("Simulation failed", call. = FALSE)

    num_nodes <- .Call(`_siminfr_num_nodes`, result_ptr)
    num_compartments <- .Call(`_siminfr_num_compartments`, result_ptr)
    tspan <- .Call(`_siminfr_tspan`, result_ptr)
    compartment_names <- .Call(`_siminfr_compartment_names`, result_ptr)

    U <- .Call(`_siminfr_trajectory`, result_ptr)
    dim(U) <- c(num_compartments, num_nodes, length(tspan))
    dimnames(U) <- list(compartment_names, NULL, tspan)

    result <- methods::new("SimInf_model",
                           ptr = result_ptr,
                           gdata = model@gdata,
                           tspan = tspan,
                           u0 = U,
                           v0 = model@v0,
                           events = model@events,
                           replicates = model@replicates,
                           num_nodes = num_nodes,
                           num_compartments = num_compartments,
                           compartment_names = compartment_names)
    methods::validObject(result)
    result
}

setMethod("show", "SimInf_model", function(object) {
    cat(sprintf("SimInf_model (%s)\n", class(object)[1]))
    cat(sprintf("  Nodes: %d\n", object@num_nodes))
    cat(sprintf("  Compartments: %s\n",
                paste(object@compartment_names, collapse = ", ")))
    cat(sprintf("  Time points: %d\n", length(object@tspan)))
    cat(sprintf("  Replicates: %d\n", object@replicates))
})
