compartments_SIR <- function() c("S", "I", "R")

select_matrix_SIR <- function() {
    matrix(c(1, 0, 0, 0,
             1, 0, 0, 0,
             1, 1, 1, 1),
           nrow = 3, ncol = 4,
           dimnames = list(compartments_SIR(), seq_len(4)))
}

SIR <- function(u0, tspan, events = NULL, beta = NULL, gamma = NULL,
                seed = 42L) {
    if (is.null(beta) || is.null(gamma))
        stop("'beta' and 'gamma' must be specified", call. = FALSE)
    if (is.data.frame(u0)) {
        num_nodes <- nrow(u0)
        u0_matrix <- as.matrix(u0)
    } else if (is.matrix(u0)) {
        num_nodes <- nrow(u0)
        u0_matrix <- u0
    } else if (is.vector(u0)) {
        num_nodes <- 1L
        u0_matrix <- matrix(u0, nrow = length(u0), ncol = 1)
    } else {
        stop("'u0' must be a data.frame, matrix, or vector", call. = FALSE)
    }

    gdata <- c(beta = beta, gamma = gamma)
    tspan <- as.numeric(tspan)
    seed <- as.integer(seed)[1]

    if (is.null(events)) {
        events <- methods::new("SimInf_events")
    }

    model_ptr <- .Call(`_siminfr_create_sir`,
          as.numeric(gdata),
          as.integer(t(u0_matrix)),
          as.numeric(tspan),
          seed)

    if (is.null(model_ptr))
        stop("Failed to create model", call. = FALSE)

    model <- methods::new("SIR",
                          ptr = model_ptr,
                          gdata = gdata,
                          tspan = tspan,
                          u0 = u0_matrix,
                          v0 = matrix(nrow = 0, ncol = 0),
                          events = events,
                          replicates = 1L,
                          num_nodes = as.integer(num_nodes),
                          num_compartments = 3L,
                          compartment_names = compartments_SIR())
    methods::validObject(model)
    model
}
