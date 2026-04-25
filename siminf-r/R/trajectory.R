trajectory <- function(model, compartments = NULL, index = NULL,
                       format = c("data.frame", "matrix")) {
    format <- match.arg(format)
    if (!is(model, "SimInf_model"))
        stop("'model' must be a SimInf_model object", call. = FALSE)

    if (is.null(index)) {
        index <- seq_along(model@tspan)
    }

    if (!is.null(compartments)) {
        c_idx <- match(compartments, model@compartment_names)
        if (any(is.na(c_idx)))
            stop("Unknown compartment(s)", call. = FALSE)
    } else {
        c_idx <- seq_len(model@num_compartments)
    }

    t_idx <- index
    tspan <- model@tspan[t_idx]
    num_nodes <- model@num_nodes
    num_compartments <- model@num_compartments

    U <- model@u0
    if (is.null(U))
        stop("Model has not been run yet. Call run() first.", call. = FALSE)

    if (format == "data.frame") {
        rows <- list()
        row_idx <- 1
        for (t in t_idx) {
            for (node in seq_len(num_nodes)) {
                vals <- U[c_idx, node, t]
                rows[[row_idx]] <- data.frame(
                    node = node,
                    time = tspan[which(t_idx == t)],
                    stringsAsFactors = FALSE
                )
                for (j in seq_along(c_idx)) {
                    rows[[row_idx]][[model@compartment_names[c_idx[j]]]] <- vals[j]
                }
                row_idx <- row_idx + 1
            }
        }
        df <- do.call("rbind", rows)
        for (j in seq_along(c_idx)) {
            df[[model@compartment_names[c_idx[j]]]] <- as.integer(
                df[[model@compartment_names[c_idx[j]]]]
            )
        }
        df
    } else {
        U
    }
}

prevalence <- function(model, ...) {
    trajectory(model, ...)
}

setMethod("plot", "SimInf_model", function(x, ...) {
    means <- rowMeans(x@u0, dims = 2)
    tspan <- x@tspan
    matplot(tspan, t(means), type = "l", xlab = "Time", ylab = "Count",
            col = seq_len(nrow(means)), lty = 1)
    legend("right", legend = x@compartment_names, col = seq_len(nrow(means)),
           lty = 1, bty = "n")
})
