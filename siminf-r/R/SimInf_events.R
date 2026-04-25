SimInf_events <- function(E = NULL, N = NULL, events = NULL, t0 = 1) {
    if (!is.null(events) && !is.data.frame(events))
        stop("'events' must be a data.frame", call. = FALSE)

    if (is.null(E)) {
        E <- Matrix::Matrix(0, nrow = 0, ncol = 0, sparse = TRUE)
    }
    if (is.null(N)) {
        N <- matrix(integer(0), nrow = 0, ncol = 0)
    }

    methods::new("SimInf_events",
                 E = E,
                 N = N,
                 event = if (!is.null(events)) events$event else integer(0),
                 time  = if (!is.null(events)) events$time else integer(0),
                 node  = if (!is.null(events)) events$node else integer(0),
                 dest  = if (!is.null(events)) events$dest else integer(0),
                 n     = if (!is.null(events)) events$n else integer(0),
                 proportion = if (!is.null(events)) events$proportion else numeric(0),
                 select = if (!is.null(events)) events$select else integer(0),
                 shift  = if (!is.null(events)) events$shift else integer(0))
}

setMethod("show", "SimInf_events", function(object) {
    cat("SimInf_events object\n")
    cat(sprintf("  Number of events: %d\n", length(object@event)))
})
