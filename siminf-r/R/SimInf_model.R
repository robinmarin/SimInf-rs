setGeneric("u0<-", function(model, value) standardGeneric("u0<-"))
setMethod("u0<-", "SimInf_model",
          function(model, value) {
              model@u0 <- as.matrix(value)
              methods::validObject(model)
              model
          })

setGeneric("shift_matrix", function(model) standardGeneric("shift_matrix"))
setGeneric("shift_matrix<-", function(model, value) standardGeneric("shift_matrix<-"))
setMethod("shift_matrix", "SimInf_model",
          function(model) model@events@N)
setMethod("shift_matrix<-", "SimInf_model",
          function(model, value) {
              model@events@N <- value
              methods::validObject(model)
              model
          })

setGeneric("select_matrix", function(model) standardGeneric("select_matrix"))
setGeneric("select_matrix<-", function(model, value) standardGeneric("select_matrix<-"))
setMethod("select_matrix", "SimInf_model",
          function(model) model@events@E)
setMethod("select_matrix<-", "SimInf_model",
          function(model, value) {
              model@events@E <- value
              methods::validObject(model)
              model
          })

setGeneric("trajectory", function(model, ...)
    standardGeneric("trajectory"))

setGeneric("gdata", function(model) standardGeneric("gdata"))
setMethod("gdata", "SimInf_model", function(model) model@gdata)

setGeneric("gdata<-", function(model, parameter, value)
    standardGeneric("gdata<-"))
setMethod("gdata<-", "SimInf_model",
          function(model, parameter, value) {
              model@gdata[parameter] <- value
              methods::validObject(model)
              model
          })

setGeneric("u0", function(model) standardGeneric("u0"))
setMethod("u0", "SimInf_model", function(model) model@u0)

setGeneric("n_compartments", function(model) standardGeneric("n_compartments"))
setMethod("n_compartments", "SimInf_model",
          function(model) model@num_compartments)

setGeneric("n_nodes", function(model) standardGeneric("n_nodes"))
setMethod("n_nodes", "SimInf_model",
          function(model) model@num_nodes)

setGeneric("events", function(model) standardGeneric("events"))
setMethod("events", "SimInf_model", function(model) model@events)

setGeneric("node_events", function(model, node)
    standardGeneric("node_events"))
setMethod("node_events", "SimInf_model",
          function(model, node) {
              ev <- model@events
              ev@event[ev@node == node]
          })
