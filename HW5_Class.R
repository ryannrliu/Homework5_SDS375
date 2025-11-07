## HW5 Class/Methods

if(!"methods" %in% loadedNamespaces()) {
  require(methods)
}

setClass(
    Class = "sparse_numeric",
    slots = c(
        value = "numeric",
        pos = "integer",
        length = "integer"
    )
)


setValidity("sparse_numeric", function(object) {
  # length must be a single positive integer
  if(length(object@length) != 1L || object@length < 0L)
    return("slot 'length' must be a single non-negative integer")
  
  # pos must be integer vector (can be length 0)
  if(!is.integer(object@pos))
    return("slot 'pos' must be integer")
  
  # value must be numeric
  if(!is.numeric(object@value))
    return("slot 'value' must be numeric")
  
  # lengths of value and pos must match
  if(length(object@value) != length(object@pos))
    return("slots 'value' and 'pos' must have the same length")
  # if pos is non-empty, entries must be within 1..length and unique and increasing
  if(length(object@pos) > 0L) {
    if(any(object@pos < 1L) || any(object@pos > object@length))
      return("entries of 'pos' must be within 1 and object@length")
    if(any(duplicated(object@pos)))
      return("entries of 'pos' must be unique")
  }
  
  TRUE
})

## numeric -> sparse_numeric
setAs("numeric", "sparse_numeric", function(from) {
  if(!is.numeric(from)) stop("can only coerce numeric vectors")
  nz <- which(from != 0)
  new("sparse_numeric",
      value = if(length(nz) > 0) as.numeric(from[nz]) else numeric(0),
      pos = if(length(nz) > 0) as.integer(nz) else integer(0),
      length = as.integer(length(from)))
})

## sparse_numeric -> numeric
setAs("sparse_numeric", "numeric", function(from) {
  out <- numeric(from@length)
  if(length(from@pos) > 0) {
    out[from@pos] <- from@value
  }
  out
})

setGeneric("sparse_add", function(x, y, ...) standardGeneric("sparse_add"))
setGeneric("sparse_sub", function(x, y, ...) standardGeneric("sparse_sub"))
setGeneric("sparse_mult", function(x, y, ...) standardGeneric("sparse_mult"))
setGeneric("sparse_crossprod", function(x, y, ...) standardGeneric("sparse_crossprod"))

## Extra generic: sparse_norm
setGeneric("sparse_norm", function(x, p = 2, ...) standardGeneric("sparse_norm"))

.check_length_equal <- function(xobj, yobj) {
  if(xobj@length != yobj@length) stop("vectors must have the same length")
}

setMethod("sparse_add", signature(x = "sparse_numeric", y = "sparse_numeric"),
          function(x, y, ...) {
            .check_length_equal(x, y)
            if(length(x@pos) == 0L && length(y@pos) == 0L) {
              return(new("sparse_numeric", value = numeric(0), pos = integer(0), length = x@length))
            }
            allpos <- sort(unique(c(x@pos, y@pos)))
            # use match to find values
            xv <- numeric(length(allpos)); yv <- numeric(length(allpos))
            if(length(x@pos) > 0) {
              m <- match(allpos, x@pos)
              xv[!is.na(m)] <- x@value[m[!is.na(m)]]
            }
            if(length(y@pos) > 0) {
              m <- match(allpos, y@pos)
              yv[!is.na(m)] <- y@value[m[!is.na(m)]]
            }
            summed <- xv + yv
            keep <- which(summed != 0)
            if(length(keep) == 0L) {
              return(new("sparse_numeric", value = numeric(0), pos = integer(0), length = x@length))
            }
            new("sparse_numeric", value = summed[keep], pos = as.integer(allpos[keep]), length = x@length)
          })

setMethod("sparse_sub", signature(x = "sparse_numeric", y = "sparse_numeric"),
          function(x, y, ...) {
            .check_length_equal(x, y)
            if(length(x@pos) == 0L && length(y@pos) == 0L) {
              return(new("sparse_numeric", value = numeric(0), pos = integer(0), length = x@length))
            }
            allpos <- sort(unique(c(x@pos, y@pos)))
            xv <- numeric(length(allpos)); yv <- numeric(length(allpos))
            if(length(x@pos) > 0) {
              m <- match(allpos, x@pos)
              xv[!is.na(m)] <- x@value[m[!is.na(m)]]
            }
            if(length(y@pos) > 0) {
              m <- match(allpos, y@pos)
              yv[!is.na(m)] <- y@value[m[!is.na(m)]]
            }
            res <- xv - yv
            keep <- which(res != 0)
            if(length(keep) == 0L) {
              return(new("sparse_numeric", value = numeric(0), pos = integer(0), length = x@length))
            }
            new("sparse_numeric", value = res[keep], pos = as.integer(allpos[keep]), length = x@length)
          })

## sparse * sparse (element-wise multiplication)
setMethod("sparse_mult", signature(x = "sparse_numeric", y = "sparse_numeric"),
          function(x, y, ...) {
            .check_length_equal(x, y)
            if(length(x@pos) == 0L || length(y@pos) == 0L) {
              return(new("sparse_numeric", value = numeric(0), pos = integer(0), length = x@length))
            }
            # Only positions in intersection can be non-zero
            inter <- intersect(x@pos, y@pos)
            if(length(inter) == 0L) {
              return(new("sparse_numeric", value = numeric(0), pos = integer(0), length = x@length))
            }
            mx <- match(inter, x@pos); my <- match(inter, y@pos)
            vals <- x@value[mx] * y@value[my]
            keep <- which(vals != 0)
            if(length(keep) == 0L) {
              return(new("sparse_numeric", value = numeric(0), pos = integer(0), length = x@length))
            }
            new("sparse_numeric", value = vals[keep], pos = as.integer(inter[keep]), length = x@length)
          })

setMethod("sparse_crossprod", signature(x = "sparse_numeric", y = "sparse_numeric"),
          function(x, y, ...) {
            .check_length_equal(x, y)
            if(length(x@pos) == 0L || length(y@pos) == 0L) return(0)
            inter <- intersect(x@pos, y@pos)
            if(length(inter) == 0L) return(0)
            mx <- match(inter, x@pos); my <- match(inter, y@pos)
            sum(x@value[mx] * y@value[my])
          })
setMethod("sparse_add", signature(x = "sparse_numeric", y = "numeric"),
          function(x, y, ...) {
            sparse_add(x, as(y, "sparse_numeric"), ...)
          })

setMethod("sparse_add", signature(x = "numeric", y = "sparse_numeric"),
          function(x, y, ...) {
            sparse_add(as(x, "sparse_numeric"), y, ...)
          })

setMethod("sparse_sub", signature(x = "sparse_numeric", y = "numeric"),
          function(x, y, ...) {
            sparse_sub(x, as(y, "sparse_numeric"), ...)
          })

setMethod("sparse_sub", signature(x = "numeric", y = "sparse_numeric"),
          function(x, y, ...) {
            sparse_sub(as(x, "sparse_numeric"), y, ...)
          })

setMethod("sparse_mult", signature(x = "sparse_numeric", y = "numeric"),
          function(x, y, ...) {
            sparse_mult(x, as(y, "sparse_numeric"), ...)
          })

setMethod("sparse_mult", signature(x = "numeric", y = "sparse_numeric"),
          function(x, y, ...) {
            sparse_mult(as(x, "sparse_numeric"), y, ...)
          })

setMethod("sparse_crossprod", signature(x = "sparse_numeric", y = "numeric"),
          function(x, y, ...) sparse_crossprod(x, as(y, "sparse_numeric"), ...))

setMethod("sparse_crossprod", signature(x = "numeric", y = "sparse_numeric"),
          function(x, y, ...) sparse_crossprod(as(x, "sparse_numeric"), y, ...))

##########################
## Map operators + - * to sparse methods for sparse x sparse
##########################

setMethod("+", signature(e1 = "sparse_numeric", e2 = "sparse_numeric"),
          function(e1, e2) sparse_add(e1, e2))

setMethod("-", signature(e1 = "sparse_numeric", e2 = "sparse_numeric"),
          function(e1, e2) sparse_sub(e1, e2))

setMethod("*", signature(e1 = "sparse_numeric", e2 = "sparse_numeric"),
          function(e1, e2) sparse_mult(e1, e2))

## Also provide methods when numeric is on one side (coerce)
setMethod("+", signature(e1 = "sparse_numeric", e2 = "numeric"),
          function(e1, e2) sparse_add(e1, as(e2, "sparse_numeric")))

setMethod("+", signature(e1 = "numeric", e2 = "sparse_numeric"),
          function(e1, e2) sparse_add(as(e1, "sparse_numeric"), e2))

setMethod("-", signature(e1 = "sparse_numeric", e2 = "numeric"),
          function(e1, e2) sparse_sub(e1, as(e2, "sparse_numeric")))

setMethod("-", signature(e1 = "numeric", e2 = "sparse_numeric"),
          function(e1, e2) sparse_sub(as(e1, "sparse_numeric"), e2))

setMethod("*", signature(e1 = "sparse_numeric", e2 = "numeric"),
          function(e1, e2) sparse_mult(e1, as(e2, "sparse_numeric")))

setMethod("*", signature(e1 = "numeric", e2 = "sparse_numeric"),
          function(e1, e2) sparse_mult(as(e1, "sparse_numeric"), e2))

##########################
## show() method (compact display)
##########################
setMethod("show", "sparse_numeric", function(object) {
  cat("An object of class \"sparse_numeric\"\n")
  cat("  length :", object@length, "\n")
  if(length(object@pos) == 0L) {
    cat("  (all zeros)\n")
  } else {
    # show up to first 10 non-zero entries
    nshow <- min(10L, length(object@pos))
    idx <- seq_len(nshow)
    cat("  non-zero entries (first", nshow, "shown):\n")
    df <- data.frame(pos = object@pos[idx], value = object@value[idx])
    print(df, row.names = FALSE)
    if(length(object@pos) > nshow) cat("  ... (", length(object@pos) - nshow, "more)\n")
  }
})

setMethod("plot", signature(x = "sparse_numeric", y = "sparse_numeric"),
          function(x, y, ...) {
            .check_length_equal(x, y)
            # positions where BOTH have non-zero entries (overlap)
            inter <- intersect(x@pos, y@pos)
            if(length(inter) == 0L) {
              plot.new()
              title(main = "No overlapping non-zero positions", ...)
              return(invisible(NULL))
            }
            mx <- match(inter, x@pos); my <- match(inter, y@pos)
            xv <- x@value[mx]; yv <- y@value[my]
            # scatterplot: x-values (values from x) vs y-values (values from y),
            # plotting at their positions (we use position for x-axis)
            plot(inter, xv, xlab = "position", ylab = "value",
                 main = "overlapping non-zero elements (x in points, y in red)",
                 pch = 19, ...)
            points(inter, yv, pch = 17)
            legend("topright", legend = c("x", "y"), pch = c(19,17))
            invisible(NULL)
          })

##########################
## Extra method: sparse_norm (L^p norm)
##########################
setMethod("sparse_norm", signature(x = "sparse_numeric"),
          function(x, p = 2, ...) {
            if(length(x@value) == 0L) return(0)
            if(p == Inf) return(max(abs(x@value)))
            sum(abs(x@value)^p)^(1/p)
          })
