#' Combine list objects
#'
#' Take a sequence of list arguments and combine element-wise.
#'
#' @param ... R lists objects of presumably the same length \code{n}.
#' Each list should have \code{n} elements, but each element of a given list can have any length.
#'
#' @param FUN the \link{function} to be used to combine list elements. Defaults to \link{c}.
# In the case of functions like \code{+}, \code{%*%}, the function name must be backquoted or quoted.
#'
#' All other arguments below are used only if \code{FUN} is missing or
#' \code{force = TRUE}. Otherwise, they are ignored. These arguments were added
#' just for the defaults function \link{c}. They do not necessary make sens for
#' user provided functions \code{FUN}. If it is certain that these arguments
#' make sense, then setting \code{force = TRUE} allows their use.
#'
#' @param unique logical, return as \code{j}th list element, the union of the \code{j}th elements of arguments \code{...}?
#' Defaults to \code{TRUE}. The alternative is to concatenate all elements, with possibly duplicates returned.
#'
#' @param keep.element.names logical, keep the names of the elements elements that are combined?
#' Only used when \code{unique = TRUE}, in which case setting \code{keep.element.names = TRUE}
#' will remove duplicates but will keep the name attributes of the list elements, if any.
#'
#' @param force logical only
#'
#' @param make.vector.names logical, name the elements of the returned list?
#' If \code{TRUE} (the default), the \code{name} attribute of the elements of the
#' \code{...} arguments with the larger length is copied to the returned result
#' (note that this \code{name} attribute can be \code{NULL}, resulting in a
#' \code{NULL} \code{name} attribute (i.e. no \code{name}!) for the result).
#'
#' @param verbose logical, should warnings be issued when the \code{...} arguments
#' have different lengths? Defaults to \code{TRUE}.
#'
#' @export lbind
#'
#' @details All arguments after \code{...} must be named.
#'
#' The elements of the \code{...} arguments are combined using \code{FUN}.
#' For the defaults function \link{c}, if all or some of these elements were
#' matrices, they will first be turned  into a vector (stacking all elements),
#' and then combined with the corresponding elements from other arguments in \code{...}.
#'
#' All arguments after \code{FUN} are used only if \code{FUN} is missing or
#' \code{force = TRUE}. Otherwise, they are ignored. These arguments were added
#' just for the defaults function \link{c}. They do not necessary make sense for
#' user provided functions \code{FUN}. Setting \code{force = TRUE} allows to
#' use these arguments with a user specified function \code{FUN}.
#'
#' If the length of any of the \code{...} arguments has a different length,
#' all missing list elements at a given position are considered \code{NULL},
#' and the warning "longer argument not a multiple of length of shorter" is
#' issued when \code{verbose = TRUE}.
#'

#
# Just calls the 'base' function mapply
lbind <- function(..., FUN = c, unique = TRUE,
                  keep.element.names = TRUE, force = FALSE,
                  make.vector.names = TRUE,
                  verbose = TRUE) {
  FUN <- match.fun(FUN)
  if (verbose)
    result <- mapply (function(...) FUN(...), ..., SIMPLIFY = FALSE)
  else
    result <- catch.conditions(mapply (function(...) FUN(...), ..., SIMPLIFY = FALSE))$value
  if (!missing(FUN) & !force) {
    return(result)
  }

  if (unique) {
    if (keep.element.names) { # If name attributes have to be kept
      result <- lapply(result, FUN = function(resj) {
        dups <- duplicated(resj)
        if (any(dups))
          resj <- resj[!dups] # Remove duplicates, if any.
        resj
      })
    }
    else { # If name attributes  can be ignored, use 'unique' which does not preserve name attributes.
      result <- lapply(result, unique)
    }
  }
  if (make.vector.names) {
    INargs <- list(...)
    largs <- sapply(INargs, FUN = length)
    names(result) <- names(INargs[[which.max(largs)[1]]])
  }
  return(result)
}
