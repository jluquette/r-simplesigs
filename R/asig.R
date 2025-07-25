#' Automatic signature detection and factorization
#'
#' Create a factor with levels in the appropriate order for a
#' character vector of unspecified mutational signature type.
#'
#' @param x Character vector of mutational signature channels from a
#'      supported signature type.
#' @param sigtype Single string corresponding to a supported signature
#'      type or `auto`. When `auto`, the appropriate mutational signature
#'      type is determined automatically. To override automatic determination,
#'      set `sigtype` to a string matching a
#'      signature type function (e.g., `"sbs96"` for [sbs96()]). There is no
#'      difference between `asig(x, sigtype="sbs96")` and running `sbs96(x)`
#'      directly, but allows convenient programmatic selection of a signature
#'      type function.
#'      IMPORTANT: if `sigtype=auto` and a type cannot be guessed, [stop()] is
#'      thrown.
#' @returns A factor vector of the guessed or requested signature type.
#' @seealso [sbs96()], [id83()], [dbs78()], [cmp257()]
#' @export
asig <- function(x=c(), sigtype=c('auto', 'sbs96', 'id83', 'dbs78')) {
    sigtype <- match.arg(sigtype)

    # Guess the type
    if (sigtype == "auto") {
        # Get the unique values once so that all comparisons below are quick.
        u <- unique(x)

        # sigfunc(.)=NA only for strings unrecognized by that signature type. By
        # calling sigfunc() rather than, say, comparing u to levels(sigfunc()),
        # we allow for more flexibility. E.g., if sigfunc() collapses
        # multiple mutation types into a single channel (e.g., collapsing A>T
        # to T>A), the levels(.) approach would fail while sum(is.na(sigfunc(.)))
        # would succeed.
        if (sum(is.na(sbs96(u))) == 0) {
            sigtype <- 'sbs96'
        } else if (sum(is.na(id83(u))) == 0) {
            sigtype <- 'id83'
        } else if (sum(is.na(dbs78(u))) == 0) {
            sigtype <- 'dbs78'
        } else if (sum(is.na(cmp257(u))) == 0) {
            sigtype <- 'cmp257'
        } else {
            stop("sigtype=auto but signature type could not be guessed")
        }
    }

    if (sigtype == 'sbs96')
        return(sbs96(x))
    if (sigtype == 'id83')
        return(id83(x))
    if (sigtype == 'dbs78')
        return(dbs78(x))
    if (sigtype == 'cmp257')
        return(cmp257(x))
}
