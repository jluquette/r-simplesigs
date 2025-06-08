#' Factorize transcribed strand from signature channel
#'
#' This function can either extract the transcribed strand tag from a
#' signature channel string (when only the `x` parameter is provided) or
#' create a new channel factor by combining strand tags passed via `x`
#' to another mutation signature channel factor vector `y`.
#'
#' @param x A character vector of the form `X:...` where X is one of the
#'      transcription strand indicators T, U, B, Q, N. For convenience,
#'      the following ... is not parsed, allowing this function to be
#'      applied to SBS, ID and DBS signatures.
#' @param y An optional factor vector of mutation signature channels.
#' @param plot If true, then only T and U factors will be created. All
#'      other possible factor levels will instead be set to NA to enable
#'      convenient filtering.
#'  
#' @return A factor vector with levels in the order: T, U, B, Q, N (when
#'      only `x` is provided) or a new factor of all (T,U,B,Q,N) x (signature
#'      channels) (in that order) when `y` is provided.
#' @export
tx <- function(x, y, plot=FALSE) {
    fct.lvls <- c('T', 'U')
    if (plot == FALSE)
        fct.lvls <- c(fct.lvls, 'B', 'Q', 'N')

    if (missing(y))
        factor(substr(x, 1, 1), levels=fct.lvls)
    else {
        mutsig.levels <- levels(y)

        # For DBS, all reverse complement-equal (AT, TA, CG, GC) and reverse-complement
        # ambiguous channels (AC, TG) cannot be assigned T or U.
        if (any(mutsig.levels == dbs78())) {
            mutsig.levels <- grep('^AC>..|^AT>..|^CG>..|^GC>..|^TA>..|^TG>..', mutsig.levels, invert=TRUE, value=TRUE)
        }

        factor(paste0(as.character(x), ':', as.character(y)),
            levels=apply(expand.grid(fct.lvls, mutsig.levels), 1, paste, collapse=':'))
    }
}

#' Create colors for transcribed strand bias plotting
#'
#' Transcribed strand colors are darker and more saturated than
#' untranscribed strand colors.
#'
#' @param colors A named character vector of colors in either hex RGB format,
#'      R color names or positive integer (i.e., any valid input to
#'      [grDevices::col2rgb]. The vector names are taken to be the factor
#'      levels of the mutation signature channels *without* transcribed
#'      strand information.
#' @param factor Scaling factor applied to saturation and value in HSV space
#'      to lighten the color.
#' @param eps Add eps to value (after converting to HSV space) before scaling
#'      by factor, or black cannot be lightened.
#' @return A character vector of hexadecimal color strings
#' @export
tx_cols <- function(colors, factor=4, eps=0.25) {
    utx.cols <- apply(t(grDevices::rgb2hsv(grDevices::col2rgb(colors))), 1,
        function(row) grDevices::hsv(row[1], row[2]/factor, min(0.8, row[3]+eps)^(1/factor)))
    stats::setNames(c(colors, utx.cols), c(paste0('T:', names(colors)), paste0('U:', names(colors))))
}
