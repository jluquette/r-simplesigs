#' Composite SBS96, ID83 and DBS78 signatures
#'
#' A composite signature that jointly represents SNVs, indels and
#' doublet base substitutions. The CMP257 signature is the concatenation
#' of the typical SBS96, ID83 and DBS78 signatures, in that order.
#'
#' @param x A character vector of recognized signature channel names from
#'      [sbs96()], [id83()] or [dbs78()].
#' @return A factor vector
#' @seealso [sbs96()], [id83()], [dbs78()].
#' @export
cmp257 <- function(x=c()) {
    # The strategy below of applying each of sbs96, id83 and dbs78 allows
    # implicit conversions. E.g., sbs96 may map AAA:A>T, which is not part
    # of the SBS96 format, to the in-format equivalent TTT:T>A. 
    s <- as.character(sbs96(x))
    i <- as.character(id83(x))
    d <- as.character(dbs78(x))
    f <- ifelse(!is.na(s), s, ifelse(!is.na(i), i, d))
    factor(f, levels=c(levels(sbs96()), levels(id83()), levels(dbs78())))
}


#' Plot CMP257 spectra with ggplot
#'
#' A simple convenience function to plot SBS257 spectra. Handles coloring, ensures
#' that channels with 0 count are not dropped, prevents a giant color guide
#' from being drawn, suppresses x-axis tick labels and enforces an aspect ratio.
#' Sadly, this plot is currently rarely useful since the SBS mutations are generally
#' an order of magnitude more common than ID and DBS mutations. As a result, the ID
#' and DBS bars are often too small to read.
#'
#' @returns A list of ggplot2 elements that can be added to a ggplot.
#' @export
geom_cmp257 <- function() {
    list(ggplot2::scale_fill_manual(values=c(sbs96_cols(), id83_cols(), dbs78_cols()), guide='none'),
        ggplot2::geom_bar(),
        ggplot2::geom_vline(xintercept=cumsum(
                c(rep(16,6),                    # SBS96
                  6,6,6,6,6,6,6,6,6,6,6,6,11,   # ID83
                  9,6,9,6,9,6,6,9,9))+0.5,
            linewidth=c(rep(0.15, 5), 0.6, rep(0.15, 12), 0.6, rep(0.15, 9))),
        ggplot2::scale_x_discrete(drop=FALSE),
        ggplot2::scale_y_continuous(expand=ggplot2::expansion(c(0, 0.05))),
        ggplot2::theme(aspect.ratio=1/12, axis.text.x=ggplot2::element_blank()))
}
