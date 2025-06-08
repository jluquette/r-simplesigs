#' Indel 83-dimensional signatures
#'
#' Create a factor with levels that represent the standard 83-dimensional
#' indel signature channels in order. The standard format is N:\{Ins|Del\}:\{C|T|R|M\}:M
#' where
#' * N is the length of the indel
#' * The indel context is: C and T represent a single base deletion of either C or T, R represents
#'   a repeat and M represents microhomology
#' * M is length of the context, which differs depending on context
#'
#' @param x Character vector of indel channels in the form specified above.
#' @return A factor vector.
#' @export
id83 <- function(x=c()) {
    id83.channel.order <- paste(
        c(rep(1,24), rep(rep(2:5,each=6),2),c(2,3,3,4,4,4,5,5,5,5,5)),
        c(rep(c('Del', 'Ins'), each=12), rep(c('Del', 'Ins'), each=24), rep('Del', 11)),
        c(rep(rep(c('C','T'), each=6), 2), rep('R',48), rep('M', 11)),
        c(rep(0:5, 12), c(1, 1,2, 1:3, 1:5)),
        sep=':')
    factor(x, levels=id83.channel.order)
}

id83_cols_map <- setNames(
    rep(c('#FBBD75', '#FC7F24',   # orange
          '#B0DB8E', '#3B9F36',   # green
          '#FBC9B6', '#F8896D', '#EE453A', '#B91C22',    # red to pink
          '#C5D4E4', '#8DBAD2', '#4D98C6', '#1D65A8',    # blue to light blue
          '#E1E1EE', '#B5B6D6', '#8684BA', '#614398'),   # purples
        c(rep(6,12), 1,2,3,5)),
    levels(id83())
)

#' Standard colors for 83-dimensional indel signatures
#'
#' Map indel signatures of the form N:\{Ins|Del\}:\{C|T|R|M\}:M to the standard
#' color, where
#' * N is the length of the indel
#' * The indel context is: C and T represent a single base deletion of either C or T, R represents
#'   a repeat and M represents microhomology
#' * M is length of the context, which differs depending on context
#'
#' @param x Character vector in the above form.
#' @return A character vector of color names.
#' @export
id83_cols <- function(x=names(id83_cols_map)) {
    id83_cols_map[x]
}

#' Plot ID83 spectra with ggplot
#'
#' A simple convenience function to quickly plot ID83 spectra. Handles coloring, ensures
#' that ID83 channels with 0 count are not dropped, prevents a 83-element (!) color guide
#' from being drawn, suppresses x-axis tick labels and enforces an aspect ratio.
#'
#' @param tx Set this to TRUE if the signatures are annotated for transcribed-strand status
#'      T, U, B, Q, N.
#' @returns A list of ggplot2 elements that can be `+`ed to a ggplot.
#' @export
geom_id83 <- function(tx=FALSE) {
    cols <- id83_cols()
    if (tx)
        cols <- tx_cols(cols)

    list(ggplot2::scale_fill_manual(values=cols, guide='none'),
        ggplot2::geom_bar(),
        ggplot2::geom_vline(xintercept=cumsum((1+tx)*c(6,6,6,6,6,6,6,6,6,6,6,6))+0.5, linewidth=0.15),
        ggplot2::scale_x_discrete(drop=FALSE),
        ggplot2::scale_y_continuous(expand=ggplot2::expansion(c(0, 0.05))),
        ggplot2::theme(aspect.ratio=1/5, axis.text.x=ggplot2::element_blank()))
}
