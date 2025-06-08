#' Raw channel set for DBS78
dbs78.channel.order <- c(
    paste('AC', c('CA', 'CG', 'CT', 'GA', 'GG', 'GT', 'TA', 'TG', 'TT'), sep='>'),
    paste('AT', c('CA', 'CC', 'CG', 'GA', 'GC', 'TA'), sep='>'),
    paste('CC', c('AA', 'AG', 'AT', 'GA', 'GG', 'GT', 'TA', 'TG', 'TT'), sep='>'),
    paste('CG', c('AT', 'GC', 'GT', 'TA', 'TC', 'TT'), sep='>'),
    paste('CT', c('AA', 'AC', 'AG', 'GA', 'GC', 'GG', 'TA', 'TC', 'TG'), sep='>'),
    paste('GC', c('AA', 'AG', 'AT', 'CA', 'CG', 'TA'), sep='>'),
    paste('TA', c('AT', 'CG', 'CT', 'GC', 'GG', 'GT'), sep='>'),
    paste('TC', c('AA', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GG', 'GT'), sep='>'),
    paste('TG', c('AA', 'AC', 'AT', 'CA', 'CC', 'CT', 'GA', 'GC', 'GT'), sep='>'),
    paste('TT', c('AA', 'AC', 'AG', 'CA', 'CC', 'CG', 'GA', 'GC', 'GG'), sep='>')
)

#' There are 144 possible dinuc mutations XY>AB where X != A and Y != B.
#' (16 choices of XY with 9 choices of AB each).
#' 6 reference dinucs XY are always redundant:
#'     AA, AG, GA, GG, CA, GT,
#' so there are 10 non-redundant reference dinucs XY.
#' Given a reference dinuc XY, there are 9 possible mutations, so 90 in total.
#' Among these 90 possible mutations, 12 are redundant due to equivalence
#' under reverse complement, e.g. AT>CC = AT>GG.
dbs78_redundancy_map <- c(
    # The 78 channels in the specification map to themselves
    setNames(dbs78.channel.order, dbs78.channel.order),
    c(
        # The 12 reverse complement redundancies
        `AT>GG`='AT>CC',
        `AT>TC`='AT>GA',
        `AT>TG`='AT>CA',
        `CG>AA`='CG>TT',
        `CG>AC`='CG>GT',
        `CG>GA`='CG>TC',
        `GC>CT`='GC>AG',
        `GC>TG`='GC>CA',
        `GC>TT`='GC>AA',
        `TA>AC`='TA>GT',
        `TA>AG`='TA>CT',
        `TA>CC`='TA>GG',
        # The 54 possible mutations using the 6 disallowed reference dinucs:
        # AA, AG, GA, GG (0 pyrimidines) and CA, GT
        # are mapped to reverse complements. Note none of the reference dinucs
        # are equal under reverse complement, so all 9 possible mutations are
        # present for each case.
        # ref=AA
        `AA>CC`='TT>GG', `AA>CG`='TT>CG', `AA>CT`='TT>AG',
        `AA>GC`='TT>GC', `AA>GG`='TT>CC', `AA>GT`='TT>AC',
        `AA>TC`='TT>GA', `AA>TG`='TT>CA', `AA>TT`='TT>AA',
        # ref=AG
        `AG>CA`='CT>TG', `AG>CC`='CT>GG', `AG>CT`='CT>AG',
        `AG>GA`='CT>TC', `AG>GC`='CT>GC', `AG>GT`='CT>AC',
        `AG>TA`='CT>TA', `AG>TC`='CT>GA', `AG>TT`='CT>AA',
        # ref=GA
        `GA>AC`='TC>GT', `GA>AG`='TC>CT', `GA>AT`='TC>AT',
        `GA>CC`='TC>GG', `GA>CG`='TC>CG', `GA>CT`='TC>AG',
        `GA>TC`='TC>GA', `GA>TG`='TC>CA', `GA>TT`='TC>AA',
        # ref=GG
        `GG>AA`='CC>TT', `GG>AC`='CC>GT', `GG>AT`='CC>AT',
        `GG>CA`='CC>TG', `GG>CC`='CC>GG', `GG>CT`='CC>AG',
        `GG>TA`='CC>TA', `GG>TC`='CC>GA', `GG>TT`='CC>AA',
        # ref=CA
        `CA>AC`='TG>GT', `CA>AG`='TG>CT', `CA>AT`='TG>AT',
        `CA>GC`='TG>GC', `CA>GG`='TG>CC', `CA>GT`='TG>AC',
        `CA>TC`='TG>GA', `CA>TG`='TG>CA', `CA>TT`='TG>AA',
        # ref=GT
        `GT>AA`='AC>TT', `GT>AC`='AC>GT', `GT>AG`='AC>CT',
        `GT>CA`='AC>TG', `GT>CC`='AC>GG', `GT>CG`='AC>CG',
        `GT>TA`='AC>TA', `GT>TC`='AC>GA', `GT>TG`='AC>CA'
    )
)

#' Dinucleotide substitution (DBS78) signature channels
#'
#' Create a factor with levels that represent the standard 78-dimensional
#' DBS78 doublet substitution channels. Dinucleotide substitutions are
#' formatted similarly to SBS, i.e., XY>AB.
#' According to Bergstrom et al. 2019, when there are multiple ways
#' to represent a DBS, the one with the most pyrimidines (Cs or Ts) across
#' all 4 bases (XY, AB) is preferred, echoing the SBS strategy. However,
#' this appears not to be true. For example, AT>GA (1 pyrimidine) is
#' selected over the equivalent AT>TC (3 pyrimidines).
#' Despite this, at least at the level of reference dinucleotides XY, the
#' 4 purine-purine combinations are not used (AA, GG, AG, GA).
#' The 10 used dinucleotide reference bases (XY) are
#' * AC, AT, CC, CG, CT, GC, TA, TC, TG, TT.
#' All 16 possible dinucleotide bases are allowed as mutant bases (AB).
#' For the 4 reference dinucleotides that are equivalent to their reverse
#' complement (AT, TA, CG, GC), there are 6 possible mutant dinucleotides.
#' For the other 6, there are 9 possible mutant dinucleotides.
#'
#' Can be tested with
#'   nucs=c('A','C','G','T')
#'   z <- dbs78(data.table(expand.grid(X=nucs, Y=nucs, A=nucs, B=nucs))\[X != A & Y != B\]\[, paste0(X,Y,'>',A,B)\])
#'   # Every DBS should appear twice except the 12 that are reverse-complement equal.
#'   table(z)
#'
#' @param x Character vector of indel channels in the form specified above.
#' @return A factor vector.
#' @export
dbs78 <- function(x=c()) {
    factor(unname(dbs78_redundancy_map[x]), levels=dbs78.channel.order)
}

dinuc_cols_map <- c(
    AC='#07B4EA', AT='#025BC4',    # blues
    CC='#96C757', CG='#045A01',    # greens
    CT='#FE8E8E', GC='#DD2421',    # pinks/reds
    TA='#FDA75B', TC='#FC7400',    # oranges
    TG='#C48DFD', TT='#42048E')    # purples
    
dbs78_cols_map <- setNames(
    rep(dinuc_cols_map, c(9, 6, 9, 6, 9, 6, 6, 9, 9, 9)),
    levels(dbs78())
)

#' Standard colors for 78-dimensional dinucleotide signatures
#'
#' Map DBS78 signatures of the form XY>AB to the standard colors.
#'
#' @param x Character vector in DBS78 format.
#' @return A character vector of color names.
#' @export
dbs78_cols <- function(x=names(dbs78_cols_map)) {
    dbs78_cols_map[x]
}

#' Plot DBS78 spectra with ggplot
#'
#' A simple convenience function to quickly plot DBS78 spectra. Handles coloring, ensures
#' that DBS78 channels with 0 count are not dropped, prevents a 78-element (!) color guide
#' from being drawn, suppresses x-axis tick labels and enforces an aspect ratio.
#'
#' @param tx Set this to TRUE if the signatures are annotated for transcribed-strand status
#'      T, U, B, Q, N.
#' @param guide Plot a color key showing the reference dinucleotides.
#' @returns A list of ggplot2 elements that can be `+`ed to a ggplot.
#' @export
geom_dbs78 <- function(guide=FALSE, tx=FALSE) {
    # Adding dinuc_cols_map here allows the user to set `fill` to the ref dinuc to
    # produce a short color key.
    if (guide) {
        guide <- ggplot2::guide_legend(title='Ref')
    } else {
        guide <- 'none'
    }

    cols <- dbs78_cols()
    if (tx)
        cols <- tx_cols(cols)

    list(
        ggplot2::scale_fill_manual(values=c(dinuc_cols_map, cols), guide=guide),
        ggplot2::geom_bar(),
        ggplot2::geom_vline(xintercept=cumsum((1+tx)*c(9,6,9,6,9,6,6,9,9))+0.5, linewidth=0.15),
        ggplot2::scale_x_discrete(drop=FALSE), #, labels=substr(levels(dbs78()), 4, 5)),
        ggplot2::scale_y_continuous(expand=ggplot2::expansion(c(0, 0.05))),
        ggplot2::theme(aspect.ratio=1/5,
            axis.text.x=ggplot2::element_text(angle=90, hjust=1, vjust=1/2))
    )
}
