# Rather than do string manipulation in R, just make a hash map
ctx_to_mutsig_map <- c(
    `AAA:A>C`="TTT:T>G", `AAC:A>C`="GTT:T>G", `AAG:A>C`="CTT:T>G", `AAT:A>C`="ATT:T>G",
    `CAA:A>C`="TTG:T>G", `CAC:A>C`="GTG:T>G", `CAG:A>C`="CTG:T>G", `CAT:A>C`="ATG:T>G",
    `GAA:A>C`="TTC:T>G", `GAC:A>C`="GTC:T>G", `GAG:A>C`="CTC:T>G", `GAT:A>C`="ATC:T>G",
    `TAA:A>C`="TTA:T>G", `TAC:A>C`="GTA:T>G", `TAG:A>C`="CTA:T>G", `TAT:A>C`="ATA:T>G",
    `AAA:A>G`="TTT:T>C", `AAC:A>G`="GTT:T>C", `AAG:A>G`="CTT:T>C", `AAT:A>G`="ATT:T>C",
    `CAA:A>G`="TTG:T>C", `CAC:A>G`="GTG:T>C", `CAG:A>G`="CTG:T>C", `CAT:A>G`="ATG:T>C",
    `GAA:A>G`="TTC:T>C", `GAC:A>G`="GTC:T>C", `GAG:A>G`="CTC:T>C", `GAT:A>G`="ATC:T>C",
    `TAA:A>G`="TTA:T>C", `TAC:A>G`="GTA:T>C", `TAG:A>G`="CTA:T>C", `TAT:A>G`="ATA:T>C",
    `AAA:A>T`="TTT:T>A", `AAC:A>T`="GTT:T>A", `AAG:A>T`="CTT:T>A", `AAT:A>T`="ATT:T>A",
    `CAA:A>T`="TTG:T>A", `CAC:A>T`="GTG:T>A", `CAG:A>T`="CTG:T>A", `CAT:A>T`="ATG:T>A",
    `GAA:A>T`="TTC:T>A", `GAC:A>T`="GTC:T>A", `GAG:A>T`="CTC:T>A", `GAT:A>T`="ATC:T>A",
    `TAA:A>T`="TTA:T>A", `TAC:A>T`="GTA:T>A", `TAG:A>T`="CTA:T>A", `TAT:A>T`="ATA:T>A",
    `ACA:C>A`="ACA:C>A", `ACC:C>A`="ACC:C>A", `ACG:C>A`="ACG:C>A", `ACT:C>A`="ACT:C>A",
    `CCA:C>A`="CCA:C>A", `CCC:C>A`="CCC:C>A", `CCG:C>A`="CCG:C>A", `CCT:C>A`="CCT:C>A",
    `GCA:C>A`="GCA:C>A", `GCC:C>A`="GCC:C>A", `GCG:C>A`="GCG:C>A", `GCT:C>A`="GCT:C>A",
    `TCA:C>A`="TCA:C>A", `TCC:C>A`="TCC:C>A", `TCG:C>A`="TCG:C>A", `TCT:C>A`="TCT:C>A",
    `ACA:C>G`="ACA:C>G", `ACC:C>G`="ACC:C>G", `ACG:C>G`="ACG:C>G", `ACT:C>G`="ACT:C>G",
    `CCA:C>G`="CCA:C>G", `CCC:C>G`="CCC:C>G", `CCG:C>G`="CCG:C>G", `CCT:C>G`="CCT:C>G",
    `GCA:C>G`="GCA:C>G", `GCC:C>G`="GCC:C>G", `GCG:C>G`="GCG:C>G", `GCT:C>G`="GCT:C>G",
    `TCA:C>G`="TCA:C>G", `TCC:C>G`="TCC:C>G", `TCG:C>G`="TCG:C>G", `TCT:C>G`="TCT:C>G",
    `ACA:C>T`="ACA:C>T", `ACC:C>T`="ACC:C>T", `ACG:C>T`="ACG:C>T", `ACT:C>T`="ACT:C>T",
    `CCA:C>T`="CCA:C>T", `CCC:C>T`="CCC:C>T", `CCG:C>T`="CCG:C>T", `CCT:C>T`="CCT:C>T",
    `GCA:C>T`="GCA:C>T", `GCC:C>T`="GCC:C>T", `GCG:C>T`="GCG:C>T", `GCT:C>T`="GCT:C>T",
    `TCA:C>T`="TCA:C>T", `TCC:C>T`="TCC:C>T", `TCG:C>T`="TCG:C>T", `TCT:C>T`="TCT:C>T",
    `AGA:G>A`="TCT:C>T", `AGC:G>A`="GCT:C>T", `AGG:G>A`="CCT:C>T", `AGT:G>A`="ACT:C>T",
    `CGA:G>A`="TCG:C>T", `CGC:G>A`="GCG:C>T", `CGG:G>A`="CCG:C>T", `CGT:G>A`="ACG:C>T",
    `GGA:G>A`="TCC:C>T", `GGC:G>A`="GCC:C>T", `GGG:G>A`="CCC:C>T", `GGT:G>A`="ACC:C>T",
    `TGA:G>A`="TCA:C>T", `TGC:G>A`="GCA:C>T", `TGG:G>A`="CCA:C>T", `TGT:G>A`="ACA:C>T",
    `AGA:G>C`="TCT:C>G", `AGC:G>C`="GCT:C>G", `AGG:G>C`="CCT:C>G", `AGT:G>C`="ACT:C>G",
    `CGA:G>C`="TCG:C>G", `CGC:G>C`="GCG:C>G", `CGG:G>C`="CCG:C>G", `CGT:G>C`="ACG:C>G",
    `GGA:G>C`="TCC:C>G", `GGC:G>C`="GCC:C>G", `GGG:G>C`="CCC:C>G", `GGT:G>C`="ACC:C>G",
    `TGA:G>C`="TCA:C>G", `TGC:G>C`="GCA:C>G", `TGG:G>C`="CCA:C>G", `TGT:G>C`="ACA:C>G",
    `AGA:G>T`="TCT:C>A", `AGC:G>T`="GCT:C>A", `AGG:G>T`="CCT:C>A", `AGT:G>T`="ACT:C>A",
    `CGA:G>T`="TCG:C>A", `CGC:G>T`="GCG:C>A", `CGG:G>T`="CCG:C>A", `CGT:G>T`="ACG:C>A",
    `GGA:G>T`="TCC:C>A", `GGC:G>T`="GCC:C>A", `GGG:G>T`="CCC:C>A", `GGT:G>T`="ACC:C>A",
    `TGA:G>T`="TCA:C>A", `TGC:G>T`="GCA:C>A", `TGG:G>T`="CCA:C>A", `TGT:G>T`="ACA:C>A",
    `ATA:T>A`="ATA:T>A", `ATC:T>A`="ATC:T>A", `ATG:T>A`="ATG:T>A", `ATT:T>A`="ATT:T>A",
    `CTA:T>A`="CTA:T>A", `CTC:T>A`="CTC:T>A", `CTG:T>A`="CTG:T>A", `CTT:T>A`="CTT:T>A",
    `GTA:T>A`="GTA:T>A", `GTC:T>A`="GTC:T>A", `GTG:T>A`="GTG:T>A", `GTT:T>A`="GTT:T>A",
    `TTA:T>A`="TTA:T>A", `TTC:T>A`="TTC:T>A", `TTG:T>A`="TTG:T>A", `TTT:T>A`="TTT:T>A",
    `ATA:T>C`="ATA:T>C", `ATC:T>C`="ATC:T>C", `ATG:T>C`="ATG:T>C", `ATT:T>C`="ATT:T>C",
    `CTA:T>C`="CTA:T>C", `CTC:T>C`="CTC:T>C", `CTG:T>C`="CTG:T>C", `CTT:T>C`="CTT:T>C",
    `GTA:T>C`="GTA:T>C", `GTC:T>C`="GTC:T>C", `GTG:T>C`="GTG:T>C", `GTT:T>C`="GTT:T>C",
    `TTA:T>C`="TTA:T>C", `TTC:T>C`="TTC:T>C", `TTG:T>C`="TTG:T>C", `TTT:T>C`="TTT:T>C",
    `ATA:T>G`="ATA:T>G", `ATC:T>G`="ATC:T>G", `ATG:T>G`="ATG:T>G", `ATT:T>G`="ATT:T>G",
    `CTA:T>G`="CTA:T>G", `CTC:T>G`="CTC:T>G", `CTG:T>G`="CTG:T>G", `CTT:T>G`="CTT:T>G",
    `GTA:T>G`="GTA:T>G", `GTC:T>G`="GTC:T>G", `GTG:T>G`="GTG:T>G", `GTT:T>G`="GTT:T>G",
    `TTA:T>G`="TTA:T>G", `TTC:T>G`="TTC:T>G", `TTG:T>G`="TTG:T>G", `TTT:T>G`="TTT:T>G"
)


#' Detect COSMIC-formatted SBS channels
#'
#' Some tools (SigProfiler*) and COSMIC represent SBS channels as N\[R>M\]N where
#' the Ns are the upstream and downstream single bases and R and A are the reference
#' and mutant bases, respectively. This function detects such formatting in a
#' vectorized way, allowing conversion.
#'
#' @param x Character vector of channel names. May be a mixture of multiple formats.
#' @return A boolean vector indicating which entries in `x` are in the COSMIC SBS
#'      format: N\[R>M\]N. Note this may not be collapsed such that R is a pyrimidine.
#' @export
is_cosmic_sbs <- function(x) {
    grepl(pattern='[ACGT]\\[[ACGT]>[ACGT]\\][ACGT]', x)
}


#' Convert COSMIC-formatted SBS channels
#'
#' Some tools (SigProfiler*) and COSMIC represent SBS channels as N\[R>M\]N where
#' the Ns are the upstream and downstream single bases and R and A are the reference
#' and mutant bases, respectively. This function converts the COSMIC format into
#' the format used in this library:
#'      N\[R>M\]N -> NRN:R>M
#'
#' @param x Character vector of SBS channel names. If `x` is already in simplesigs
#'      format, it is returned unchanged.
#' @return A character vector (not a factor!) in simplesigs format. Factor conversion
#'      is not performed in case, e.g., 192-to-96 dimensional collapse is not desired.
#' @export
cosmic_sbs_to_simplesigs <- function(x) {
    ifelse(is_cosmic_sbs(x),
        paste0(substr(x,1,1), substr(x,3,3), substr(x,7,7), ':', substr(x,3,5)),
        x)
}



#' Collapse 192 and 96-dimensional SBS spectra to 6-dimensional SBS6
#'
#' The SBS6 format only encodes the reference and mutated base with no local
#' nucleotide context. I.e., NRN:R>M becomes simply R>M.
#'
#' @param x Character vector of any SBS channel with 1 upstream and 1 downstream
#'      base (i.e., trinucleotide context).
#' @param convert_cosmic If `TRUE`, convert COSMIC-style SBS channel names via
#'      [cosmic_sbs_to_simplesigs()]. This comes with a small performance cost, so
#'      when `x` is very large and known to be in simplesigs format, setting to
#'      `FALSE` to skip this conversion may be useful.
#' @return A factor vector in SBS6 order.
#' @export
sbs6 <- function (x=c(), convert_cosmic=TRUE) {
    if (convert_cosmic)
        x <- cosmic_sbs_to_simplesigs(x)

    sbs6.channel.order <- c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G')
    factor(substr(x, 5, 7), levels = sbs6.channel.order)
}


#' Collapse 192-dimensional spectra to 96
#'
#' There are 192 combinations of the 16 trinucleotide contexts and 12 single base
#' changes. However, these are usually collapsed into 96 dimensions by standardizing
#' on the pyrimidine reference base. I.e., an A>C mutation is also a T:G mutation on
#' the opposing strand. This function maps 192-dimensional spectra to the usual
#' 96-dimensional representation using a lookup table, specifically avoiding string
#' parsing.  N.B. the 192-dimensional description DOES NOT refer to the commonly
#' used transcribed strand-aware mutational signature.
#'
#' @param x Character vector of trinucleotide contexts plus base changes. E.g.,
#'      "ACA:C>T".
#' @return A factor with the ref base standardized to the pyrimidine. The factor
#'      levels are ordered to match the standard SBS96 plots.
#' @seealso [sbs96()]
#' @export
ctx_to_mutsig <- function(x) {
    sbs96(ctx_to_mutsig_map[x])
}


#' Trinucleotide + base change 96-dimensional signatures
#'
#' Create a factor with levels that represent the standard 96-dimensional
#' single base substitution (SBS) mutational signature in order.
#'
#' @param x Character vector with format `trinucleotide context`:`base change`
#'      where base change is `reference base`>`alt base` and is
#'      from the perspective of the pyrimidine at the site. For example, an A
#'      to C mutation would be represented as T>G. A full example is: ATA:T>G.
#' @param convert_cosmic If `TRUE`, convert COSMIC-style SBS channel names via
#'      [cosmic_sbs_to_simplesigs()]. This comes with a small performance cost, so
#'      when `x` is very large and known to be in simplesigs format, setting to
#'      `FALSE` to skip this conversion may be useful.
#' @return A factor vector.
#' @seealso [ctx_to_mutsig] to collapse from a 192-dimensional SBS signature to
#'      a 96-dimensional signature and factorize in one step.
#' @export
sbs96 <- function (x=c(), convert_cosmic=TRUE) {
    # allow nested calls - e.g., sbs96(sbs96(x))
    if (is.factor(x))
        x <- as.character(x)

    if (convert_cosmic)
        x <- cosmic_sbs_to_simplesigs(x)

    sbs96.channel.order <- paste0(rep(c("A", "C", "G", "T"),
        each = 4), rep(c("C", "T"), each = 48), rep(c("A", "C",
        "G", "T"), times = 4), ":", rep(c("C", "T"), each = 48),
        ">", c(rep(c("A", "G", "T"), each = 16), rep(c("A", "C",
            "G"), each = 16)))
    factor(x, levels = sbs96.channel.order)
}

sbs_cols_map <- c(
    `C>A`='#00BFFF', `C>G`='#000000', `C>T`='#EE2C2C', `T>A`='#B6B6B6', `T>C`='#29CA00', `T>G`='#F99BAF',
    `G>T`='#00BFFF', `G>C`='#000000', `G>A`='#EE2C2C', `A>T`='#B6B6B6', `A>G`='#29CA00', `A>C`='#F99BAF')
    #`C>A`='deepskyblue', `C>G`='black', `C>T`='firebrick2', `T>A`='grey', `T>C`='chartreuse3', `T>G`='pink2',
    #`G>T`='deepskyblue', `G>C`='black', `G>A`='firebrick2', `A>T`='grey', `A>G`='chartreuse3', `A>C`='pink2')

#' Standard colors for single base substitutions
#'
#' Map base change strings of the form `ref base`>`alt base` to standard colors.
#'
#' @param x Character vector with the form "A>C". Purine ref bases are converted
#'      to pyrimidines. If not specified, will return a named vector of base change to
#'      color mappings.
#' @return A character vector of color names.
#' @export
sbs_cols <- function(x=names(sbs_cols_map)) {
    sbs_cols_map[x]
}

sbs96_cols_map <- setNames(
    sbs_cols(substr(levels(sbs96()), 5, 7)),  # AAA:C>A substr(5,7) is the muttype
    levels(sbs96())
)

#' Standard colors for 96-dimensional single base substitution signatures
#'
#' Map trinucleotide context + base change strings of the form "N`ref`N:`ref base`>`alt base`
#' to standard colors. E.g., "ACA:C>A".
#'
#' @param x Character vector with the form "N`ref`N:`ref`>`alt`". These character vectors
#'      must be expressed such that ref is a pyrimidine. This requires reverse complementing
#'      the trinucleotide context. If left unspecified, a named vector mapping each of the
#'      possible 96 trinucleotide context+base change combinations to the standard color is
#'      returned.
#' @return A character vector of color names.
#' @export
sbs96_cols <- function(x=names(sbs96_cols_map)) {
    sbs96_cols_map[x]
}


#' Plot SBS96 spectra with ggplot
#'
#' A simple convenience function to quickly plot SBS96 spectra. Handles coloring, ensures
#' that SBS96 channels with 0 count are not dropped, prevents a 96-element (!) color guide
#' from being drawn, suppresses x-axis tick labels and enforces an aspect ratio.
#'
#' @param tx Set this to TRUE if the signatures are annotated for transcribed-strand status
#'      T, U, B, Q, N.
#' @returns A list of ggplot2 elements that can be `+`ed to a ggplot.
#' @export
geom_sbs96 <- function(tx=FALSE) {
    cols <- sbs96_cols()
    if (tx)
        cols <- tx_cols(cols)
    list(ggplot2::scale_fill_manual(values=cols, guide='none'),
        ggplot2::geom_vline(xintercept=cumsum((1+tx)*rep(16,5))+0.5, linewidth=0.15),
        ggplot2::geom_bar(),
        ggplot2::scale_x_discrete(drop=FALSE),
        ggplot2::scale_y_continuous(expand=ggplot2::expansion(c(0, 0.05))),
        ggplot2::theme(aspect.ratio=1/4, axis.text.x=ggplot2::element_blank()))
}


#' Fancy SBS96 spectrum plot
#'
#' Make a fancier SBS96 spectrum plot with a header row of colored boxes describing the
#' base changes and monospaced font x-axis labels describing the trinucleotide contexts.
#' This plot hijacks [facet_grid()] to create the top-strip of colored boxes, so is not
#' as flexible as [geom_sbs96()] and is generally only useful to plot one or a very
#' small number of mutation signatures. This is a great choice for plotting a single
#' signature in a publication figure. Requires the ggh4x library.
#' IMPORTANT: this function hijacks facet_grid to create the colorful strip of mutation
#' types along the top edge. So the usual flexibility and utility of faceting is lost
#' when using this style of plot. 
#'
#' @param data data.frame containing a character column of SBS96 channels and a
#'      numeric weight column. IMPORTANT! Unlike other plotting functions, `data`
#'      must contain explicit entries for all SBS96 channels (e.g., 0s)!
#' @param mutsig Unquoted expression yielding a character vector suitable for conversion
#'      by [sbs96()]. Usually the name of the column in `data` containing SBS96 channels. Do
#'      not specify this as a string.
#' @param weight Unquoted expression to compute channel weights (bar heights). Usually
#'      the name of the column in `data`.
#' @param rows Unquoted expression to define rows in plot since facet_grid is hijacked.
#' @param scale Passed to facet_grid2(scale). Set to "free" to allow y-axis limits to change
#'      between `rows` facets.
#' @param aspect.ratio Numeric value describing the aspect ratio of all 6 panels combined.
#' @return A ggplot object.
#' @export
plot_fancy_sbs96 <- function(data, mutsig, weight, rows=NULL, aspect.ratio=1/4, scale='free_x') {
    color_strip <- ggh4x::strip_themed(
        background_x=ggh4x::elem_list_rect(linewidth=rep(0, 6), fill=unique(sbs96_cols())),
        text_x=ggh4x::elem_list_text(face=rep('bold',6), colour=c(rep(c('white'), each=6)))
    )
    ggplot2::ggplot(data, ggplot2::aes(x=sbs96({{ mutsig }}), fill=sbs96({{ mutsig }}), weight={{ weight }})) +
        ggplot2::scale_fill_manual(values=sbs96_cols(), guide='none') +
        ggplot2::geom_bar() + #width=1) +
        ggplot2::geom_vline(xintercept=16 + 0.9/2) +    # geom_bar width is 0.9 by default
        ggplot2::theme_classic() +
        ggplot2::scale_y_continuous(expand=ggplot2::expansion(c(0, 0.05))) +
        ggplot2::xlab('Trinucleotide context') +
        ggplot2::ylab(ggplot2::element_blank()) +
        ggplot2::scale_x_discrete(
            drop=TRUE,
            expand=ggplot2::expansion(c(0,0)),
            breaks=levels(sbs96()),
            labels=substr(levels(sbs96()), 1, 3)) + 
        ggplot2::theme(aspect.ratio=6 * aspect.ratio,
            panel.spacing=ggplot2::unit(0.1, units='mm'),
            axis.text.x=ggplot2::element_text(angle=90, size=5, hjust=1, vjust=1/2, family='mono')) +
        ggh4x::facet_grid2(rows=ggplot2::vars({{ rows }}),
            cols=ggplot2::vars(sbs6( {{ mutsig }})), scale=scale, strip=color_strip)
}


#' Order trinucleotide contexts
#'
#' @param x Character vector for which the first 3 letters in each element are
#'      trinucleotides.
#' @param convert_cosmic If `TRUE`, convert COSMIC-style SBS channel names via
#'      [cosmic_sbs_to_simplesigs()]. This comes with a small performance cost, so
#'      when `x` is very large and known to be in simplesigs format, setting to
#'      `FALSE` to skip this conversion may be useful.
#' @return A factor vector of flanking nucleotides in the SBS96 order.
#' @export
ctx3 <- function(x=c(), convert_cosmic=TRUE) {
    if (convert_cosmic)
        x <- cosmic_sbs_to_simplesigs(x)

    ctx3.channel.order <- c(
        'A_A', 'A_C', 'A_G', 'A_T',
        'C_A', 'C_C', 'C_G', 'C_T',
        'G_A', 'G_C', 'G_G', 'G_T',
        'T_A', 'T_C', 'T_G', 'T_T')
    factor(paste0(substr(x,1,1), '_', substr(x,3,3)), levels=ctx3.channel.order)
}


#' Unsuccessful attempt to make plot_fancy_sbs96 work for non-full spectra
#'
#' plot_fancy_sbs96 requires values for all 96 SBS96 channels which can be inconvenient.
#' Tried to make a function that would automatically handle the case where not all values
#' are present *without manipulating the data*. Looks like this isn't possible, but this
#' is the closest I've gotten. The main problem is using facet_grid to layout the colored
#' strip along the top of the plot. Using faceting in this way causes so many unusual
#' problems.
#' @param data don't use
#' @param mutsig don't use
#' @param weight don't use
#' @param aspect.ratio don't use
#' @return don't use
#'
plot_fancy_sbs96_2 <- function(data, mutsig, weight=1, aspect.ratio=1/4) {
    color_strip <- ggh4x::strip_themed(
        background_x=ggh4x::elem_list_rect(linewidth=rep(0, 6), fill=unique(sbs96_cols())),
        text_x=ggh4x::elem_list_text(face=rep('bold',6), colour=c(rep(c('white'), each=6)))
    )
    ggplot2::ggplot(data, ggplot2::aes(x=ctx3({{ mutsig }}), fill=sbs96({{ mutsig }}), weight={{ weight }})) +
        ggplot2::scale_fill_manual(values=sbs96_cols(), guide='none') +
        ggplot2::geom_bar() + #width=1) +
        ggplot2::geom_vline(xintercept=16 + 0.9/2) +    # geom_bar width is 0.9 by default
        ggplot2::theme_classic() +
        ggplot2::scale_y_continuous(expand=ggplot2::expansion(c(0, 0.05))) +
        ggplot2::xlab('Trinucleotide context') +
        ggplot2::ylab(ggplot2::element_blank()) +
        ggplot2::scale_x_discrete(
            drop=FALSE,
            expand=ggplot2::expansion(c(0,0))) +
            #labels=substr(levels(sbs96()), 1, 3)) + 
        ggplot2::theme(aspect.ratio=6 * aspect.ratio,
                       panel.spacing=ggplot2::unit(0.1, units='mm'),
                       axis.text.x=ggplot2::element_text(angle=90, size=5, hjust=1, vjust=1/2, family='mono')) +
        ggh4x::facet_grid2(cols=ggplot2::vars(sbs6( {{ mutsig }})), scale='free_x', strip=color_strip, drop=FALSE)
}
