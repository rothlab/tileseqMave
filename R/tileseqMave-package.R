# Package-level imports for tileseqMave.
#
# All NAMESPACE import directives are generated from the roxygen tags in this
# file (run devtools::document() after changing them). 

# --- Whole-package imports ---------------------------------------------------

#' @import yogitools
#' @import hgvsParseR

# --- Function-level imports --------------------------------------------------

#' @importFrom hash hash has.key keys
#' @importFrom pbmcapply pbmclapply
#' @importFrom optimization optim_nm
#' @importFrom RJSONIO toJSON fromJSON
#' @importFrom mclust Mclust
#' @importFrom beeswarm beeswarm bxplot
#' @importFrom stats setNames rnorm quantile na.omit median dpois rbeta qnorm
#'   predict.lm lm cor as.formula coefficients cov dnorm predict rpois sd var
#' @importFrom utils read.csv write.csv packageVersion data write.table head
#'   read.delim combn
#' @importFrom graphics par text abline mtext grid rect barplot axis legend
#'   layout hist plot.new pairs lines points arrows curve polygon segments
#'   strwidth
#' @importFrom grDevices pdf dev.off rainbow colorRampPalette gray.colors png
NULL

# --- Packages used only by inst/scripts/ -------------------------------------
# inst/scripts/ use them, but no code below R/ does. The references below only
# make R CMD check see them as used, without loading them when tileseqMave is
# loaded. The function is never called and never exported.
ignore_unused_imports <- function() {
  argparser::arg_parser
  bitops::bitAnd
  parallel::mclapply
  yogilog::new.logger
  # inst/scripts/popcodeOligos.R still calls Tm_NN/Tm_GC/Tm_Wallace, which
  # current TmCalculator releases no longer export (they are tm_nn etc.).
  # The script needs updating before TmCalculator is genuinely usable here.
  TmCalculator::tm_nn
}