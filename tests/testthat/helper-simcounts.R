# Writes minimal marginalCounts.csv and sampleDepths.csv for a parsed params
# object. Not a realistic simulation: just enough structure and variation
# for calcEnrichment to run deterministically.

simulateCountFiles <- function(params, outDir, seed = 1, tiles = NULL,
                               perTile = 25, depth = 2e5) {
  set.seed(seed)
  dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
  samples <- params$samples
  if (is.null(tiles)) tiles <- head(params$tiles[, "Tile Number"], 2)

  # sampleDepths.csv: sample table plus depth columns
  depths <- samples
  depths$alignedreads <- depth
  depths$wtreads <- depth
  depths$mutreads <- depth
  depths$depth <- depth
  write.csv(depths, file.path(outDir, "sampleDepths.csv"), row.names = FALSE)

  # variants: a mix of synonymous, missense and nonsense changes per tile
  aas <- c("A", "G", "L", "K", "D", "V")
  variants <- do.call(rbind, lapply(tiles, function(tl) {
    row <- which(params$tiles[, "Tile Number"] == tl)
    pos <- sample(
      params$tiles[row, "Start AA"]:params$tiles[row, "End AA"], perTile
    )
    from <- sample(aas, perTile, replace = TRUE)
    kind <- sample(c("syn", "mis", "non"), perTile, replace = TRUE,
                   prob = c(.2, .6, .2))
    to <- ifelse(kind == "syn", from,
          ifelse(kind == "non", "*",
                 sapply(from, function(f) sample(setdiff(aas, f), 1))))
    data.frame(pos = pos, from = from, to = to, kind = kind)
  }))
  n <- nrow(variants)
  aaChange <- paste0(variants$from, variants$pos, variants$to)

  # label columns, in the order calcEnrichment expects ([,1:4])
  out <- data.frame(
    hgvsc       = sprintf("c.%d_%d", variants$pos * 3, seq_len(n)), # unique
    hgvsp       = paste0("p.", aaChange),
    codonChange = paste0("AAA", variants$pos, "GGG"),
    aaChange    = aaChange
  )

  # "true" baseline frequency and fitness per variant
  f0 <- 10^runif(n, -4.5, -3)
  phi <- ifelse(variants$kind == "syn", 1,
         ifelse(variants$kind == "non", 0.02, runif(n, 0.05, 1.2)))

  selConds <- getSelects(params)
  wtConds <- unique(unlist(lapply(
    c(selConds, getNonselects(params)), getWTControlFor, params
  )))
  eff <- depth * 0.95  # depth drop 0.05, below the default maxDrop of 0.2

  for (cond in unique(samples$Condition)) {
    mu <- if (cond %in% wtConds) {
      rep(2e-5, n)
    } else if (cond %in% selConds) {
      f0 * phi
    } else {
      f0
    }
    condRows <- samples$Condition == cond
    for (tp in unique(samples$`Time point`[condRows])) {
      reps <- unique(samples$Replicate[condRows & samples$`Time point` == tp])
      for (repi in reps) {
        cnt <- rpois(n, mu * eff * exp(rnorm(n, 0, 0.1)))
        pre <- sprintf("%s.t%s.rep%d.", cond, tp, repi)
        out[[paste0(pre, "count")]] <- cnt
        out[[paste0(pre, "frequency")]] <- cnt / eff
        out[[paste0(pre, "effectiveDepth")]] <- eff
      }
    }
  }

  con <- file(file.path(outDir, "marginalCounts.csv"), open = "w")
  writeLines("# MARGINAL VARIANT COUNTS #", con)
  write.csv(out, con, row.names = FALSE)
  close(con)
  invisible(out)
}