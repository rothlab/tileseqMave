toMavedbCountTable <- function(data) {
  cn <- colnames(data)
  cn[cn == "hgvsc"] <- "hgvs_nt"
  cn[cn == "hgvsp"] <- "hgvs_pro"
  colnames(data) <- cn

  internalCols <- c(
    "codonChanges", "codonHGVS", "aaChanges", "aaChangeHGVS",
    "codonChange", "aaChange"
  )

  data[, !(cn %in% internalCols), drop = FALSE]
}
