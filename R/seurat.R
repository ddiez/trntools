#' trn_fit
#'
#' Fit a TRN model.
#'
#' @param x dataset with gene expression levels.
#' @param db list of transcription factors.
#' @param assay assay for Seurat objects.
#' @param slot slot for Seurat objects.
#' @param scale whether to scale the data (only when slot is not scale.data).
#' @param filter whether to filter out genes with zero standard deviation.
#' @param ... arguments passed down to methods.
#'
#' @export
#'
trn_fit <- function(x, db, ...) {
  UseMethod("trn_fit")
}

#' @rdname trn_fit
#' @export
trn_fit.Seurat <- function(x, db, assay = NULL, slot = "data", scale = TRUE, filter = TRUE, ...) {
  x <- SeuratObject::GetAssayData(x, slot = slot, assay = assay) %>% as.matrix()

  if (filter)
    x <- x[! apply(x, 1, sd) == 0, ]

  if (slot != "scale.data" && scale) {
    x <- t(scale(t(x)))
  }

  tfs <- intersect(db, rownames(x))

  tfa <- x[tfs, ]
  x <- x[setdiff(rownames(x), tfs), ]

  trn_lm(t(tfa), t(x))
}
