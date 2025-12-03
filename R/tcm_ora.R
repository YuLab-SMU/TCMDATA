#' A wrapper function of `clusterProfiler::enricher()` to do herb-target enrichment analysis.
#'
#' @param genes Character vector. Query genes such as DEGs or disease targets.
#' @param bg_universe Character vector or NULL. Background gene set. If NULL, defaults to all unique targets in TCMDATA.
#' @param type Character. Specifies which column of `tcm_data` to use as herb names. Default is `"Herb_pinyin_name"`.
#' @param pvalueCutoff Numeric. P-value threshold for significance (default 0.05).
#' @param qvalueCutoff Numeric. Q-value (FDR) threshold for significance (default 0.2).
#' @param pAdjustMethod Character. Multiple testing correction method (default `"BH"`).
#' @param minGSSize Integer. Minimum herb gene-set size to include (default 10).
#' @param maxGSSize Integer. Maximum herb gene-set size to include (default 5000).
#' @importFrom clusterProfiler enricher
#' @importFrom stats na.omit
#' @return A \code{enrichResult} instance
#' @export

herb_enricher <- function(genes,
                          bg_universe = NULL,
                          type = c("Herb_pinyin_name","Herb_cn_name","Herb_en_name"),
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2,
                          pAdjustMethod = "BH",
                          minGSSize = 10,
                          maxGSSize = 5000){
  type <- match.arg(type)

  if (is.null(bg_universe)){
    bg_universe <- unique(na.omit(as.character(tcm_data$target)))
  }
  df_map <- tcm_data[, c(type, "target")]

  TERM2GENE <- unique(
    na.omit(
      data.frame(
        term = df_map[[type]],
        gene = df_map$target,
        stringsAsFactors = FALSE)))

  query_genes <- unique(na.omit(as.character(genes)))

  herb_enrich <- clusterProfiler::enricher(
    gene = query_genes,
    universe = bg_universe,
    TERM2GENE = TERM2GENE,
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff,
    pAdjustMethod = pAdjustMethod,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize)

  return(herb_enrich)

}
