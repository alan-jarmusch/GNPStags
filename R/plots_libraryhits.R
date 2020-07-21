#' Produce Plots from Library Hits
#'
#' This function works to generate plots which summarize the 4 different types of tags, including the GNPS Tags, NCBI taxonomy tags, UBERON tags, and lifestyle tags contribtued by users.
#' @param path Path to the downloaded table of "View all Library Hits" from molecular networking (classic), library search workflows, or feature-based molecular networking in GNPS.
#' @keywords
#' @export
#' @examples plots_libraryhit(path = "METABOLOMICS-SNETS-V2-f1a1f3a6-view_all_annotations_DB-main.tsv")
#' @author Alan K. Jarmusch 2020
#' @export
#'
plots_libraryhits <- function(path) {
  summary_TAGS_libraryhits()
  summary_UBERON_libraryhits()
  summary_NCBI_libraryhits()
  summary_Lifestyle_libraryhits()
}