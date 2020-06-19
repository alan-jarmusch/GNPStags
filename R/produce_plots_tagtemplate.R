#' Produce Plots from Tag Template
#'
#' This function works to generate plots from user-completed tag template files (providing the most complete information about the data) which summarize the 4 different types of tags, including the GNPS Tags, NCBI taxonomy tags, UBERON tags, and lifestyle tags contribtued by users.
#' @param path Path to the downloaded GNPS tag template (Google Sheets).
#' @keywords
#' @export
#' @examples plots_tagtemplate(path = "TEST GNPS Tag - Batch Tag Template.tsv")
#' @author Alan K. Jarmusch 2020
#' @export
#'
plots_tagtemplate <- function(path) {
  summary_TAGS_tagtemplate()
  summary_UBERON_tagtemplate()
  summary_NCBI_tagtemplate()
  summary_Lifestyle_tagtemplate()
}
