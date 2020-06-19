#' Summary of UBERON Tags in GNPS - Library Annotations
#'
#' This function works to generate a summary of the UBERON tags associated with annotated chemicals resulting from spectral library matching in GNPS. UBERON tags are associated with planar InChIKeys.
#' @param library_annotations Path to the downloaded table of "View all Library Hits" from molecular networking (classic) or library search workflows in GNPS.
#' @keywords
#' @export
#' @examples 
#' @author Alan K. Jarmusch 2020
#' @export
#'
summary_UBERON_libraryhits <- function(){
  
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # INPUT: Master Tag Sheet (.tsv)
  df_master <- fread("https://docs.google.com/spreadsheets/d/e/2PACX-1vTNNhRZYDQ9IzS-B5MoulYi3Fmpr5H1STPGYIUzjrpSw2TyoJ0Yp9wKwMmz7T_j5zPjVeYU9rJIG_pM/pub?gid=830441611&single=true&output=tsv",sep="\t", header=TRUE, stringsAsFactors = FALSE)
  df_master_CCMSLIB <- fread("https://docs.google.com/spreadsheets/d/e/2PACX-1vTNNhRZYDQ9IzS-B5MoulYi3Fmpr5H1STPGYIUzjrpSw2TyoJ0Yp9wKwMmz7T_j5zPjVeYU9rJIG_pM/pub?gid=748450796&single=true&output=tsv",sep="\t", header=TRUE, stringsAsFactors = FALSE)

  # Assign Colors to Tag Types
  group.colors <- c(`InChIKey-Planar` = "#d8b365", CCMSLIB = "#5ab4ac")

  df_libraryhits <- fread(input = path, sep="\t", header=TRUE, stringsAsFactors = FALSE)

  df_libraryhits_tags <- separate(df_libraryhits, col = `InChIKey-Planar`, sep = "-", into = c("term1_match")) %>% select(term1_match)
  df_libraryhits_tags <- merge(df_libraryhits_tags, df_master, by.x="term1_match", by.y="InChI_Key_Planar")
  
  if (nrow(df_libraryhits_tags) > 0) {
    df_libraryhits_tags <- df_libraryhits_tags %>% distinct()
    colnames(df_libraryhits_tags)[1] <- "InChIKey-Planar"
    tags <- cbind(df_libraryhits_tags[,"InChIKey-Planar"],separate(df_libraryhits_tags[,"UBERONBodyPartName"], UBERONBodyPartName, paste0("X",1:max(sapply(strsplit(df_libraryhits_tags$UBERONBodyPartName,"\\|"),length))), sep="\\|"))
    colnames(tags)[2] <- "terms"
    try(tags <- gather(tags, key="tags", value="terms", 2:length(tags)), silent = TRUE)
    tag_summary <- tags %>% group_by(terms) %>% tally
    tag_summary <- subset(tag_summary, tag_summary$terms != "NA" & tag_summary$terms != "")
    tag_summary$type <- "InChIKey-Planar"
  }  else {
    print("no tags")
    tag_summary <- data.frame()
  }

  df_CCMSLIB_tags <- subset(df_libraryhits, df_libraryhits$`InChIKey-Planar` == "N/A")
  df_CCMSLIB_tags <- merge(df_CCMSLIB_tags, df_master_CCMSLIB, by.x="SpectrumID", by.y="CCMSLIB")
  if (nrow(df_CCMSLIB_tags) > 0) {
    df_CCMSLIB_tags <- df_CCMSLIB_tags %>% distinct()
    colnames(df_CCMSLIB_tags)[1] <- "CCMSLIB"
    tags_CCMS <- cbind(df_CCMSLIB_tags[,"CCMSLIB"],separate(df_CCMSLIB_tags[,"UBERONBodyPartName"], UBERONBodyPartName, paste0("X",1:max(sapply(strsplit(df_CCMSLIB_tags$UBERONBodyPartName,"\\|"),length))), sep="\\|"))
    colnames(tags_CCMS)[2] <- "terms"
    tags_CCMS <- subset(tags_CCMS, tags_CCMS$terms != "NA" & tags_CCMS$terms != "")
    if (nrow(tags_CCMS) > 0) {
      try(tags_CCMS <- gather(tags_CCMS, key="tags", value="terms", 2:length(tags)), silent = TRUE)
      tags_CCMS_summary <- tags_CCMS %>% group_by(terms) %>% tally
      tags_CCMS_summary <- subset(tags_CCMS_summary, tags_CCMS_summary$terms != "NA" & tags_CCMS_summary$terms != "")
      tags_CCMS_summary$type <- "CCMSLIB" }
    else {
      print("no tags")
      tags_CCMS_summary <- data.frame()
    }
  }  else {
    print("no tags")
    tags_CCMS_summary <- data.frame()
  }

  hit_table <- rbind(tag_summary,tags_CCMS_summary)

  # Plot
  plot_summary <- ggplot(data=hit_table, aes(x= reorder(as.factor(terms),n), y=as.numeric(n)))+
    geom_point(aes(group=type, colour=type), pch=16, alpha=0.7, cex=2.0)+
    geom_text(aes(label=n), size=2.5)+
    scale_y_continuous(limits=c(0,max(as.numeric(hit_table$n))+max(as.numeric(hit_table$n))*0.05))+
    scale_colour_manual(values=group.colors)+
    coord_flip()+
    theme_minimal()+
    theme(panel.grid.major=element_line(colour ="grey75",size=0.25, linetype="dotted"),
          panel.grid.minor=element_blank(),
          axis.ticks=element_line(colour ="black",size=0.25, linetype="solid"),
          axis.text=element_text(colour="black",size=6),
          axis.text.x=element_text(angle=0, hjust=0.5, vjust=1.25),
          axis.line=element_line(colour="black",size=0.25, linetype="solid"),
          axis.title=element_text(colour="black",size=6),
          strip.text.x = element_text(colour="black",size=6),
          #aspect.ratio=1.0,
          legend.title = element_text(colour="black", size=6),
          legend.text = element_text(colour="black", size=6),
          legend.justification = c(1,0),
          legend.position = c(1,0))+
    labs(x="UBERON Tags", y="Number of Tags")
  print(plot_summary)
  ggsave(file="GNPS_LibraryHitTable_Summary_UBERON.pdf", width = 3.25, height = 3.25, units = "in", dpi=300, useDingbats=FALSE)
  ggsave(file="GNPS_LibraryHitTable_Summary_UBERON.png", width = 3.25, height = 3.25, units = "in")
}
