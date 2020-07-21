#' Summary of GNPS Lifestyle Tags - Library Annotations
#'
#' This function works to generate a summary of the GNPS Lifestyle tags associated with annotated chemicals resulting from spectral library matching in GNPS. GNPS Lifestyle tags are associated with planar InChIKeys.
#' @param library_annotations Path to the downloaded table of "View all Library Hits" from molecular networking (classic) or library search workflows in GNPS.
#' @keywords
#' @export
#' @examples 
#' @author Alan K. Jarmusch 2020
#' @export

summary_Lifestyle_libraryhits <- function(){
  
  suppressMessages(library(data.table))
  suppressMessages(library(dplyr))
  suppressMessages(library(tidyr))
  suppressMessages(library(ggplot2))
  suppressMessages(library(zip))
  
  # INPUT: Master Tag Sheet (.tsv)
  df_master <- fread("https://docs.google.com/spreadsheets/d/e/2PACX-1vTNNhRZYDQ9IzS-B5MoulYi3Fmpr5H1STPGYIUzjrpSw2TyoJ0Yp9wKwMmz7T_j5zPjVeYU9rJIG_pM/pub?gid=830441611&single=true&output=tsv",sep="\t", header=TRUE, stringsAsFactors = FALSE)
  df_master_CCMSLIB <- fread("https://docs.google.com/spreadsheets/d/e/2PACX-1vTNNhRZYDQ9IzS-B5MoulYi3Fmpr5H1STPGYIUzjrpSw2TyoJ0Yp9wKwMmz7T_j5zPjVeYU9rJIG_pM/pub?gid=748450796&single=true&output=tsv",sep="\t", header=TRUE, stringsAsFactors = FALSE)

  # Assign Colors to Tag Types
  group.colors <- c(`InChIKey-Planar` = "#d8b365", CCMSLIB = "#5ab4ac")
  
  df_libraryhits <- fread(input = path, sep="\t", header=TRUE, stringsAsFactors = FALSE)
  studyinfo_n_row_annotations <- as.numeric(length(df_libraryhits$`#Scan#`)) #numbers out
  studyinfo_n_rows_with_InChI <- subset(df_libraryhits, df_libraryhits$`InChIKey-Planar` != "N/A")
  studyinfo_n_unique_InChI <- as.numeric(length(unique(studyinfo_n_rows_with_InChI$`InChIKey-Planar`))) #numbers out
  studyinfo_n_rows_with_InChI <- as.numeric(length(studyinfo_n_rows_with_InChI$`InChIKey-Planar`)) #numbers out
  
  df_libraryhits_tags <- separate(df_libraryhits, col = `InChIKey-Planar`, sep = "-", into = c("term1_match")) %>% select(term1_match)
  df_libraryhits_tags <- merge(df_libraryhits_tags, df_master, by.x="term1_match", by.y="InChI_Key_Planar")
  colnames(df_libraryhits_tags)[1] <- "InChI_Key_Planar"
  studyinfo_n_unique_InChI_w_tag <- as.numeric(length(unique(df_libraryhits_tags$InChI_Key_Planar))) #numbers out
  
  if (nrow(df_libraryhits_tags) > 0) {
    df_libraryhits_tags <- df_libraryhits_tags %>% distinct()
    colnames(df_libraryhits_tags)[1] <- "InChIKey-Planar"
    tags <- cbind(df_libraryhits_tags[,"InChIKey-Planar"],separate(df_libraryhits_tags[,"Lifestyle_Tag"], Lifestyle_Tag, paste0("X",1:max(sapply(strsplit(df_libraryhits_tags$Lifestyle_Tag,"\\|"),length))), sep="\\|"))
    tags <- gather(tags, key="tags", value="terms", 2:length(tags))
    tag_summary <- tags %>% group_by(terms) %>% tally
    tag_summary <- subset(tag_summary, tag_summary$terms != "NA" & tag_summary$terms != "")
    tag_summary$type <- "InChIKey-Planar"
  }  else {
    print("no tags")
    tag_summary <- data.frame()
  }
  
  tag_summary <- tag_summary %>% mutate(percentchemicalstagged = 100*(n / studyinfo_n_unique_InChI_w_tag))
  
  df_CCMSLIB_tags <- subset(df_libraryhits, df_libraryhits$`InChIKey-Planar` == "N/A")
  studyinfo_CCMS_n_rows <- as.numeric(length(df_CCMSLIB_tags$`#Scan#`)) #numbers out
  studyinfo_CCMS_n_rows_with_uniqueCCMS <- as.numeric(length(unique(df_CCMSLIB_tags$SpectrumID))) #numbers out
  df_CCMSLIB_tags <- merge(df_CCMSLIB_tags, df_master_CCMSLIB, by.x="SpectrumID", by.y="CCMSLIB")
  colnames(df_CCMSLIB_tags)[1] <- "CCMSLIB"
  studyinfo_CCMS_n_unique_CCMS_w_tag <- as.numeric(length(unique(df_CCMSLIB_tags$CCMSLIB))) #numbers out
  
  if (nrow(df_CCMSLIB_tags) > 0) {
    df_CCMSLIB_tags <- df_CCMSLIB_tags %>% distinct()
    colnames(df_CCMSLIB_tags)[1] <- "CCMSLIB"
    tags_CCMS <- cbind(df_CCMSLIB_tags[,"CCMSLIB"],separate(df_CCMSLIB_tags[,"Lifestyle_Tag"], Lifestyle_Tag, paste0("X",1:max(sapply(strsplit(df_CCMSLIB_tags$Lifestyle_Tag,"\\|"),length))), sep="\\|"))
    colnames(tags_CCMS)[2] <- "terms"
    tags_CCMS <- subset(tags_CCMS, tags_CCMS$terms != "NA" & tags_CCMS$terms != "")
    if (nrow(tags_CCMS) > 0) {
      try(tags_CCMS <- gather(tags_CCMS, key="tags", value="terms", 2:length(tags)), silent = TRUE)
      tags_CCMS_summary <- tags_CCMS %>% group_by(terms) %>% tally
      tags_CCMS_summary <- subset(tags_CCMS_summary, tags_CCMS_summary$terms != "NA" & tags_CCMS_summary$terms != "")
      tags_CCMS_summary$type <- "CCMSLIB" 
      tags_CCMS_summary <- tags_CCMS_summary %>% mutate(percentchemicalstagged = 100*(n / studyinfo_CCMS_n_unique_CCMS_w_tag))}
    else {
      print("no tags")
      tags_CCMS_summary <- data.frame()
    }
  }  else {
    print("no tags")
    tags_CCMS_summary <- data.frame()
  }
  
  hit_table <- rbind(tag_summary,tags_CCMS_summary)
  write.csv(hit_table, "GNPS_hittable_LifestyleTag.csv", row.names=FALSE)
  
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
    labs(x="Lifestyle Tags", y="Number of Tags")
  print(plot_summary)
  ggsave(file="GNPS_LibraryHitTable_Summary_Lifestyle.pdf", width = 3.25, height = 3.25, units = "in", dpi=300, useDingbats=FALSE)
  ggsave(file="GNPS_LibraryHitTable_Summary_Lifestyle.png", width = 3.25, height = 3.25, units = "in")
  
  # Plot
  plot_summary <- ggplot(data=hit_table, aes(x= reorder(as.factor(terms),percentchemicalstagged), y=as.numeric(percentchemicalstagged)))+
    geom_point(aes(group=type, colour=type), pch=16, alpha=0.7, cex=2.0)+
    geom_text(aes(label=round(percentchemicalstagged,2)), size=2.5)+
    scale_y_continuous(limits=c(0,max(as.numeric(hit_table$percentchemicalstagged))+max(as.numeric(hit_table$percentchemicalstagged))*0.05))+
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
    labs(x="Lifestyle Tags", y="% of Unique Planar InChIKeys or CCMSLIB IDs")
  print(plot_summary)
  ggsave(file="GNPS_LibraryHitTable_SummaryPercentage_Lifestyle.pdf", width = 3.25, height = 3.25, units = "in", dpi=300, useDingbats=FALSE)
  ggsave(file="GNPS_LibraryHitTable_SummaryPercentage_Lifestyle.png", width = 3.25, height = 3.25, units = "in")
  
  #numbers
  tag_numbers <- rbind(studyinfo_n_row_annotations, studyinfo_n_rows_with_InChI, 
                       studyinfo_n_unique_InChI,studyinfo_n_unique_InChI_w_tag, studyinfo_CCMS_n_rows,
                       studyinfo_CCMS_n_rows_with_uniqueCCMS, studyinfo_CCMS_n_unique_CCMS_w_tag)
  rownames(tag_numbers) <- c("n_GNPSannotations", "n_GNPSannotations_with_InChIKeyPlanar", 
                             "n_unique_InChIKeyPlanar", "n_InChIKeyPlanar_with1moretag", 
                             "n_GNPSannotations_with_CCMS", "n_unique_CCMS", "n_CCMS_with1moretag")
  colnames(tag_numbers) <- "LifestyleTag"
  write.csv(tag_numbers, "number_LifestyleTag.csv", row.names = TRUE)
  
  zipr("zip_Lifestyle.zip", c("GNPS_hittable_LifestyleTag.csv",
                           "GNPS_LibraryHitTable_Summary_Lifestyle.pdf", 
                           "GNPS_LibraryHitTable_Summary_Lifestyle.png",
                           "GNPS_LibraryHitTable_SummaryPercentage_Lifestyle.pdf", 
                           "GNPS_LibraryHitTable_SummaryPercentage_Lifestyle.png",
                           "number_LifestyleTag.csv"))
  }
