#' Summary of UBERON Tags in GNPS - Tag Template
#'
#' This function works to generate a summary of the UBERON tags associated with annotated chemicals resulting from spectral library matching in GNPS. Users can choose to add new tags using the GNPS tag template. The GNPS tag template can be downloaded and used to produce the summary plots. UBERON tags are associated with planar InChIKeys.
#' @param tag_template Path to the downloaded GNPS tag template (Google Sheets).
#' @keywords
#' @export
#' @examples
#' @author Alan K. Jarmusch 2020
#' @export
#'

summary_UBERON_tagtemplate <- function(){

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

  df_tagtemplate <- fread(input = path, sep="\t", header=TRUE, stringsAsFactors = FALSE)
    studyinfo_n_row_annotations <- as.numeric(length(df_tagtemplate$`INPUT-#Scan#`)) #numbers out
  tags <- cbind(df_master$InChI_Key_Planar,separate(df_master[,"UBERONBodyPartName"], UBERONBodyPartName, paste0("X",1:max(sapply(strsplit(df_master$UBERONBodyPartName,"\\|"),length))), sep="\\|"))
  colnames(tags)[1] <- "InChI_Key_Planar"
  add_tags <- df_tagtemplate %>% select(`INPUT-InChIKey-Planar`,`Add_InCHIKey-Planar`,`Contribute-UBERONTag1`,`Contribute-UBERONTag2`,`Contribute-UBERONTag3`)
  add_tags$`INPUT-InChIKey-Planar` <- gsub("N/A", "",x = add_tags$`INPUT-InChIKey-Planar`)
  add_tags <- add_tags %>% mutate_all(na_if,"") %>% mutate(InChI_Key_Planar = coalesce(`INPUT-InChIKey-Planar`, `Add_InCHIKey-Planar`))
  add_tags <- add_tags[,-c(1:2)]
  study_InChI <- as.data.frame(unique(add_tags$InChI_Key_Planar))
  colnames(study_InChI)[1] <- "InChI_Key_Planar"
    studyinfo_n_rows_with_InChI <- as.data.frame(add_tags$InChI_Key_Planar) #numbers out
    studyinfo_n_rows_with_InChI <- subset(studyinfo_n_rows_with_InChI, studyinfo_n_rows_with_InChI$`add_tags$InChI_Key_Planar` != "NA") #numbers out
    studyinfo_n_unique_InChI <- unique(studyinfo_n_rows_with_InChI$`add_tags$InChI_Key_Planar`) #numbers out
    studyinfo_n_rows_with_InChI <- as.numeric(length(studyinfo_n_rows_with_InChI$`add_tags$InChI_Key_Planar`)) #numbers out
    studyinfo_n_unique_InChI <- as.numeric(length(studyinfo_n_unique_InChI)) #numbers out
  add_tags <- merge(tags, add_tags, by = "InChI_Key_Planar", all=TRUE)
  add_tags <- add_tags %>% unite(TAGS, 2:length(add_tags), sep = "|") %>% mutate(TAGS = gsub("(?<![a-zA-Z])NA\\||\\|NA(?![a-zA-Z])|\\|NA$", '', TAGS, perl = T))
  add_tags$TAGS <- sapply(strsplit(add_tags$TAGS, "\\|"), function(x) TAGS = paste(unique(x), collapse = "|"))
  add_tags <- subset(add_tags, add_tags$InChI_Key_Planar != "NA" & add_tags$TAGS != "NA" &  add_tags$TAGS != "")

  df_libraryhits_tags <- merge(study_InChI, add_tags, by="InChI_Key_Planar")
    studyinfo_n_unique_InChI_w_tag <- as.numeric(length(unique(df_libraryhits_tags$InChI_Key_Planar))) #numbers out
  
  if (nrow(df_libraryhits_tags) > 0) {
    df_libraryhits_tags <- df_libraryhits_tags %>% distinct()
    colnames(df_libraryhits_tags)[1] <- "InChIKey-Planar"
    df_libraryhits_tags <- as.data.table(df_libraryhits_tags)
    tags <- cbind(df_libraryhits_tags[,"InChIKey-Planar"],separate(df_libraryhits_tags[,"TAGS"], TAGS, paste0("X",1:max(sapply(strsplit(df_libraryhits_tags$TAGS,"\\|"),length))), sep="\\|"))
    tags <- gather(tags, key="tags", value="terms", 2:length(tags))
    tag_summary <- tags %>% group_by(terms) %>% tally
    tag_summary <- subset(tag_summary, tag_summary$terms != "NA" & tag_summary$terms != "")
    tag_summary$type <- "InChIKey-Planar"
  }  else {
    print("no tags")
    tag_summary <- data.frame()
  }
  
  tag_summary <- tag_summary %>% mutate(percentchemicalstagged = 100*(n / studyinfo_n_unique_InChI_w_tag))
    
  if (sum(df_master_CCMSLIB$UBERONBodyPartName != "") > 0) {
  tags <- cbind(df_master_CCMSLIB$CCMSLIB,separate(df_master_CCMSLIB[,"UBERONBodyPartName"], UBERONBodyPartName, paste0("X",1:max(sapply(strsplit(df_master_CCMSLIB$UBERONBodyPartName,"\\|"),length))), sep="\\|"))
  colnames(tags)[1] <- "CCMSLIB"
  } else {
    print("no tags")
  }

  df_CCMSLIB_tags <- subset(df_tagtemplate, df_tagtemplate$`INPUT-InChIKey-Planar` == "N/A" & df_tagtemplate$`Cannot_Add_InChIKey-Planar` == "No InChIKey-Planar")
    studyinfo_CCMS_n_rows <- as.numeric(length(df_CCMSLIB_tags$`INPUT-#Scan#`)) #numbers out
  add_tags <- df_CCMSLIB_tags %>% select(`INPUT-SpectrumID`,`Contribute-UBERONTag1`,`Contribute-UBERONTag2`,`Contribute-UBERONTag3`)
  study_SpectrumID <- as.data.frame(unique(add_tags$`INPUT-SpectrumID`))
  colnames(study_SpectrumID)[1] <- "CCMSLIB"
    studyinfo_CCMS_n_rows_with_uniqueCCMS <- as.numeric(length(study_SpectrumID$CCMSLIB)) #numbers out
  add_tags <- merge(tags, add_tags, by.x ="CCMSLIB", by.y= "INPUT-SpectrumID", all=TRUE)
  add_tags <- add_tags %>% mutate_all(na_if,"")
  add_tags <- add_tags %>% unite(TAGS, 2:length(add_tags), sep = "|") %>% mutate(TAGS = gsub("(?<![a-zA-Z])NA\\||\\|NA(?![a-zA-Z])|\\|NA$", '', TAGS, perl = T))
  add_tags$TAGS <- sapply(strsplit(add_tags$TAGS, "\\|"), function(x) TAGS = paste(unique(x), collapse = "|"))
  add_tags <- subset(add_tags, add_tags$CCMSLIB != "NA" & add_tags$TAGS != "NA" &  add_tags$TAGS != "")

  df_CCMSLIB_tags <- as.data.table(merge(study_SpectrumID, add_tags, by="CCMSLIB"))
    studyinfo_CCMS_n_unique_CCMS_w_tag <- as.numeric(length(unique(df_CCMSLIB_tags$CCMSLIB))) #numbers out
  
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
  write.csv(hit_table, "GNPS_hittable_UBERONTag.csv", row.names=FALSE)
  
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
  ggsave(file="GNPS_TagTemplate_Summary_UBERON.pdf", width = 3.25, height = 3.25, units = "in", dpi=300, useDingbats=FALSE)
  ggsave(file="GNPS_TagTemplate_Summary_UBERON.png", width = 3.25, height = 3.25, units = "in")
  
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
    labs(x="UBERON Tags", y="% of Unique Planar InChIKeys or CCMSLIB IDs")
  print(plot_summary)
  ggsave(file="GNPS_TagTemplate_SummaryPercentage_UBERON.pdf", width = 3.25, height = 3.25, units = "in", dpi=300, useDingbats=FALSE)
  ggsave(file="GNPS_TagTemplate_SummaryPercentage_UBERON.png", width = 3.25, height = 3.25, units = "in")
  
  #numbers
  tag_numbers <- rbind(studyinfo_n_row_annotations, studyinfo_n_rows_with_InChI, 
                       studyinfo_n_unique_InChI,studyinfo_n_unique_InChI_w_tag, studyinfo_CCMS_n_rows,
                       studyinfo_CCMS_n_rows_with_uniqueCCMS, studyinfo_CCMS_n_unique_CCMS_w_tag)
  rownames(tag_numbers) <- c("n_GNPSannotations", "n_GNPSannotations_with_InChIKeyPlanar", 
                             "n_unique_InChIKeyPlanar", "n_InChIKeyPlanar_with1moretag", 
                             "n_GNPSannotations_with_CCMS", "n_unique_CCMS", "n_CCMS_with1moretag")
  colnames(tag_numbers) <- "UBERONTag"
  write.csv(tag_numbers, "numbers_UBERONTag.csv", row.names = TRUE)
  
  zipr("UBERON_zip.zip", c("GNPS_hittable_UBERONTag.csv",
                         "GNPS_TagTemplate_Summary_UBERON.pdf", 
                         "GNPS_TagTemplate_Summary_UBERON.png",
                         "GNPS_TagTemplate_SummaryPercentage_UBERON.pdf", 
                         "GNPS_TagTemplate_SummaryPercentage_UBERON.png",
                         "numbers_UBERONTag.csv"))
}

