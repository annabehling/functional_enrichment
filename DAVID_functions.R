## load libraries

library(stringr)
library(stringi)
library(dplyr)
library(ggplot2)

## load functions

#filter a raw DAVID GO enrichment table
filter_david_go <- function(david_out){
  #david_out : raw DAVID GO enrichment table
  raw_tab <- read.table(david_out, header = TRUE, sep = "\t", quote = '"', stringsAsFactors = FALSE) #read in DAVID output
  red_tab <- raw_tab[,c(2,8,3,5)] #get only columns we want, and order them
  red_tab$GO_ID <- sapply(str_split(red_tab[,1], "~"), "[", 1) #split first column on the "~", and make first part GO_ID column
  red_tab$GO_term <- sapply(str_split(red_tab[,1], "~"), "[", 2) #split first column on the "~" and make second part GO term column
  fin_tab <- red_tab[,c(5,6,2,3,4)] #reorder columns
  colnames(fin_tab) <- c("GO ID", "GO term", "background count", "subset count", "p value") #rename columns
  filt_tab <- fin_tab[fin_tab[, 5] < 0.05, ] #only print rows where p value < 0.05
  filt_tab[, 5] <- format(filt_tab[, 5], scientific = TRUE) #display p values in scientific notation
  filt_tab
}

#filter all raw DAVID GO enrichment tables in a directory
read_david_go <- function(results_dir){
  #results_dir : directory containing raw DAVID GO enrichment tables
  fnames <- list.files(results_dir, "goterm", full.names=TRUE) #list all files in directory containing "goterm" character string
  lapply(fnames, filter_david_go) #apply the filtering function to to all GO enrichment files
}

#plot enrichment results for a given GO ontology and level
GO_scatter <- function(allo_f_david_out, hh_f_david_out, allo_p_david_out, hh_p_david_out, allo_a_david_out, hh_a_david_out, go_cat, y_label){
  #allo_f_david_out : allopolyploid fungi output of read_david_go()
  #hh_f_david_out : homoploid hybrid fungi output of read_david_go()
  #allo_p_david_out : allopolyploid plants output of read_david_go()
  #hh_p_david_out : homoploid hybrid plants output of read_david_go()
  #allo_a_david_out : allopolyploid animals output of read_david_go()
  #hh_a_david_out : homoploid hybrid animals output of read_david_go()
  #go_cat : number corresponding to GO ontology and level (1-5 : BP level 1-5, 6-10 : CC level 1-5. 11-15 : MF level 1-5)
  #y_label : y axis label
  if(nrow(allo_f_david_out[[go_cat]]) > 0) allo_f_david_out[[go_cat]]$system <- "allo fungi"
  if(nrow(hh_f_david_out[[go_cat]]) > 0) hh_f_david_out[[go_cat]]$system <- "HH fungi"
  if(nrow(allo_p_david_out[[go_cat]]) > 0) allo_p_david_out[[go_cat]]$system <- "allo plants"
  if(nrow(hh_p_david_out[[go_cat]]) > 0) hh_p_david_out[[go_cat]]$system <- "HH plants"
  if(nrow(allo_a_david_out[[go_cat]]) > 0) allo_a_david_out[[go_cat]]$system <- "allo animals"
  if(nrow(hh_a_david_out[[go_cat]]) > 0) hh_a_david_out[[go_cat]]$system <- "HH animals"
  
  scatter_tab <- rbind(allo_f_david_out[[go_cat]], hh_f_david_out[[go_cat]],
                       allo_p_david_out[[go_cat]], hh_p_david_out[[go_cat]],
                       allo_a_david_out[[go_cat]], hh_a_david_out[[go_cat]])
  
  scatter_tab$neglog10_p <- -log10(as.numeric(as.character(scatter_tab$`p value`))) #important to do as.character for as.numeric
  
  ggplot(data = scatter_tab,
         mapping = aes(y=`GO term`,
                       x=factor(system, levels=c("allo fungi", "HH fungi", "allo plants", "HH plants", "allo animals", "HH animals")),
                       color=`neglog10_p`,
                       size=`subset count`)) +
    scale_x_discrete(drop = FALSE) +
    geom_point() +
    scale_color_gradient(low="blue", high="red") +
    labs(x="Biological system", y=y_label, size="Annotated", color="-log10(p value)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), #formats background
          panel.background = element_rect(colour = "black", size=0.5), #adds a plot border
          axis.title = element_text(face="bold"), #formats axis titles
          axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) #formats angle and position of x axis labels
}

#find common GO terms across all representative species 
common_terms <- function(df1, df2, df3, df4, df5, df6){
  #df1 - df6 : output of read_david_go()
  join_1 <- inner_join(df1, df2, by="GO term")
  join_2 <- inner_join(join_1, df3, by="GO term")
  join_3 <- inner_join(join_2, df4, by="GO term")
  join_4 <- inner_join(join_3, df5, by="GO term")
  join_5 <- inner_join(join_4, df6, by="GO term")
  join_5
}