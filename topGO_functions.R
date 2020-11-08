## load libraries

library(topGO)
library(stringr)
library(ggplot2)

## load functions

#get the list of interesting genes from the classes_df dataframe that are in the HEBl, HEBi and HER categories
get_int_genes <- function(classes_df, gene_names){
  #classes_df : dataframe with one gene ID per row and expression category classifications for each gene
  #gene_names : character vector of all genes IDs with a non-NA expression category classification
  HEBl_geneList <- factor(as.integer(classes_df$classification == "HEBl")) #get HEBl genes
  names(HEBl_geneList) <- gene_names
  
  HEBi_geneList <- factor(as.integer(classes_df$classification == "HEBi")) #get HEBi genes
  names(HEBi_geneList) <- gene_names
  
  HER_geneList <- factor(as.integer(classes_df$classification == "HER")) #get HER genes
  names(HER_geneList) <- gene_names
  
  list(HEBl_geneList, HEBi_geneList, HER_geneList)
}

#make the topGO data object
make_topGO_DO <- function(int_gene_list, ontology, gene2GO_list){
  #int_gene_list : named list of interesting genes
  #ontology : two letter GO ontology code (BP / CC / MF)
  #gene2GO_list : list containing all genes with mapped GO IDs
  topGO_data <- new("topGOdata", ontology = ontology, allGenes = int_gene_list,
                    annot = annFUN.gene2GO, gene2GO = gene2GO_list)
  fishers_result <- runTest(topGO_data, statistic = "fisher")
  fishers_table <- GenTable(topGO_data, Fishers = fishers_result, useLevels = TRUE, topNodes = 100) #get 100 most significant terms
  fishers_table$Ontology <- ontology
  fishers_table
}

#run enrichment analyses for the BP, CC, MF ontologies and filter the results
all_enriched <- function(int_genes, all_golist){
  #int_genes : indexed output from get_int_genes() where [[1]] = HEBl, [[2]] = HEBi, [[3]] = HER
  #all_golist : gene2GO list containing all genes (that have a non-NA expression classification), with mapped GO IDs 
  BP <- make_topGO_DO(int_genes, "BP", all_golist) #BP topGO data object
  CC <- make_topGO_DO(int_genes, "CC", all_golist) #CC topGO data object
  MF <- make_topGO_DO(int_genes, "MF", all_golist) #MF topGO data object
  
  all_tab <- rbind(BP, CC, MF) #bind ontology tables together
  all_tab <- all_tab[, c(1:3, 8, 4, 5, 7)] #get specific columns and reorder
  colnames(all_tab) <- c("GO ID", "GO term", "level", "ontology", "background count", "subset count", "p value") #rename columns
  str_replace(all_tab$`p value`, "< 1e-30", "1e-31") #for reordering, temporarily remove "<", which for topGO appears at "< 1e-30"
  all_tab <- all_tab[order(as.numeric(all_tab$`p value`)), ] #order rows by p value, as.numeric() needed as some have scientific notation
  str_replace(all_tab$`p value`, "1e-31", "< 1e-30") #re-add in "<"
  all_tab[as.numeric(all_tab$`p value`) < 0.05, ] #print out only results where p < 0.05
}

#plot enrichment results for a particular expression classification
topGO_scatter <- function(allo_f_topgo_out, hh_f_topgo_out, allo_p_topgo_out, hh_p_topgo_out, allo_a_topgo_out, hh_a_topgo_out){
  #allo_f_topgo_out : allopolyploid fungi output of all_enriched()
  #hh_f_topgo_out : homoploid hybrid fungi output of all_enriched()
  #allo_p_topgo_out : allopolyploid plants output of all_enriched()
  #hh_p_topgo_out : homoploid hybrid plants ouput of all_enriched()
  #allo_a_topgo_out : allopolyploid animals output of all_enriched()
  #hh_a_topgo_out : homoploid hybrid animals output of all_enriched()
  if(nrow(allo_f_topgo_out) > 0) allo_f_topgo_out$system <- "allo fungi"
  if(nrow(hh_f_topgo_out) > 0) hh_f_topgo_out$system <- "HH fungi"
  if(nrow(allo_p_topgo_out) > 0) allo_p_topgo_out$system <- "allo plants"
  if(nrow(hh_p_topgo_out) > 0) hh_p_topgo_out$system <- "HH plants"
  if(nrow(allo_a_topgo_out) > 0) allo_a_topgo_out$system <- "allo animals"
  if(nrow(hh_a_topgo_out) > 0) hh_a_topgo_out$system <- "HH animals"
  
  scatter_tab <- rbind(allo_f_topgo_out, hh_f_topgo_out,
                       allo_p_topgo_out, hh_p_topgo_out,
                       allo_a_topgo_out, hh_a_topgo_out)
  
  scatter_tab$neglog10_p <- -log10(as.numeric(as.character(scatter_tab$`p value`))) #important to do as.character for as.numeric
  
  ggplot(data = scatter_tab,
         mapping = aes(y=`GO term`,
                       x=factor(system, levels=c("allo fungi", "HH fungi", "allo plants", "HH plants", "allo animals", "HH animals")),
                       color=`neglog10_p`,
                       size=`subset count`)) +
    scale_x_discrete(drop=F) +
    geom_point() +
    scale_color_gradient(low="blue", high="red") +
    labs(x="Biological system", y="GO term", size="Annotated", color="-log10(p value)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), #formats background
          panel.background = element_rect(colour = "black", size=0.5), #adds a plot border
          axis.title = element_text(face="bold")) #formats axis titles
}