# functional_enrichment

Functional enrichment analyses can be used to identify an over-representation of biological functions in a subset of interesting genes, when compared with the entire background gene set.

In this instance, we are investigating the evolution of hybrid gene expression, so the subset of interesting genes are those with a particular expression pattern, relative to their parental species.

Briefly, the four possible categories for hybrid gene expression are:

1. **Parental expression inheritance** (PEI): any parental bias (or lack thereof) is maintained in the hybrid.
2. **Homeolog expression blending** (HEBl): a parental expression bias is lost in the hybrid.
3. **Homeolog expression bias** (HEBi): a hybrid expression bias has arisen from no parental bias.
4. **Homeolog expression reversal** (HER): an expression bias in the parents is reversed in the hybrid.

More information on these classes can be found in the following publication:

Cox, M.P., T. Dong, G. Shen, Y. Dalvi, D.B. Scott and A.R.D. Ganley. 2014. An interspecific fungal hybrid reveals cross-kingdom rules for allopolyploid gene expression patterns. *PLoS Genetics* 10: e1004180. [https://doi.org/10.1371/journal.pgen.1004180](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004180)

## Description

The following code utilises two routes of gene ontology (GO) functional enrichment analysis - [DAVID bioinformatics resource](https://david.ncifcrf.gov/summary.jsp "DAVID") and [topGO](https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf "topGO").

**DAVID** in an online tool that can perform functional annotation and functional enrichment analyses. It requires the user to submit a gene list with particular sequence identifiers, and gives the user the possibility of submitting a background gene list, or using the species background provided by the DAVID knowledgebase. Species must be present in the DAVID knowledgebase for usage of the tool.

Results of the enrichment analyses can then be downloaded as `.txt` files.  
DAVID enables GO enrichment analyses at specific GO levels (1-5). 

**topGO** is an R package that can perform GO term enrichment analysis on functionally-annotated genes. Unlike in DAVID, the topGO GO enrichment analyses are level-independent.

## Installation

To install the required scripts, first clone the **functional_enrichment** repository.
```
git clone https://github.com/annabehling/functional_enrichment
```

## Processing of raw DAVID functional enrichment analysis results

To run the functions found in the file `DAVID_functions.R`, you will need a directory containing raw GO enrichment tables, downloaded from the DAVID bioinformatics resource.

Example files for DAVID GO term enrichment among homeolog expression bias (HEBi) genes can be found in `files/DAVID_HEBi`.

First load the functions:
```{r}
source("DAVID_functions.R")
```

To read in and filter all of the example enrichment tables, applying a *p* value threshold of 0.05, run:
```{r}
allo_f_hebi_david <- read_david_go("./files/DAVID_HEBi") #allopolyploid fungi
hh_f_hebi_david <- read_david_go("./files/DAVID_HEBi") #homoploid hybrid fungi
allo_p_hebi_david <- read_david_go("./files/DAVID_HEBi") #allopolyploid plants
hh_p_hebi_david <- read_david_go("./files/DAVID_HEBi") #homoploid hybrid plants
allo_a_hebi_david <- read_david_go("./files/DAVID_HEBi") #allopolyploid animals
hh_a_hebi_david <- read_david_go("./files/DAVID_HEBi") #homoploid hybrid animals
```

We then can check that all data files are present; we expect there to be 15 tables x 6 systems = 90 data tables total.  
Even if some are empty, it is important that the total number of files are present so that the `go_cat` argument of the following plotting argument is accurate.
```{r}
length(c(allo_f_hebi_david, hh_f_hebi_david, allo_p_hebi_david, hh_p_hebi_david, allo_a_hebi_david, hh_a_hebi_david)) #90
```

The plotting code also assumes that the raw enrichment file names meet a number of criteria as in the example file names, in order for the plotting function to be accurate.
Specifically, file names must contain:
* the character string **goterm**
* a reference to the GO ontology, either abbreviated (e.g. **bp**), or in full (e.g. **biological process**) - as a primary means of indexing the file names
* a reference to the GO level (e.g. **1 - 5**) - as a secondary means of indexing the file names

To plot the cross-kingdom results of the functional enrichment analysis at a given GO ontology level, run:
```{r}
GO_scatter(allo_f_hebi_david, hh_f_hebi_david, allo_p_hebi_david, hh_p_hebi_david, allo_a_hebi_david, hh_a_hebi_david, go_cat = 2, 
           y_label = "GO term (Biological Process level 2)")
```

In this example code, we have plotted from biological process level 2.  
The output should of these functions should match the following plot, which can also be found as an example output file in `files/DAVID_BP_2.png`.

![Image of DAVID output plot](files/DAVID_BP_2.png "DAVID output plot")  
The plot shows no common enriched GO terms at biological process level 2.

Alternatively, to see if there are any common GO terms across any GO ontology level, run:
```{r}
mapply(FUN = common_terms, allo_f_hebi_david, hh_f_hebi_david, allo_p_hebi_david, hh_p_hebi_david, allo_a_hebi_david, hh_a_hebi_david, SIMPLIFY = FALSE)
```

## topGO functional enrichment analysis and processing

To run the functions found in the file `topGO_functions.R`, you will need a list of background genes with GO annotations, in **gene2GO** format. This means that every gene is present in the list, with zero, one or many GO IDs annotated to each gene. 

You will also need a dataframe with at least two columns; one containing gene IDs (that match the gene2GO list IDs) and one containing the expression category classification for each gene.  
More information on how this dataframe can be made from raw read count data can be found at [this repository](https://github.com/annabehling/DEA_and_fit "github.com/annabehling/DEA_and_fit").