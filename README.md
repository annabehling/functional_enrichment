# functional_enrichment

Functional enrichment analyses can be used to identify an over-representation of biological functions in a subset of interesting genes, when compared with the entire background gene set.

In this instance, the subset of interesting genes is those with a particular expression pattern of interest.

## Description

The following code utilises two routes of gene ontology (GO) functional enrichment analysis - [DAVID bioinformatics resource](https://david.ncifcrf.gov/summary.jsp "DAVID") and [topGO](https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf "topGO").

DAVID - what it does (annotation, enrichment), what it takes as input (list, background(or not)), entry point to this code ...
topGO - what it does (enrichment), what it takes as input, entry point to this code ...

## Installation

To install the required scripts, first clone the **DEA_and_fit** repository.
```
git clone https://github.com/annabehling/functional_enrichment
```

## Usage - DAVID functional enrichment

To run the functions found in the file `DAVID_functions.R`, you will need a directory containing raw GO enrichment tables, downloaded from the [DAVID bioinformatics resource](https://david.ncifcrf.gov/).

Example files for DAVID GO term enrichment among Homeolog Expression Bias (HEBi) genes can be found in `files/DAVID_HEBi`.

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

## Usage - topGO functional enrichment analysis