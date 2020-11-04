# functional_enrichment

Functional enrichment analyses can be used to identify an over-representation of biological functions in a subset of interesting genes, when compared with the entire background gene set.

In this instance, the subset of interesting genes are those that have certain expression patterns.

## Description

This code utilises two means of ... DAVID and topGO. both GO related, one does annotations, one needs annotations first

DAVID entry steps
-
-

topGO entry steps

## Installation

To install the required scripts, first clone the **DEA_and_fit** repository.
```
git clone https://github.com/annabehling/functional_enrichment
```

## Usage - DAVID

To run the functions found in the file `DAVID_functions.R`, you will need a directory containing raw GO enrichment tables, downloaded from the [DAVID bioinformatics resource](https://david.ncifcrf.gov/).

Example files for Homeolog Expression Bias (HEBi) enrichment can be found in ...

First load the functions:
```{r}
source("DAVID_functions.R")
```

Next, we will perform functiona

Example files have been provided

## Output - DAVID

## Usage - topGO

## Output - topGO

