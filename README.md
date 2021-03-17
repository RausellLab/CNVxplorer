
# CNVxplorer: A webtool for the clinical interpretation of CNVs in rare disease patients

----------------------------------------

Welcome to the official Github repository of the **CNVxplorer** webserver presented at the BioRxiv preprint [CNVxplorer: a web tool to assist clinical interpretation of CNVs in rare disease patients](https://www.biorxiv.org/)


# Table of contents

  - [Overview](#Overview)
  - [Features](#Features)
  - [Access](#Access)
  - [Installation](#Installation)

## Overview

<p align="center">

<img  src="https://github.com/frequena/cnvxplorer/blob/master/doc/overview.jpg">

</p>

## Feature highlights

  - Query a genomic interval, cytoband, or upload a file with multiple
    genomic regions.
  - Explore CNV syndromes, pathogenic and non-pathogenic CNVs, disease
    variants (ClinVar, GWAS), and denovo variants (denovo.db) that
    overlap with the query CNV(s).
  - Explore protein-coding genes targetted by the query CNV(s), including
    disease annotations from five reference sources (OMIM, ORPHANET, 
    ClinGen, DECIPHER, and Genomics England PanelApp).
  - Find regulatory elements overlapping within the genomic intervals
    of the query CNV(s) as well as their target genes (within or outside
    the query CNV(s) genomic intervals.
  - Identify the genomic elements associated to the query CNV(s) with
    the closest phenotypic similarity to the patient's clinical signs.
  - Check the phenotype of KO experiments of the orthologous genes in mouse
    models.
  - Perform gene-set level analyses covering functional and pathway enrichment,
    tissue specific expression and gene/protein interaction network analyses.
  - Find PubMed articles related to the query CNV(s), analyse their keywork
    co-presence network, and filter articles by association
    with Mendelian diseases and genes.

### Access

You can find an instance of CNVxplorer running at this address:

<http://cnvxplorer.com>

## Installation

## Docker installation (3 lines)

``` bash

# Note: the first session after the deployment is slower since the application loads all the data required

git clone https://github.com/frequena/cnvxplorer.git

docker build -t cnvxplorer . # The tag "cnvxplorer" is optional

docker run -d -p 3838:3838 cnvxplorer # -p (specify port) -d (detached mode)

# The port 3838 is optional. Please make sure you set a port not blocked by firewalls.
# If you change the port number (3838) by any other, make sure to set it in the Dockerfile 
# (EXPOSE instruction)
```

### Local installation

``` bash

git clone https://github.com/frequena/cnvxplorer.git

cd cnvxplorer

gunzip -c local_data.RData.gz > local_data.RData

R -e 'load('local_data.RData')'

# Once you loaded the R environment with the files and functions
# comment the same two lines above in global.R so the app doesn't need
# to load the data every session but only once.

# Make sure you have all the packages installed
R -e 'shiny:runApp()'
```
