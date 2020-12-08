
# CNVxplorer: A webtool for the clinical interpretation of CNVs

# Table of contents

  - [Overview](#Overview)
  - [Features](#Features)
  - [Access](#Access)
  - [Installation](#Installation)
  - [Contributions](#Contributions)

## Overview

<p align="center">

<img  src="https://github.com/frequena/cnvxplorer/blob/master/doc/overview.jpg">

</p>

## Features

  - Query a genomic interval, cytoband, or upload a file with multiple
    genomic regions and analyze it simultaneously.
  - Compare the genetic evidence of disease genes from five different
    databases (OMIM, ORPHANET, ClinGen, DECIPHER, and Genomics England
    PanelApp).
  - Explore CNV syndromes, pathogenic and non-pathogenic CNVs, disease
    variants (ClinVar, GWAS), and denovo variants (denovo.db) that
    overlap with the CNVs.
  - Find regulatory elements disrupted by the variant and their target
    genes.
  - Consider those target genes that do not map the input variant(s) in
    the subsequent analysis (optional).
  - Identify relevant diseases and genes based on phenotypic similarity
    with the clinical characteristics of the patient.
  - Check the phenotype of the orthologous genes in experiments with
    mouse models.
  - Find PubMed articles related to your variants and filter the
    relevant ones by keyword, number of citations, or by association
    with the OMIM database.

### Access

You can find an instance of CNVxplorer running at this address:

<http://compute.users-imagine.fr:3838/>

## Installation

### Local installation

``` bash

git clone https://github.com/frequena/cnvxplorer.git

cd cnvxplorer

gunzip -c local_data.RData.gz > local_data.RData

R -e 'load('local_data.RData')'

# Once you loaded the R environment with the files and functions
# comment the same two lines above in global.R so the app doesn't need
# to load the data every session but only once.

R -e 'shiny:runApp()'
```

## Docker installation (3 lines)

``` bash

git clone https://github.com/frequena/cnvxplorer.git

# The tag "cnvxplorer" is optional, you can choose your own tag
docker build -t cnvxplorer .

# -p Specify the port. The default port is 3838. If you select a different port, you
# need to change it before in the Dockerfile (command: EXPOSE)
# Please, make sure the port is listening and not blocked by firewalls
# -d Detached mode
docker run -d -p 3838:3838 cnvxplorer

# Note: the first session after the deployment is slower since the application loads all the data required
```
