
# CNVxplorer: a web tool to assist clinical interpretation of CNVs in rare disease patients

-----

Welcome to the official Github repository of the **CNVxplorer**
webserver available at <http://cnvxplorer.com>

The CNVxplorer manuscript is available on:

Requena, Francisco, et al. [“CNVxplorer: a web tool to assist clinical
interpretation of CNVs in rare disease patients.” Nucleic acids
research:
gkab347](https://academic.oup.com/nar/article/49/W1/W93/6279834).

# Table of contents

  - [Overview](#Overview)
  - [Features highlight](#features-highlight)
  - [Availability](#Availability)
  - [Local installation](#local-installation)
  - [Authors and contact](#authors-and-contact)
  - [License](#License)
  - [Disclaimer](#Disclaimer)
  - [References](#References)
  - [News](#News)

## Overview

CNVxplorer, a web server suited for the functional assessment of CNVs in
a clinical diagnostic setting. CNVxplorer mines a comprehensive set of
clinical, genomic, and epigenomic features associated with CNVs. It
provides sequence constraint metrics, impact on regulatory elements and
topologically associating domains, as well as expression patterns.
Analyses offered cover (a) agreement with patient phenotypes; (b)
visualizations of associations among genes, including those mediated by
regulatory elements, and transcription factors; (c) enrichment on
functional and pathway annotations; and (d) co-occurrence of terms
across PubMed publications related to the query CNVs. A flexible
evaluation workflow allows dynamic re-interrogation in clinical
sessions

<p align="center">

<img  src="https://github.com/RausellLab/CNVxplorer/blob/master/doc/Overview.svg">

</p>

## Features highlight

  - Query a genomic interval, cytoband, or upload a file with multiple
    genomic regions.
  - Explore CNV syndromes, pathogenic and non-pathogenic CNVs, disease
    variants (ClinVar, GWAS), and denovo variants (denovo.db) that
    overlap with the query CNV(s).
  - Explore protein-coding genes targetted by the query CNV(s),
    including disease annotations from five reference sources (OMIM,
    ORPHANET, ClinGen, DECIPHER, and Genomics England PanelApp).
  - Find regulatory elements overlapping within the genomic intervals of
    the query CNV(s) as well as their target genes (within or outside
    the query CNV(s) genomic intervals.
  - Identify the genomic elements associated to the query CNV(s) with
    the closest phenotypic similarity to the patient’s clinical signs.
  - Check the phenotype of KO experiments of the orthologous genes in
    mouse models.
  - Perform gene-set level analyses covering functional and pathway
    enrichment, tissue specific expression and gene/protein interaction
    network analyses.
  - Find PubMed articles related to the query CNV(s), analyse their
    keywork co-presence network, and filter articles by association with
    Mendelian diseases and genes.

## Availability

CNVxplorer is publicly available at <http://cnvxplorer.com>

Detailed tutorials, comprehensive documentation, and a Frequently Asked
Questions section are provided.

In addition, a stand-alone open-source R implementation with a shiny
interface is offered in this repository, allowing its deployment as a
private server through a Docker image without external dependencies.
Instructions to locally deploy the application are provided in the next
section.

## Local installation

### Docker installation

``` bash

# Note: the first session after the deployment is slower since the application loads all the data required

git clone https://github.com/RausellLab/CNVxplorer.git

mv CNVxplorer/Dockerfile .

docker build -t cnvxplorer . # The tag "cnvxplorer" is optional

docker run -d -p 3838:3838 cnvxplorer # -p (specify port) -d (detached mode)

# The port 3838 is optional. Please make sure you set a port not blocked by firewalls.
# If you change the port number (3838) by any other, make sure to set it in the Dockerfile 
# (EXPOSE instruction)
```

### Installation of CNVxplorer as a local private server

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

## Authors and contact

CNVxplorer has been developed by Francisco Requena and Antonio Rausell,
at the [Clinical Bioinformatics
Laboratory](https://www.institutimagine.org/en/antonio-rausell-161) of
the [Imagine Institute](https://www.institutimagine.org/en/) in Paris,
France.

CNVxplorer is the result of a close collaboration with Prof. Valérie
Malan and Prof. Serge Romana, from the [Cytogenetics Service of the
Necker Hospital for Sick
Children](http://hopital-necker.aphp.fr/histologie) (APHP) and the
[Bioinformatics
Platform](https://www.institutimagine.org/en/patrick-nitschke-201) of
the Imagine Institute headed by Patrick Nitschké.

Please address comments and questions about CNVxplorer to: \*
**Francisco Requena** -
[francisco.requena@institutimagine.org](francisco.requena@institutimagine.org)
\* **Antonio Rausell** -
[antonio.rausell@institutimagine.org](antonio.rausell@institutimagine.org)

## License

This project is licensed under the GNU General Public License 3 - see
the [LICENSE](LICENSE) file for details

See the License for the specific language governing permissions and
limitations under the License.

Copyright 2021 Clinical BioInformatics Laboratory - Institut Imagine

## Disclaimer

CNVxplorer or any document available from this server are distributed on
an “AS IS” BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
express, implied, or statutory, including, but not limited to, any
implied warranties of merchantability, fitness for a particular purpose
and freedom from infringement, or that CNVxplorer or any documents
available from this server will be error-free.

In no event will the Imagine Institute, the Clinical Bioinformatics lab,
or any of its members be liable for any damages, including but not
limited to direct, indirect, special, or consequential damages, arising
out of, resulting from, or in any way connected with the use of
CNVxplorer or documents available from it.

## References

Requena, Francisco, et al. “CNVxplorer: a web tool to assist clinical
interpretation of CNVs in rare disease patients.” Nucleic acids
research: gkab347.
[(Access)](https://academic.oup.com/nar/article/49/W1/W93/6279834)

## News

You may follow us in Twitter for regular news and updates:
<https://twitter.com/AntonioRausell>
