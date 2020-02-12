FROM rocker/shiny-verse:3.6.1

RUN apt-get update && apt-get install libcurl4-openssl-dev \ 
libv8-3.14-dev \ 
libbz2-dev \
zlib1g-dev \
libbz2-1.0 \
libbz2-ocaml \
libbz2-ocaml-dev \
libjpeg-dev \
libncurses5-dev \ 
libncursesw5-dev \
liblzma-dev  -y &&\
mkdir -p /var/lib/shiny-server/bookmarks/shiny
	


RUN R -e "devtools::install_github('thomasp85/patchwork')"
RUN R -e "devtools::install_github('glin/reactable')"


RUN R -e "install.packages(c('shinydashboard', \
'shinyjs', \
'prettydoc' \
'ontologySimilarity' \
'ontologyIndex' \
'formattable', \
'hrbrthemes', \
'import', \
'BiocManager', \
'tablerDash', \
'shinyEffects', \
'echarts4r', \
'shinyWidgets', \
'DT', \
'shinymanager', \
'gghighlight', \
'shinycssloaders', \
'plotly', \
'waiter', \
'DescTools', \
'rentrez', \
'reactable', \
'ggridges', \
'UpSetR', \
'randomForest', \
'XML', \
'shiny'))"


RUN Rscript -e "BiocManager::install('rols')"
RUN Rscript -e "BiocManager::install('ReactomePA')"
RUN Rscript -e "BiocManager::install('Gviz')"

RUN Rscript -e "BiocManager::install(c('karyoploteR', 'Rhtslib', 'ReactomePA', 'clusterProfiler', 'org.Hs.eg.db', 'DOSE', 'enrichplo', 'chromPlot'))"



# Copy the app to the image
COPY cnvxplorer /srv/shiny-server/


# Make all app files readable
RUN chmod -R +r /srv/shiny-server/

EXPOSE 3838

CMD ["/usr/bin/shiny-server.sh"] 
