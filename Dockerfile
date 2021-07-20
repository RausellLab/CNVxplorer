FROM rocker/shiny-verse:3.6.1

# install linux libraries
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
	
# install github packages
RUN R -e "devtools::install_github('thomasp85/patchwork')"
RUN R -e "devtools::install_github('glin/reactable')"
RUN R -e "devtools::install_github('rnabioco/valr')"
RUN R -e "devtools::install_github('frequena/bioloupe')"

# install CRAN packages
RUN R -e "install.packages(c('shinydashboard', \
'shinyjs', \
'dplyr', \
'highcharter', \
'ggraph', \
'ggrepel', \
'shinyhelper', \
'igraph', \
'widyr', \
'prettydoc', \
'ontologySimilarity', \
'ontologyIndex', \
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
'DescTools', \
'rentrez', \
'reactable', \
'ggridges', \
'UpSetR', \
'shinyalert', \
'networkD3', \
'XML', \
'shiny'))"


# install Bioconductor packages
RUN Rscript -e "BiocManager::install(c('karyoploteR','Gviz','TissueEnrich', 'Rhtslib', 'ReactomePA', 'clusterProfiler', 'org.Hs.eg.db', 'DOSE', 'chromPlot'))"

# Copy the app to the image
COPY CNVxplorer /srv/shiny-server/

# Updated configuration file 
RUN cat /srv/shiny-server/shiny-server.txt > /etc/shiny-server/shiny-server.conf

RUN echo 'allow_app_override;' >> /etc/shiny-server/shiny-server.conf

# Make all app files readable
RUN chmod -R +r /srv/shiny-server/

EXPOSE 3838

CMD ["/usr/bin/shiny-server.sh"] 
