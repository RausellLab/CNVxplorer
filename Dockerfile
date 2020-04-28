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
'highcharter', \
'arulesCBA', \
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
'plotly', \
'waiter', \
'DescTools', \
'rentrez', \
'reactable', \
'ggridges', \
'UpSetR', \
'randomForest', \
'shinyalert', \
'XML', \
'shiny'))"


RUN R -e "devtools::install_github('rnabioco/valr')"



RUN Rscript -e "BiocManager::install('ReactomePA')"
RUN Rscript -e "BiocManager::install('Gviz')"
RUN Rscript -e "BiocManager::install('TissueEnrich')"
RUN Rscript -e "BiocManager::install(c('karyoploteR', 'Rhtslib', 'ReactomePA', 'clusterProfiler', 'org.Hs.eg.db', 'DOSE', 'chromPlot'))"



# Copy the app to the image
COPY cnvxplorer /srv/shiny-server/

# Unzip local data

RUN gunzip /srv/shiny-server/local_data.RData.gz 

RUN echo 'allow_app_override;' >> /etc/shiny-server/shiny-server.conf


# Make all app files readable
RUN chmod -R +r /srv/shiny-server/

EXPOSE 3838

CMD ["/usr/bin/shiny-server.sh"] 
