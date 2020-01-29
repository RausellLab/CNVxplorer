FROM rocker/shiny-verse:3.6.1

RUN apt-get update && apt-get install libcurl4-openssl-dev libv8-3.14-dev  libbz2-dev zlib1g-dev libncurses5-dev libncursesw5-dev liblzma-dev  -y &&\
  mkdir -p /var/lib/shiny-server/bookmarks/shiny
	

RUN R -e "devtools::install_github('thomasp85/patchwork')"
RUN R -e "devtools::install_github('glin/reactable')"

RUN R -e "install.packages('prettydoc')"

RUN R -e "install.packages(c('shinydashboard', \
'shinyjs', \
'formattable' \
'hrbrthemes' \
'import', \
'BiocManager', \
'tablerDash', \
'shinyEffects', \
'echarts4r', \
'shinyWidgets', \
'DT', \
'shinymanager', \
'gghighlight', \
'ReactomePA', \
'shinycssloaders', \
'plotly', \
'waiter', \
'DescTools', \
'clusterProfiler', \
'org.Hs.eg.db', \
'DOSE', \
'enrichplot', \
'rentrez', \
'reactable', \
'ggridges', \
'UpSetR', \
'randomForest', \
'chromPlot', \
'Gviz', \
'XML', \
'shiny'))"

RUN Rscript -e "BiocManager::install('rols')"

RUN Rscript -e "BiocManager::install(c('karyoploteR', 'Rhtslib', 'ReactomePA', 'clusterProfiler', 'org.Hs.eg.db', 'DOSE', 'enrichplot','rentrez', 'chromPlot', 'Gviz'))"



# Copy the app to the image
COPY cnvxplore /srv/shiny-server/


# Make all app files readable
RUN chmod -R +r /srv/shiny-server/

EXPOSE 3838

CMD ["/usr/bin/shiny-server.sh"] 
