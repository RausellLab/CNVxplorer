# CNVxplorer: A webtool for the clinical interpretation of CNVs


## Overview

<p align="center">
  <img  src="https://github.com/frequena/cnvxplorer/blob/master/www/overview.jpg">
</p>

### Access

You can find an instance of CNVxplorer running at this address:

[http://compute.users-imagine.fr:3838/](http://compute.users-imagine.fr:3838/)



## Installation

### Local installation


```{bash}

git clone https://github.com/frequena/cnvxplorer.git

cd cnvxplorer

gunzip -c local_data.RData.gz > local_data.RData

R -e 'load('local_data.RData')'

# Once you loaded the R environment with the files and functions
# comment the same two lines above in global.R so the app doesn't need
# to load the data every session but only once.

R -e 'shiny:runApp()'

```

## Docker installation


```{bash}

git clone https://github.com/frequena/cnvxplorer.git

# The tag "cnvxplorer" is optional, you can choose your own tag
docker build -t cnvxplorer .

# The default port is 3838. If you specify a different port, you
# need to change it before in the Dockerfile (command: EXPOSE)
docker run -p 3838:3838 cnvxplorer

```

