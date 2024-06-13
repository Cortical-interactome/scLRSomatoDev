# Base image https://hub.docker.com/u/rocker/
FROM rocker/shiny:4.3.2

# system libraries of general use
## install debian packages
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    libxml2-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libpq-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libhdf5-dev \
    libmagick++-dev \
    libglpk-dev \
    patch \
    libudunits2-dev \
    libproj-dev \
    libglu1-mesa

## update system libraries
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean

# install renv & restore packages
RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"
WORKDIR /project
COPY renv.lock renv.lock
ENV RENV_PATHS_LIBRARY renv/library
RUN R -e "renv::restore()"
RUN R -e "install.packages('markdown', dependencies=TRUE, repos='http://cran.rstudio.com/')"

# copy necessary files
## app folder
COPY /scLRSomatoDev ./app

# expose port
EXPOSE 3838

# run app on container start
CMD ["R", "-e", "shiny::runApp('/app', host = '0.0.0.0', port = 3838)"]
