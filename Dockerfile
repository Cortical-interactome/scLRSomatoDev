# Base image https://hub.docker.com/u/rocker/
FROM rocker/shiny:4.3.2

# system libraries of general use
## install debian packages and update system libraries in a single layer
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
    libglu1-mesa \
    && apt-get upgrade -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# install renv & restore packages in a single layer
WORKDIR /project
COPY renv.lock renv.lock
ENV RENV_PATHS_LIBRARY=renv/library
RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))" && \
    R -e "renv::restore()" && \
    R -e "install.packages(c('markdown', 'rbibutils'), dependencies=TRUE, repos='http://cran.rstudio.com/')"

# copy necessary files into the /app directory
WORKDIR /app
COPY scLRSomatoDev/server.R .
COPY scLRSomatoDev/ui.R .
COPY scLRSomatoDev/utils.R .
COPY scLRSomatoDev/www ./www

# expose port
EXPOSE 3838

# run app on container start with fixed syntax
CMD ["R", "-e", "shiny::runApp('/app', host = '0.0.0.0', port = 3838)"]
