FROM rocker/rstudio:4.5.1

# Install ALL common system dependencies for R packages
RUN apt-get update && apt-get install -y \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libglpk-dev \
    libgmp-dev \
    libharfbuzz-dev \
    libicu-dev \
    libjpeg-dev \
    libpng-dev \
    libssl-dev \
    libtiff5-dev \
    libx11-dev \
    libxml2-dev \
    libhdf5-dev \
    libuv1-dev \
    patch \
    make \
    pandoc \
    zlib1g-dev \
    cmake \
    gdal-bin \ 
    libgdal-dev \
    libudunits2-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /home/rstudio/

CMD ["/init"]
