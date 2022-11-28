FROM centos:centos7

# Setup EPEL
RUN yum -y install epel-release && yum -y update

# set locales
RUN echo "LANG=en_US.utf8" >> /etc/locale.conf \
    && localedef -c -f UTF-8 -i en_US en_US.UTF-8 \
    && export LC_ALL=en_US.UTF-8

RUN yum makecache && yum -y install \
    which \
    libgfortran \
    gcc-gfortran \
    gcc \
    gcc-c++ \
    libatlas-base-dev \
    bzip2-devel \
    cairo \
    libcurl4-openssl-dev \
    libcurl-devel \
    libgsl-dev \
    libicu-dev \
    xz-devel \
    nlopt-devel \
    nlopt \
    pango \
    pcre2-devel \
    pcre-devel \
    libssl-dev \
    libtcl8.6 \
    libtiff-devel \
    tk-8.6.8 \
    libxml2-devel \
    locales \
    tzdata \
    zlib-devel \
    wget \
    make \
    pandoc \
    openssl-devel \
    java \
    perl \
    hdf5-devel \
    # BLAS libraries
    openblas-devel \
    lapack-devel \
    scalapack-common \
    openmpi-devel \
    openmpi \
    scalapack-openmpi-devel \
    scalapack-openmpi \
    scalapack-openmpi-static 


# install latest version of cmake
RUN wget https://github.com/Kitware/CMake/releases/download/v3.23.4/cmake-3.23.4.tar.gz \
    && tar -zxvf cmake-3.23.4.tar.gz \
    && cd cmake-3.23.4/ \
    && ./bootstrap \
    && make \
    && make install \
    && cd /

# install flexiblas
RUN curl -O https://csc.mpi-magdeburg.mpg.de/mpcsc/software/flexiblas/flexiblas-3.2.1.tar.gz \
    && tar -xzvf flexiblas-3.2.1.tar.gz \
    && cd flexiblas-3.2.1 \
    && mkdir build \
    && cd build \
    && cmake ../ \
    && make \
    && make install

RUN echo -e "\nLD_LIBRARY_PATH=/usr/local/lib64 \nexport LD_LIBRARY_PATH" >> ~/.bashrc && source ~/.bashrc

# install latest version of R
RUN curl -O https://cran.rstudio.com/src/base/R-4/R-4.2.1.tar.gz \
    && tar -xzvf R-4.2.1.tar.gz \
    && cd R-4.2.1 \
    && ./configure \
    --prefix=/opt/R/4.2.1 \
    --enable-memory-profiling \
    --enable-R-shlib \
    --enable-BLAS-shlib \
    --with-blas=flexiblas \
    --with-lapack \
    --with-x=no \
    --with-readline=no \
    && make \
    && make install \
    && ln -s /opt/R/4.2.1/bin/R /usr/local/bin/R \
    && ln -s /opt/R/4.2.1/bin/Rscript /usr/local/bin/Rscript \
    && cd /

RUN echo "options(repos=c('http://cran.us.r-project.org'))" >> ~/.Rprofile
RUN Rscript -e "install.packages('flexiblas')"

RUN mkdir R && mkdir data && mkdir tests
COPY ./tests /tests
COPY ./data /data
COPY ./R /R
CMD ["bash"]