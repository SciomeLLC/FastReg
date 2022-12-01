FROM ubuntu

RUN apt-get update && apt-get install -y build-essential 


ENV DEBIAN_FRONTEND=noninteractive 
RUN apt-get -y install \
    gfortran \
    bzip2 \
    libcairo2-dev \
    libssl-dev \
    libgsl-dev \
    libicu-dev \
    libnlopt-dev \
    libpcre3 \
    libpcre3-dev \
    libpcre2-dev \
    libssl-dev \
    tk \
    pandoc \
    wget \
    make \
    perl \
    libhdf5-serial-dev \
    libopenblas-dev \
    libblas-dev \
    liblapack-dev \
    openmpi-bin \
    libopenmpi-dev \
    cmake \
    fakeroot \
    default-jdk

RUN apt-get -y install curl linux-tools-common libopenblas-openmp-dev libopenblas-serial-dev libopenblas64-openmp-dev libopenblas64-pthread-dev libopenblas64-serial-dev libblis-openmp-dev libblis-pthread-dev libblis-serial-dev libblis64-openmp-dev libblis64-pthread-dev libblis64-serial-dev libblas64-dev liblapack64-dev

# install flexiblas
RUN curl -O https://csc.mpi-magdeburg.mpg.de/mpcsc/software/flexiblas/flexiblas-3.2.1.tar.gz \
    && tar -xzvf flexiblas-3.2.1.tar.gz \
    && cd flexiblas-3.2.1 \
    && fakeroot dpkg-buildpackage -us -uc \
    && dpkg -i ../libflexiblas-*.deb

# && mkdir build \
# && cd build \
# && cmake ../ \
# && make \
# && make install

# RUN echo -e "\nLD_LIBRARY_PATH=/usr/local/lib64 \nexport LD_LIBRARY_PATH" >> ~/.bashrc && source ~/.bashrc

# # install latest version of R
# RUN curl -O https://cran.rstudio.com/src/base/R-4/R-4.2.1.tar.gz \
#     && tar -xzvf R-4.2.1.tar.gz \
#     && cd R-4.2.1 \
#     && ./configure \
#     --prefix=/opt/R/4.2.1 \
#     --enable-memory-profiling \
#     --enable-R-shlib \
#     --enable-BLAS-shlib \
#     --with-blas=flexiblas \
#     --with-lapack \
#     --with-x=no \
#     --with-readline=no \
#     && make \
#     && make install \
#     && ln -s /opt/R/4.2.1/bin/R /usr/local/bin/R \
#     && ln -s /opt/R/4.2.1/bin/Rscript /usr/local/bin/Rscript \
#     && cd /