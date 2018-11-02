FROM nvidia/cuda:9.0-devel

RUN apt-get update && apt-get -y install wget git gromacs doxygen gfortran libfftw3-dev && \
    wget https://repo.continuum.io/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p "/opt/miniconda" && \
    rm -f miniconda.sh && \
    mkdir /.conda && chmod -R 777 /.conda && chmod -R 777 /opt/miniconda

ENV PATH="/opt/miniconda/bin:${PATH}"