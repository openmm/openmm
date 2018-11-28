FROM ubuntu:xenial

ENV PATH="/opt/miniconda/bin:${PATH}"

RUN apt-get update && \
    apt-get install -y wget git gromacs doxygen gfortran libfftw3-dev gcc g++ bzip2 automake make lsb-core && \
    wget https://repo.continuum.io/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p "/opt/miniconda" && \
    rm -f miniconda.sh && \
    mkdir /.conda &&  chmod -R 777 /.conda &&  chmod -R 777 /opt/miniconda && \
    conda install -y cmake numpy scipy pytest swig cython
