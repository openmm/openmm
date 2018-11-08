FROM ubuntu:xenial

ENV PATH="/opt/miniconda/bin:${PATH}"

RUN sudo update-ca-certificates && \
    sudo apt-get update && \
    sudo apt-get install --no-install-recommends -y \
        wget git gromacs doxygen gfortran libfftw3-dev gcc g++ && \
    sudo wget https://repo.continuum.io/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh -O miniconda.sh && \
    sudo bash miniconda.sh -b -p "/opt/miniconda" && \
    sudo rm -f miniconda.sh && \
    sudo mkdir /.conda && sudo chmod -R 777 /.conda && sudo chmod -R 777 /opt/miniconda && \
    conda install -y cmake numpy scipy pytest swig
