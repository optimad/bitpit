FROM ubuntu:20.04 AS bitpit-ubuntu-base
LABEL stage=bitpit-ubuntu-base
ENV DEBIAN_FRONTEND noninteractive
RUN sed -i "s/htt[p|ps]:\/\/archive.ubuntu.com\/ubuntu\//mirror:\/\/mirrors.ubuntu.com\/mirrors.txt/g" /etc/apt/sources.list && \
    apt-get update && \
    apt-get install -y \
        build-essential \
        cmake \
        clang \
        libblas-dev  \
        liblapack-dev  \
        liblapacke-dev \
        mpich \
        libxml2-dev \
        libmetis-dev \
        git \
        wget \
        python3 \
        rapidjson-dev

FROM bitpit-ubuntu-base AS bitpit-ubuntu-mpi
LABEL stage=bitpit-ubuntu-mpi
RUN cd /opt && \
    wget https://www.mpich.org/static/downloads/4.0.3/mpich-4.0.3.tar.gz && \
    tar xzf mpich-4.0.3.tar.gz && \
    cd mpich-4.0.3 && \
    ./configure --prefix=/usr/local && \
    make -j10 && \
    make install && \
    cd .. && \
    rm -rf mpich-4.0.3*

FROM bitpit-ubuntu-base AS bitpit-ubuntu-boost
LABEL stage=bitpit-ubuntu-boost
RUN cd /opt && \
    wget --no-verbose https://boostorg.jfrog.io/artifactory/main/release/1.71.0/source/boost_1_71_0.tar.gz && \
    tar xzf boost_1_71_0.tar.gz && \
    cd boost_1_71_0 && \
    ./bootstrap.sh --prefix=/opt/boost-install && ./b2 -j 10 install && \
    cd .. && \
    rm -rf boost_1_71_0

FROM bitpit-ubuntu-mpi AS bitpit-ubuntu-petsc
LABEL stage=bitpit-ubuntu-petsc
RUN cd /opt && \
    git clone --depth 1 --single-branch -b release https://gitlab.com/petsc/petsc.git petsc-code && \
    cd petsc-code && \
    for MPI in "OFF" "ON"; do \
        if [ "${MPI}" = "ON" ]; then \
            PETSC_DIR="/opt/petsc-install"; \
            PETSC_MPI_OPTIONS="--with-mpi=1 --with-mpi-dir=/usr/local/"; \
        else \
            PETSC_DIR="/opt/petsc-serial-install"; \
            PETSC_MPI_OPTIONS="--with-mpi=0"; \
        fi && \
        ./configure \
            ${PETSC_MPI_OPTIONS} \
            --prefix=${PETSC_DIR} \
            --with-debug=0 \
            --with-shared-libraries=1 \
            --with-scalar-type=real && \
            COPTFLAGS='-O3 -march=native -mtune=native' \
            CXXOPTFLAGS='-O3 -march=native -mtune=native' \
            FOPTFLAGS='-O3 -march=native -mtune=native' \
        make -j10 && \
        make install; \
    done && \
    cd .. && \
    rm -rf petsc-code

FROM bitpit-ubuntu-mpi AS bitpit-environment
LABEL stage=bitpit-environment
COPY --from=bitpit-ubuntu-boost /opt/boost-install/ /usr/local/
COPY --from=bitpit-ubuntu-petsc /opt/petsc-install /opt/petsc-install
COPY --from=bitpit-ubuntu-petsc /opt/petsc-serial-install /opt/petsc-serial-install

FROM bitpit-environment AS bitpit-build
LABEL stage=bitpit-build
COPY . /opt/bitpit
ARG MPI
ARG DEBUG
ARG COMPILER
RUN cd /opt/bitpit && \
    mkdir build && \
    cd build && \
    if [ "$DEBUG" = "ON" ]; then \
      BUILD_TYPE="Debug"; \
    else \
      BUILD_TYPE="Release"; \
    fi && \
    if [ "${MPI}" = "ON" ]; then \
        PETSC_DIR="/opt/petsc-install"; \
    else \
        PETSC_DIR="/opt/petsc-serial-install"; \
    fi && \
    if [ "${COMPILER}" = "clang" ]; then \
        C_COMPILER="clang"; \
        CXX_COMPILER="clang++"; \
        Fortran_COMPILER="gfortran"; \
    else \
        C_COMPILER="gcc"; \
        CXX_COMPILER="g++"; \
        Fortran_COMPILER="gfortran"; \
    fi && \
    cmake \
        -DBITPIT_BUILD_EXAMPLES=ON \
        -DBITPIT_BUILD_SHARED_LIBS=ON \
        -DBITPIT_ENABLE_INTEGRATION_TESTS=ON \
        -DBITPIT_ENABLE_MPI=${MPI} \
        -DBITPIT_ENABLE_UNIT_TESTS=ON \
        -DBUILD_TESTING=ON \
        -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
        -DCMAKE_C_COMPILER=${C_COMPILER} \
        -DCMAKE_CXX_COMPILER=${CXX_COMPILER} \
        -DCMAKE_Fortran_COMPILER=${Fortran_COMPILER} \
        -DPETSC_DIR=${PETSC_DIR} \
        ../ && \
    make -j10

FROM bitpit-build AS bitpit-test
LABEL stage=bitpit-test
RUN cd /opt/bitpit/build && \
    ctest --output-on-failure
