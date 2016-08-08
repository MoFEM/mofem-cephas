FROM ubuntu:precise

MAINTAINER Lukasz.Kaczmarczyk@glasgow.ac.uk

ADD mofem /mofem

# Set environment
ENV BUILD_DIR /build
ENV MOFEM_SRC_DIR /mofem/
ENV MOFEM_BUILD_DIR $BUILD_DIR/lib
ENV MOFEM_INSTALL_DIR $BUILD_DIR/um

ENV PETSC_VERSION 3.7.2
ENV PETSC_DIR=/opt/petsc
ENV PETSC_ARCH=arch-linux2-c-opt

ENV BUILD_TYPE=MinSizeRel
ENV CMAKE /opt/bin/cmake

RUN apt-get update \
&& apt-get install -y \
build-essential \
gcc \
g++ \
libopenmpi-dev \
openmpi-bin \
openssh-server \
libblas-dev \
liblapack-dev \
&& rm -rf /var/lib/apt/lists/*

RUN wget https://www.dropbox.com/s/axa48op5mevukfj/mofem_env.tar.gz \
&& tar -xzf mofem_env.tar.gz \
&& rm -f mofem_env.tar.gz \

RUN echo 'root:mofem' | chpasswd \
&& adduser -q developer \
&& mkdir $BUILD_DIR \
&& chown developer $BUILD_DIR

USER developer

RUN mkdir $MOFEM_BUILD_DIR \
&& cd $MOFEM_BUILD_DIR \
&& $CMAKE \
-DCMAKE_BUILD_TYPE=MinSizeRel \
-DBUILD_SHARED_LIBS=1 \
-DCMAKE_CXX_FLAGS="-Wall" \
-DPETSC_DIR=$PETSC_DIR \
-DPETSC_ARCH=$PETSC_ARCH \
-DMOAB_DIR=$PETSC_DIR/$PETSC_ARCH \
-DWITH_ADOL-C=1 \
-DCMAKE_INSTALL_PREFIX=$MOFEM_INSTALL_DIR \
$MOFEM_SRC_DIR \
&& cd $MOFEM_BUILD_DIR \
&& make install \
&& ctest --output-on-failure

RUN cd $MOFEM_INSTALL_DIR  \
&& $CMAKE -DCMAKE_CXX_FLAGS="-Wall" users_modules \
&& cd $MOFEM_INSTALL_DIR/basic_finite_elements/atom_tests \
&& make \
&& ctest --output-on-failure
