FROM likask/ubuntu_mofem:v0.7

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

RUN adduser -q developer
RUN mkdir $BUILD_DIR
RUN chown developer $BUILD_DIR
USER developer
RUN mkdir $MOFEM_BUILD_DIR

RUN cd $MOFEM_BUILD_DIR \
&& cmake \
-DCMAKE_BUILD_TYPE=$MinSizeRel \
-DCMAKE_CXX_FLAGS="-Wall -std=c++11" \
-DPETSC_DIR=$PETSC_DIR \
-DPETSC_ARCH=$PETSC_ARCH \
-DMOAB_DIR=$PETSC_DIR/$PETSC_ARCH \
-DADOL-C_DIR=/usr \
-DWITH_TETGEN=1 \
-DBUILD_SHARED_LIBS=yes \
-DCMAKE_INSTALL_PREFIX=$MOFEM_INSTALL_DIR \
$MOFEM_SRC_DIR

RUN cd $MOFEM_BUILD_DIR && make install
RUN cd $MOFEM_BUILD_DIR && ctest --output-on-failure

RUN cd $MOFEM_INSTALL_DIR && cmake -DCMAKE_BUILD_TYPE=$MinSizeRel -DWITH_METAIO=1 -DCMAKE_CXX_FLAGS="-Wall -std=c++11" users_modules
RUN cd $MOFEM_INSTALL_DIR && make
RUN cd $MOFEM_INSTALL_DIR && ctest --output-on-failure
