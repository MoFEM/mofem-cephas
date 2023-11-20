#!/bin/bash

# Determine number of cores to compile code
NBCORES=$(cat /proc/cpuinfo | grep processor | wc -l)
NB_MAX=12
if (($NBCORES > $NB_MAX)); then
  NB=$NB_MAX
else
  NB=$NBCORES
fi
echo Nb. of cores $NBCORES and nb. of cores used to compilation $NB

set -e

echo "Configure users modules"
cd $MOFEM_INSTALL_DIR

if [ -e $MOFEM_SRC_DIR/users_modules/*bone_remodelling* ]
then
/opt/local/bin/cmake \
-DCMAKE_BUILD_TYPE=Release \
-DBUILD_SHARED_LIBS=yes \
-DCMAKE_CXX_FLAGS="-Wall -Wno-sign-compare" \
-DMPI_RUN_FLAGS="--allow-run-as-root" \
-DWITH_METAIO=1 \
users_modules;
else
/opt/local/bin/cmake \
-DCMAKE_BUILD_TYPE=Release \
-DBUILD_SHARED_LIBS=yes \
-DCMAKE_CXX_FLAGS="-Wall -Wno-sign-compare" \
-DMPI_RUN_FLAGS="--allow-run-as-root" \
users_modules
fi

echo "Run basic tests"
cd $MOFEM_INSTALL_DIR
make -k -j $NB; /bin/true
ctest --output-on-failure -D Experimental
make clean
echo "All done"
