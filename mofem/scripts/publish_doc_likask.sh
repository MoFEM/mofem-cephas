#!/bin/sh

RSYNC=/usr/bin/rsync
PROJECT_SOURCE_DIR=/Users/xiaoyi/mofem_installation/mofem-cephas/mofem
BINDIR=/Users/xiaoyi/mofem_installation/lib

$RSYNC -avz html "campus\lk58p"@userweb.eng.gla.ac.uk:/home/CAMPUS/lk58p/html/lukasz.kaczmarczyk/MoFem
