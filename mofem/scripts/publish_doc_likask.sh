#!/bin/sh

RSYNC=/usr/bin/rsync
PROJECT_SOURCE_DIR=/Users/karollewandowski/moFEM/mofem-cephas/mofem
BINDIR=/Users/karollewandowski/moFEM/lib

$RSYNC -avz --delete html "campus\lk58p"@userweb.eng.gla.ac.uk:/home/CAMPUS/lk58p/html/lukasz.kaczmarczyk/MoFem
