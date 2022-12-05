#!/bin/sh

RSYNC=/usr/bin/rsync
PROJECT_SOURCE_DIR=/Users/karollewandowski/moFEM/mofem-cephas/mofem
BINDIR=/Users/karollewandowski/moFEM/lib

$RSYNC -avz --delete html lukasz@cdash.eng.gla.ac.uk:/var/www/html/mofem
