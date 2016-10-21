#!/bin/sh

RSYNC=/usr/bin/rsync
PROJECT_SOURCE_DIR=/Users/xiaoyi/mofem_installation/mofem-cephas/mofem
BINDIR=/Users/xiaoyi/mofem_installation/lib

$RSYNC -avz html lukasz@cdash.eng.gla.ac.uk:/var/www/html/mofem
