#!/bin/sh

SFTP=/usr/bin/sftp
PROJECT_SOURCE_DIR=/Users/xiaoyi/mofem_installation/mofem-cephas/mofem

$SFTP -C -b $PROJECT_SOURCE_DIR/scripts/publish_doc "campus\lk58p"@userweb.eng.gla.ac.uk
