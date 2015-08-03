#!/bin/sh

SFTP=/usr/bin/sftp
PROJECT_SOURCE_DIR=/mnt/home/MyBuild/mofem-cephas/mofem_v0.2

$SFTP -C -b $PROJECT_SOURCE_DIR/scripts/publish_doc "campus\lk58p"@userweb.eng.gla.ac.uk
