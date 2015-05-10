#!/bin/sh

SFTP=/usr/bin/sftp
CMAKE_CURRENT_SOURCE_DIR=/Users/likask/MyBuild/mofem-cephas/mofem_v0.2/scripts

$SFTP -C -b $CMAKE_CURRENT_SOURCE_DIR/publish_doc "campus\lk58p"@userweb.eng.gla.ac.uk
