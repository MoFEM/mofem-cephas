#!/bin/sh

GIT=/usr/bin/git
PROJECT_SOURCE_DIR=/mnt/home/MyBuild/mofem-cephas/mofem_v0.2

$GIT shortlog -s -e -n  | \
  sed 's/\@/ at /g' | \
  tee $PROJECT_SOURCE_DIR/doc/contributors_list
