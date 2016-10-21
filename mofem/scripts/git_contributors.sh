#!/bin/sh

GIT=/usr/bin/git
PROJECT_SOURCE_DIR=/Users/xiaoyi/mofem_installation/mofem-cephas/mofem

$GIT shortlog -s -e -n  | \
  sed 's/\@/ at /g' | \
  tee $PROJECT_SOURCE_DIR/doc/contributors_list
