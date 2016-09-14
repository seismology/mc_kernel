#!/bin/bash
TRAVIS_ROOT="$1"

if [ ! -f "$TRAVIS_ROOT/bin/git-lfs" ]; then
  wget https://github.com/github/git-lfs/releases/download/v1.2.1/git-lfs-linux-amd64-1.2.1.tar.gz
  tar xvfz git-lfs-linux-amd64-1.2.1.tar.gz 
  mv git-lfs-1.2.1/git-lfs $TRAVIS_ROOT/bin/
fi
git lfs install
