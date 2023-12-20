#!/usr/bin/sh

INSTALL_DIR=${PWD}/prefix
emcmake cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -S . -B lhbuild &&\
cd lhbuild/;\
emmake make LDFLAGS=-all-static -j16;\
cd ..