#!/usr/bin/sh
source ./EMSCRIPTEN_CONFIG
INSTALL_DIR=${PWD}/prefix
emcmake cmake -DMEMORY64=${MEMORY64} -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -S . -B lhbuild &&\
cd lhbuild/;\
emmake make LDFLAGS=-all-static -j16;\
cd ..