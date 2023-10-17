#!/usr/bin/sh

INSTALL_DIR=${PWD}/prefix
emcmake cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -S . -B lhbuild &&\
cd lhbuild/;\
emmake make LDFLAGS=-all-static -j16;\
cd ..
emcc -lembind -O2 -o liblhasa.js lhbuild/liblhasa.a
cp liblhasa.{js,wasm} lhrs/
cd lhrs/
cargo make build
cd ..