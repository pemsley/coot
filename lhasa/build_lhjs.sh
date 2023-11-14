#!/usr/bin/sh

./build_lhasa.sh
cd lhjs/
cp ../lhbuild/lhasa.{worker.js,js,wasm} ./
# nothing for now
cd ..