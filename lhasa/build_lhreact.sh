#!/usr/bin/sh

./build_lhasa.sh
cd lhreact/
cp ../lhbuild/lhasa.{worker.js,js,wasm} ./public/
cp ../lhbuild/lhasa.d.ts ./src/
# nothing for now
# todo: distribute
cd ..