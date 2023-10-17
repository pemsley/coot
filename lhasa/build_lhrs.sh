#!/usr/bin/sh

./build_lhasa.sh
cd lhrs/
cp ../lhbuild/lhasa.{worker.js,js,wasm} ./
cargo make build_release
cd ..