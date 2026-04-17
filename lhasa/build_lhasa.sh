#!/usr/bin/sh

# This is the final build step: Produce WASM+JS Lhasa module, linking against dependencies.

if [ "x$LHASA_MAIN_DIR" = "x" ]; then
    # Based on '$0', i.e. the command being executed, finds the absolute path of where this script is located
    if command -v greadlink > /dev/null 2>&1; then
        LHASA_MAIN_DIR=`dirname -- "$( greadlink -f -- "$0"; )"`
    else
        LHASA_MAIN_DIR=`dirname -- "$( readlink -f -- "$0"; )"`
    fi
    echo "Using LHASA_MAIN_DIR: $LHASA_MAIN_DIR"
else
    echo "Using LHASA_MAIN_DIR from environment: $LHASA_MAIN_DIR"
fi

if [ -e "$LHASA_MAIN_DIR/lhasa_build_functions" ]; then
    . "$LHASA_MAIN_DIR/lhasa_build_functions"
else
    # The fail function is not available here
    echo "Cannot find lhasa_build_functions. Exiting."
    exit 1
fi

setcolor cyan
echo "Running Emscripten CMake to create build dir at $LHASA_CMAKE_BUILD_DIR..."
setcolor reset

emcmake cmake -DMEMORY64=${MEMORY64} -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -S "$LHASA_MAIN_DIR" -B "$LHASA_CMAKE_BUILD_DIR" ||\
    fail "Failed to configure the build with EmscriptenCMake."

setcolor cyan
echo "Building Lhasa with Emscripten Make..."
setcolor reset

cd "$LHASA_CMAKE_BUILD_DIR" &&\
emmake make LDFLAGS=-all-static -j${NUMPROCS} ||\
    fail "Failed to build Lhasa with Emscripten Make."


cd "$LHASA_MAIN_DIR"

setcolor cyan
echo "Done building Lhasa."
setcolor reset
