#!/bin/sh

if command -v greadlink > /dev/null 2>&1; then
    SOURCE_DIR=`dirname -- "$( greadlink -f -- "$0"; )"`
else
    SOURCE_DIR=`dirname -- "$( readlink -f -- "$0"; )"`
fi

BUILD_DIR=${PWD}/build
INSTALL_DIR=${PWD}/prefix

NUMPROCS=`nproc --all`

echo "Sources are in ${SOURCE_DIR}"
echo "Building in ${BUILD_DIR}"
echo "Installing in ${INSTALL_DIR}"

mkdir -p ${INSTALL_DIR}
mkdir -p ${BUILD_DIR}


#Boost (has to be built in source tree as far as I am aware)
cd ${SOURCE_DIR}/deps/boost
./bootstrap.sh --with-libraries=serialization,regex,chrono,date_time,filesystem,iostreams,program_options,thread,math,random,system &&\
./b2 toolset=emscripten link=static variant=release threading=single runtime-link=static thread system filesystem regex serialization chrono date_time program_options random -j ${NUMPROCS} &&\
./b2 toolset=emscripten link=static variant=release threading=single runtime-link=static install --prefix=${INSTALL_DIR}
cd ${BUILD_DIR}

emar q ${INSTALL_DIR}/lib/libboost_chrono.a ${INSTALL_DIR}/lib/libboost_chrono.bc
emar q ${INSTALL_DIR}/lib/libboost_date_time.a ${INSTALL_DIR}/lib/libboost_date_time.bc
emar q ${INSTALL_DIR}/lib/libboost_filesystem.a ${INSTALL_DIR}/lib/libboost_filesystem.bc
emar q ${INSTALL_DIR}/lib/libboost_iostreams.a ${INSTALL_DIR}/lib/libboost_iostreams.bc
emar q ${INSTALL_DIR}/lib/libboost_math_c99.a ${INSTALL_DIR}/lib/libboost_math_c99.bc
emar q ${INSTALL_DIR}/lib/libboost_math_c99f.a ${INSTALL_DIR}/lib/libboost_math_c99f.bc
emar q ${INSTALL_DIR}/lib/libboost_math_c99l.a ${INSTALL_DIR}/lib/libboost_math_c99l.bc
emar q ${INSTALL_DIR}/lib/libboost_math_tr1.a ${INSTALL_DIR}/lib/libboost_math_tr1.bc
emar q ${INSTALL_DIR}/lib/libboost_math_tr1f.a ${INSTALL_DIR}/lib/libboost_math_tr1f.bc
emar q ${INSTALL_DIR}/lib/libboost_math_tr1l.a ${INSTALL_DIR}/lib/libboost_math_tr1l.bc
emar q ${INSTALL_DIR}/lib/libboost_program_options.a ${INSTALL_DIR}/lib/libboost_program_options.bc
emar q ${INSTALL_DIR}/lib/libboost_random.a ${INSTALL_DIR}/lib/libboost_random.bc
emar q ${INSTALL_DIR}/lib/libboost_regex.a ${INSTALL_DIR}/lib/libboost_regex.bc
emar q ${INSTALL_DIR}/lib/libboost_serialization.a ${INSTALL_DIR}/lib/libboost_serialization.bc
emar q ${INSTALL_DIR}/lib/libboost_system.a ${INSTALL_DIR}/lib/libboost_system.bc
emar q ${INSTALL_DIR}/lib/libboost_wserialization.a ${INSTALL_DIR}/lib/libboost_wserialization.bc

#RDKit
mkdir -p ${BUILD_DIR}/rdkit_build
cd ${BUILD_DIR}/rdkit_build
emcmake cmake -Dboost_iostreams_DIR="${INSTALL_DIR}/lib/cmake/boost_iostreams-1.83.0"  \
              -Dboost_system_DIR="${INSTALL_DIR}/lib/cmake/boost_system-1.83.0"        \
              -Dboost_headers_DIR="${INSTALL_DIR}/lib/cmake/boost_headers-1.83.0"      \
              -Dboost_serialization_DIR="${INSTALL_DIR}/lib/cmake/boost_serialization-1.83.0"      \
              -DBoost_DIR="${INSTALL_DIR}/lib/cmake/Boost-1.83.0"                      \
              -DRDK_BUILD_PYTHON_WRAPPERS=OFF                                          \
              -DRDK_INSTALL_STATIC_LIBS=ON                                             \
              -DRDK_INSTALL_INTREE=OFF                                                 \
              -DRDK_BUILD_SLN_SUPPORT=OFF                                              \
              -DRDK_TEST_MMFF_COMPLIANCE=OFF                                           \
              -DRDK_BUILD_CPP_TESTS=OFF                                                \
              -DRDK_USE_BOOST_SERIALIZATION=ON                                         \
              -DRDK_BUILD_THREADSAFE_SSS=OFF                                           \
              -DBoost_INCLUDE_DIR="${INSTALL_DIR}/include"                             \
              -DBoost_USE_STATIC_LIBS=ON                                               \
              -DBoost_USE_STATIC_RUNTIME=ON                                            \
              -DBoost_DEBUG=TRUE                                                       \
              -DCMAKE_CXX_FLAGS="-s USE_PTHREADS=1 -pthread -Wno-enum-constexpr-conversion -D_HAS_AUTO_PTR_ETC=0" \
              -DCMAKE_INSTALL_PREFIX="${INSTALL_DIR}" "${SOURCE_DIR}/deps/rdkit"       \
              -DRDK_OPTIMIZE_POPCNT=OFF
emmake make -j ${NUMPROCS}
emmake make install
cd ${BUILD_DIR}


export CHOST="wasm32-unknown-linux"
export ax_cv_c_float_words_bigendian=no

export MESON_CROSS="${BUILD_DIR}/emscripten-crossfile.meson"

cat > "${BUILD_DIR}/emscripten-crossfile.meson" <<END
[binaries]
c = 'emcc'
cpp = 'em++'
ld = 'wasm-ld'
ar = 'emar'
ranlib = 'emranlib'
pkgconfig = ['emconfigure', 'pkg-config']

# https://docs.gtk.org/glib/cross-compiling.html#cross-properties
[properties]
growing_stack = true
have_c99_vsnprintf = true
have_c99_snprintf = true
have_unix98_printf = true

# Ensure that '-s PTHREAD_POOL_SIZE=*' is not injected into .pc files
[built-in options]
c_thread_count = 0
cpp_thread_count = 0

[host_machine]
system = 'emscripten'
cpu_family = 'wasm32'
cpu = 'wasm32'
endian = 'little'
END

export EM_PKG_CONFIG_PATH=${INSTALL_DIR}/lib/pkgconfig/
export PKG_CONFIG_PATH=${INSTALL_DIR}/lib/pkgconfig/
export EM_PKG_CONFIG_LIBDIR=${INSTALL_DIR}/lib/
export PKG_CONFIG_LIBDIR=${INSTALL_DIR}/lib/

# Graphene
cd ${SOURCE_DIR}/deps/graphene-1.10.8/
CFLAGS="-s USE_PTHREADS" LDFLAGS=" -lpthread" meson setup ${BUILD_DIR}/graphene_build \
    --prefix=${INSTALL_DIR} \
    --cross-file=$MESON_CROSS \
    --default-library=static \
    --buildtype=release \
    -Dtests=false && \
    meson install -C ${BUILD_DIR}/graphene_build

# Libsigc++
cd ${SOURCE_DIR}/deps/libsigcplusplus-3.6.0/
meson setup ${BUILD_DIR}/libsigcplusplus_build \
    --prefix=${INSTALL_DIR} \
    --libdir=lib \
    --cross-file=$MESON_CROSS \
    --default-library=static \
    -Dc_link_args='-pthread' \
    -Dcpp_args='-s USE_PTHREADS=1 -pthread' \
    --buildtype=release && \
    meson install -C ${BUILD_DIR}/libsigcplusplus_build
