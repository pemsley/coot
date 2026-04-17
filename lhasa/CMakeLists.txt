cmake_minimum_required(VERSION 3.25)
project(Lhasa)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()


# -Wno-unused-command-line-argument 

if(MEMORY64 EQUAL "1")
  set(CMAKE_CXX_FLAGS "-Wall -sMEMORY64=1 --emit-tsd lhasa.d.ts -s USE_PTHREADS=1 -pthread -fwasm-exceptions")
else()
  set(CMAKE_CXX_FLAGS "-Wall --emit-tsd lhasa.d.ts -s USE_PTHREADS=1 -pthread -fwasm-exceptions")
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
set(CMAKE_C_FLAGS_RELEASE "-O2")
set(CMAKE_CXX_STANDARD 23)



include_directories(${CMAKE_INSTALL_PREFIX}/include)

# Needed so that find_package() searches the prefix where dependencies are installed.
# Stored in the cache so it propagates into find_package() calls inside dependency config files.
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_INSTALL_PREFIX})

set(RDKit_DIR ${CMAKE_INSTALL_PREFIX}/lib/cmake/rdkit)
set(Eigen3_DIR ${CMAKE_INSTALL_PREFIX}/share/eigen3/cmake)

# Resolving Boost is a fucking shit-show

# Pre-find Boost in config mode; required because rdkit-config.cmake calls find_dependency(Boost)
# and CMake 4.x removed the legacy FindBoost module (CMP0167)
# Collect all versioned Boost cmake config dirs (e.g. lib/cmake/Boost-1.89.0) into a list
file(GLOB _boost_cmake_dirs "${CMAKE_INSTALL_PREFIX}/lib/cmake/Boost-*")
# Pick the first (and normally only) entry as Boost_DIR so find_package(Boost CONFIG) can locate BoostConfig.cmake
list(GET _boost_cmake_dirs 0 Boost_DIR)
# BoostConfig.cmake internally calls find_package(boost_<comp>) for each component, using HINTS to
# locate versioned per-component dirs (e.g. boost_headers-1.89.0). CMake 4.x no longer resolves
# these from HINTS alone, so we pre-populate each boost_<comp>_DIR from the installed layout.
file(GLOB _boost_comp_dirs "${CMAKE_INSTALL_PREFIX}/lib/cmake/boost_*")
foreach(_dir IN LISTS _boost_comp_dirs)
    # Extract the directory name (e.g. "boost_headers-1.89.0") from the full path
    get_filename_component(_comp_ver "${_dir}" NAME)
    # Strip the version suffix to get the cmake package name (e.g. "boost_headers")
    string(REGEX REPLACE "-[0-9]+\\..*" "" _comp "${_comp_ver}")
    # Set boost_<comp>_DIR so find_package(boost_<comp>) resolves without searching
    set(${_comp}_DIR "${_dir}" CACHE PATH "" FORCE)
endforeach()

find_package(Boost CONFIG REQUIRED COMPONENTS iostreams serialization)
find_package(RDKit REQUIRED)
find_package(Eigen3 REQUIRED)

set(ENV{EM_PKG_CONFIG_PATH} ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig/)
set(ENV{PKG_CONFIG_PATH} ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig/)
set(ENV{EM_PKG_CONFIG_LIBDIR} ${CMAKE_INSTALL_PREFIX}/lib/)
set(ENV{PKG_CONFIG_LIBDIR} ${CMAKE_INSTALL_PREFIX}/lib/)

if("${EMSCRIPTEN_VERSION}" VERSION_GREATER 3.1.26)
  add_link_options(-sDEFAULT_PTHREAD_STACK_SIZE=2MB -sSTACK_SIZE=5MB)
endif()

if(MEMORY64 EQUAL "1")
  add_link_options(-sMEMORY64=1)
endif()

add_link_options(-lembind -fwasm-exceptions)

include(FindPkgConfig)
pkg_check_modules(DEPS REQUIRED sigc++-3.0 graphene-1.0)

link_directories(${CMAKE_INSTALL_PREFIX}/lib)
include_directories(${DEPS_INCLUDE_DIRS})

add_library(lh STATIC lhasa.hpp lhasa.cpp)
target_link_libraries(lh PRIVATE RDKit::rdkit_base Eigen3::Eigen)


target_sources(lh PRIVATE
  glog_replacement.hpp
  glog_replacement.cpp
  ../layla/qed.hpp
  ../layla/qed.cpp
  ../layla/utils.hpp
  ../layla/utils.cpp
  ../layla/ligand_editor_canvas.hpp
  ../layla/ligand_editor_canvas.cpp
  ../layla/ligand_editor_canvas/core.hpp
  ../layla/ligand_editor_canvas/core.cpp
  ../layla/ligand_editor_canvas/model.hpp
  ../layla/ligand_editor_canvas/model.cpp
  ../layla/ligand_editor_canvas/tools.hpp
  ../layla/ligand_editor_canvas/tools.cpp
  ../layla/ligand_editor_canvas/render.hpp
  ../layla/ligand_editor_canvas/render.cpp
  ../utils/base64-encode-decode.hh
  ../utils/base64-encode-decode.cc
)

add_executable(lhasa embind.cpp embind.hpp)

set_target_properties(lhasa PROPERTIES LINK_FLAGS " -s ENVIRONMENT=web,worker -s MODULARIZE=1 -s EXPORT_NAME='Lhasa' -s ALLOW_MEMORY_GROWTH=1 -s ASSERTIONS")
#set_target_properties(lhasa PROPERTIES LINK_FLAGS " -s ENVIRONMENT=web,worker -s EXPORT_ES6=1 -s MODULARIZE=1 -s EXPORT_NAME='Lhasa' -s ALLOW_MEMORY_GROWTH=1 -s ASSERTIONS")

target_link_libraries(lhasa PRIVATE
    RDKit::Depictor
    RDKit::coordgen
    RDKit::MolTransforms
    RDKit::SmilesParse
    RDKit::MolDraw2D
    RDKit::GraphMol
    RDKit::RDInchiLib
    RDKit::Inchi
    RDKit::RDGeneral
    RDKit::RDGeometryLib
    RDKit::FileParsers
    RDKit::RDChemDrawLib
    RDKit::ChemDraw
    RDKit::CIPLabeler
    RDKit::expat
    RDKit::Descriptors
    RDKit::ChemTransforms
    RDKit::SubstructMatch
    RDKit::GenericGroups
    RDKit::DataStructs
    Boost::iostreams
    Boost::serialization
    gemmi_cpp
    RDKit::maeparser
    ${DEPS_LIBRARIES}
    lh
)
