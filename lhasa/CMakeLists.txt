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

# todo: use find_package() for rdkit
include_directories(${CMAKE_INSTALL_PREFIX}/include/rdkit)

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
    RDKitDepictor 
    RDKitcoordgen
    RDKitMolTransforms
    RDKitSmilesParse 
    RDKitMolDraw2D 
    RDKitGraphMol 
    RDKitRDInchiLib
    RDKitInchi
    RDKitRDGeneral 
    RDKitRDGeometryLib
    RDKitFileParsers 
    RDKitDescriptors
    RDKitChemTransforms
    RDKitSubstructMatch
    RDKitGenericGroups
    RDKitDataStructs
    boost_iostreams 
    boost_serialization
    ${DEPS_LIBRARIES}
    lh
)
