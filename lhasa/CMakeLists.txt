cmake_minimum_required(VERSION 3.25)
project(Lhasa)
set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wno-enum-constexpr-conversion -s USE_PTHREADS=1 -pthread")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
set(CMAKE_C_FLAGS_RELEASE "-O2")
set(CMAKE_CXX_STANDARD 17)


link_directories(${CMAKE_INSTALL_PREFIX}/lib)
include_directories(${CMAKE_INSTALL_PREFIX}/include)
include_directories(${CMAKE_INSTALL_PREFIX}/include/graphene-1.0)
include_directories(${CMAKE_INSTALL_PREFIX}/lib/graphene-1.0/include)
include_directories(${CMAKE_INSTALL_PREFIX}/include/rdkit)

if("${EMSCRIPTEN_VERSION}" VERSION_GREATER 3.1.26)
add_link_options(-sDEFAULT_PTHREAD_STACK_SIZE=2MB -sSTACK_SIZE=5MB -lembind)
endif()

add_library(lhasa STATIC lhasa.hpp lhasa.cpp)

target_sources(lhasa PRIVATE
  ../layla/utils.hpp
  ../layla/utils.cpp
  # ../layla/ligand_editor_canvas.hpp
  # ../layla/ligand_editor_canvas.cpp
  # ../layla/ligand_editor_canvas/core.hpp
  # ../layla/ligand_editor_canvas/core.cpp
  ../layla/ligand_editor_canvas/core.hpp
  # ../layla/ligand_editor_canvas/core.cpp
  # ../layla/ligand_editor_canvas/model.hpp
  # ../layla/ligand_editor_canvas/model.cpp
  ../layla/ligand_editor_canvas/tools.hpp
  ../layla/ligand_editor_canvas/tools.cpp
)

target_link_libraries(lhasa PRIVATE 
    RDKitDepictor 
    RDKitSmilesParse 
    RDKitMolDraw2D 
    RDKitGraphMol 
    RDKitRDGeneral 
    RDKitFileParsers 
    RDKitSubstructMatch 
    boost_iostreams 
    boost_serialization
    graphene-1.0
)

# target_link_libraries(lhasa PUBLIC ...)