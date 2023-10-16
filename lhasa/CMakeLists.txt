cmake_minimum_required(VERSION 3.13)
project(Lhasa)
set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
   
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -Wno-enum-constexpr-conversion -I${CMAKE_INSTALL_PREFIX}/include -s USE_PTHREADS=1 -pthread")
set(CMAKE_C_FLAGS_RELEASE "-O2")
set(CMAKE_CXX_STANDARD 17)
  


link_directories(${CMAKE_INSTALL_PREFIX}/lib)
include_directories(${CMAKE_INSTALL_PREFIX}/include/rdkit)

if("${EMSCRIPTEN_VERSION}" VERSION_GREATER 3.1.26)
add_link_options(-sDEFAULT_PTHREAD_STACK_SIZE=2MB -sSTACK_SIZE=5MB)
endif()

set(LHASA_SOURCES lhasa.cpp)

add_library(lhasa STATIC ${LHASA_SOURCES})

# target_link_libraries(lhasa PUBLIC ...)