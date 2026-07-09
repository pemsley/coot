cmake_minimum_required(VERSION 3.25)
project(layla CXX)

# Standalone build of Layla (Coot's ligand editor), translated from Makefile.am
# (and the hand-written Makefile.hgonomeg). It links against an already-installed
# Coot prefix for libcoot-geometry / -lidia-core / -utils, and uses system GTK4,
# RDKit, Boost, Python, libcurl and libpng.

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Install prefix that holds the built Coot libraries and their headers
# (headers under <prefix>/include/coot, e.g. utils/coot-utils.hh, geometry/..).
set(COOT_PREFIX "$ENV{HOME}/freshcoot/experimental_systemlibsbuild/prefix2"
    CACHE PATH "Prefix containing libcoot-geometry/-lidia-core/-utils and coot headers")

# --- Dependencies -----------------------------------------------------------
find_package(PkgConfig REQUIRED)
pkg_check_modules(GTK4 REQUIRED IMPORTED_TARGET gtk4)
pkg_check_modules(CURL REQUIRED IMPORTED_TARGET libcurl)
# mmdb2 is header-only for Layla (pulled in via the coot geometry/lidia-core
# headers); linking is carried by the coot .so's, so we only need its includes.
pkg_check_modules(MMDB2 IMPORTED_TARGET mmdb2)

find_package(PNG REQUIRED)          # used by similar_ligands.cpp
find_package(Eigen3 REQUIRED)       # pulled in by RDKit's MolTransforms.h
find_package(Threads REQUIRED)

# RDKit: system install, headers under /usr/include/rdkit. Link the components
# Layla uses, mirroring Makefile.hgonomeg (+ FMCS for similar-ligand search).
# They are shared libraries, so transitive RDKit deps resolve automatically.
find_path(RDKIT_INCLUDE_DIR GraphMol/RDKitBase.h PATH_SUFFIXES rdkit REQUIRED)
set(RDKIT_COMPONENTS
    Depictor MolTransforms SmilesParse MolDraw2D GraphMol RDInchiLib RDGeneral
    RDGeometryLib FileParsers Descriptors ChemTransforms SubstructMatch
    GenericGroups DataStructs FMCS RDChemDrawLib)
set(RDKIT_LIBRARIES "")
foreach(_c IN LISTS RDKIT_COMPONENTS)
  find_library(RDKIT_${_c}_LIB "RDKit${_c}" REQUIRED)
  list(APPEND RDKIT_LIBRARIES "${RDKIT_${_c}_LIB}")
endforeach()

find_library(BOOST_IOSTREAMS_LIB boost_iostreams REQUIRED)
find_library(BOOST_SERIALIZATION_LIB boost_serialization REQUIRED)

find_library(GEMMI_LIB gemmi_cpp)   # optional; normally pulled in via coot libs

# Coot libraries from the install prefix.
find_library(COOT_GEOMETRY_LIB   coot-geometry   PATHS "${COOT_PREFIX}/lib" REQUIRED)
find_library(COOT_LIDIA_CORE_LIB coot-lidia-core PATHS "${COOT_PREFIX}/lib" REQUIRED)
find_library(COOT_UTILS_LIB      coot-utils      PATHS "${COOT_PREFIX}/lib" REQUIRED)

# --- libcoot-layla ----------------------------------------------------------
add_library(coot-layla SHARED
    layla_embedded.cpp
    generators.cpp
    inchi_key_database.cpp
    notifier.cpp
    signals.cpp
    state.cpp
    ui.cpp
    utils.cpp
    file_io.cpp
    qed.cpp
    ligand_editor_canvas.cpp
    ligand_editor_canvas/core.cpp
    ligand_editor_canvas/model.cpp
    ligand_editor_canvas/tools.cpp
    ligand_editor_canvas/render.cpp
    similar_ligands.cpp
    network_utils.cpp)

target_include_directories(coot-layla PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${COOT_PREFIX}/include
    ${COOT_PREFIX}/include/coot
    ${RDKIT_INCLUDE_DIR})

target_compile_definitions(coot-layla PUBLIC
    USE_LIBCURL=1
    VERSION="standalone"
    PKGDATADIR="${CMAKE_INSTALL_PREFIX}/share/coot")

target_link_libraries(coot-layla PUBLIC
    ${COOT_GEOMETRY_LIB} ${COOT_LIDIA_CORE_LIB} ${COOT_UTILS_LIB}
    PkgConfig::GTK4 PkgConfig::CURL
    ${RDKIT_LIBRARIES}
    ${BOOST_IOSTREAMS_LIB} ${BOOST_SERIALIZATION_LIB}
    PNG::PNG
    Eigen3::Eigen
    Threads::Threads
    ${CMAKE_DL_LIBS})
if(TARGET PkgConfig::MMDB2)
  target_link_libraries(coot-layla PUBLIC PkgConfig::MMDB2)
endif()
if(GEMMI_LIB)
  target_link_libraries(coot-layla PUBLIC ${GEMMI_LIB})
endif()

# --- layla-bin --------------------------------------------------------------
add_executable(layla-bin main.cpp)
target_link_libraries(layla-bin PRIVATE coot-layla)
# Export symbols so the embedded Python interpreter can resolve them; bake an
# rpath to the coot prefix so its .so's are found at runtime.
target_link_options(layla-bin PRIVATE -Wl,--export-dynamic)
set_target_properties(layla-bin PROPERTIES
    BUILD_RPATH   "${COOT_PREFIX}/lib"
    INSTALL_RPATH "${COOT_PREFIX}/lib")

# --- launcher script (layla.in -> layla) ------------------------------------
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(exec_prefix "\${prefix}")
set(libexecdir "${CMAKE_INSTALL_PREFIX}/libexec")
configure_file(layla.in "${CMAKE_CURRENT_BINARY_DIR}/layla" @ONLY)

# --- install ----------------------------------------------------------------
install(TARGETS coot-layla LIBRARY DESTINATION lib)
install(TARGETS layla-bin  RUNTIME DESTINATION libexec)
install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/layla" DESTINATION bin)
install(FILES layla.ui DESTINATION share/coot/ui)
