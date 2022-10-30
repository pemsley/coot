project(coot)
cmake_minimum_required(VERSION 3.13)

set(CMAKE_CXX_STANDARD 11)

# set(CMAKE_EXECUTABLE_SUFFIX ".js")

set(coot_src ${CMAKE_CURRENT_SOURCE_DIR}/../coot)
set(dependency_install_dir $ENV{HOME}/autobuild/$ENV{OS}-$ENV{HOST}-gtk3)

add_definitions(-DANALYSE_REFINEMENT_TIMING -DVERSION="1.0.06-pre" -DPACKAGE_VERSION="1" -DHAVE_CXX_THREAD -DHAVE_BOOST -DHAVE_GSL -DPKGDATADIR="${coot_src}/pixmaps" -DDATADIR="${coot_src}/data")

#  pyrogen, compat, auxil, goograph, src ignored.

add_library(cootapi SHARED
            ${coot_src}/analysis/daca.cc
            ${coot_src}/analysis/kolmogorov.cc
            ${coot_src}/analysis/stats.cc
            ${coot_src}/analysis/bfkurt.cc
            ${coot_src}/analysis/mogul.cc
            ${coot_src}/analysis/b-factor-histogram.cc
            ${coot_src}/analysis/typed-distances.cc
            ${coot_src}/analysis/cablam.cc
            ${coot_src}/density-contour/density-contour-triangles.cc
            ${coot_src}/density-contour/Vectors.cpp
            ${coot_src}/density-contour/CIsoSurface.cpp
            ${coot_src}/density-contour/occlusion.cc
            ${coot_src}/ccp4mg-utils/ppmutil.cc
            ${coot_src}/ccp4mg-utils/quat.cc
            ${coot_src}/ccp4mg-utils/mgtree.cc
            ${coot_src}/ccp4mg-utils/cartesian.cc
            ${coot_src}/ccp4mg-utils/geomutil.cc
            ${coot_src}/ccp4mg-utils/matrix.cc
            ${coot_src}/skeleton/dm-utils.cc
            ${coot_src}/skeleton/BuildCas.cc
            ${coot_src}/skeleton/graphical_skel.cc
            ${coot_src}/ideal/pepflip.cc
            ${coot_src}/ideal/parallel-planes.cc
            ${coot_src}/ideal/simple-restraint.cc
            ${coot_src}/ideal/ng.cc
            ${coot_src}/ideal/torsion-bonds.cc
            ${coot_src}/ideal/chi-squareds.cc
            ${coot_src}/ideal/trans-peptide.cc
            ${coot_src}/ideal/model-bond-deltas.cc
            ${coot_src}/ideal/make-restraints.cc
            ${coot_src}/ideal/pull-restraint.cc
            ${coot_src}/ideal/neutron.cc
            ${coot_src}/ideal/coot-compare-residues.cc
            ${coot_src}/ideal/link-restraints.cc
            ${coot_src}/ideal/mods.cc
            ${coot_src}/ideal/extra-restraints-kk.cc
            ${coot_src}/ideal/chirals.cc
            ${coot_src}/ideal/extra-restraints.cc
            ${coot_src}/ideal/gradients.cc
            ${coot_src}/ideal/regularize-minimol.cc
            ${coot_src}/ideal/flanking.cc
            ${coot_src}/ideal/distortion.cc
            ${coot_src}/ideal/simple-restraint-kk.cc
            ${coot_src}/ideal/crankshaft.cc
            ${coot_src}/ideal/phi-psi.cc
            ${coot_src}/ideal/lograma.cpp
            ${coot_src}/ideal/process-dfs-in-range.cc
            ${coot_src}/ideal/write-png.cc
            ${coot_src}/ideal/zo-rama.cc
            ${coot_src}/tw/tw.cc
            ${coot_src}/utils/radix.cc
            ${coot_src}/utils/coot-utils.cc
            ${coot_src}/utils/hsv-rgb.cc
            ${coot_src}/utils/backward.cpp
            ${coot_src}/utils/align-utils.cc
            ${coot_src}/utils/threads.cc
            ${coot_src}/utils/coot-fasta.cc
            ${coot_src}/utils/win-compat.cc
            ${coot_src}/utils/atomic-numbers.cc
            ${coot_src}/utils/split-indices.cc
            ${coot_src}/utils/gdk-keysyms.cc
            ${coot_src}/utils/dodec.cc
            ${coot_src}/utils/logging.cc
            ${coot_src}/utils/pir-alignment.cc
            ${coot_src}/db-main/db-strands.cc
            ${coot_src}/high-res/residue-distortions.cc
            ${coot_src}/high-res/coot-atom-graph.cc
            ${coot_src}/high-res/sequence-assignment.cc
            ${coot_src}/high-res/high-res.cc
            ${coot_src}/protein_db/protein_db.cpp
            ${coot_src}/protein_db/protein_db_utils.cpp
            ${coot_src}/protein_db/rna_db.cpp
            ${coot_src}/coot-utils/merge-atom-selections.cc
            ${coot_src}/coot-utils/coot-least-squares.cc
            ${coot_src}/coot-utils/shiftfield.cpp
            ${coot_src}/coot-utils/stack-and-pair.cc
            ${coot_src}/coot-utils/spherical-harmonics.cc
            ${coot_src}/coot-utils/peak-search-from-list.cc
            ${coot_src}/coot-utils/helix-analysis.cc
            ${coot_src}/coot-utils/exp-fit.cc
            ${coot_src}/coot-utils/segmap.cc
            ${coot_src}/coot-utils/lidia-core-functions.cc
            ${coot_src}/coot-utils/cablam-markup.cc
            ${coot_src}/coot-utils/fib-sphere.cc
            ${coot_src}/coot-utils/atom-selection-container.cc
            ${coot_src}/coot-utils/coot-coord-utils-glyco.cc
            ${coot_src}/coot-utils/coot-coord-utils.cc
            ${coot_src}/coot-utils/merge-molecules.cc
            ${coot_src}/coot-utils/coot-coord-extras.cc
            ${coot_src}/coot-utils/water-coordination.cc
            ${coot_src}/coot-utils/bonded-pairs.cc
            ${coot_src}/coot-utils/c-beta-deviations.cc
            ${coot_src}/coot-utils/atom-tools.cc
            ${coot_src}/coot-utils/blob-line.cc
            ${coot_src}/coot-utils/coot-rama.cc
            ${coot_src}/coot-utils/secondary-structure-headers.cc
            ${coot_src}/coot-utils/helix-like.cc
            ${coot_src}/coot-utils/polar-atoms.cc
            ${coot_src}/coot-utils/plane-utils.cc
            ${coot_src}/coot-utils/coot-tree-extras.cc
            ${coot_src}/coot-utils/coot-map-heavy.cc
            ${coot_src}/coot-utils/jed-flip.cc
            ${coot_src}/coot-utils/pepflip-using-difference-map.cc
            ${coot_src}/coot-utils/atom-overlaps.cc
            ${coot_src}/coot-utils/bricks.cc
            ${coot_src}/coot-utils/mutate.cc
            ${coot_src}/coot-utils/elastic.cc
            ${coot_src}/coot-utils/coot-map-utils.cc
            ${coot_src}/coot-utils/coot-h-bonds.cc
            ${coot_src}/coot-utils/coot-coord-utils-nucleotides.cc
            ${coot_src}/coot-utils/glyco-torsions.cc
            ${coot_src}/coot-utils/fast-eigens.cc
            ${coot_src}/coot-utils/coot-shelx-ins.cc
            ${coot_src}/coot-utils/contacts-by-bricks.cc
            ${coot_src}/coot-utils/lsq-improve.cc
            ${coot_src}/coot-utils/trim.cc
            ${coot_src}/coot-utils/coot_shiftfield.cpp
            ${coot_src}/coot-utils/coot-fffear.cc
            ${coot_src}/coot-utils/coot-coord-lsq.cc
            ${coot_src}/coot-utils/strand-fragments.cc
            ${coot_src}/coot-utils/bonded-atoms.cc
            ${coot_src}/coot-utils/fragment-container.cc
            ${coot_src}/coot-utils/maps-spherical.cc
            ${coot_src}/coot-utils/peak-search.cc
            ${coot_src}/coot-utils/xmap-stats.cc
            ${coot_src}/coot-utils/sfcalc-genmap.cc
            ${coot_src}/coot-utils/reduce.cc
            ${coot_src}/coot-utils/contact-info.cc
            ${coot_src}/coot-utils/hole.cc
            ${coot_src}/coot-utils/slurp-map.cc
            ${coot_src}/coot-utils/read-sm-cif.cc
            ${coot_src}/coords/phenix-geo-bonds.cc
            ${coot_src}/coords/Bond_lines_ext.cc
            ${coot_src}/coords/mmdb-crystal.cc
            ${coot_src}/coords/utils.cc
            ${coot_src}/coords/phenix-geo.cc
            ${coot_src}/coords/mmdb.cc
            ${coot_src}/coords/Cartesian.cc
            ${coot_src}/coords/loop-path.cc
            ${coot_src}/coords/coot-close.cc
            ${coot_src}/coords/cos-sin.cc
            ${coot_src}/coords/rotamer-markup.cc
            ${coot_src}/coords/mmdb-extras.cc
            ${coot_src}/coords/Bond_lines.cc
            ${coot_src}/geometry/cif-parse-link.cc
            ${coot_src}/geometry/cif-parse.cc
            ${coot_src}/geometry/main-chain.cc
            ${coot_src}/geometry/cif-parse-mon-lib.cc
            ${coot_src}/geometry/protein-geometry.cc
            ${coot_src}/geometry/dictionary-residue.cc
            ${coot_src}/geometry/dict-utils.cc
            ${coot_src}/geometry/read-sbase.cc
            ${coot_src}/geometry/hydrophobic.cc
            ${coot_src}/geometry/energy-lib.cc
            ${coot_src}/geometry/srs-interface.cc
            ${coot_src}/geometry/available-comp-id.cc
            ${coot_src}/geometry/residue-and-atom-specs.cc
            ${coot_src}/geometry/protein-geometry-mmdb.cc
            ${coot_src}/geometry/mol-utils.cc
            ${coot_src}/geometry/dreiding.cc
            ${coot_src}/geometry/chem_mods.cc
            ${coot_src}/geometry/link.cc
            ${coot_src}/geometry/lbg-graph.cc
            ${coot_src}/geometry/glyco-tools.cc
            ${coot_src}/geometry/mol-utils-2.cc
            ${coot_src}/geometry/protein-donor-acceptors.cc
            ${coot_src}/build/CalphaBuild.cc
            ${coot_src}/mini-mol/mini-mol-utils.cc
            ${coot_src}/mini-mol/mini-mol.cc
            ${coot_src}/mini-mol/atom-quads.cc
            ${coot_src}/cootilus/nautilus-tools.cpp
            ${coot_src}/cootilus/cootilus-build.cpp
            ${coot_src}/cootilus/nucleicacid_db.cpp
            ${coot_src}/cootilus/nautilus-join.cpp
            ${coot_src}/cootilus/nautilus-target.cpp
            ${coot_src}/cootilus/nautilus-ss-find.cpp
            ${coot_src}/lidia-core/lbg-shared.cc
            ${coot_src}/lidia-core/cod-types.cc
            ${coot_src}/lidia-core/cod-atom-types.cc
            ${coot_src}/lidia-core/chemical-feature-clusters.cc
            ${coot_src}/lidia-core/bond-record-container-t.cc
            ${coot_src}/lidia-core/atom-types.cc
            ${coot_src}/lidia-core/cod-atom-type-t.cc
            ${coot_src}/lidia-core/lbg-molfile.cc
            ${coot_src}/lidia-core/lig-build.cc
            ${coot_src}/lidia-core/rdkit-interface.cc
            ${coot_src}/lidia-core/get-residue.cc
            ${coot_src}/lidia-core/bond-table-record-t.cc
            ${coot_src}/angles/AngleInfo-angle-torsions-autogen.cc
            ${coot_src}/angles/AngleInfo-theta_2-autogen.cc
            ${coot_src}/angles/AngleInfo.cc
            ${coot_src}/ligand/dipole.cc
            ${coot_src}/ligand/chi-angles.cc
            ${coot_src}/ligand/trace-2.cc
            ${coot_src}/ligand/monomer-utils.cc
            ${coot_src}/ligand/ligand-score-card.cc
            ${coot_src}/ligand/multi-peptide.cc
            ${coot_src}/ligand/chi-angles-autogen.cc
            ${coot_src}/ligand/ligand.cc
            ${coot_src}/ligand/helix-placement.cc
            ${coot_src}/ligand/backrub-rotamer.cc
            ${coot_src}/ligand/dunbrack.cc
            ${coot_src}/ligand/primitive-chi-angles.cc
            ${coot_src}/ligand/rotamer-tables.cc
            ${coot_src}/ligand/side-chain-densities.cc
            ${coot_src}/ligand/richardson-rotamers.cc
            ${coot_src}/ligand/ligands-db.cc
            ${coot_src}/ligand/ligand-extras.cc
            ${coot_src}/ligand/residue_by_phi_psi.cc
            ${coot_src}/ligand/base-pairing.cc
            ${coot_src}/ligand/ideal-rna.cc
            ${coot_src}/ligand/wiggly-ligand.cc
            ${coot_src}/ligand/torsion-general.cc
            ${coot_src}/ligand/joining-residues.cc
            ${coot_src}/ligand/richardson-rotamer.cc
            ${coot_src}/ligand/side-chain.cc
            ${coot_src}/ligand/rigid-body.cc
            ${coot_src}/ligand/trace.cc
            ${coot_src}/ligand/smith-waterman.cc
            ${coot_src}/ligand/tube-finder.cc
            ${coot_src}/ligand/fast-ss-search.cc
            ${coot_src}/ligand/rotamer.cc
            ${coot_src}/cootaneer/buccaneer-sequence.cpp
            ${coot_src}/cootaneer/cootaneer-sequence.cpp
            ${coot_src}/cootaneer/buccaneer-lib.cpp
            ${coot_src}/cootaneer/buccaneer-prot.cpp
            ${coot_src}/src/coot-nomenclature.cc
            ${coot_src}/db-main/db-main.cc
            ${coot_src}/compat/getopt.cc
            ${coot_src}/compat/getopt1.cc
            ${coot_src}/api/prideout-octasphere.cc
            ${coot_src}/api/oct.cc
            ${coot_src}/api/interfaces.cc
            ${coot_src}/api/colour-functions.cc
            ${coot_src}/api/coot-colour.cc
            ${coot_src}/api/coot_molecule.cc
            ${coot_src}/api/coot_molecule_maps.cc
            ${coot_src}/api/coot_molecule_bonds.cc
            ${coot_src}/api/molecules_container.cc
            ${coot_src}/api/molecules_container_modelling.cc
            ${coot_src}/api/rama-plot-phi-psi.cc
            ${coot_src}/api/simple-mesh.cc
            ${coot_src}/api/make-a-dodec.cc
            ${coot_src}/api/cylinder.cc
            ${coot_src}/api/add-terminal-residue.cc

)

# target_include_directories(coot PRIVATE ${coot_src} ${coot_src}/coords ${coot_src}/utils ../include ../example ../checkout/ssm ../checkout/mmdb2 ../checkout/clipper ../checkout/fftw-2.1.5/fftw ../checkout/fftw-2.1.5/rfftw ../checkout/libccp4/ccp4 ../boost)

target_include_directories(cootapi PRIVATE ${coot_src} ${coot_src}/coords ${coot_src}/utils ${dependency_install_dir}/include ${dependency_install_dir}/include/python3.9)
target_link_directories(cootapi PUBLIC ${dependency_install_dir}/lib)
target_link_libraries(cootapi PUBLIC clipper-minimol clipper-phs clipper-mmdb clipper-cif clipper-core clipper-contrib clipper-ccp4 fftw rfftw mmdb2 ccp4c png gsl gslcblas)

# add in -DCLIPPER_HAS_TOP8000=1
set_target_properties(cootapi PROPERTIES COMPILE_FLAGS "-DHAVE_BOOST -DUSE_LIBPNG=1 -DEMSCRIPTEN=1")


add_executable(test-molecules-container ${coot_src}/api/test_molecules_container.cc)
# # target_link_options(test-molecules-container PRIVATE -sALLOW_MEMORY_GROWTH=1 -sNODERAWFS=1 --pre-js coot_env.js)
target_link_directories(test-molecules-container PUBLIC ${dependency_install_dir}/lib .)

# how to use the above: set(dependency_install_dir $ENV{HOME}/autobuild/$ENV{OS}-$ENV{HOST}-gtk3)
target_link_libraries(test-molecules-container cootapi)
target_include_directories(test-molecules-container PRIVATE ${coot_src} ${coot_src}/src
  ../checkout/clipper/ ../checkout/mmdb2 ../checkout/mmdb2/mmdb2 ../checkout/ssm
  ../checkout/libccp4/ccp4 ../fftw ../rfftw ${dependency_install_dir}/include)
set_target_properties(test-molecules-container PROPERTIES COMPILE_FLAGS "-DHAVE_BOOST -DUSE_LIBPNG=1 -DEMSCRIPTEN=1")

# target_compile_options(test-molecules-container PUBLIC -DFFTW_ENABLE_FLOAT -D_CCP4_NODE_EXAMPLE_ -DNODERAWFS)



FILE(WRITE ${CMAKE_BINARY_DIR}/coot/coot_env.js
Module.preRun = [function() {ENV.COOT_DATA_DIR="${coot_src}"}, function() {ENV.CCP4_LIB="${CMAKE_BINARY_DIR}/ccp4_lib"}]
)

FILE(MAKE_DIRECTORY "data/coot/data/metal/")
FILE(MAKE_DIRECTORY "data/ccp4_lib/data/monomers/a")
FILE(MAKE_DIRECTORY "data/ccp4_lib/data/monomers/c")
FILE(MAKE_DIRECTORY "data/ccp4_lib/data/monomers/g")
FILE(MAKE_DIRECTORY "data/ccp4_lib/data/monomers/p")
FILE(MAKE_DIRECTORY "data/ccp4_lib/data/monomers/h")
FILE(MAKE_DIRECTORY "data/ccp4_lib/data/monomers/i")
FILE(MAKE_DIRECTORY "data/ccp4_lib/data/monomers/l")
FILE(MAKE_DIRECTORY "data/ccp4_lib/data/monomers/m")
FILE(MAKE_DIRECTORY "data/ccp4_lib/data/monomers/s")
FILE(MAKE_DIRECTORY "data/ccp4_lib/data/monomers/t")
FILE(MAKE_DIRECTORY "data/ccp4_lib/data/monomers/v")
FILE(MAKE_DIRECTORY "data/ccp4_lib/data/monomers/e")
FILE(MAKE_DIRECTORY "data/ccp4_lib/data/monomers/u")
FILE(MAKE_DIRECTORY "data/ccp4_lib/data/monomers/d")
FILE(MAKE_DIRECTORY "data/ccp4_lib/data/monomers/h")
FILE(MAKE_DIRECTORY "data/ccp4_lib/data/monomers/n")
FILE(MAKE_DIRECTORY "data/ccp4_lib/data/monomers/list")
FILE(COPY ${coot_src}/data/metal/metal-N-distance.table DESTINATION data/coot/data/metal )
FILE(COPY ${coot_src}/data/metal/metal-O-distance.table DESTINATION data/coot/data/metal )
FILE(COPY ${coot_src}/data/metal/metal-S-distance.table DESTINATION data/coot/data/metal )

# install(TARGETS mini-rsr-web DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../react-app/public/wasm)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/mini-rsr-web.wasm
              ${CMAKE_CURRENT_BINARY_DIR}/mini-rsr-web.worker.js
              ${CMAKE_CURRENT_BINARY_DIR}/mini-rsr-web.js
              ${CMAKE_CURRENT_BINARY_DIR}/mini-rsr-web.data
              ${CMAKE_CURRENT_BINARY_DIR}/coot_env_web.js
              DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../react-app/public/wasm)

