cmake_minimum_required(VERSION 3.13)
cmake_policy(SET CMP0048 NEW)
project(coot)
find_package(PkgConfig)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# C++20 is only used because newest versions of RDKit (March 2025 and newer) seems to require it
# Without RDKit, Coot currently needs C++17
set(CMAKE_CXX_STANDARD 20)

if (APPLE)
else()
   find_package(CCP4 REQUIRED mmdb2 ssm clipper-core clipper-mmdb clipper-ccp4 clipper-contrib clipper-minimol clipper-cif)
endif(APPLE)

find_package(BLAS)
find_package(GSL REQUIRED)
find_package(PNG REQUIRED)
find_package(GEMMI REQUIRED)
find_package(Threads REQUIRED)
find_package(Python REQUIRED)

find_package(Boost COMPONENTS iostreams regex thread serialization REQUIRED)
find_package(glm REQUIRED)

set(THREADS_PREFER_PTHREAD_FLAG ON)

find_package(Python 3
  REQUIRED COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule)

# This is flaky.
# If it doesn't work, first try setting nanobind_DIR to wherever_your_python_is/lib/python3.$your_python_version/site-packages/nanobind/cmake
# If that doesn't work, make sure that CMAKE_INSTALL_PREFIX is set
find_package(nanobind CONFIG REQUIRED)

# only via homebrew on macOS (not tried to compile them myself) - see below also though.
if (APPLE)
else()
   find_library(libelf elf)
   find_library(libdwarf dw)
endif(APPLE)

# From Charles
if (WIN32)
   string(REPLACE "\\" "/" Python_SITEARCH "${Python_SITEARCH}")
   string(REPLACE "\\" "/" Python_SITELIB  "${Python_SITELIB}")
endif(WIN32)


set(coot_src ${CMAKE_CURRENT_SOURCE_DIR})

add_definitions(-DANALYSE_REFINEMENT_TIMING -DPACKAGE_VERSION="1.1.20"
  -DHAVE_CXX_THREAD -DHAVE_SSMLIB -DHAVE_BOOST -DUSE_LIBPNG=1 -DCLIPPER_HAS_TOP8000 -DUSE_GEMMI -DHAVE_BOOST_THREAD -DHAVE_GSL
  -DPKGDATADIR="${CMAKE_INSTALL_PREFIX}/share/coot" -DXDATADIR="${CMAKE_INSTALL_PREFIX}/share")

# 20250314-PE: Charles, I guess that you will need to comment out this block
#
add_definitions(-DMAKE_ENHANCED_LIGAND_TOOLS="1")
message("-- Found ENHANCED_LIGAND_TOOLS setting - now requiring the RDKit")
find_package(Cairo REQUIRED)
find_package(RDKit CONFIG COMPONENTS RDGeneral REQUIRED)
set(RDKit_libs_for_libcootapi RDKit::rdkit_base RDKit::Descriptors RDKit::ForceFieldHelpers RDKit::ForceField RDKit::Optimizer RDKit::MolDraw2D  RDKit::MolChemicalFeatures)

# we can be on Apple, and also have access to homebrew dwarf and elf - but this typically not the case
# should this actually be ENABLE_BACKWARD? (util/backward.cpp is not conditionally compiled).
# This doesn't seem to wort at the moment though - maybe a  problem with libdwarf.
if (ENABLE_DWARF_ELF)
   find_library(libelf elf)
   find_library(libdwarf dw)
   add_definitions(-DBACKWARD_HAS_DW)
   message("-- Enabled dwarf and elf")
endif()

#  pyrogen, auxil, goograph, src ignored.

add_library(cootmoleculestotriangles SHARED
            ${coot_src}/MoleculesToTriangles/CXXSurface/CXXBall.cpp
            ${coot_src}/MoleculesToTriangles/CXXSurface/CXXChargeTable.cpp
            ${coot_src}/MoleculesToTriangles/CXXSurface/CXXCircle.cpp
            ${coot_src}/MoleculesToTriangles/CXXSurface/CXXCircleNode.cpp
            ${coot_src}/MoleculesToTriangles/CXXSurface/CXXCoord.cpp
            ${coot_src}/MoleculesToTriangles/CXXSurface/CXXCreator.cpp
            ${coot_src}/MoleculesToTriangles/CXXSurface/CXXException.cpp
            ${coot_src}/MoleculesToTriangles/CXXSurface/CXXFFTSolventMap.cpp
            ${coot_src}/MoleculesToTriangles/CXXSurface/CXXFortranFile.cpp
            ${coot_src}/MoleculesToTriangles/CXXSurface/CXXNewHood.cpp
            ${coot_src}/MoleculesToTriangles/CXXSurface/CXXQADSurface.cpp
            ${coot_src}/MoleculesToTriangles/CXXSurface/CXXQADSurfaceStaticData.cpp
            ${coot_src}/MoleculesToTriangles/CXXSurface/CXXSpace.cpp
            ${coot_src}/MoleculesToTriangles/CXXSurface/CXXSphereElement.cpp
            ${coot_src}/MoleculesToTriangles/CXXSurface/CXXSphereNode.cpp
            ${coot_src}/MoleculesToTriangles/CXXSurface/CXXSphereTriangle.cpp
            ${coot_src}/MoleculesToTriangles/CXXSurface/CXXSphereTriangleEdge.cpp
            ${coot_src}/MoleculesToTriangles/CXXSurface/CXXSurface.cpp
            ${coot_src}/MoleculesToTriangles/CXXSurface/CXXSurfaceMaker.cpp
            ${coot_src}/MoleculesToTriangles/CXXSurface/CXXSurfaceVertex.cpp
            ${coot_src}/MoleculesToTriangles/CXXSurface/CXXTorusElement.cpp
            ${coot_src}/MoleculesToTriangles/CXXSurface/CXXTorusNode.cpp
            ${coot_src}/MoleculesToTriangles/CXXSurface/CXXTriangle.cpp
            ${coot_src}/MoleculesToTriangles/CXXSurface/CXXUtils.cpp
            ${coot_src}/MoleculesToTriangles/CXXSurface/CXXVCN.cpp
            #${coot_src}/MoleculesToTriangles/CXXSurface/ElectroSurface.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/AtomPropertyRampColorRule.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/BallsPrimitive.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/BondsPrimitive.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/BoxSectionPrimitive.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/ColorRule.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/ColorScheme.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/CompoundSelection.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/CXXMatrix.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/CylinderPoint.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/CylindersPrimitive.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/DirectedPositionedEntity.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/DiscreteSegment.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/DishyBase.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/DisplayPrimitive.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/ElectrostaticColorScheme.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/FlatFanPrimitive.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/Light.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/LinesPrimitive.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/MolecularRepresentation.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/MolecularRepresentationInstance.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/MyMolecule.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/NRStuff.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/oglPolyhedron.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/Polyhedron.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/PositionedEntity.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/Quaternion.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/Renderer.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/Representation.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/RepresentationInstance.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/RotatedTranslatedScaledEntity.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/SecondaryColorScheme.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/SolidColorRule.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/SpherePrimitive.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/SticksPrimitive.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/SurfacePrimitive.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/Vector2.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/Vector3.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/VertexColorNormalPrimitive.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/lfit.c
            #  ${coot_src}/MoleculesToTriangles/CXXClasses/tubes.cc - not ready yet
)


add_library(cootapi SHARED
            ${coot_src}/lidia-core/svg-molecule.cc
            ${coot_src}/lidia-core/cod-types.cc
            ${coot_src}/lidia-core/cod-atom-types.cc
            ${coot_src}/lidia-core/chemical-feature-clusters.cc
            ${coot_src}/lidia-core/bond-record-container-t.cc
            ${coot_src}/lidia-core/atom-types.cc
            ${coot_src}/lidia-core/cod-atom-type-t.cc
            ${coot_src}/lidia-core/rdkit-interface.cc
            ${coot_src}/lidia-core/bond-table-record-t.cc
            ${coot_src}/api/rigid-body-fit.cc
            ${coot_src}/api/blender-mesh.cc
            ${coot_src}/api/add-terminal-residue.cc
            ${coot_src}/api/interfaces.cc
            ${coot_src}/api/instancing.cc
            ${coot_src}/api/colour-functions.cc
            ${coot_src}/api/coot-colour.cc
            ${coot_src}/api/coot-molecule.cc
            ${coot_src}/api/coot-molecule-backup.cc
            ${coot_src}/api/coot-molecule-blender.cc
            ${coot_src}/api/coot-molecule-bonds.cc
            ${coot_src}/api/coot-molecule-bonds-instanced.cc
            ${coot_src}/api/coot-molecule-json.cc
            ${coot_src}/api/coot-molecule-maps.cc
            ${coot_src}/api/coot-molecule-change-chain-id.cc
            ${coot_src}/api/coot-molecule-moltris.cc
            ${coot_src}/api/coot-molecule-validation.cc
            ${coot_src}/api/coot-molecule-merge-molecules.cc
            ${coot_src}/api/coot-molecule-modelling.cc
            ${coot_src}/api/coot-molecule-refine.cc
            ${coot_src}/api/coot-molecule-replace-fragment.cc
            ${coot_src}/api/ghost-molecule-display.cc
            ${coot_src}/api/model-analysis.cc
            ${coot_src}/api/molecules-container.cc
            ${coot_src}/api/molecules-container-blender.cc
            ${coot_src}/api/molecules-container-bucca-ml-growing.cc
            ${coot_src}/api/molecules-container-gltf.cc
            ${coot_src}/api/molecules-container-maps.cc
            ${coot_src}/api/molecules-container-mmrrcc.cc
            ${coot_src}/api/molecules-container-superpose.cc
            ${coot_src}/api/molecules-container-modelling.cc
            ${coot_src}/api/molecules-container-ligand-fitting.cc
            ${coot_src}/api/mmcif-tests.cc
            ${coot_src}/api/phi-psi-prob.cc
            ${coot_src}/api/moorhen-h-bonds.cc
            ${coot_src}/api/rama-plot-phi-psi.cc
            ${coot_src}/analysis/daca.cc
            ${coot_src}/analysis/kolmogorov.cc
            ${coot_src}/analysis/stats.cc
            ${coot_src}/analysis/bfkurt.cc
            ${coot_src}/analysis/mogul.cc
            ${coot_src}/analysis/b-factor-histogram.cc
            ${coot_src}/analysis/typed-distances.cc
            ${coot_src}/analysis/cablam.cc
            ${coot_src}/analysis/chi-squared.cc
            ${coot_src}/density-contour/density-contour-triangles.cc
            ${coot_src}/density-contour/Vectors.cpp
            ${coot_src}/density-contour/CIsoSurface.cpp
            ${coot_src}/density-contour/occlusion.cc
            ${coot_src}/density-contour/gaussian-surface.cc
            ${coot_src}/ccp4mg-utils/quat.cc
            ${coot_src}/ccp4mg-utils/mgtree.cc
            ${coot_src}/ccp4mg-utils/cartesian.cc
            ${coot_src}/ccp4mg-utils/geomutil.cc
            ${coot_src}/ccp4mg-utils/matrix.cc
            ${coot_src}/skeleton/dm-utils.cc
            ${coot_src}/skeleton/BuildCas.cc
            ${coot_src}/skeleton/graphical_skel.cc
            ${coot_src}/ideal/add-linked-cho.cc
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
            ${coot_src}/ideal/zo-rama.cc
            ${coot_src}/utils/align-utils.cc
            ${coot_src}/utils/base64-encode-decode.cc
            ${coot_src}/utils/atomic-numbers.cc
            ${coot_src}/utils/backward.cpp
            ${coot_src}/utils/colour-holder.cc
            ${coot_src}/utils/coot-fasta.cc
            ${coot_src}/utils/coot-utils.cc
            ${coot_src}/utils/dodec.cc
            ${coot_src}/utils/radix.cc
            ${coot_src}/utils/gdk-keysyms.cc
            ${coot_src}/utils/hsv-rgb.cc
            ${coot_src}/utils/pir-alignment.cc
            ${coot_src}/utils/setup-syminfo.cc
            ${coot_src}/utils/split-indices.cc
            ${coot_src}/utils/logging.cc
            ${coot_src}/utils/static-logging.cc
            ${coot_src}/utils/threads.cc
            ${coot_src}/utils/win-compat.cc
            ${coot_src}/db-main/db-strands.cc
            ${coot_src}/high-res/residue-distortions.cc
            ${coot_src}/high-res/coot-atom-graph.cc
            ${coot_src}/high-res/sequence-assignment.cc
            ${coot_src}/high-res/high-res.cc
            ${coot_src}/protein_db/protein_db.cpp
            ${coot_src}/protein_db/protein_db_utils.cpp
            ${coot_src}/protein_db/rna_db.cpp
            ${coot_src}/pli/flev-annotations.cc
            ${coot_src}/pli/flev.cc
            ${coot_src}/pli/flev-annotations.cc
            ${coot_src}/pli/flev-attached-hydrogens.cc
            ${coot_src}/pli/optimise-residue-circles.cc
            ${coot_src}/pli/residue-circles.cc
            ${coot_src}/pli/dots-representation-info.cc
            ${coot_src}/pli/pi-stacking.cc
            ${coot_src}/pli/protein-ligand-interactions.cc
            ${coot_src}/pli/pi-stacking.cc
            ${coot_src}/pli/sdf-interface.cc
            ${coot_src}/pli/specs.cc
            ${coot_src}/coot-utils/cfc.cc
            ${coot_src}/coot-utils/coot-align.cc
            ${coot_src}/coot-utils/dict-link-info.cc
            ${coot_src}/coot-utils/cmtz-interface.cc
            ${coot_src}/coot-utils/dict-link-info.cc
            ${coot_src}/coot-utils/lidia-core-functions.cc
            ${coot_src}/coot-utils/atom-selection-container.cc
            ${coot_src}/coot-utils/merge-atom-selections.cc
            ${coot_src}/coot-utils/coot-least-squares.cc
            ${coot_src}/coot-utils/diff-diff-map-peaks.cc
            ${coot_src}/coot-utils/gltf-export.cc
            ${coot_src}/coot-utils/shiftfield.cpp
            ${coot_src}/coot-utils/stack-and-pair.cc
            ${coot_src}/coot-utils/spherical-harmonics.cc
            ${coot_src}/coot-utils/peak-search-from-list.cc
            ${coot_src}/coot-utils/helix-analysis.cc
            ${coot_src}/coot-utils/segmap.cc
            ${coot_src}/coot-utils/lidia-core-functions.cc
            ${coot_src}/coot-utils/cablam-markup.cc
            ${coot_src}/coot-utils/fib-sphere.cc
            ${coot_src}/coot-utils/find-water-baddies.cc
            ${coot_src}/coot-utils/atom-selection-container.cc
            ${coot_src}/coot-utils/blob-line.cc
            ${coot_src}/coot-utils/coot-coord-utils-glyco.cc
            ${coot_src}/coot-utils/coot-coord-utils.cc
            ${coot_src}/coot-utils/merge-molecules.cc
            ${coot_src}/coot-utils/texture-as-floats.cc
            ${coot_src}/coot-utils/coot-coord-extras.cc
            ${coot_src}/coot-utils/water-coordination.cc
            ${coot_src}/coot-utils/bonded-pairs.cc
            ${coot_src}/coot-utils/c-beta-deviations.cc
            ${coot_src}/coot-utils/atom-tools.cc
            ${coot_src}/coot-utils/coot-rama.cc
            ${coot_src}/coot-utils/secondary-structure-headers.cc
            ${coot_src}/coot-utils/helix-like.cc
            ${coot_src}/coot-utils/polar-atoms.cc
            ${coot_src}/coot-utils/plane-utils.cc
            ${coot_src}/coot-utils/pae.cc
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
            ${coot_src}/coot-utils/merge-C-and-N-terminii.cc
            ${coot_src}/coot-utils/trim.cc
            ${coot_src}/coot-utils/coot_shiftfield.cpp
            ${coot_src}/coot-utils/coot-fffear.cc
            ${coot_src}/coot-utils/coot-coord-lsq.cc
            ${coot_src}/coot-utils/strand-fragments.cc
            ${coot_src}/coot-utils/bonded-atoms.cc
            ${coot_src}/coot-utils/fragment-container.cc
            ${coot_src}/coot-utils/positron.cc
            ${coot_src}/coot-utils/maps-spherical.cc
            ${coot_src}/coot-utils/peak-search.cc
            ${coot_src}/coot-utils/pugixml.cpp
            ${coot_src}/coot-utils/xmap-stats.cc
            ${coot_src}/coot-utils/sfcalc-genmap.cc
            ${coot_src}/coot-utils/reduce.cc
            ${coot_src}/coot-utils/contact-info.cc
            ${coot_src}/coot-utils/hole.cc
            ${coot_src}/coot-utils/slurp-map.cc
            ${coot_src}/coot-utils/read-sm-cif.cc
            ${coot_src}/coot-utils/exp-fit.cc
            ${coot_src}/coot-utils/cylinder.cc
            ${coot_src}/coot-utils/cylinder-with-rotation-translation.cc
            ${coot_src}/coot-utils/prideout-octasphere.cc
            ${coot_src}/coot-utils/oct.cc
            ${coot_src}/coot-utils/shapes.cc
            ${coot_src}/coot-utils/simple-mesh.cc
            ${coot_src}/coot-utils/gl-matrix.cc
            ${coot_src}/coot-utils/gl-matrix-eigensystem.cc
            ${coot_src}/coords/Bond_lines.cc
            ${coot_src}/coords/Bond_lines_ext.cc
            ${coot_src}/coords/coot-close.cc
            ${coot_src}/coords/cos-sin.cc
            ${coot_src}/coords/Cartesian.cc
            ${coot_src}/coords/graphics-line.cc
            ${coot_src}/coords/mmdb-crystal.cc
            ${coot_src}/coords/loop-path.cc
            ${coot_src}/coords/mmdb.cc
            ${coot_src}/coords/mmdb-extras.cc
            ${coot_src}/coords/phenix-geo.cc
            ${coot_src}/coords/phenix-geo-bonds.cc
            ${coot_src}/coords/rotamer-markup.cc
            ${coot_src}/coords/utils.cc
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
            ${coot_src}/geometry/intercept-3GP.cc
            ${coot_src}/build/CalphaBuild.cc
            ${coot_src}/mini-mol/mini-mol-utils.cc
            ${coot_src}/mini-mol/mini-mol.cc
            ${coot_src}/mini-mol/atom-quads.cc
            ${coot_src}/cootilus/cootilus-build.cpp
            ${coot_src}/cootilus/nautilus-join.cpp
            ${coot_src}/cootilus/nautilus-ss-find.cpp
            ${coot_src}/cootilus/nautilus-target.cpp
            ${coot_src}/cootilus/nautilus-tools.cpp
            ${coot_src}/cootilus/nucleicacid_db.cpp
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
            ${coot_src}/ligand/dSFMT.c
            ${coot_src}/ligand/chi-angles.cc
            ${coot_src}/ligand/monomer-utils.cc
            ${coot_src}/ligand/ligand-score-card.cc
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
            ${coot_src}/ligand/libres-tracer.cc
            ${coot_src}/ligand/new-residue-by-3-phi-psi.cc
            ${coot_src}/ligand/rama-rsr-extend-fragments.cc
            ${coot_src}/ligand/residue_by_phi_psi.cc
            ${coot_src}/ligand/base-pairing.cc
            ${coot_src}/ligand/ideal-rna.cc
            ${coot_src}/ligand/wiggly-ligand.cc
            ${coot_src}/ligand/torsion-general.cc
            ${coot_src}/ligand/joining-residues.cc
            ${coot_src}/ligand/richardson-rotamer.cc
            ${coot_src}/ligand/side-chain.cc
            ${coot_src}/ligand/rigid-body.cc
            ${coot_src}/ligand/smith-waterman.cc
            ${coot_src}/ligand/tube-finder.cc
            ${coot_src}/ligand/fast-ss-search.cc
            ${coot_src}/ligand/rotamer.cc
            #${coot_src}/cootaneer/buccaneer-sequence.cpp
            #${coot_src}/cootaneer/cootaneer-sequence.cpp
            #${coot_src}/cootaneer/buccaneer-lib.cpp
            #${coot_src}/cootaneer/buccaneer-prot.cpp
            ${coot_src}/src/coot-nomenclature.cc
            ${coot_src}/db-main/db-main.cc
            ${coot_src}/compat/getopt.cc
            ${coot_src}/compat/getopt1.cc
            ${coot_src}/ideal/write-png.cc
            )

# wretched mae
# find_library(MAEPARSER_LIBS maeparser)
# target_link_libraries(cootapi PRIVATE ${MAEPARSER_LIBS})

nanobind_add_module(coot_headless_api ${coot_src}/api/molecules-container-nanobind.cc
                                      ${coot_src}/api/molecules-container-nanobind-functions.cc)

# on macOS only?
# target_compile_options(cootapi PRIVATE -fsanitize=address -fno-omit-frame-pointer)
# target_link_libraries(cootapi PRIVATE -fsanitize=address) # Linker flags

set(MAKE_COOT_HEADLESS_API_PYI 1)

if (MAKE_COOT_HEADLESS_API_PYI)
  # thanks Jordan (install coot_headless_api.pyi in the right place)
  nanobind_add_stub(coot_headless_api_ext_stub
    MODULE coot_headless_api
    OUTPUT ${Python_SITELIB}/coot_headless_api.pyi
    PYTHON_PATH $<TARGET_FILE_DIR:coot_headless_api>
    DEPENDS coot_headless_api
  )
  # MARKER_FILE py.typed
endif()

# macOS doesn't have libelf
if (ENABLE_DWARF)
  # find_library(libdwarf dw)
  #  add_definitions(-DBACKWARD_HAS_DW)
  #  message("-- Enabled dwarf (macOS)")
  find_library(LIBDWARF dwarf PATHS /opt/homebrew/lib)
     if (NOT LIBDWARF)
        message(FATAL_ERROR "libdwarf not found. Please install with 'brew install libdwarf'")
     endif()
  target_link_libraries(cootapi PRIVATE ${LIBDWARF})
endif()

target_include_directories(cootmoleculestotriangles PRIVATE ${MMDB2_INCLUDE_DIR}
   ${CLIPPER-CORE_INCLUDE_DIR} ${coot_src} ${CLIPPER-MMDB_INCLUDE_DIR} ${CLIPPER-CCP4_INCLUDE_DIR} ${FFTW2_INCLUDE_DIRS})

# use MMDB2_INCLUDE_DIR to find the RDKIT for now.
target_include_directories(cootapi PRIVATE ${coot_src} ${coot_src}/coords ${coot_src}/utils ${coot_src}/ccp4mg-utils ${coot_src}/ccp4mg-utils/util ${coot_src}/ccp4mg-utils/pygl ${MMDB2_INCLUDE_DIR} ${SSM_INCLUDE_DIR} ${CLIPPER-CORE_INCLUDE_DIR} ${CLIPPER-MMDB_INCLUDE_DIR} ${CLIPPER-CCP4_INCLUDE_DIR} ${FFTW2_INCLUDE_DIRS} ${RDKit_INCLUDE_DIRS} ${GLM_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS} ${PNG_PNG_INCLUDE_DIR} ${RDKit_INCLUDE_DIR})

target_link_libraries(cootmoleculestotriangles PRIVATE ${SSM_LIBRARY} ${CLIPPER-MINIMOL_LIBRARY} ${CLIPPER-MMDB_LIBRARY} ${CLIPPER-CIF_LIBRARY} ${CLIPPER-CORE_LIBRARY} ${CLIPPER-CCP4_LIBRARY} ${CLIPPER-CONTRIB_LIBRARY} ${FFTW2_LIBRARY} ${RFFTW2_LIBRARY} ${MMDB2_LIBRARY} ${CCP4C_LIBRARY})

# Charles adds Boost::iostreams Boost::filesystem Boost::system. I don't know why Boost::filesystem was added - I have removed it.
# 2025-09-20 14:15 PE: Boost::system is no longer a thing - I have removed it
target_link_libraries(cootapi PRIVATE cootmoleculestotriangles ${SSM_LIBRARY} ${CLIPPER-MINIMOL_LIBRARY} ${CLIPPER-MMDB_LIBRARY} ${CLIPPER-CIF_LIBRARY} ${CLIPPER-CORE_LIBRARY} ${CLIPPER-CCP4_LIBRARY} ${CLIPPER-CONTRIB_LIBRARY} ${FFTW2_LIBRARY} ${RFFTW2_LIBRARY} ${MMDB2_LIBRARY} ${CCP4C_LIBRARY} ${ZLIB_LIBRARIES} ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARY} ${RDKit_libs_for_libcootapi} ${PNG_LIBRARY} Boost::thread Boost::iostreams gemmi::gemmi_cpp ${libdwarf} ${libelf})

target_compile_features(cootapi PRIVATE cxx_std_17)

# Charles adds this MSVC CXX compiler test
if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  target_compile_options(cootapi PRIVATE
	/D_CRT_SECURE_NO_WARNINGS
	/D_SILENCE_ALL_CXX17_DEPRECATION_WARNINGS
	/bigobj
  )
  target_link_libraries(cootapi PRIVATE netapi32.lib)
endif()

target_include_directories(coot_headless_api PRIVATE ${coot_src} ${coot_src}/coords ${coot_src}/utils ${coot_src}/ccp4mg-utils
                                                     ${coot_src}/ccp4mg-utils/util ${coot_src}/ccp4mg-utils/pygl
                                                     ${MMDB2_INCLUDE_DIR} ${SSM_INCLUDE_DIR} ${CLIPPER-CORE_INCLUDE_DIR}
                                                     ${CLIPPER-MMDB_INCLUDE_DIR} ${CLIPPER-CCP4_INCLUDE_DIR}
                                                     ${FFTW2_INCLUDE_DIRS}
                                                     ${RDKit_INCLUDE_DIRS} ${GLM_INCLUDE_DIRS}
                                                     ${GSL_INCLUDE_DIRS} ${PNG_PNG_INCLUDE_DIR})

add_executable(coot-findligand ${coot_src}/ligand/find-ligand.cc)

add_executable(coot-findwaters ${coot_src}/ligand/find-waters.cc)

add_executable(coot-mini-rsr ${coot_src}/ideal/mini-rsr.cc)

add_executable(coot-res-tracer-libcootapi ${coot_src}/ligand/res-tracer.cc)

add_executable(test-molecules-container ${coot_src}/api/test-molecules-container.cc
                                        ${coot_src}/api/filo-tests.cc
                                        ${coot_src}/api/lucrezia-tests.cc)

add_executable(mmcif-tests ${coot_src}/api/mmcif-tests-main.cc ${coot_src}/api/mmcif-tests.cc)

target_compile_features(coot-mini-rsr PRIVATE cxx_std_17)

target_compile_features(coot-res-tracer-libcootapi PRIVATE cxx_std_17)

target_link_libraries(coot_headless_api PUBLIC cootapi  ${CLIPPER-CORE_LIBRARY} ${CLIPPER-CCP4_LIBRARY} ${MMDB2_LIBRARY}
                      ${RDKit_libs_for_libcootapi} ${libdwarf} ${libelf})

target_include_directories(coot-findligand PRIVATE ${MMDB2_INCLUDE_DIR} ${CLIPPER-CORE_INCLUDE_DIR} ${coot_src})

target_include_directories(coot-findwaters PRIVATE ${MMDB2_INCLUDE_DIR} ${CLIPPER-CORE_INCLUDE_DIR} ${coot_src})

target_include_directories(coot-mini-rsr PRIVATE ${MMDB2_INCLUDE_DIR} ${coot_src} ${PNG_INCLUDE_DIR})

target_include_directories(coot-res-tracer-libcootapi PRIVATE ${MMDB2_INCLUDE_DIR} ${coot_src} ${PNG_INCLUDE_DIR})

target_link_libraries(coot-findligand PUBLIC cootapi ${CLIPPER-CORE_LIBRARY} ${CLIPPER-CCP4_LIBRARY} ${MMDB2_LIBRARY} Boost::thread Threads::Threads)

target_link_libraries(coot-findwaters PUBLIC cootapi ${CLIPPER-CORE_LIBRARY} ${CLIPPER-CCP4_LIBRARY} ${MMDB2_LIBRARY})

target_link_libraries(coot-mini-rsr PUBLIC cootapi ${CLIPPER-CORE_LIBRARY} ${CLIPPER-CCP4_LIBRARY} ${MMDB2_LIBRARY} Threads::Threads)

target_link_libraries(coot-res-tracer-libcootapi PUBLIC cootapi ${CLIPPER-CORE_LIBRARY} ${CLIPPER-CCP4_LIBRARY} ${MMDB2_LIBRARY} Threads::Threads)

target_include_directories(test-molecules-container PRIVATE ${coot_src} ${coot_src}/coords ${coot_src}/utils ${coot_src}/ccp4mg-utils ${coot_src}/ccp4mg-utils/util ${coot_src}/ccp4mg-utils/pygl ${MMDB2_INCLUDE_DIR}/include ${SSM_INCLUDE_DIR} ${CLIPPER-CORE_INCLUDE_DIR} ${CLIPPER-MMDB_INCLUDE_DIR} ${CLIPPER-CCP4_INCLUDE_DIR} ${FFTW2_INCLUDE_DIRS} ${RDKit_INCLUDE_DIRS} ${PNG_INCLUDE_DIR})

target_link_libraries(test-molecules-container PUBLIC cootapi ${SSM_LIBRARY} ${CLIPPER-MINIMOL_LIBRARY} ${CLIPPER-MMDB_LIBRARY} ${CLIPPER-CIF_LIBRARY} ${CLIPPER-CORE_LIBRARY} ${CLIPPER-CCP4_LIBRARY} ${CLIPPER-CONTRIB_LIBRARY} ${FFTW2_LIBRARY} ${RFFTW2_LIBRARY} ${MMDB2_LIBRARY} ${CCP4C_LIBRARY} ${ZLIB_LIBRARIES} ${GSL_LIBRARIES} ${RDKit_libs_for_libcootapi} ${PNG_LIBRARY} gemmi::gemmi_cpp ${libdwarf} ${libelf} Threads::Threads)

target_include_directories(mmcif-tests PRIVATE ${MMDB2_INCLUDE_DIR} ${CLIPPER-MMDB_INCLUDE_DIR})

target_link_libraries(mmcif-tests PUBLIC ${MMDB2_LIBRARY})

set_target_properties(cootapi PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

set_target_properties(cootapi PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

set_target_properties(coot_headless_api PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

set_target_properties(coot_headless_api PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

install(TARGETS cootmoleculestotriangles cootapi)

install(FILES standard-residues.pdb  DESTINATION share/coot)

install(FILES data/syminfo.lib       DESTINATION share/coot/data)

install(DIRECTORY rama-data          DESTINATION share/coot)

install(DIRECTORY data/metal         DESTINATION share/coot/data)

install(DIRECTORY data/cho-acedrg    DESTINATION share/coot/data)

install(DIRECTORY data/cho-models    DESTINATION share/coot/data)

install(DIRECTORY data/cho-links     DESTINATION share/coot/data)

install(DIRECTORY data/pdb-templates DESTINATION share/coot/data)

install(DIRECTORY ${coot_src}/analysis     DESTINATION include/coot/analysis     FILES_MATCHING PATTERN "*.hh" PATTERN "*.h")
install(DIRECTORY ${coot_src}/angles       DESTINATION include/coot/angles       FILES_MATCHING PATTERN "*.hh" PATTERN "*.h")
install(DIRECTORY ${coot_src}/api          DESTINATION include/coot/api          FILES_MATCHING PATTERN "*.hh" PATTERN "*.h")
install(DIRECTORY ${coot_src}/auxil        DESTINATION include/coot/auxil        FILES_MATCHING PATTERN "*.hh" PATTERN "*.h")
install(DIRECTORY ${coot_src}/build        DESTINATION include/coot/build        FILES_MATCHING PATTERN "*.hh" PATTERN "*.h")
install(DIRECTORY ${coot_src}/ccp4mg-utils DESTINATION include/coot/ccp4mg-utils FILES_MATCHING PATTERN "*.hh" PATTERN "*.h")
install(DIRECTORY ${coot_src}/compat       DESTINATION include/coot/compat       FILES_MATCHING PATTERN "*.hh" PATTERN "*.h")
install(DIRECTORY ${coot_src}/coot-utils   DESTINATION include/coot/coot-utils   FILES_MATCHING PATTERN "*.hh" PATTERN "*.h")
install(DIRECTORY ${coot_src}/cootaneer    DESTINATION include/coot/cootaneer    FILES_MATCHING PATTERN "*.hh" PATTERN "*.h")
install(DIRECTORY ${coot_src}/cootilus     DESTINATION include/coot/cootilus     FILES_MATCHING PATTERN "*.hh" PATTERN "*.h")
install(DIRECTORY ${coot_src}/coords       DESTINATION include/coot/coords       FILES_MATCHING PATTERN "*.hh" PATTERN "*.h")
install(DIRECTORY ${coot_src}/db-main      DESTINATION include/coot/db-main      FILES_MATCHING PATTERN "*.hh" PATTERN "*.h")
install(DIRECTORY ${coot_src}/geometry     DESTINATION include/coot/geometry     FILES_MATCHING PATTERN "*.hh" PATTERN "*.h")
install(DIRECTORY ${coot_src}/high-res     DESTINATION include/coot/high-res     FILES_MATCHING PATTERN "*.hh" PATTERN "*.h")
install(DIRECTORY ${coot_src}/ideal        DESTINATION include/coot/ideal        FILES_MATCHING PATTERN "*.hh" PATTERN "*.h")
install(DIRECTORY ${coot_src}/lidia-core   DESTINATION include/coot/lidia-core   FILES_MATCHING PATTERN "*.hh" PATTERN "*.h")
install(DIRECTORY ${coot_src}/ligand       DESTINATION include/coot/ligand       FILES_MATCHING PATTERN "*.hh" PATTERN "*.h")
install(DIRECTORY ${coot_src}/mini-mol     DESTINATION include/coot/mini-mol     FILES_MATCHING PATTERN "*.hh" PATTERN "*.h")
install(DIRECTORY ${coot_src}/pli          DESTINATION include/coot/pli          FILES_MATCHING PATTERN "*.hh" PATTERN "*.h")
install(DIRECTORY ${coot_src}/protein_db   DESTINATION include/coot/protein_db   FILES_MATCHING PATTERN "*.hh" PATTERN "*.h")
install(DIRECTORY ${coot_src}/pyrogen      DESTINATION include/coot/pyrogen      FILES_MATCHING PATTERN "*.hh" PATTERN "*.h")
install(DIRECTORY ${coot_src}/skeleton     DESTINATION include/coot/skeleton     FILES_MATCHING PATTERN "*.hh" PATTERN "*.h")
install(DIRECTORY ${coot_src}/utils        DESTINATION include/coot/utils        FILES_MATCHING PATTERN "*.hh" PATTERN "*.h")
install(DIRECTORY ${coot_src}/density-contour      DESTINATION include/coot/density-contour      FILES_MATCHING PATTERN "*.hh" PATTERN "*.h")
install(DIRECTORY ${coot_src}/MoleculesToTriangles DESTINATION include/coot/MoleculesToTriangles FILES_MATCHING PATTERN "*.hh" PATTERN "*.h")

install(DIRECTORY ${CMAKE_SOURCE_DIR}/monomers DESTINATION share/coot/lib/data
        FILES_MATCHING PATTERN "*" PATTERN "*Makefile.*" EXCLUDE)

# Thanks Charles B.
install(TARGETS coot_headless_api DESTINATION ${Python_SITELIB})

install(TARGETS coot_headless_api DESTINATION lib)

# thanks merkys
set_target_properties(cootapi PROPERTIES SOVERSION 1.1)
set_target_properties(cootmoleculestotriangles PROPERTIES SOVERSION 1.1)
