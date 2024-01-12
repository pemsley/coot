cmake_policy(SET CMP0048 NEW)
project(coot)
cmake_minimum_required(VERSION 3.13)
find_package(PkgConfig)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(CMAKE_CXX_STANDARD 17)

find_package(CCP4 REQUIRED mmdb2 ssm clipper-core clipper-mmdb clipper-ccp4 clipper-contrib clipper-minimol clipper-cif)

find_package(BLAS)
find_package(GSL REQUIRED)
find_package(PNG REQUIRED)

# from Charles B.

# The problem with the RDKit package is that it requires Eigen, which is not in the Coot build system
# OK, well I will add it now and then have to change the build script later.
# 20231023-PE I have added Eigen3 to the build-it-3-3 script.
find_package(RDKit CONFIG COMPONENTS RDGeneral REQUIRED)
find_package(Boost COMPONENTS iostreams system python regex thread serialization REQUIRED)
find_package(glm REQUIRED)


# Charles wants 3.9 here, but if I change to 3.9 then it will not build for me
# Hmm...
find_package(Python 3.11
  REQUIRED COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule)

find_package(nanobind CONFIG REQUIRED)

set(coot_src ${CMAKE_CURRENT_SOURCE_DIR})

add_definitions(-DANALYSE_REFINEMENT_TIMING -DVERSION="1.x-pre" -DPACKAGE_VERSION="1"
  -DHAVE_CXX_THREAD -DHAVE_BOOST -DHAVE_BOOST_THREAD -DHAVE_GSL
  -DPKGDATADIR="${CMAKE_INSTALL_PREFIX}/share/coot" -DXDATADIR="${CMAKE_INSTALL_PREFIX}/share")

#  pyrogen, compat, auxil, goograph, src ignored.

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
            ${coot_src}/pli/pi-stacking.cc
            ${coot_src}/coot-utils/cmtz-interface.cc
            ${coot_src}/coot-utils/lidia-core-functions.cc
            ${coot_src}/coot-utils/atom-selection-container.cc
            #${coot_src}/coot-utils/dictionary-atom-types.cc
            #${coot_src}/coot-utils/cod-db-access.cc
            ${coot_src}/src/coot-hydrogens.cc
            ${coot_src}/src/cfc.cc
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
            # ${coot_src}/MoleculesToTriangles/CXXClasses/RendererGL.cpp
            # ${coot_src}/MoleculesToTriangles/CXXClasses/RendererGLSL.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/Representation.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/RepresentationInstance.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/RotatedTranslatedScaledEntity.cpp
            # ${coot_src}/MoleculesToTriangles/CXXClasses/SceneSetup.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/SecondaryColorScheme.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/SolidColorRule.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/SpherePrimitive.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/SticksPrimitive.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/SurfacePrimitive.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/Vector2.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/Vector3.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/VertexColorNormalPrimitive.cpp
            ${coot_src}/MoleculesToTriangles/CXXClasses/lfit.c
            ${coot_src}/api/rigid-body-fit.cc
            ${coot_src}/api//molecules_container_bucca_ml_growing.cc
            ${coot_src}/api/molecules_container_mmrrcc.cc
            ${coot_src}/api/molecules_container_superpose.cc
            ${coot_src}/api/coot_molecule_validation.cc
            ${coot_src}/api/molecules_container_modelling.cc
            ${coot_src}/api/add-terminal-residue.cc
            ${coot_src}/api/interfaces.cc
            ${coot_src}/api/colour-functions.cc
            ${coot_src}/api/coot-colour.cc
            ${coot_src}/api/coot_molecule.cc
            ${coot_src}/api/coot_molecule_backup.cc
            ${coot_src}/api/coot_molecule_bonds.cc
            ${coot_src}/api/coot_molecule_bonds_instanced.cc
            ${coot_src}/api/coot_molecule_maps.cc
            ${coot_src}/api/moorhen-h-bonds.cc
            ${coot_src}/api/coot_molecule_moltris.cc
            ${coot_src}/api/phi-psi-prob.cc
            ${coot_src}/api/molecules_container_ligand_fitting.cc
            ${coot_src}/api/coot_molecule_merge_molecules.cc
            ${coot_src}/api/coot_molecule_refine.cc
            ${coot_src}/api/coot_molecule_replace_fragment.cc
            ${coot_src}/api/add-terminal-residue.cc
            ${coot_src}/api/molecules_container.cc
            ${coot_src}/api/molecules_container_bucca_ml_growing.cc
            ${coot_src}/api/rama-plot-phi-psi.cc
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
            ${coot_src}/density-contour/gaussian-surface.cc
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
            #${coot_src}/ideal/mini-rsr.cc
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
            ${coot_src}/tw/tw.cc
            ${coot_src}/utils/colour-holder.cc
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
            #${coot_src}/protein_db/protein_db.cpp
            #${coot_src}/protein_db/protein_db_utils.cpp
            #${coot_src}/protein_db/rna_db.cpp
            #${coot_src}/pli/specs.cc
            ${coot_src}/pli/flev-annotations.cc
            ${coot_src}/pli/protein-ligand-interactions.cc
            ${coot_src}/pli/pi-stacking.cc
            ${coot_src}/pli/sdf-interface.cc
            ${coot_src}/coot-utils/merge-atom-selections.cc
            ${coot_src}/coot-utils/coot-least-squares.cc
            ${coot_src}/coot-utils/diff-diff-map-peaks.cc
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
            ${coot_src}/coot-utils/exp-fit.cc
            ${coot_src}/coot-utils/cylinder.cc
            ${coot_src}/coot-utils/cylinder-with-rotation-translation.cc
            ${coot_src}/coot-utils/prideout-octasphere.cc
            ${coot_src}/coot-utils/oct.cc
            ${coot_src}/coot-utils/shapes.cc
            ${coot_src}/coot-utils/simple-mesh.cc
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
            #${coot_src}/build/testcabuild.cc
            ${coot_src}/mini-mol/mini-mol-utils.cc
            ${coot_src}/mini-mol/mini-mol.cc
            ${coot_src}/mini-mol/atom-quads.cc
            ${coot_src}/cootilus/nautilus-tools.cpp
            ${coot_src}/cootilus/cootilus-build.cpp
            ${coot_src}/cootilus/nucleicacid_db.cpp
            ${coot_src}/cootilus/nautilus-join.cpp
            ${coot_src}/cootilus/nautilus-target.cpp
            #${coot_src}/cootilus/nautilus-ss-find.cpp
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
            #${coot_src}/ligand/trace-2.cc
            ${coot_src}/ligand/monomer-utils.cc
            ${coot_src}/ligand/ligand-score-card.cc
            #${coot_src}/ligand/multi-peptide.cc
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
            #${coot_src}/ligand/trace.cc
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
            ${coot_src}/src/gl-matrix.cc
            ${coot_src}/ideal/write-png.cc
            ${coot_src}/buccaneer_ml_growing/buccaneer-NN-Select.cpp
            ${coot_src}/buccaneer_ml_growing/buccaneer-NN-build-cyc.cpp
            ${coot_src}/buccaneer_ml_growing/buccaneer-NN-build.cpp
            ${coot_src}/buccaneer_ml_growing/buccaneer-NN-features.cpp
            ${coot_src}/buccaneer_ml_growing/buccaneer-NN-join.cpp
            ${coot_src}/buccaneer_ml_growing/buccaneer-NN-model.cpp
            ${coot_src}/buccaneer_ml_growing/buccaneer-build.cpp
            ${coot_src}/buccaneer_ml_growing/buccaneer-correct.cpp
            ${coot_src}/buccaneer_ml_growing/buccaneer-filter.cpp
            ${coot_src}/buccaneer_ml_growing/buccaneer-find.cpp
            ${coot_src}/buccaneer_ml_growing/buccaneer-grow.cpp
            #${coot_src}/buccaneer_ml_growing/buccaneer-join.cpp
            ${coot_src}/buccaneer_ml_growing/buccaneer-known.cpp
            ${coot_src}/buccaneer_ml_growing/buccaneer-lib.cpp
            ${coot_src}/buccaneer_ml_growing/buccaneer-link.cpp
            ${coot_src}/buccaneer_ml_growing/buccaneer-merge.cpp
            ${coot_src}/buccaneer_ml_growing/buccaneer-ncsbuild.cpp
            ${coot_src}/buccaneer_ml_growing/buccaneer-prep.cpp
            ${coot_src}/buccaneer_ml_growing/buccaneer-prot.cpp
            ${coot_src}/buccaneer_ml_growing/buccaneer-prune.cpp
            ${coot_src}/buccaneer_ml_growing/buccaneer-sequence.cpp
            ${coot_src}/buccaneer_ml_growing/buccaneer-tidy.cpp
            ${coot_src}/buccaneer_ml_growing/buccaneer-util.cpp
            #${coot_src}/buccaneer_ml_growing/cbuccaneer.cpp
            #${coot_src}/buccaneer_ml_growing/cbuccaneerNN.cpp
            ${coot_src}/buccaneer_ml_growing/csequins.cpp
            ${coot_src}/buccaneer_ml_growing/k2c_activations.cpp
            ${coot_src}/buccaneer_ml_growing/k2c_convolution_layers.cpp
            ${coot_src}/buccaneer_ml_growing/k2c_core_layers.cpp
            ${coot_src}/buccaneer_ml_growing/k2c_embedding_layers.cpp
            ${coot_src}/buccaneer_ml_growing/k2c_helper_functions.cpp
            ${coot_src}/buccaneer_ml_growing/k2c_merge_layers.cpp
            ${coot_src}/buccaneer_ml_growing/k2c_normalization_layers.cpp
            ${coot_src}/buccaneer_ml_growing/k2c_pooling_layers.cpp
            ${coot_src}/buccaneer_ml_growing/k2c_recurrent_layers.cpp
            ${coot_src}/buccaneer_ml_growing/ml-grow-next.cpp
            ${coot_src}/buccaneer_ml_growing/ml-grow-prev.cpp
            ${coot_src}/buccaneer_ml_growing/ml-grow.cpp
            ${coot_src}/buccaneer_ml_growing/protein_db.cpp
            ${coot_src}/buccaneer_ml_growing/simplex-lib.cpp
            ${coot_src}/buccaneer_ml_growing/simulate-lib.cpp
)

# use MMDB2_INCLUDE_DIR to find the RDKIT for now.
target_include_directories(cootapi PRIVATE ${coot_src} ${coot_src}/coords ${coot_src}/utils ${coot_src}/ccp4mg-utils ${coot_src}/ccp4mg-utils/util ${coot_src}/ccp4mg-utils/pygl ${MMDB2_INCLUDE_DIR} ${SSM_INCLUDE_DIR} ${CLIPPER-CORE_INCLUDE_DIR} ${CLIPPER-MMDB_INCLUDE_DIR} ${CLIPPER-CCP4_INCLUDE_DIR} ${FFTW2_INCLUDE_DIRS} ${RDKit_INCLUDE_DIRS} ${GLM_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS} ${PNG_PNG_INCLUDE_DIR})

# Charles adds Boost::iostreams Boost::filesystem Boost::system. I don't know why Boost::filesystem was added - I have removed it.
target_link_libraries(cootapi PUBLIC ${SSM_LIBRARY} ${CLIPPER-MINIMOL_LIBRARY} ${CLIPPER-MMDB_LIBRARY} ${CLIPPER-CIF_LIBRARY} ${CLIPPER-CORE_LIBRARY} ${CLIPPER-CCP4_LIBRARY} ${CLIPPER-CONTRIB_LIBRARY} ${FFTW2_LIBRARY} ${RFFTW2_LIBRARY} ${MMDB2_LIBRARY} ${CCP4C_LIBRARY} ${ZLIB_LIBRARIES} ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARY} ${PNG_LIBRARY} Boost::iostreams Boost::system)

set_target_properties(cootapi PROPERTIES COMPILE_FLAGS "-DHAVE_SSMLIB -DHAVE_BOOST -DUSE_LIBPNG=1 -DCLIPPER_HAS_TOP8000")

# Charles adds this MSVC CXX compiler test
if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  target_compile_options(cootapi PRIVATE
	/D_CRT_SECURE_NO_WARNINGS
	/D_SILENCE_ALL_CXX17_DEPRECATION_WARNINGS
	/bigobj
  )
  target_link_libraries(cootapi PRIVATE netapi32.lib)
endif()

nanobind_add_module(chapi ${coot_src}/api/molecules_container_nanobind.cc)

target_include_directories(chapi PRIVATE ${coot_src} ${coot_src}/coords ${coot_src}/utils ${coot_src}/ccp4mg-utils ${coot_src}/ccp4mg-utils/util ${coot_src}/ccp4mg-utils/pygl ${MMDB2_INCLUDE_DIR} ${SSM_INCLUDE_DIR} ${CLIPPER-CORE_INCLUDE_DIR} ${CLIPPER-MMDB_INCLUDE_DIR} ${CLIPPER-CCP4_INCLUDE_DIR} ${FFTW2_INCLUDE_DIRS} ${RDKit_INCLUDE_DIRS} ${GLM_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS} ${PNG_PNG_INCLUDE_DIR})

target_link_libraries(chapi PUBLIC cootapi ${SSM_LIBRARY} ${CLIPPER-MINIMOL_LIBRARY} ${CLIPPER-MMDB_LIBRARY} ${CLIPPER-CIF_LIBRARY} ${CLIPPER-CORE_LIBRARY} ${CLIPPER-CCP4_LIBRARY} ${CLIPPER-CONTRIB_LIBRARY} ${FFTW2_LIBRARY} ${RFFTW2_LIBRARY} ${MMDB2_LIBRARY} ${CCP4C_LIBRARY} ${ZLIB_LIBRARIES} ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARY} ${PNG_LIBRARY})

set_target_properties(chapi PROPERTIES COMPILE_FLAGS "-DHAVE_SSMLIB -DHAVE_BOOST -DUSE_LIBPNG=1 -DCLIPPER_HAS_TOP8000")

# I think Charles is right. *Coot* should have test-molecules-container, but this CMakeLists.txt file
# is for python chapi, it shouldn't be making executables.
# add_executable(test-molecules-container ${coot_src}/api/test_molecules_container.cc ${coot_src}/api/filo-tests.cc)

# target_include_directories(test-molecules-container PRIVATE ${coot_src} ${coot_src}/coords ${coot_src}/utils ${coot_src}/ccp4mg-utils ${coot_src}/ccp4mg-utils/util ${coot_src}/ccp4mg-utils/pygl ${MMDB2_INCLUDE_DIR}/include ${SSM_INCLUDE_DIR} ${CLIPPER-CORE_INCLUDE_DIR} ${CLIPPER-MMDB_INCLUDE_DIR} ${CLIPPER-CCP4_INCLUDE_DIR} ${FFTW2_INCLUDE_DIRS})

# target_link_libraries(test-molecules-container PUBLIC cootapi ${SSM_LIBRARY} ${CLIPPER-MINIMOL_LIBRARY} ${CLIPPER-MMDB_LIBRARY} ${CLIPPER-CIF_LIBRARY} ${CLIPPER-CORE_LIBRARY} ${CLIPPER-CCP4_LIBRARY} ${CLIPPER-CONTRIB_LIBRARY} ${FFTW2_LIBRARY} ${RFFTW2_LIBRARY} ${MMDB2_LIBRARY} ${CCP4C_LIBRARY} ${ZLIB_LIBRARIES} ${GSL_LIBRARIES} ${GSLBLAS_LIBRARY} ${PNG_LIBRARY})

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/syminfo.lib
DESTINATION ${CMAKE_INSTALL_PREFIX}/data)

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/standard-residues.pdb
DESTINATION ${CMAKE_INSTALL_PREFIX}/share/coot)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/pdb-templates/
DESTINATION ${CMAKE_INSTALL_PREFIX}/share/coot)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/rama-data/
DESTINATION ${CMAKE_INSTALL_PREFIX}/share/coot/rama-data)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data/metal/
DESTINATION ${CMAKE_INSTALL_PREFIX}/share/coot/data/metal)
