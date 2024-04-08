#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/vector.h>
#include <nanobind/ndarray.h>
#include <nanobind/operators.h>

#include "molecules-container.hh"
#include "mini-mol/mini-mol-utils.hh"

#include "clipper/core/ramachandran.h"
#include "clipper/clipper-ccp4.h"

#include "coot-utils/g_triangle.hh"
#include "coords/mmdb-crystal.h"

namespace nb = nanobind;

struct RamachandranInfo {
    std::string chainId;
    int seqNum;
    std::string insCode;
    std::string restype;
    double phi;
    double psi;
    bool isOutlier;
    bool is_pre_pro;
};

struct ResiduePropertyInfo {
    std::string chainId;
    int seqNum;
    std::string insCode;
    std::string restype;
    double property;
};

class molecules_container_js : public molecules_container_t {
    public:
        explicit molecules_container_js(bool verbose=true) : molecules_container_t(verbose) {
        }

        int writePDBASCII(int imol, const std::string &file_name) {
            const char *fname_cp = file_name.c_str();
            return get_mol(imol)->WritePDBASCII(fname_cp);
        }
        int writeCIFASCII(int imol, const std::string &file_name) {
            const char *fname_cp = file_name.c_str();
            return get_mol(imol)->WriteCIFASCII(fname_cp);
        }
        int writeCCP4Map(int imol, const std::string &file_name) {
            auto xMap = (*this)[imol].xmap;
            auto clipperMap = clipper::CCP4MAPfile();
            clipperMap.open_write(file_name);
            clipperMap.export_xmap(xMap);
            return 0;
        }
};

NB_MODULE(chapi, m) {
    nb::class_<clipper::Coord_orth>(m,"Coord_orth")
    .def(nb::init<const clipper::ftype&, const clipper::ftype&, const clipper::ftype&>())
    .def("x", &clipper::Coord_orth::x)
    .def("y", &clipper::Coord_orth::y)
    .def("z", &clipper::Coord_orth::z)
    ;
    nb::class_<coot::util::map_molecule_centre_info_t>(m,"map_molecule_centre_info_t")
    .def_ro("success", &coot::util::map_molecule_centre_info_t::success)
    .def_ro("updated_centre", &coot::util::map_molecule_centre_info_t::updated_centre)
    .def_ro("suggested_contour_level", &coot::util::map_molecule_centre_info_t::suggested_contour_level)
    ;
    nb::class_<clipper::Cell_descr>(m,"Cell_descr")
    .def(nb::init<const clipper::ftype&, const clipper::ftype&, const clipper::ftype&, const clipper::ftype&, const clipper::ftype&, const clipper::ftype&>())
    .def("a",         &clipper::Cell_descr::a)
    .def("b",         &clipper::Cell_descr::b)
    .def("c",         &clipper::Cell_descr::c)
    .def("alpha",     &clipper::Cell_descr::alpha)
    .def("beta",      &clipper::Cell_descr::beta)
    .def("gamma",     &clipper::Cell_descr::gamma)
    .def("alpha_deg", &clipper::Cell_descr::alpha_deg)
    .def("beta_deg",  &clipper::Cell_descr::beta_deg)
    .def("gamma_deg", &clipper::Cell_descr::gamma_deg)
    .def("format",    &clipper::Cell_descr::format)
    ;
    nb::class_<clipper::Cell, clipper::Cell_descr>(m,"Cell")
    .def(nb::init<>())
    .def(nb::init<const clipper::Cell_descr &>())
    .def("a_star", &clipper::Cell::a_star)
    .def("b_star", &clipper::Cell::b_star)
    .def("c_star", &clipper::Cell::c_star)
    .def("alpha_star", &clipper::Cell::alpha_star)
    .def("beta_star", &clipper::Cell::beta_star)
    .def("gamma_star", &clipper::Cell::gamma_star)
    .def("descr", &clipper::Cell::descr)
    .def("is_null", &clipper::Cell::is_null)
    .def("init", &clipper::Cell::init)
    ;
    nb::class_<clipper::Xmap_base>(m,"Xmap_base")
    .def("cell", &clipper::Xmap_base::cell)
    ;
    nb::class_<clipper::String>(m,"Clipper_String")
    .def(nb::init<>())
    .def(nb::init<const std::string>())
    ;
    nb::class_<clipper::Xmap<float>, clipper::Xmap_base>(m,"Xmap_float")
    .def(nb::init<>())
    ;
    nb::class_<clipper::CCP4MAPfile>(m,"CCP4MAPfile")
    .def(nb::init<>())
    .def("open_read",&clipper::CCP4MAPfile::open_read)
    .def("open_write",&clipper::CCP4MAPfile::open_write)
    .def("close_read",&clipper::CCP4MAPfile::close_read)
    .def("close_write",&clipper::CCP4MAPfile::close_write)
    ;
    nb::class_<mmdb::Atom>(m,"Atom")
    .def(nb::init<>())
    .def_prop_rw("x",[](mmdb::Atom &t) { return t.x ; },[](mmdb::Atom &t, float value) { t.x = value; })
    .def_prop_rw("y",[](mmdb::Atom &t) { return t.y ; },[](mmdb::Atom &t, float value) { t.y = value; })
    .def_prop_rw("z",[](mmdb::Atom &t) { return t.z ; },[](mmdb::Atom &t, float value) { t.z = value; })
    .def_prop_rw("serNum",[](mmdb::Atom &t) { return t.serNum ; },[](mmdb::Atom &t, float value) { t.serNum = value; })
    .def_prop_rw("occupancy",[](mmdb::Atom &t) { return t.occupancy ; },[](mmdb::Atom &t, float value) { t.occupancy = value; })
    .def_prop_rw("tempFactor",[](mmdb::Atom &t) { return t.tempFactor ; },[](mmdb::Atom &t, float value) { t.tempFactor = value; })
    .def_prop_rw("charge",[](mmdb::Atom &t) { return t.charge ; },[](mmdb::Atom &t, float value) { t.charge = value; })
    .def_prop_rw("sigX",[](mmdb::Atom &t) { return t.sigX ; },[](mmdb::Atom &t, float value) { t.sigX = value; })
    .def_prop_rw("sigY",[](mmdb::Atom &t) { return t.sigY ; },[](mmdb::Atom &t, float value) { t.sigY = value; })
    .def_prop_rw("sigZ",[](mmdb::Atom &t) { return t.sigZ ; },[](mmdb::Atom &t, float value) { t.sigZ = value; })
    .def_prop_rw("sigOcc",[](mmdb::Atom &t) { return t.sigOcc ; },[](mmdb::Atom &t, float value) { t.sigOcc = value; })
    .def_prop_rw("sigTemp",[](mmdb::Atom &t) { return t.sigTemp ; },[](mmdb::Atom &t, float value) { t.sigTemp = value; })
    .def_prop_rw("u11",[](mmdb::Atom &t) { return t.u11 ; },[](mmdb::Atom &t, float value) { t.u11 = value; })
    .def_prop_rw("u22",[](mmdb::Atom &t) { return t.u22 ; },[](mmdb::Atom &t, float value) { t.u22 = value; })
    .def_prop_rw("u33",[](mmdb::Atom &t) { return t.u33 ; },[](mmdb::Atom &t, float value) { t.u33 = value; })
    .def_prop_rw("u13",[](mmdb::Atom &t) { return t.u13 ; },[](mmdb::Atom &t, float value) { t.u13 = value; })
    .def_prop_rw("u23",[](mmdb::Atom &t) { return t.u23 ; },[](mmdb::Atom &t, float value) { t.u23 = value; })
    .def_prop_rw("Het",[](mmdb::Atom &t) { return t.Het ; },[](mmdb::Atom &t, bool value) { t.Het = value; })
    .def_prop_rw("Ter",[](mmdb::Atom &t) { return t.Ter ; },[](mmdb::Atom &t, bool value) { t.Ter = value; })
    .def("GetNBonds",&mmdb::Atom::GetNBonds)
    .def("GetModelNum",&mmdb::Atom::GetModelNum)
    .def("GetSeqNum",&mmdb::Atom::GetSeqNum)
    .def("GetLabelSeqID",&mmdb::Atom::GetLabelSeqID)
    .def("GetLabelEntityID",&mmdb::Atom::GetLabelEntityID)
    .def("GetSSEType",&mmdb::Atom::GetSSEType)
    .def("isTer",&mmdb::Atom::isTer)
    .def("isMetal",&mmdb::Atom::isMetal)
    .def("isSolvent",&mmdb::Atom::isSolvent)
    .def("isInSelection",&mmdb::Atom::isInSelection)
    .def("isNTerminus",&mmdb::Atom::isNTerminus)
    .def("isCTerminus",&mmdb::Atom::isCTerminus)
    .def("GetResidueNo",&mmdb::Atom::GetResidueNo)
    .def("GetIndex",&mmdb::Atom::GetIndex)
    .def("GetAtomName",&mmdb::Atom::GetAtomName)
    .def("SetAtomName",nb::overload_cast<const char*>(&mmdb::Atom::SetAtomName))
    .def("GetChainID",&mmdb::Atom::GetChainID)
    .def("GetLabelAsymID",&mmdb::Atom::GetLabelAsymID)
    .def("GetLabelCompID",&mmdb::Atom::GetLabelCompID)
    .def("GetInsCode",&mmdb::Atom::GetInsCode)
    ;
    nb::class_<mmdb::Residue>(m,"Residue")
    .def(nb::init<>())
    .def_prop_rw("seqNum",[](mmdb::Residue &t) { return t.seqNum ; },[](mmdb::Residue &t, int value) { t.seqNum = value; })
    .def_prop_rw("label_seq_id",[](mmdb::Residue &t) { return t.label_seq_id ; },[](mmdb::Residue &t, int value) { t.label_seq_id = value; })
    .def_prop_rw("label_entity_id",[](mmdb::Residue &t) { return t.label_entity_id ; },[](mmdb::Residue &t, int value) { t.label_entity_id = value; })
    .def_prop_rw("index",[](mmdb::Residue &t) { return t.index ; },[](mmdb::Residue &t, int value) { t.index = value; })
    .def_prop_rw("nAtoms",[](mmdb::Residue &t) { return t.nAtoms ; },[](mmdb::Residue &t, int value) { t.nAtoms = value; })
    .def("GetModelNum",&mmdb::Residue::GetModelNum)
    .def("GetSeqNum",&mmdb::Residue::GetSeqNum)
    .def("GetLabelSeqID",&mmdb::Residue::GetLabelSeqID)
    .def("GetLabelEntityID",&mmdb::Residue::GetLabelEntityID)
    .def("GetResidueNo",&mmdb::Residue::GetResidueNo)
    .def("GetNofAltLocations",&mmdb::Residue::GetNofAltLocations)
    .def("isAminoacid",&mmdb::Residue::isAminoacid)
    .def("isNucleotide",&mmdb::Residue::isNucleotide)
    .def("isDNARNA",&mmdb::Residue::isDNARNA)
    .def("isSugar",&mmdb::Residue::isSugar)
    .def("isSolvent",&mmdb::Residue::isSolvent)
    .def("isModRes",&mmdb::Residue::isModRes)
    .def("isInSelection",&mmdb::Residue::isInSelection)
    .def("isNTerminus",&mmdb::Residue::isNTerminus)
    .def("isCTerminus",&mmdb::Residue::isCTerminus)
    .def("GetResName",&mmdb::Residue::GetResName)
    .def("GetChainID",&mmdb::Residue::GetChainID)
    .def("GetLabelAsymID",&mmdb::Residue::GetLabelAsymID)
    .def("GetLabelCompID",&mmdb::Residue::GetLabelCompID)
    .def("GetInsCode",&mmdb::Residue::GetInsCode)
    .def("GetAtom", nb::overload_cast<int>(&mmdb::Residue::GetAtom), nb::rv_policy::reference)
    .def("GetNumberOfAtoms", nb::overload_cast<>(&mmdb::Residue::GetNumberOfAtoms))
    .def("GetNumberOfAtoms_countTers", nb::overload_cast<bool>(&mmdb::Residue::GetNumberOfAtoms))
    ;
    nb::class_<molecules_container_t>(m,"molecules_container_t")
    .def(nb::init<bool>())
    .def("M2T_updateFloatParameter",&molecules_container_t::M2T_updateFloatParameter)
    .def("SSM_superpose",&molecules_container_t::SSM_superpose)
    .def("add_alternative_conformation",&molecules_container_t::add_alternative_conformation)
    .def("add_colour_rule",&molecules_container_t::add_colour_rule)
    .def("add_colour_rules_multi",&molecules_container_t::add_colour_rules_multi)
    .def("add_hydrogen_atoms",&molecules_container_t::add_hydrogen_atoms)
    .def("add_target_position_restraint",&molecules_container_t::add_target_position_restraint)
    .def("add_target_position_restraint_and_refine",&molecules_container_t::add_target_position_restraint_and_refine)
    .def("add_terminal_residue_directly_using_cid",&molecules_container_t::add_terminal_residue_directly_using_cid)
    .def("add_to_non_drawn_bonds",&molecules_container_t::add_to_non_drawn_bonds)
    .def("add_waters",&molecules_container_t::add_waters)
    .def("all_molecule_contact_dots",&molecules_container_t::all_molecule_contact_dots)
    .def("apply_transformation_to_atom_selection",&molecules_container_t::apply_transformation_to_atom_selection)
    .def("average_map",&molecules_container_t::average_map)
    .def("associate_data_mtz_file_with_map",&molecules_container_t::associate_data_mtz_file_with_map)
    .def("auto_fit_rotamer",&molecules_container_t::auto_fit_rotamer)
    .def("auto_read_mtz",&molecules_container_t::auto_read_mtz)
    .def("calculate_new_rail_points",&molecules_container_t::calculate_new_rail_points)
    .def("change_to_first_rotamer",&molecules_container_t::change_to_first_rotamer)
    .def("change_to_next_rotamer",&molecules_container_t::change_to_next_rotamer)
    .def("change_to_previous_rotamer",&molecules_container_t::change_to_previous_rotamer)
    .def("cis_trans_convert",&molecules_container_t::cis_trans_convert)
    .def("clear_extra_restraints",&molecules_container_t::clear_extra_restraints)
    .def("clear_non_drawn_bonds",&molecules_container_t::clear_non_drawn_bonds)
    .def("clear_refinement",&molecules_container_t::clear_refinement)
    .def("clear_target_position_restraints",&molecules_container_t::clear_target_position_restraints)
    .def("close_molecule",&molecules_container_t::close_molecule)
    .def("connect_updating_maps",&molecules_container_t::connect_updating_maps)
    .def("contact_dots_for_ligand",&molecules_container_t::contact_dots_for_ligand)
    .def("copy_fragment_for_refinement_using_cid",&molecules_container_t::copy_fragment_for_refinement_using_cid)
    .def("copy_fragment_using_cid",&molecules_container_t::copy_fragment_using_cid)
    .def("copy_fragment_using_residue_range",&molecules_container_t::copy_fragment_using_residue_range)
    .def("delete_atom",&molecules_container_t::delete_atom)
    .def("delete_atom_using_cid",&molecules_container_t::delete_atom_using_cid)
    .def("delete_chain_using_cid",&molecules_container_t::delete_chain_using_cid)
    .def("delete_colour_rules",&molecules_container_t::delete_colour_rules)
    .def("delete_hydrogen_atoms",&molecules_container_t::delete_hydrogen_atoms)
    .def("delete_residue",&molecules_container_t::delete_residue)
    .def("delete_residue_atoms_using_cid",&molecules_container_t::delete_residue_atoms_using_cid)
    .def("delete_residue_atoms_with_alt_conf",&molecules_container_t::delete_residue_atoms_with_alt_conf)
    .def("delete_residue_using_cid",&molecules_container_t::delete_residue_using_cid)
    .def("delete_side_chain",&molecules_container_t::delete_side_chain)
    .def("delete_using_cid",&molecules_container_t::delete_using_cid)
    .def("density_correlation_analysis",&molecules_container_t::density_correlation_analysis)
    .def("density_fit_analysis",&molecules_container_t::density_fit_analysis)
    .def("difference_map_peaks",&molecules_container_t::difference_map_peaks)
    .def("eigen_flip_ligand", nb::overload_cast<int, const std::string&>                        (&molecules_container_t::eigen_flip_ligand_using_cid))
    .def("file_name_to_string",&molecules_container_t::file_name_to_string)
    .def("fill_partial_residue",&molecules_container_t::fill_partial_residue)
    .def("fill_rotamer_probability_tables",&molecules_container_t::fill_rotamer_probability_tables)
    .def("fit_ligand",&molecules_container_t::fit_ligand)
    .def("fit_ligand_right_here",&molecules_container_t::fit_ligand_right_here)
    .def("fit_to_map_by_random_jiggle",&molecules_container_t::fit_to_map_by_random_jiggle)
    .def("fit_to_map_by_random_jiggle_using_cid",&molecules_container_t::fit_to_map_by_random_jiggle_using_cid)
    .def("fit_to_map_by_random_jiggle_using_cid",&molecules_container_t::fit_to_map_by_random_jiggle_using_cid)
    .def("flipPeptide",       nb::overload_cast<int, const coot::atom_spec_t&,const std::string&>(&molecules_container_t::flip_peptide))
    .def("flipPeptide_cid",   nb::overload_cast<int, const std::string&,      const std::string&>(&molecules_container_t::flip_peptide_using_cid))
    .def("flip_hand",&molecules_container_t::flip_hand)
    .def("fourier_shell_correlation",&molecules_container_t::fourier_shell_correlation)
    .def("generate_self_restraints",&molecules_container_t::generate_self_restraints)
    .def("geometry_init_standard",&molecules_container_t::geometry_init_standard)
    .def("get_active_atom",&molecules_container_t::get_active_atom)
    .def("get_atom",&molecules_container_t::get_atom, nb::rv_policy::reference)
    .def("get_bonds_mesh",&molecules_container_t::get_bonds_mesh)
    .def("get_bonds_mesh_for_selection_instanced",&molecules_container_t::get_bonds_mesh_for_selection_instanced)
    .def("get_bonds_mesh_instanced",&molecules_container_t::get_bonds_mesh_instanced)
    .def("get_cell",&molecules_container_t::get_cell)
    .def("get_chains_in_model",&molecules_container_t::get_chains_in_model)
    .def("get_chemical_features_mesh",&molecules_container_t::get_chemical_features_mesh)
    .def("get_colour_rules",&molecules_container_t::get_colour_rules)
    .def("get_colour_table_for_blender", &molecules_container_t::get_colour_table_for_blender)
    .def("get_density_at_position", &molecules_container_t::get_density_at_position)
    .def("get_gaussian_surface",&molecules_container_t::get_gaussian_surface)
    .def("get_goodsell_style_mesh_instanced",&molecules_container_t::get_goodsell_style_mesh_instanced)
    .def("get_gphl_chem_comp_info",&molecules_container_t::get_gphl_chem_comp_info)
    .def("get_group_for_monomer",&molecules_container_t::get_group_for_monomer)
    .def("get_groups_for_monomers",&molecules_container_t::get_groups_for_monomers)
    .def("get_hb_type",&molecules_container_t::get_hb_type)
    // .def("get_h_bonds",&molecules_container_t::get_h_bonds)
    // .def("get_interesting_places",&molecules_container_t::get_interesting_places)
    .def("get_map_contours_mesh",&molecules_container_t::get_map_contours_mesh)
    .def("get_map_molecule_centre",&molecules_container_t::get_map_molecule_centre)
    .def("get_map_rmsd_approx",&molecules_container_t::get_map_rmsd_approx)
    .def("get_map_weight",&molecules_container_t::get_map_weight)
    .def("get_molecular_representation_mesh",&molecules_container_t::get_molecular_representation_mesh)
    .def("get_molecule_centre",&molecules_container_t::get_molecule_centre)
    .def("get_molecule_name",&molecules_container_t::get_molecule_name)
    .def("get_monomer",&molecules_container_t::get_monomer)
    .def("get_monomer_and_position_at",&molecules_container_t::get_monomer_and_position_at)
    .def("get_monomer_from_dictionary",&molecules_container_t::get_monomer_from_dictionary)
    .def("get_number_of_molecules",&molecules_container_t::get_number_of_molecules)
    .def("get_r_factor_stats",&molecules_container_t::get_r_factor_stats)
    .def("get_rama_plot_restraints_weight",&molecules_container_t::get_rama_plot_restraints_weight)
    .def("get_ramachandran_validation_markup_mesh",&molecules_container_t::get_ramachandran_validation_markup_mesh)
    //Using allow_raw_pointers(). Perhaps suggests we need to do something different from exposing mmdb pointers to JS.
    .def("get_residue",&molecules_container_t::get_residue, nb::rv_policy::reference)
    .def("get_residue_name",&molecules_container_t::get_residue_name)
    .def("get_residue_names_with_no_dictionary",&molecules_container_t::get_residue_names_with_no_dictionary)
    .def("get_residue_using_cid",&molecules_container_t::get_residue_using_cid)
    .def("get_residues_near_residue",&molecules_container_t::get_residues_near_residue)
    .def("get_rotamer_dodecs",&molecules_container_t::get_rotamer_dodecs)
    .def("get_rotamer_dodecs_instanced",&molecules_container_t::get_rotamer_dodecs_instanced)
    .def("get_single_letter_codes_for_chain",&molecules_container_t::get_single_letter_codes_for_chain)
    .def("get_svg_for_residue_type",&molecules_container_t::get_svg_for_residue_type)
    .def("get_symmetry",&molecules_container_t::get_symmetry)
    .def("get_torsion_restraints_weight",&molecules_container_t::get_torsion_restraints_weight)
    .def("get_triangles_for_blender", &molecules_container_t::get_triangles_for_blender)
    .def("get_use_gemmi",&molecules_container_t::get_use_gemmi)
    .def("get_use_rama_plot_restraints",&molecules_container_t::get_use_rama_plot_restraints)
    .def("get_use_torsion_restraints",&molecules_container_t::get_use_torsion_restraints)
    .def("get_vertices_for_blender", &molecules_container_t::get_vertices_for_blender)
    .def("go_to_blob",&molecules_container_t::go_to_blob)
    .def("import_cif_dictionary",&molecules_container_t::import_cif_dictionary)
    .def("init_refinement_of_molecule_as_fragment_based_on_reference",&molecules_container_t::init_refinement_of_molecule_as_fragment_based_on_reference)
    .def("is_a_difference_map",&molecules_container_t::is_a_difference_map)
    .def("is_valid_map_molecule",&molecules_container_t::is_valid_map_molecule)
    .def("is_valid_model_molecule",&molecules_container_t::is_valid_model_molecule)
    .def("jed_flip",          nb::overload_cast<int, const std::string&, bool>           (&molecules_container_t::jed_flip))
    .def("make_mesh_for_bonds_for_blender", &molecules_container_t::make_mesh_for_bonds_for_blender)
    .def("make_mesh_for_gaussian_surface_for_blender", &molecules_container_t::make_mesh_for_gaussian_surface_for_blender)
    .def("make_mesh_for_goodsell_style_for_blender", &molecules_container_t::make_mesh_for_goodsell_style_for_blender)
    .def("make_mesh_for_map_contours_for_blender", &molecules_container_t::make_mesh_for_map_contours_for_blender)
    .def("make_mesh_for_molecular_representation_for_blender", &molecules_container_t::make_mesh_for_molecular_representation_for_blender)
    .def("mask_map_by_atom_selection",&molecules_container_t::mask_map_by_atom_selection)
    .def("merge_molecules", nb::overload_cast<int,const std::string &>(&molecules_container_t::merge_molecules))
    .def("minimize_energy",&molecules_container_t::minimize_energy)
    .def("mmrrcc",&molecules_container_t::mmrrcc)
    .def("move_molecule_to_new_centre",&molecules_container_t::move_molecule_to_new_centre)
    .def("multiply_residue_temperature_factors",&molecules_container_t::multiply_residue_temperature_factors)
    .def("mutate",&molecules_container_t::mutate)
    .def("new_positions_for_atoms_in_residues",&molecules_container_t::new_positions_for_atoms_in_residues)
    .def("new_positions_for_residue_atoms",&molecules_container_t::new_positions_for_residue_atoms)
    .def("non_standard_residue_types_in_model",&molecules_container_t::non_standard_residue_types_in_model)
    .def("pepflips_using_difference_map",&molecules_container_t::pepflips_using_difference_map)
    .def("peptide_omega_analysis",&molecules_container_t::peptide_omega_analysis)
    .def("print_secondary_structure_info",&molecules_container_t::print_secondary_structure_info)
    .def("rail_points_total",&molecules_container_t::rail_points_total)
    .def("ramachandran_analysis",&molecules_container_t::ramachandran_analysis)
    .def("ramachandran_validation",&molecules_container_t::ramachandran_validation)
    .def("read_ccp4_map",&molecules_container_t::read_ccp4_map)
    .def("read_mtz",&molecules_container_t::read_mtz)
    .def("read_pdb",&molecules_container_t::read_pdb)
    .def("redo",&molecules_container_t::redo)
    .def("refine",&molecules_container_t::refine)
    .def("refine_residue_range",&molecules_container_t::refine_residue_range)
    .def("refine_residues_using_atom_cid",&molecules_container_t::refine_residues_using_atom_cid)
    .def("regen_map",&molecules_container_t::regen_map)
    .def("replace_fragment",&molecules_container_t::replace_fragment)
    .def("replace_map_by_mtz_from_file",&molecules_container_t::replace_map_by_mtz_from_file)
    .def("replace_molecule_by_model_from_file",&molecules_container_t::replace_molecule_by_model_from_file)
    .def("residues_with_missing_atoms",&molecules_container_t::residues_with_missing_atoms)
    .def("rigid_body_fit",&molecules_container_t::rigid_body_fit)
    .def("rotamer_analysis",&molecules_container_t::rotamer_analysis)
    .def("set_colour_wheel_rotation_base",&molecules_container_t::set_colour_wheel_rotation_base)
    .def("set_draw_missing_residue_loops",&molecules_container_t::set_draw_missing_residue_loops)
    .def("set_draw_missing_residue_loops",&molecules_container_t::set_draw_missing_residue_loops)
    .def("set_imol_refinement_map",&molecules_container_t::set_imol_refinement_map)
    .def("set_make_backups",&molecules_container_t::set_make_backups)
    .def("set_map_sampling_rate",&molecules_container_t::set_map_sampling_rate)
    .def("set_map_weight",&molecules_container_t::set_map_weight)
    .def("set_molecule_name",&molecules_container_t::set_molecule_name)
    .def("set_rama_plot_restraints_weight",&molecules_container_t::set_rama_plot_restraints_weight)
    .def("set_refinement_is_verbose",&molecules_container_t::set_refinement_is_verbose)
    .def("set_refinement_is_verbose",&molecules_container_t::set_refinement_is_verbose)
    .def("set_show_timings",&molecules_container_t::set_show_timings)
    .def("set_torsion_restraints_weight",&molecules_container_t::set_torsion_restraints_weight)
    .def("set_use_gemmi",&molecules_container_t::set_use_gemmi)
    .def("set_use_rama_plot_restraints",&molecules_container_t::set_use_rama_plot_restraints)
    .def("set_use_torsion_restraints",&molecules_container_t::set_use_torsion_restraints)
    .def("set_user_defined_atom_colour_by_selection",&molecules_container_t::set_user_defined_atom_colour_by_selection)
    .def("set_user_defined_bond_colours",&molecules_container_t::set_user_defined_bond_colours)
    .def("sfcalc_genmap",&molecules_container_t::sfcalc_genmap)
    .def("sfcalc_genmaps_using_bulk_solvent",&molecules_container_t::sfcalc_genmaps_using_bulk_solvent)
    .def("sharpen_blur_map",&molecules_container_t::sharpen_blur_map)
    .def("sharpen_blur_map_with_resample",&molecules_container_t::sharpen_blur_map_with_resample)
    .def("side_chain_180",    nb::overload_cast<int, const std::string&>                         (&molecules_container_t::side_chain_180))
    .def("split_multi_model_molecule",&molecules_container_t::split_multi_model_molecule)
    .def("test_origin_cube",&molecules_container_t::test_origin_cube)
    .def("undo",&molecules_container_t::undo)
    .def("unmodelled_blobs",&molecules_container_t::unmodelled_blobs)
    .def("write_coordinates",&molecules_container_t::write_coordinates)
    .def("write_map",&molecules_container_t::write_map)
    .def("write_png",&molecules_container_t::write_png)
    ;
    nb::class_<molecules_container_js, molecules_container_t>(m,"molecules_container_py")
    .def(nb::init<bool>())
    .def("writePDBASCII",&molecules_container_js::writePDBASCII)
    .def("writeCIFASCII",&molecules_container_js::writeCIFASCII)
    .def("writeCCP4Map",&molecules_container_js::writeCCP4Map)
    ;
    nb::class_<coot::simple_rotamer>(m,"simple_rotamer")
    .def("P_r1234",&coot::simple_rotamer::P_r1234)
    .def("Probability_rich",&coot::simple_rotamer::Probability_rich)
    .def("get_chi",&coot::simple_rotamer::get_chi)
    ;
    nb::class_<merge_molecule_results_info_t>(m,"merge_molecule_results_info_t")
    .def_ro("chain_id", &merge_molecule_results_info_t::chain_id)
    .def_prop_ro("spec",[](merge_molecule_results_info_t &t) { return t.spec ; })
    .def_ro("is_chain", &merge_molecule_results_info_t::is_chain)
    ;
    nb::class_<coot::residue_validation_information_t>(m,"residue_validation_information_t")
    .def_ro("function_value", &coot::residue_validation_information_t::function_value)
    .def_ro("label", &coot::residue_validation_information_t::label)
    .def_prop_ro("residue_spec",[](coot::residue_validation_information_t &t) { return t.residue_spec ; })
    .def_prop_ro("atom_spec",[](coot::residue_validation_information_t &t) { return t.atom_spec ; })
    ;
    nb::class_<coot::chain_validation_information_t>(m,"chain_validation_information_t")
    .def_ro("chain_id", &coot::chain_validation_information_t::chain_id)
    .def_ro("rviv", &coot::chain_validation_information_t::rviv)
    ;
    nb::class_<coot::validation_information_t>(m,"validation_information_t")
    .def_ro("name", &coot::validation_information_t::name)
    .def_ro("type", &coot::validation_information_t::type)
    .def_ro("cviv", &coot::validation_information_t::cviv)
    .def("get_index_for_chain",&coot::validation_information_t::get_index_for_chain)
    ;
    nb::class_<molecules_container_t::fit_ligand_info_t>(m, "fit_ligand_info_t")
    .def_ro("imol", &molecules_container_t::fit_ligand_info_t::imol)
    .def_ro("cluster_idx", &molecules_container_t::fit_ligand_info_t::cluster_idx)
    .def_ro("ligand_idx", &molecules_container_t::fit_ligand_info_t::ligand_idx)
    .def("get_fitting_score", &molecules_container_t::fit_ligand_info_t::get_fitting_score)
    .def("get_cluster_volume", &molecules_container_t::fit_ligand_info_t::get_cluster_volume)
    ;
    nb::class_<coot::residue_spec_t>(m,"residue_spec_t")
    .def(nb::init<const std::string &, int, const std::string &>())
    .def_rw("model_number",&coot::residue_spec_t::model_number)
    .def_rw("chain_id",&coot::residue_spec_t::chain_id)
    .def_rw("res_no",&coot::residue_spec_t::res_no)
    .def_rw("ins_code",&coot::residue_spec_t::ins_code)
    .def_rw("int_user_data",&coot::residue_spec_t::int_user_data)
    ;
    nb::class_<coot::atom_spec_t>(m,"atom_spec_t")
    .def(nb::init<const std::string &, int, const std::string &, const std::string &, const std::string &>())
    .def_rw("chain_id",&coot::atom_spec_t::chain_id)
    .def_rw("res_no",&coot::atom_spec_t::res_no)
    .def_rw("ins_code",&coot::atom_spec_t::ins_code)
    .def_rw("atom_name",&coot::atom_spec_t::atom_name)
    .def_rw("alt_conf",&coot::atom_spec_t::alt_conf)
    .def_rw("int_user_data",&coot::atom_spec_t::int_user_data)
    .def_rw("float_user_data",&coot::atom_spec_t::float_user_data)
    .def_rw("string_user_data",&coot::atom_spec_t::string_user_data)
    .def_rw("model_number",&coot::atom_spec_t::model_number)
    ;
    nb::class_<generic_3d_lines_bonds_box_t>(m,"generic_3d_lines_bonds_box_t")
    .def_ro("line_segments", &generic_3d_lines_bonds_box_t::line_segments)
    ;
    nb::class_<coot::CartesianPair>(m,"CartesianPair")
    .def("getStart", &coot::CartesianPair::getStart)
    .def("getFinish", &coot::CartesianPair::getFinish)
    .def("amplitude", &coot::CartesianPair::amplitude)
    ;
    nb::class_<RamachandranInfo>(m,"RamachandranInfo")
    .def(nb::init<>())
    .def_rw("chainId", &RamachandranInfo::chainId)
    .def_rw("seqNum", &RamachandranInfo::seqNum)
    .def_rw("insCode", &RamachandranInfo::insCode)
    .def_rw("restype", &RamachandranInfo::restype)
    .def_rw("phi", &RamachandranInfo::phi)
    .def_rw("psi", &RamachandranInfo::psi)
    .def_rw("isOutlier", &RamachandranInfo::isOutlier)
    .def_rw("is_pre_pro", &RamachandranInfo::is_pre_pro)
    ;
    nb::class_<ResiduePropertyInfo>(m,"ResiduePropertyInfo")
    .def(nb::init<>())
    .def_rw("chainId", &ResiduePropertyInfo::chainId)
    .def_rw("seqNum", &ResiduePropertyInfo::seqNum)
    .def_rw("insCode", &ResiduePropertyInfo::insCode)
    .def_rw("restype", &ResiduePropertyInfo::restype)
    .def_rw("property", &ResiduePropertyInfo::property)
    ;
    //TODO = spped up the return of these meshes
    nb::class_<coot::instancing_data_type_A_t>(m,"instancing_data_type_A_t")
    .def_prop_ro("colour", [](coot::instancing_data_type_A_t &m) {
        float data[4] = {m.colour[0], m.colour[1], m.colour[2], m.colour[3]};
        std::vector<float> dest;
        dest.insert(dest.begin(), std::begin(data), std::end(data));
        return dest;
    })
    .def_prop_ro("size", [](coot::instancing_data_type_A_t &m) {
        float data[3] = {m.size[0], m.size[1], m.size[2]};
        std::vector<float> dest;
        dest.insert(dest.begin(), std::begin(data), std::end(data));
        return dest;
    })
    .def_prop_ro("position", [](coot::instancing_data_type_A_t &m) {
        float data[3] = {m.position[0], m.position[1], m.position[2]};
        std::vector<float> dest;
        dest.insert(dest.begin(), std::begin(data), std::end(data));
        return dest;
    })
    ;
    nb::class_<coot::instancing_data_type_B_t>(m,"instancing_data_type_B_t")
    .def_prop_ro("colour", [](coot::instancing_data_type_B_t &m) {
        float data[4] = {m.colour[0], m.colour[1], m.colour[2], m.colour[3]};
        std::vector<float> dest;
        dest.insert(dest.begin(), std::begin(data), std::end(data));
        return dest;
    })
    .def_prop_ro("size", [](coot::instancing_data_type_B_t &m) {
        float data[3] = {m.size[0], m.size[1], m.size[2]};
        std::vector<float> dest;
        dest.insert(dest.begin(), std::begin(data), std::end(data));
        return dest;
    })
    .def_prop_ro("position", [](coot::instancing_data_type_B_t &m) {
        float data[3] = {m.position[0], m.position[1], m.position[2]};
        std::vector<float> dest;
        dest.insert(dest.begin(), std::begin(data), std::end(data));
        return dest;
    })
    .def_prop_ro("orientation", [](coot::instancing_data_type_B_t &m) {
        float data[16] = {
            m.orientation[0][0], m.orientation[0][1], m.orientation[0][2], m.orientation[0][3],
            m.orientation[1][0], m.orientation[1][1], m.orientation[1][2], m.orientation[1][3],
            m.orientation[2][0], m.orientation[2][1], m.orientation[2][2], m.orientation[2][3],
            m.orientation[3][0], m.orientation[3][1], m.orientation[3][2], m.orientation[3][3],
        };
        std::vector<float> dest;
        dest.insert(dest.begin(), std::begin(data), std::end(data));
        return dest;
    })
    ;
    nb::class_<coot::instanced_geometry_t>(m,"instanced_geometry_t")
    .def_ro("vertices",          &coot::instanced_geometry_t::vertices)
    .def_ro("triangles",         &coot::instanced_geometry_t::triangles)
    .def_ro("instancing_data_A", &coot::instanced_geometry_t::instancing_data_A)
    .def_ro("instancing_data_B", &coot::instanced_geometry_t::instancing_data_B)
    .def_ro("name",              &coot::instanced_geometry_t::name)
    ;
    nb::class_<coot::instanced_mesh_t>(m,"instanced_mesh_t")
    .def_ro("geom",   &coot::instanced_mesh_t::geom)
    .def_ro("markup", &coot::instanced_mesh_t::markup)
    ;
    nb::class_<coot::util::phi_psi_t>(m,"phi_psi_t")
    .def("phi",               &coot::util::phi_psi_t::phi)
    .def("psi",               &coot::util::phi_psi_t::psi)
    .def("label",             &coot::util::phi_psi_t::label)
    .def("residue_name",      &coot::util::phi_psi_t::residue_name)
    .def("is_filled",         &coot::util::phi_psi_t::is_filled)
    .def("is_pre_pro",        &coot::util::phi_psi_t::is_pre_pro)
    .def_ro("ins_code",       &coot::util::phi_psi_t::ins_code)
    .def_ro("chain_id",       &coot::util::phi_psi_t::chain_id)
    .def_ro("residue_number", &coot::util::phi_psi_t::residue_number)
    ;
    nb::class_<coot::Cartesian>(m,"Cartesian")
    .def("x", &coot::Cartesian::x)
    .def("y", &coot::Cartesian::y)
    .def("z", &coot::Cartesian::z)
    ;
    nb::class_<coot::api::vnc_vertex>(m,"vnc_vertex")
    .def(nb::init<const glm::vec3 &, const glm::vec3 &, const glm::vec4 &>())
    .def_prop_ro("pos", [](coot::api::vnc_vertex &m) {
        const float data[3] = {m.pos[0], m.pos[1], m.pos[2]};
        std::vector<float> dest;
        dest.insert(dest.begin(), std::begin(data), std::end(data));
        return dest;
    })
    .def_prop_ro("normal", [](coot::api::vnc_vertex &m) {
        float data[3] = {m.normal[0], m.normal[1], m.normal[2]};
        std::vector<float> dest;
        dest.insert(dest.begin(), std::begin(data), std::end(data));
        return dest;
    })
    .def_prop_ro("color", [](coot::api::vnc_vertex &m) {
        float data[4] = {m.color[0], m.color[1], m.color[2], m.color[3]};
        std::vector<float> dest;
        dest.insert(dest.begin(), std::begin(data), std::end(data));
        return dest;
    })
    ;
    nb::class_<coot::api::vn_vertex>(m,"vn_vertex")
    .def(nb::init<const glm::vec3 &, const glm::vec3 &>())
    .def_prop_ro("pos", [](coot::api::vn_vertex &m) {
        const float data[3] = {m.pos[0], m.pos[1], m.pos[2]};
        std::vector<float> dest;
        dest.insert(dest.begin(), std::begin(data), std::end(data));
        return dest;
    })
    .def_prop_ro("normal", [](coot::api::vn_vertex &m) {
        float data[3] = {m.normal[0], m.normal[1], m.normal[2]};
        std::vector<float> dest;
        dest.insert(dest.begin(), std::begin(data), std::end(data));
        return dest;
    })
    ;
    nb::class_<coot::molecule_t::rotamer_change_info_t>(m,"rotamer_change_info_t")
    .def_ro("rank",                   &coot::molecule_t::rotamer_change_info_t::rank)
    .def_ro("name",                   &coot::molecule_t::rotamer_change_info_t::name)
    .def_ro("richardson_probability", &coot::molecule_t::rotamer_change_info_t::richardson_probability)
    .def_ro("status",                 &coot::molecule_t::rotamer_change_info_t::status)
    ;
    nb::class_<g_triangle>(m,"g_triangle")
    .def(nb::init<unsigned int,unsigned int,unsigned int>())
    //member way
    .def_prop_ro("point_id", [](g_triangle &m) {
        const unsigned data[3] = {m.point_id[0], m.point_id[1], m.point_id[2]};
        std::vector<unsigned> dest;
        dest.insert(dest.begin(), std::begin(data), std::end(data));
        return dest;
    })
    //accessor way
    .def("set_point_id", [](g_triangle &m,std::vector<unsigned> d_in) {
        m.point_id[0] = d_in[0];
        m.point_id[1] = d_in[1];
        m.point_id[2] = d_in[2];
    })
    .def("get_point_id", [](g_triangle &m) {
        std::vector<unsigned> ret;
        ret.push_back(m.point_id[0]);
        ret.push_back(m.point_id[1]);
        ret.push_back(m.point_id[2]);
        return ret;
    });
    ;
    nb::class_<Cell_Translation>(m,"Cell_Translation")
    .def(nb::init<>())
    .def(nb::init<int,int,int>())
    .def_ro("us", &Cell_Translation::us)
    .def_ro("ws", &Cell_Translation::ws)
    .def_ro("vs", &Cell_Translation::vs)
    ;
    nb::class_<symm_trans_t>(m,"symm_trans_t")
    .def_ro("symm_as_string",&symm_trans_t::symm_as_string)
    .def("is_identity",&symm_trans_t::is_identity)
    .def("add_shift",&symm_trans_t::add_shift)
    .def("isym",&symm_trans_t::isym)
    .def("x",&symm_trans_t::x)
    .def("y",&symm_trans_t::y)
    .def("z",&symm_trans_t::z)
    ;
    nb::class_<coot::simple_mesh_t>(m,"simple_mesh_t")
    .def_ro("vertices",  &coot::simple_mesh_t::vertices)
    .def_ro("triangles", &coot::simple_mesh_t::triangles)
    .def_ro("status",    &coot::simple_mesh_t::status)
    .def_ro("name",      &coot::simple_mesh_t::name)
    ;
    // nb::class_<coot::blender_mesh_t>(m,"blender_mesh_t")
    //    .def_ro("vertices",  &coot::blender_mesh_t::vertices)
    //    .def_ro("normals",   &coot::blender_mesh_t::normals)
    //    .def_ro("triangles", &coot::blender_mesh_t::triangles)
    // ;

    nb::class_<coot::util::density_correlation_stats_info_t>(m,"density_correlation_stats_info_t")
    .def_ro("n",          &coot::util::density_correlation_stats_info_t::n)
    .def_ro("sum_xy",     &coot::util::density_correlation_stats_info_t::sum_xy)
    .def_ro("sum_sqrd_x", &coot::util::density_correlation_stats_info_t::sum_sqrd_x)
    .def_ro("sum_sqrd_y", &coot::util::density_correlation_stats_info_t::sum_sqrd_y)
    .def_ro("sum_x",      &coot::util::density_correlation_stats_info_t::sum_x)
    .def_ro("sum_y",      &coot::util::density_correlation_stats_info_t::sum_y)
    .def("var_x",         &coot::util::density_correlation_stats_info_t::var_x)
    .def("var_y",         &coot::util::density_correlation_stats_info_t::var_y)
    .def("correlation",   &coot::util::density_correlation_stats_info_t::correlation)
    ;

    nb::class_<superpose_results_t>(m,"superpose_results_t")
       .def_ro("superpose_info",     &superpose_results_t::superpose_info) // a json file (string)
       .def_ro("alignment",          &superpose_results_t::alignment)
       .def_ro("alignment_info_vec", &superpose_results_t::alignment_info_vec)
       .def_ro("aligned_pairs",      &superpose_results_t::aligned_pairs)
    ;

    nb::class_<moorhen::h_bond>(m,"h_bond")
        .def_ro("hb_hydrogen",&moorhen::h_bond::hb_hydrogen)
        .def_ro("donor",&moorhen::h_bond::donor)
        .def_ro("acceptor",&moorhen::h_bond::acceptor)
        .def_ro("donor_neigh",&moorhen::h_bond::donor_neigh)
        .def_ro("acceptor_neigh",&moorhen::h_bond::acceptor_neigh)
        .def_ro("angle_1",&moorhen::h_bond::angle_1)
        .def_ro("angle_2",&moorhen::h_bond::angle_2)
        .def_ro("angle_3",&moorhen::h_bond::angle_3)
        .def_ro("dist",&moorhen::h_bond::dist)
        .def_ro("ligand_atom_is_donor",&moorhen::h_bond::ligand_atom_is_donor)
        .def_ro("hydrogen_is_ligand_atom",&moorhen::h_bond::hydrogen_is_ligand_atom)
        .def_ro("bond_has_hydrogen_flag",&moorhen::h_bond::bond_has_hydrogen_flag)
    ;

    nb::class_<moorhen::h_bond_atom>(m,"h_bond_atom")
        .def_ro("serial",&moorhen::h_bond_atom::serial)
        .def_ro("x",&moorhen::h_bond_atom::x)
        .def_ro("y",&moorhen::h_bond_atom::y)
        .def_ro("z",&moorhen::h_bond_atom::z)
        .def_ro("charge",&moorhen::h_bond_atom::charge)
        .def_ro("occ",&moorhen::h_bond_atom::occ)
        .def_ro("b_iso",&moorhen::h_bond_atom::b_iso)
        .def_ro("element",&moorhen::h_bond_atom::element)
        .def_ro("name",&moorhen::h_bond_atom::name)
        .def_ro("model",&moorhen::h_bond_atom::model)
        .def_ro("chain",&moorhen::h_bond_atom::chain)
        .def_ro("res_no",&moorhen::h_bond_atom::res_no)
        .def_ro("residue_name",&moorhen::h_bond_atom::residue_name)
        .def_ro("altLoc",&moorhen::h_bond_atom::altLoc)
    ;
    nb::class_<coot::phi_psi_prob_t>(m,"phi_psi_prob_t")
    .def_ro("phi_psi", &coot::phi_psi_prob_t::phi_psi)
    .def_ro("position", &coot::phi_psi_prob_t::position)
    .def_ro("is_allowed_flag", &coot::phi_psi_prob_t::is_allowed_flag)
    .def("residue_name", &coot::phi_psi_prob_t::residue_name)
    .def("is_allowed", &coot::phi_psi_prob_t::is_allowed)
    ;
    nb::class_<coot::molecule_t::moved_atom_t>(m,"moved_atom_t")
    .def(nb::init<const std::string&, const std::string&, float, float, float, int>())
    .def_ro("atom_name", &coot::molecule_t::moved_atom_t::atom_name)
    .def_ro("alt_conf", &coot::molecule_t::moved_atom_t::alt_conf)
    .def_ro("x", &coot::molecule_t::moved_atom_t::x)
    .def_ro("y", &coot::molecule_t::moved_atom_t::y)
    .def_ro("z", &coot::molecule_t::moved_atom_t::z)
    .def_ro("index", &coot::molecule_t::moved_atom_t::index)
    ;
    nb::class_<coot::molecule_t::interesting_place_t>(m,"interesting_place_t")
    .def(nb::init<const std::string &, const coot::residue_spec_t &, const clipper::Coord_orth &, const std::string &>())
    .def(nb::init<const std::string &, const clipper::Coord_orth &, const std::string &>())
    .def_ro("feature_type", &coot::molecule_t::interesting_place_t::feature_type)
    .def_ro("residue_spec", &coot::molecule_t::interesting_place_t::residue_spec)
    .def_ro("button_label", &coot::molecule_t::interesting_place_t::button_label)
    .def_ro("feature_value", &coot::molecule_t::interesting_place_t::feature_value)
    .def_ro("badness", &coot::molecule_t::interesting_place_t::badness)
    .def_ro("x", &coot::molecule_t::interesting_place_t::x)
    .def_ro("y", &coot::molecule_t::interesting_place_t::y)
    .def_ro("z", &coot::molecule_t::interesting_place_t::z)
    ;
    nb::class_<coot::molecule_t::moved_residue_t>(m,"moved_residue_t")
    .def(nb::init<const std::string&, int, const std::string&>())
    .def_ro("chain_id", &coot::molecule_t::moved_residue_t::chain_id)
    .def_ro("res_no", &coot::molecule_t::moved_residue_t::res_no)
    .def_ro("ins_code", &coot::molecule_t::moved_residue_t::ins_code)
    .def_ro("moved_atoms", &coot::molecule_t::moved_residue_t::moved_atoms)
    .def("add_atom",&coot::molecule_t::moved_residue_t::add_atom)
    ;
}
