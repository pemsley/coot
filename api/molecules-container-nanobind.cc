
#include <filesystem>
#include <unordered_map>
#include <sstream>

#include <stdlib.h> // for getenv()
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/vector.h>
#include <nanobind/ndarray.h>
#include <nanobind/operators.h>

#include <clipper/core/ramachandran.h>
#include <clipper/clipper-ccp4.h>

#include "coot-utils/pugixml.hpp"
#include "coords/mmdb-crystal.hh"
#include "coot-utils/acedrg-types-for-residue.hh"
#include "coot-utils/g_triangle.hh"
#include "mini-mol/mini-mol-utils.hh"

#if NB_VERSION_MAJOR // for flychecking
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>        // For std::string conversion
#include <nanobind/stl/shared_ptr.h>    // <-- CRITICAL: Provides std::shared_ptr bindings
#include <nanobind/stl/vector.h>        // Typically useful for RDKit
#endif


#include "molecules-container.hh"

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

// Helper to cache and retrieve docstrings from XML
std::unordered_map<std::string, std::string> docstring_cache;

std::string get_docstring_from_xml(const std::string& func_name) {

   // this is the relative path for standard out-of-tree build
   std::string api_doxygen_xml_file_name_1 =
      "../coot/api/doxy-sphinx/xml/classmolecules__container__t.xml";
   // for a build inside the source:
   std::string api_doxygen_xml_file_name_2 =
      "../api/doxy-sphinx/xml/classmolecules__container__t.xml";
   // fill this:
   std::string api_doxygen_xml_file_name;
   if (std::filesystem::exists(std::filesystem::path(api_doxygen_xml_file_name_1)))
      api_doxygen_xml_file_name = api_doxygen_xml_file_name_1;
   if (std::filesystem::exists(std::filesystem::path(api_doxygen_xml_file_name_2)))
      api_doxygen_xml_file_name = api_doxygen_xml_file_name_2;
   // try to find the xml file using CONDA_PREFIX - idea from eunos-1128
   const char *e = getenv("CONDA_PREFIX");
   if (e) {
      std::filesystem::path conda_prefix(e);
      std::filesystem::path xml_dir = conda_prefix / "share" / "doxy-sphinx" / "xml";
      std::filesystem::path full_path = xml_dir / "classmolecules__container__t.xml";
      if (std::filesystem::exists(full_path))
         api_doxygen_xml_file_name = full_path.string();
   }

   auto convert_type = [] (const std::string &s_in) {
      std::string s = s_in;
      if (s_in == "const std::string &") s = "str";
      if (s_in == "std::string")         s = "str";
      if (s_in == "void")                s = "None";
      if (s_in == "std::vector<")        s = "list";
      if (s_in == "std::vector< std::pair< double, double > >") s = "list";
      if (s_in == "std::vector< std::pair< std::string, std::string > >") s = "list";
      if (s_in.find("std::vector<") != std::string::npos) s = "list";
      if (s_in.compare(0,10, "std::pair<") == 0) s= "tuple"; // needs checking. Use starts_with() in C++20
      return s;
   };

   class arg_info_t {
   public:
      arg_info_t(const std::string &n, const std::string &t) : name(n), type(t) {}
      std::string name;
      std::string type;
      std::string description;
   };

   auto update_arg_in_args = [] (const std::string &arg_name, const std::string &descr,
                                 std::vector<arg_info_t> &args) {

      for(auto &arg : args) {
         if (arg.name == arg_name) {
            arg.description = descr;
         }
      }
   };

   if (docstring_cache.empty()) {
      pugi::xml_document doc;
      if (!doc.load_file(api_doxygen_xml_file_name.c_str())) {
         std::cout << "WARNING:: doxygen file " << api_doxygen_xml_file_name
                   << " not found - nanobind API docummentation will not be generated"
                   << std::endl;
         return "";
      }
      auto compounddef = doc.child("doxygen").child("compounddef");
      for (auto sectiondef : compounddef.children("sectiondef")) {
         for (auto member : sectiondef.children("memberdef")) {
            auto name_elem = member.child("name");
            if (!name_elem) continue;
            std::string name = name_elem.child_value();
            std::ostringstream oss;
            std::vector<arg_info_t> args;

            // Collect all <para> from <briefdescription>
            auto brief = member.child("briefdescription");
            if (brief) {
               for (auto para : brief.children("para")) {
                  std::string para_text = para.text().get();
                  if (!para_text.empty())
                     oss << para_text << "\n";
               }
            }

            auto type = member.child("type");
            std::string type_string;
            if (type)
               type_string = convert_type(type.text().get());

            // can be many params
            for (auto param : member.children("param")) {
               auto p_type = param.child("type");
               auto p_declname = param.child("declname");
               if (p_type) {
                  if (p_declname) {
                     std::string tt = convert_type(p_type.text().get());
                     arg_info_t ai(p_declname.text().get(), tt);
                     args.push_back(ai);
                  }
               }
            }

	    std::string return_type_docs;

            // Collect all <para> from <detaileddescription>
            auto detailed = member.child("detaileddescription");
            if (detailed) {
               unsigned int n_para = 0;
               for (auto para : detailed.children("para")) {
                  n_para++;
                  std::string para_text = para.text().get();
                  if (!para_text.empty()) {
		     if (n_para > 1)
			oss << "\n    ";
		     if (para_text[0] == '\n')
			para_text.erase(0,1); // remove first char
                     oss << para_text << "\n";
                  }
                  for (auto parameterlist : para.children("parameterlist")) {
                     for (auto parameteritem : parameterlist.children("parameteritem")) {
                        std::string parameter_name_text;
                        for (auto parameternamelist : parameteritem.children("parameternamelist")) {
                           for (auto parametername : parameternamelist.children("parametername")) {
                              parameter_name_text = parametername.text().get();
                           }
                        }
                        for (auto parameterdescription : parameteritem.children("parameterdescription")) {
                           for (auto d_para : parameterdescription.children("para")) {
                              std::string t =  d_para.text().get();
                              if (! parameter_name_text.empty()) {
                                 if (! t.empty()) {
                                    update_arg_in_args(parameter_name_text, t, args); // modify an arg in args
                                 }
                              }
                           }
                        }
                     }
                  }
                  for (auto simplesect : para.children("simplesect")) {
		     if (std::string(simplesect.attribute("kind").value()) == "return") {
			for(auto ss_para : simplesect.children("para")) {
			   std::string return_type_doc = ss_para.text().get();
			   return_type_docs += return_type_doc;
			}
		     }
		  }
               }
               if (! args.empty()) {
                  oss << "\n";
                  oss << "    Args:\n";
                  for (const auto &arg : args) {
                     oss << "        " << arg.name << " (" << arg.type << "): " << arg.description << "\n";
                  }
               }
               if (! type_string.empty()) {
                  oss << "\n";
                  oss << "    Returns:\n";
		  if (return_type_docs.empty())
		     oss << "        " << type_string << "\n";
		  else
		     oss << "        " << type_string << ": " << return_type_docs << "\n";
               }
            }
            docstring_cache[name] = oss.str();
         }
      }
   }
   auto it = docstring_cache.find(func_name);
   if (it != docstring_cache.end()) {
      if (false) { // debugging - this can be quiet now
	 std::cout << "function:" << func_name << "()" << std::endl;
	 std::cout << it->second << std::endl;
      }
      return it->second;
   } else {
      std::cout << "::function " << func_name << " not found"
		<< " - out of " << docstring_cache.size() << " docstrings" << std::endl;
   }
   return "";
}

NB_MODULE(coot_headless_api, m) {
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
    .def(nb::init<bool>(), nb::arg("be_verbose_when_reading_dictionary"), "molecules container Documentation")
    .def("M2T_updateFloatParameter",
         &molecules_container_t::M2T_updateFloatParameter,
         nb::arg("imol"), nb::arg("param_name"), nb::arg("value"),
         get_docstring_from_xml("M2T_updateFloatParameter").c_str())
    .def("M2T_updateIntParameter",
         &molecules_container_t::M2T_updateIntParameter,
         nb::arg("imol"), nb::arg("param_name"), nb::arg("value"),
         get_docstring_from_xml("M2T_updateIntParameter").c_str())
    .def("SSM_superpose",
         &molecules_container_t::SSM_superpose,
         nb::arg("imol_ref"), nb::arg("chain_id_ref"),
         nb::arg("imol_mov"), nb::arg("chain_id_mov"),
         get_docstring_from_xml("SSM_superpose").c_str())
    .def("add_alternative_conformation",
         &molecules_container_t::add_alternative_conformation,
         nb::arg("imol_model"), nb::arg("cid"),
         get_docstring_from_xml("add_alternative_conformation").c_str())
    .def("add_colour_rule",
         &molecules_container_t::add_colour_rule,
         nb::arg("imol"), nb::arg("selection_cid"), nb::arg("colour"),
         get_docstring_from_xml("add_colour_rule").c_str())
    .def("add_colour_rules_multi",
         &molecules_container_t::add_colour_rules_multi,
         nb::arg("imol"), nb::arg("selections_and_colours_combo_string"),
         get_docstring_from_xml("add_colour_rules_multi").c_str())
    .def("add_hydrogen_atoms",
         &molecules_container_t::add_hydrogen_atoms,
         nb::arg("imol_model"),
         get_docstring_from_xml("add_hydrogen_atoms").c_str())
    .def("add_lsq_superpose_match",
         &molecules_container_t::add_lsq_superpose_match,
         nb::arg("chain_id_ref"), nb::arg("res_no_ref_start"), nb::arg("res_no_ref_end"), nb::arg("chain_id_mov"), nb::arg("res_no_mov_start"), nb::arg("res_no_mov_end"), nb::arg("match_type"),
         get_docstring_from_xml("add_lsq_superpose_match").c_str())
    .def("add_lsq_superpose_atom_match",
         &molecules_container_t::add_lsq_superpose_atom_match,
         nb::arg("chain_id_ref"), nb::arg("res_no_ref"), nb::arg("atom_name_ref"), nb::arg("chain_id_mov"), nb::arg("res_no_mov"), nb::arg("atom_name_mov"),
         get_docstring_from_xml("add_lsq_superpose_atom_match").c_str())
    .def("add_named_glyco_tree",
         &molecules_container_t::add_named_glyco_tree,
         nb::arg("imol_model"), nb::arg("imol_map"), nb::arg("glycosylation_name"), nb::arg("asn_chain_id"), nb::arg("asn_res_no"),
         get_docstring_from_xml("add_named_glyco_tree").c_str())
    .def("add_target_position_restraint",
         &molecules_container_t::add_target_position_restraint,
         nb::arg("imol"), nb::arg("atom_cid"), nb::arg("pos_x"), nb::arg("pos_y"), nb::arg("pos_z"),
         get_docstring_from_xml("add_target_position_restraint").c_str())
    .def("add_target_position_restraint_and_refine",
         &molecules_container_t::add_target_position_restraint_and_refine,
         nb::arg("imol"), nb::arg("atom_cid"), nb::arg("pos_x"), nb::arg("pos_y"), nb::arg("pos_z"), nb::arg("n_cycles"),
         get_docstring_from_xml("add_target_position_restraint_and_refine").c_str())
    .def("add_terminal_residue_directly",
         &molecules_container_t::add_terminal_residue_directly,
         nb::arg("imol"), nb::arg("chain_id"), nb::arg("res_no"), nb::arg("ins_code"),
         get_docstring_from_xml("add_terminal_residue_directly").c_str())
    .def("add_terminal_residue_directly_using_cid",
         &molecules_container_t::add_terminal_residue_directly_using_cid,
         nb::arg("imol"), nb::arg("cid"),
         get_docstring_from_xml("add_terminal_residue_directly_using_cid").c_str())
    .def("add_to_non_drawn_bonds",
         &molecules_container_t::add_to_non_drawn_bonds,
         nb::arg("imol"), nb::arg("atom_selection_cid"),
         get_docstring_from_xml("add_to_non_drawn_bonds").c_str())
    .def("add_waters",
         &molecules_container_t::add_waters,
         nb::arg("imol_model"), nb::arg("imol_map"),
         get_docstring_from_xml("add_waters").c_str())
    .def("all_molecule_contact_dots",
         &molecules_container_t::all_molecule_contact_dots,
         nb::arg("imol"), nb::arg("smoothness_factor"),
         get_docstring_from_xml("all_molecule_contact_dots").c_str())
    .def("apply_transformation_to_atom_selection",
         &molecules_container_t::apply_transformation_to_atom_selection,
         nb::arg("imol"), nb::arg("atoms_selection_cid"), nb::arg("n_atoms"), nb::arg("m00"), nb::arg("m01"), nb::arg("m02"), nb::arg("m10"), nb::arg("m11"), nb::arg("m12"), nb::arg("m20"), nb::arg("m21"), nb::arg("m22"), nb::arg("c0"), nb::arg("c1"), nb::arg("c2"), nb::arg("t0"), nb::arg("t1"), nb::arg("t2"),
         get_docstring_from_xml("apply_transformation_to_atom_selection").c_str())
    .def("assign_sequence",
         &molecules_container_t::assign_sequence,
         nb::arg("imol_model"), nb::arg("imol_map"),
         get_docstring_from_xml("assign_sequence").c_str())
    .def("associate_data_mtz_file_with_map",
         &molecules_container_t::associate_data_mtz_file_with_map,
         nb::arg("imol"), nb::arg("data_mtz_file_name"), nb::arg("f_col"), nb::arg("sigf_col"), nb::arg("free_r_col"),
         get_docstring_from_xml("associate_data_mtz_file_with_map").c_str())
    .def("associate_sequence",
         &molecules_container_t::associate_sequence,
         nb::arg("imol"), nb::arg("name_or_chain_id"), nb::arg("sequence"),
         get_docstring_from_xml("associate_sequence").c_str())
    .def("auto_fit_rotamer",
         &molecules_container_t::auto_fit_rotamer,
         nb::arg("imol"), nb::arg("chain_id"), nb::arg("res_no"), nb::arg("ins_code"), nb::arg("alt_conf"), nb::arg("imol_map"),
         get_docstring_from_xml("auto_fit_rotamer").c_str())
    .def("auto_read_mtz",
         &molecules_container_t::auto_read_mtz,
         nb::arg("file_name"),
         get_docstring_from_xml("auto_read_mtz").c_str())
    .def("average_map",
         &molecules_container_t::average_map,
         nb::arg("imol_maps"), nb::arg("scales"),
         get_docstring_from_xml("average_map").c_str())
    .def("calculate_new_rail_points",
         &molecules_container_t::calculate_new_rail_points,
         get_docstring_from_xml("calculate_new_rail_points").c_str())
    .def("change_to_first_rotamer",
         &molecules_container_t::change_to_first_rotamer,
         nb::arg("imol"), nb::arg("residue_cid"), nb::arg("alt_conf"),
         get_docstring_from_xml("change_to_first_rotamer").c_str())
    .def("change_to_next_rotamer",
         &molecules_container_t::change_to_next_rotamer,
         nb::arg("imol"), nb::arg("residue_cid"), nb::arg("alt_conf"),
         get_docstring_from_xml("change_to_next_rotamer").c_str())
    .def("change_to_previous_rotamer",
         &molecules_container_t::change_to_previous_rotamer,
         nb::arg("imol"), nb::arg("residue_cid"), nb::arg("alt_conf"),
         get_docstring_from_xml("change_to_previous_rotamer").c_str())
    .def("cis_trans_convert",
         &molecules_container_t::cis_trans_convert,
         nb::arg("imol"), nb::arg("atom_cid"),
         get_docstring_from_xml("cis_trans_convert").c_str())
    .def("clear_extra_restraints",
         &molecules_container_t::clear_extra_restraints,
         nb::arg("imol"),
         get_docstring_from_xml("clear_extra_restraints").c_str())
    .def("clear_lsq_matches",
         &molecules_container_t::clear_lsq_matches,
         get_docstring_from_xml("clear_lsq_matches").c_str())
    .def("clear_non_drawn_bonds",
         &molecules_container_t::clear_non_drawn_bonds,
         nb::arg("imol"),
         get_docstring_from_xml("clear_non_drawn_bonds").c_str())
    .def("clear_refinement",
         &molecules_container_t::clear_refinement,
         nb::arg("imol"),
         get_docstring_from_xml("clear_refinement").c_str())
    .def("clear_target_position_restraints",
         &molecules_container_t::clear_target_position_restraints,
         nb::arg("imol"),
         get_docstring_from_xml("clear_target_position_restraints").c_str())
    .def("change_alt_locs",
         &molecules_container_t::change_alt_locs,
         nb::arg("imol"), nb::arg("cid"), nb::arg("change_mode"),
         get_docstring_from_xml("change_alt_locs").c_str())
    .def("close_molecule",
         &molecules_container_t::close_molecule,
         nb::arg("imol"),
         get_docstring_from_xml("close_molecule").c_str())
    .def("connect_updating_maps",
         &molecules_container_t::connect_updating_maps,
         nb::arg("imol_model"), nb::arg("imol_with_data_info_attached"), nb::arg("imol_map_2fofc"), nb::arg("imol_map_fofc"),
         get_docstring_from_xml("connect_updating_maps").c_str())
    .def("contact_dots_for_ligand",
         &molecules_container_t::contact_dots_for_ligand,
         nb::arg("imol"), nb::arg("cid"), nb::arg("smoothness_factor"),
         get_docstring_from_xml("contact_dots_for_ligand").c_str())
    .def("copy_fragment_for_refinement_using_cid",
         &molecules_container_t::copy_fragment_for_refinement_using_cid,
         nb::arg("imol"), nb::arg("multi_cid"),
         get_docstring_from_xml("copy_fragment_for_refinement_using_cid").c_str())
    .def("copy_fragment_using_cid",
         &molecules_container_t::copy_fragment_using_cid,
         nb::arg("imol"), nb::arg("multi_cid"),
         get_docstring_from_xml("copy_fragment_using_cid").c_str())
    .def("copy_fragment_using_residue_range",
         &molecules_container_t::copy_fragment_using_residue_range,
         nb::arg("imol"), nb::arg("chain_id"), nb::arg("res_no_start"), nb::arg("res_no_end"),
         get_docstring_from_xml("copy_fragment_using_residue_range").c_str())
    .def("copy_molecule",
         &molecules_container_t::copy_molecule,
         nb::arg("imol"),
         get_docstring_from_xml("copy_molecule").c_str())
    .def("dedust_map",
          &molecules_container_t::dedust_map,
          nb::arg("imol"),
          get_docstring_from_xml("dedust_map").c_str())
    .def("delete_all_carbohydrate",
         &molecules_container_t::delete_all_carbohydrate,
         nb::arg("imol"),
         get_docstring_from_xml("delete_all_carbohydrate").c_str())
    .def("delete_atom",
         &molecules_container_t::delete_atom,
         nb::arg("imol"), nb::arg("chain_id"), nb::arg("res_no"), nb::arg("ins_code"), nb::arg("atom_name"), nb::arg("alt_conf"),
         get_docstring_from_xml("delete_atom").c_str())
    .def("delete_atom_using_cid",
         &molecules_container_t::delete_atom_using_cid,
         nb::arg("imol"), nb::arg("cid"),
         get_docstring_from_xml("delete_atom_using_cid").c_str())
    .def("delete_chain_using_cid",
         &molecules_container_t::delete_chain_using_cid,
         nb::arg("imol"), nb::arg("cid"),
         get_docstring_from_xml("delete_chain_using_cid").c_str())
    .def("delete_colour_rules",
         &molecules_container_t::delete_colour_rules,
         nb::arg("imol"),
         get_docstring_from_xml("delete_colour_rules").c_str())
    .def("delete_hydrogen_atoms",
         &molecules_container_t::delete_hydrogen_atoms,
         nb::arg("imol_model"),
         get_docstring_from_xml("delete_hydrogen_atoms").c_str())
    .def("delete_residue",
         &molecules_container_t::delete_residue,
         nb::arg("imol"), nb::arg("chain_id"), nb::arg("res_no"), nb::arg("ins_code"),
         get_docstring_from_xml("delete_residue").c_str())
    .def("delete_residue_atoms_using_cid",
         &molecules_container_t::delete_residue_atoms_using_cid,
         nb::arg("imol"), nb::arg("cid"),
         get_docstring_from_xml("delete_residue_atoms_using_cid").c_str())
    .def("delete_residue_atoms_with_alt_conf",
         &molecules_container_t::delete_residue_atoms_with_alt_conf,
         nb::arg("imol"), nb::arg("chain_id"), nb::arg("res_no"), nb::arg("ins_code"), nb::arg("alt_conf"),
         get_docstring_from_xml("delete_residue_atoms_with_alt_conf").c_str())
    .def("delete_residue_using_cid",
         &molecules_container_t::delete_residue_using_cid,
         nb::arg("imol"), nb::arg("cid"),
         get_docstring_from_xml("delete_residue_using_cid").c_str())
    .def("delete_side_chain",
         &molecules_container_t::delete_side_chain,
         nb::arg("imol"), nb::arg("chain_id"), nb::arg("res_no"), nb::arg("ins_code"),
         get_docstring_from_xml("delete_side_chain").c_str())
    .def("delete_using_cid",
         &molecules_container_t::delete_using_cid,
         nb::arg("imol"), nb::arg("cid"), nb::arg("scope"),
         get_docstring_from_xml("delete_using_cid").c_str())
    .def("density_correlation_analysis",
         &molecules_container_t::density_correlation_analysis,
         nb::arg("imol_model"), nb::arg("imol_map"),
         get_docstring_from_xml("density_correlation_analysis").c_str())
    .def("density_fit_analysis",
         &molecules_container_t::density_fit_analysis,
         nb::arg("imol_model"), nb::arg("imol_map"),
         get_docstring_from_xml("density_fit_analysis").c_str())
    .def("dictionary_atom_name_map",
         &molecules_container_t::dictionary_atom_name_map,
         nb::arg("comp_id_1"), nb::arg("imol_1"), nb::arg("comp_id_2"), nb::arg("imol_2"),
         get_docstring_from_xml("dictionary_atom_name_map").c_str())
    .def("difference_map_peaks",
         &molecules_container_t::difference_map_peaks,
         nb::arg("imol_map"), nb::arg("imol_protein"), nb::arg("n_rmsd"),
         get_docstring_from_xml("difference_map_peaks").c_str())
    .def("eigen_flip_ligand",
         &molecules_container_t::eigen_flip_ligand,
         nb::arg("imol"), nb::arg("chain_id"), nb::arg("res_no"), nb::arg("ins_code"),
         get_docstring_from_xml("eigen_flip_ligand").c_str())
    .def("eigen_flip_ligand_using_cid",
         &molecules_container_t::eigen_flip_ligand_using_cid,
         nb::arg("imol"), nb::arg("residue_id"),
         get_docstring_from_xml("eigen_flip_ligand_using_cid").c_str())
    .def("export_chemical_features_as_gltf",
         &molecules_container_t::export_chemical_features_as_gltf,
         nb::arg("imol"), nb::arg("cid"), nb::arg("file_name"),
         get_docstring_from_xml("export_chemical_features_as_gltf").c_str())
    .def("export_molecular_representation_as_gltf",
         &molecules_container_t::export_molecular_representation_as_gltf,
         nb::arg("imol"), nb::arg("atom_selection_cid"), nb::arg("colour_scheme"), nb::arg("style"), nb::arg("secondary_structure_usage_flag"), nb::arg("file_name"),
         get_docstring_from_xml("export_molecular_representation_as_gltf").c_str())
    .def("export_model_molecule_as_gltf",
         &molecules_container_t::export_model_molecule_as_gltf,
         nb::arg("imol"), nb::arg("selection_cid"), nb::arg("mode"), nb::arg("against_a_dark_background"), nb::arg("bonds_width"), nb::arg("atom_radius_to_bond_width_ratio"), nb::arg("smoothness_factor"), nb::arg("draw_hydrogen_atoms_flag"), nb::arg("draw_missing_residue_loops"), nb::arg("file_name"),
         get_docstring_from_xml("export_model_molecule_as_gltf").c_str())
    .def("export_map_molecule_as_gltf",
         &molecules_container_t::export_map_molecule_as_gltf,
         nb::arg("imol"), nb::arg("pos_x"), nb::arg("pos_y"), nb::arg("pos_z"), nb::arg("radius"), nb::arg("contour_level"), nb::arg("file_name"),
         get_docstring_from_xml("export_map_molecule_as_gltf").c_str())
    .def("find_water_baddies",
         &molecules_container_t::find_water_baddies,
         nb::arg("imol_model"), nb::arg("imol_map"), nb::arg("b_factor_lim"), nb::arg("outlier_sigma_level"), nb::arg("min_dist"), nb::arg("max_dist"), nb::arg("ignore_part_occ_contact_flag"), nb::arg("ignore_zero_occ_flag"),
         get_docstring_from_xml("find_water_baddies").c_str())
    .def("file_name_to_string",
         &molecules_container_t::file_name_to_string,
         nb::arg("file_name"),
         get_docstring_from_xml("file_name_to_string").c_str())
    .def("fill_partial_residue",
         &molecules_container_t::fill_partial_residue,
         nb::arg("imol"), nb::arg("chain_id"), nb::arg("res_no"), nb::arg("ins_code"),
         get_docstring_from_xml("fill_partial_residue").c_str())
    .def("fill_partial_residues",
         &molecules_container_t::fill_partial_residues,
         nb::arg("imol"),
         get_docstring_from_xml("fill_partial_residues").c_str())
    .def("fill_rotamer_probability_tables",
         &molecules_container_t::fill_rotamer_probability_tables,
         get_docstring_from_xml("fill_rotamer_probability_tables").c_str())
    .def("fit_ligand",
         &molecules_container_t::fit_ligand,
         nb::arg("imol_protein"), nb::arg("imol_map"), nb::arg("imol_ligand"), nb::arg("n_rmsd"), nb::arg("use_conformers"), nb::arg("n_conformers"),
         get_docstring_from_xml("fit_ligand").c_str())
    .def("fit_ligand_right_here",
         &molecules_container_t::fit_ligand_right_here,
         nb::arg("imol_protein"), nb::arg("imol_map"), nb::arg("imol_ligand"), nb::arg("x"), nb::arg("y"), nb::arg("z"), nb::arg("n_rmsd"), nb::arg("use_conformers"), nb::arg("n_conformers"),
         get_docstring_from_xml("fit_ligand_right_here").c_str())
    .def("fit_to_map_by_random_jiggle",
         &molecules_container_t::fit_to_map_by_random_jiggle,
         nb::arg("imol"), nb::arg("res_spec"), nb::arg("n_trials"), nb::arg("translation_scale_factor"),
         get_docstring_from_xml("fit_to_map_by_random_jiggle").c_str())
    .def("fit_to_map_by_random_jiggle_using_cid",
         &molecules_container_t::fit_to_map_by_random_jiggle_using_cid,
         nb::arg("imol"), nb::arg("cid"), nb::arg("n_trials"), nb::arg("translation_scale_factor"),
         get_docstring_from_xml("fit_to_map_by_random_jiggle_using_cid").c_str())
    .def("fit_to_map_by_random_jiggle_using_cid",
         &molecules_container_t::fit_to_map_by_random_jiggle_using_cid,
         nb::arg("imol"), nb::arg("cid"), nb::arg("n_trials"), nb::arg("translation_scale_factor"),
         get_docstring_from_xml("fit_to_map_by_random_jiggle_using_cid").c_str())
    .def("flip_peptide_using_cid",
         nb::overload_cast<int, const std::string&, const std::string&>(&molecules_container_t::flip_peptide_using_cid),
         get_docstring_from_xml("flip_peptide_using_cid").c_str())
    .def("flip_hand",
         &molecules_container_t::flip_hand,
         nb::arg("imol_map"),
         get_docstring_from_xml("flip_hand").c_str())
    .def("flood",
         &molecules_container_t::flood,
         nb::arg("imol_model"), nb::arg("imol_map"), nb::arg("n_rmsd"),
         get_docstring_from_xml("flood").c_str())
    .def("fourier_shell_correlation",
         &molecules_container_t::fourier_shell_correlation,
         nb::arg("imol_map_1"), nb::arg("imol_map_2"),
         get_docstring_from_xml("fourier_shell_correlation").c_str())
    .def("generate_self_restraints",
         &molecules_container_t::generate_self_restraints,
         nb::arg("imol"), nb::arg("local_dist_max"),
         get_docstring_from_xml("generate_self_restraints").c_str())
    .def("geometry_init_standard",
         &molecules_container_t::geometry_init_standard,
         get_docstring_from_xml("geometry_init_standard").c_str())
    .def("get_active_atom",
         &molecules_container_t::get_active_atom,
         nb::arg("x"), nb::arg("y"), nb::arg("z"), nb::arg("displayed_model_molecules_list"),
         get_docstring_from_xml("get_active_atom").c_str())
    .def("get_acedrg_atom_types",
         &molecules_container_t::get_acedrg_atom_types,
         nb::arg("compound_id"), nb::arg("imol_enc"),
         get_docstring_from_xml("get_acedrg_atom_types").c_str())
    .def("get_acedrg_atom_types_for_ligand",
         &molecules_container_t::get_acedrg_atom_types_for_ligand,
         nb::arg("imol"), nb::arg("residue_cid"),
         get_docstring_from_xml("get_acedrg_atom_types_for_ligand").c_str())
    .def("get_atom_differences",
         &molecules_container_t::get_atom_differences,
         nb::arg("imol1"), nb::arg("imol2"),
         get_docstring_from_xml("get_atom_differences").c_str())
    .def("get_atom_using_cid",
         &molecules_container_t::get_atom_using_cid,
         nb::arg("imol"), nb::arg("atom_cid"),
         get_docstring_from_xml("get_atom_using_cid").c_str())
    .def("get_atom_overlaps",
         &molecules_container_t::get_atom_overlaps,
         nb::arg("imol"),
         get_docstring_from_xml("get_atom_overlaps").c_str())
    .def("get_atom_overlap_score",
         &molecules_container_t::get_atom_overlap_score,
	 nb::arg("imol"),
         get_docstring_from_xml("get_atom_overlap_score").c_str())
    .def("get_bonds_mesh",
         &molecules_container_t::get_bonds_mesh,
         nb::arg("imol"), nb::arg("mode"), nb::arg("against_a_dark_background"), nb::arg("bond_width"), nb::arg("atom_radius_to_bond_width_ratio"), nb::arg("smoothness_factor"),
         get_docstring_from_xml("get_bonds_mesh").c_str())
    .def("get_bonds_mesh_for_selection_instanced",
         &molecules_container_t::get_bonds_mesh_for_selection_instanced,
         nb::arg("imol"), nb::arg("atom_selection_cid"), nb::arg("mode"), nb::arg("against_a_dark_background"), nb::arg("bond_width"), nb::arg("atom_radius_to_bond_width_ratio"), nb::arg("show_atoms_as_aniso_flag"), nb::arg("show_aniso_atoms_as_ortep_flag"), nb::arg("draw_hydrogen_atoms_flag"), nb::arg("smoothness_factor"),
         get_docstring_from_xml("get_bonds_mesh_for_selection_instanced").c_str())
    .def("get_bonds_mesh_instanced",
         &molecules_container_t::get_bonds_mesh_instanced,
         nb::arg("imol"), nb::arg("mode"), nb::arg("against_a_dark_background"), nb::arg("bond_width"), nb::arg("atom_radius_to_bond_width_ratio"), nb::arg("show_atoms_as_aniso_flag"), nb::arg("show_aniso_atoms_as_ortep_flag"), nb::arg("draw_hydrogen_atoms_flag"), nb::arg("smoothness_factor"),
         get_docstring_from_xml("get_bonds_mesh_instanced").c_str())
    .def("get_cell",
         &molecules_container_t::get_cell,
         nb::arg("imol"),
         get_docstring_from_xml("get_cell").c_str())
    .def("get_chains_in_model",
         &molecules_container_t::get_chains_in_model,
         nb::arg("imol"),
         get_docstring_from_xml("get_chains_in_model").c_str())
    .def("get_chemical_features_mesh",
         &molecules_container_t::get_chemical_features_mesh,
         nb::arg("imol"), nb::arg("cid"),
         get_docstring_from_xml("get_chemical_features_mesh").c_str())
    .def("get_colour_rules",
         &molecules_container_t::get_colour_rules,
         nb::arg("imol"),
         get_docstring_from_xml("get_colour_rules").c_str())
    .def("get_colour_table_for_blender",
         &molecules_container_t::get_colour_table_for_blender,
         nb::arg("imol"),
         get_docstring_from_xml("get_colour_table_for_blender").c_str())
    .def("get_density_at_position",
         &molecules_container_t::get_density_at_position,
         nb::arg("imol_map"), nb::arg("x"), nb::arg("y"), nb::arg("z"),
         get_docstring_from_xml("get_density_at_position").c_str())
    .def("get_dictionary_conformers",
         &molecules_container_t::get_dictionary_conformers,
         nb::arg("comp_id"), nb::arg("imol_enc"), nb::arg("remove_internal_clash_conformers"),
         get_docstring_from_xml("get_dictionary_conformers").c_str())
    .def("get_distances_between_atoms_of_residues",
         &molecules_container_t::get_distances_between_atoms_of_residues,
         nb::arg("imol"), nb::arg("cid_res_1"), nb::arg("cid_res_2"), nb::arg("dist_max"),
         get_docstring_from_xml("get_distances_between_atoms_of_residues").c_str())
    .def("get_gaussian_surface",
         &molecules_container_t::get_gaussian_surface,
         nb::arg("imol"), nb::arg("sigma"), nb::arg("contour_level"), nb::arg("box_radius"), nb::arg("grid_scale"), nb::arg("b_factor"),
         get_docstring_from_xml("get_gaussian_surface").c_str())
    .def("get_goodsell_style_mesh_instanced",
         &molecules_container_t::get_goodsell_style_mesh_instanced,
         nb::arg("imol"), nb::arg("colour_wheel_rotation_step"), nb::arg("saturation"), nb::arg("goodselliness"),
         get_docstring_from_xml("get_goodsell_style_mesh_instanced").c_str())
    .def("get_gphl_chem_comp_info",
         &molecules_container_t::get_gphl_chem_comp_info,
         nb::arg("compound_id"), nb::arg("imol_enc"),
         get_docstring_from_xml("get_gphl_chem_comp_info").c_str())
    .def("get_group_for_monomer",
         &molecules_container_t::get_group_for_monomer,
         nb::arg("residue_name"),
         get_docstring_from_xml("get_group_for_monomer").c_str())
    .def("get_groups_for_monomers",
         &molecules_container_t::get_groups_for_monomers,
         nb::arg("residue_names"),
         get_docstring_from_xml("get_groups_for_monomers").c_str())
    .def("get_hb_type",
         &molecules_container_t::get_hb_type,
         nb::arg("compound_id"), nb::arg("imol_enc"), nb::arg("atom_name"),
         get_docstring_from_xml("get_hb_type").c_str())
    .def("get_header_info",
         &molecules_container_t::get_header_info,
         nb::arg("imol"),
         get_docstring_from_xml("get_header_info").c_str())
    .def("get_h_bonds",&molecules_container_t::get_h_bonds,
         nb::arg("imol"), nb::arg("cid_str"), nb::arg("mcdonald_and_thornton_mode"))
    .def("get_HOLE",
         &molecules_container_t::get_HOLE,
         nb::arg("imol"), nb::arg("start_pos_x"), nb::arg("start_pos_y"), nb::arg("start_pos_z"), nb::arg("end_pos_x"), nb::arg("end_pos_y"), nb::arg("end_pos_z"),
         get_docstring_from_xml("get_HOLE").c_str())
    .def("get_imol_enc_any",
         &molecules_container_t::get_imol_enc_any,
         get_docstring_from_xml("get_imol_enc_any").c_str())
    .def("get_ligand_validation_vs_dictionary",
         &molecules_container_t::get_ligand_validation_vs_dictionary,
         nb::arg("imol"), nb::arg("ligand_cid"), nb::arg("include_non_bonded_contacts"),
         get_docstring_from_xml("get_ligand_validation_vs_dictionary").c_str())
    .def("get_ligand_distortion",
         &molecules_container_t::get_ligand_distortion,
         nb::arg("imol"), nb::arg("ligand_cid"), nb::arg("include_non_bonded_contacts"),
         get_docstring_from_xml("get_ligand_distortion").c_str())
    .def("get_lsq_matrix",
         &molecules_container_t::get_lsq_matrix,
         nb::arg("imol_ref"), nb::arg("imol_mov"), nb::arg("summary_to_screen"),
         get_docstring_from_xml("get_lsq_matrix").c_str())
    .def("get_map_contours_mesh",
         &molecules_container_t::get_map_contours_mesh,
         nb::arg("imol"), nb::arg("position_x"), nb::arg("position_y"), nb::arg("position_z"), nb::arg("radius"), nb::arg("contour_level"),
         get_docstring_from_xml("get_map_contours_mesh").c_str())
    .def("get_map_contours_mesh_using_other_map_for_colours",
         &molecules_container_t::get_map_contours_mesh_using_other_map_for_colours,
         nb::arg("imol_ref"), nb::arg("imol_map_for_colouring"), nb::arg("position_x"), nb::arg("position_y"), nb::arg("position_z"), nb::arg("radius"), nb::arg("contour_level"), nb::arg("other_map_for_colouring_min_value"), nb::arg("other_map_for_colouring_max_value"), nb::arg("invert_colour_ramp"),
         get_docstring_from_xml("get_map_contours_mesh_using_other_map_for_colours").c_str())
    .def("get_map_molecule_centre",
         &molecules_container_t::get_map_molecule_centre,
         nb::arg("imol"),
         get_docstring_from_xml("get_map_molecule_centre").c_str())
    .def("get_map_rmsd_approx",
         &molecules_container_t::get_map_rmsd_approx,
         nb::arg("imol_map"),
         get_docstring_from_xml("get_map_rmsd_approx").c_str())
    .def("get_map_weight",
         &molecules_container_t::get_map_weight,
         get_docstring_from_xml("get_map_weight").c_str())
    .def("get_mean_and_variance_of_density_for_non_water_atoms",
         &molecules_container_t::get_mean_and_variance_of_density_for_non_water_atoms,
         nb::arg("imol_coords"), nb::arg("imol_map"),
         get_docstring_from_xml("get_mean_and_variance_of_density_for_non_water_atoms").c_str())
    .def("get_median_temperature_factor",
         &molecules_container_t::get_median_temperature_factor,
         nb::arg("imol"),
         get_docstring_from_xml("get_median_temperature_factor").c_str())
    .def("get_missing_residue_ranges",
         &molecules_container_t::get_missing_residue_ranges,
         nb::arg("imol"),
         get_docstring_from_xml("get_missing_residue_ranges").c_str())
    .def("get_molecular_representation_mesh",
         &molecules_container_t::get_molecular_representation_mesh,
         nb::arg("imol"), nb::arg("cid"), nb::arg("colour_scheme"), nb::arg("style"), nb::arg("secondary_structure_usage_flag"),
         get_docstring_from_xml("get_molecular_representation_mesh").c_str())
    .def("get_molecule_centre",
         &molecules_container_t::get_molecule_centre,
         nb::arg("imol"),
         get_docstring_from_xml("get_molecule_centre").c_str())
    .def("get_molecule_name",
         &molecules_container_t::get_molecule_name,
         nb::arg("imol"),
         get_docstring_from_xml("get_molecule_name").c_str())
    .def("get_molecule_selection_as_json",
         &molecules_container_t::get_molecule_selection_as_json,
         nb::arg("imol"), nb::arg("cid"),
         get_docstring_from_xml("get_molecule_selection_as_json").c_str())
    .def("get_monomer",
         &molecules_container_t::get_monomer,
         nb::arg("monomer_name"),
         get_docstring_from_xml("get_monomer").c_str())
    .def("get_monomer_and_position_at",
         &molecules_container_t::get_monomer_and_position_at,
         nb::arg("comp_id"), nb::arg("imol"), nb::arg("x"), nb::arg("y"), nb::arg("z"),
         get_docstring_from_xml("get_monomer_and_position_at").c_str())
    .def("get_monomer_from_dictionary",
         &molecules_container_t::get_monomer_from_dictionary,
         nb::arg("comp_id"), nb::arg("imol"), nb::arg("idealised_flag"),
         get_docstring_from_xml("get_monomer_from_dictionary").c_str())
    .def("get_number_of_atoms",
         &molecules_container_t::get_number_of_atoms,
         nb::arg("imol"),
         get_docstring_from_xml("get_number_of_atoms").c_str())
    .def("get_number_of_atoms_in_residue",
         &molecules_container_t::get_number_of_atoms_in_residue,
         nb::arg("imol"), nb::arg("residue_cid"),
         get_docstring_from_xml("get_number_of_atoms_in_residue").c_str())
    .def("get_number_of_hydrogen_atoms",
         &molecules_container_t::get_number_of_hydrogen_atoms,
         nb::arg("imol"),
         get_docstring_from_xml("get_number_of_hydrogen_atoms").c_str())
    .def("get_number_of_molecules",
         &molecules_container_t::get_number_of_molecules,
         get_docstring_from_xml("get_number_of_molecules").c_str())
    .def("get_number_of_map_sections",
         &molecules_container_t::get_number_of_map_sections,
         nb::arg("imol_map"), nb::arg("axis_id"),
         get_docstring_from_xml("get_number_of_map_sections").c_str())
    .def("get_octahemisphere",
         &molecules_container_t::get_octahemisphere,
         nb::arg("n_divisions"),
         get_docstring_from_xml("get_octahemisphere").c_str())
    .def("get_overlaps_for_ligand",
         &molecules_container_t::get_overlaps_for_ligand,
         nb::arg("imol"), nb::arg("ligand_cid"),
         get_docstring_from_xml("get_overlaps_for_ligand").c_str())
    .def("get_pucker_analysis_info",
         &molecules_container_t::get_pucker_analysis_info,
         nb::arg("imol"),
         get_docstring_from_xml("get_pucker_analysis_info").c_str())
    .def("get_q_score",
         &molecules_container_t::get_q_score,
         nb::arg("imol_model"), nb::arg("imol_map"),
         get_docstring_from_xml("get_q_score").c_str())
    .def("get_q_score_for_cid",
         &molecules_container_t::get_q_score_for_cid,
         nb::arg("imol_model"), nb::arg("cid"), nb::arg("imol_map"),
         get_docstring_from_xml("get_q_score_for_cid").c_str())
    .def("get_r_factor_stats",
         &molecules_container_t::get_r_factor_stats,
         get_docstring_from_xml("get_r_factor_stats").c_str())
    .def("get_radius_of_gyration",
         &molecules_container_t::get_radius_of_gyration,
         get_docstring_from_xml("get_radius_of_gyration").c_str())
    .def("get_rama_plot_restraints_weight",
         &molecules_container_t::get_rama_plot_restraints_weight,
         get_docstring_from_xml("get_rama_plot_restraints_weight").c_str())
    // maybe these will work in future - or maybe just delete them
    // .def("get_rdkit_mol",
    //      &molecules_container_t::get_rdkit_mol,
    //      nb::arg("res_name"), nb::arg("imol_enc"),
    //      get_docstring_from_xml("get_rdkit_mol").c_str())
    // .def("get_rdkit_mol_shared",
    //      &molecules_container_t::get_rdkit_mol_shared,
    //      nb::arg("res_name"), nb::arg("imol_enc"),
    //      get_docstring_from_xml("get_rdkit_mol_shared").c_str())
    .def("get_rdkit_mol_pickle_base64",
         &molecules_container_t::get_rdkit_mol_pickle_base64,
         nb::arg("res_name"), nb::arg("imol_enc"),
         get_docstring_from_xml("get_rdkit_mol_pickle_base64").c_str())
    .def("get_ramachandran_validation_markup_mesh",
         &molecules_container_t::get_ramachandran_validation_markup_mesh,
         nb::arg("imol"),
         get_docstring_from_xml("get_ramachandran_validation_markup_mesh").c_str())
    //Using allow_raw_pointers(). Perhaps suggests we need to do something different from exposing mmdb pointers to JS.
    .def("get_residue_average_position",
         &molecules_container_t::get_residue_average_position,
         nb::arg("imol"), nb::arg("cid"),
         get_docstring_from_xml("get_residue_average_position").c_str())
    .def("get_residue_CA_position",
         &molecules_container_t::get_residue_CA_position,
         nb::arg("imol"), nb::arg("cid"),
         get_docstring_from_xml("get_residue_CA_position").c_str())
    .def("get_residue_name",
         &molecules_container_t::get_residue_name,
         nb::arg("imol"), nb::arg("chain_id"), nb::arg("res_no"), nb::arg("ins_code"),
         get_docstring_from_xml("get_residue_name").c_str())
    .def("get_residue_names_with_no_dictionary",
         &molecules_container_t::get_residue_names_with_no_dictionary,
         nb::arg("imol"),
         get_docstring_from_xml("get_residue_names_with_no_dictionary").c_str())
    .def("get_residue_sidechain_average_position",
         &molecules_container_t::get_residue_sidechain_average_position,
         nb::arg("imol"), nb::arg("cid"),
         get_docstring_from_xml("get_residue_sidechain_average_position").c_str())
    .def("get_residue_using_cid",
         &molecules_container_t::get_residue_using_cid,
         nb::arg("imol"), nb::arg("cid"),
         get_docstring_from_xml("get_residue_using_cid").c_str())
    .def("get_residues_near_residue",
         &molecules_container_t::get_residues_near_residue,
         nb::arg("imol"), nb::arg("residue_cid"), nb::arg("dist"),
         get_docstring_from_xml("get_residues_near_residue").c_str())
    .def("get_rotamer_dodecs",
         &molecules_container_t::get_rotamer_dodecs,
         nb::arg("imol"),
         get_docstring_from_xml("get_rotamer_dodecs").c_str())
    .def("get_rotamer_dodecs_instanced",
         &molecules_container_t::get_rotamer_dodecs_instanced,
         nb::arg("imol"),
         get_docstring_from_xml("get_rotamer_dodecs_instanced").c_str())
    .def("get_single_letter_codes_for_chain",
         &molecules_container_t::get_single_letter_codes_for_chain,
         nb::arg("imol"), nb::arg("chain_id"),
         get_docstring_from_xml("get_single_letter_codes_for_chain").c_str())
    .def("get_spherical_variance",
         &molecules_container_t::get_spherical_variance,
         nb::arg("imol_map"), nb::arg("imol_model"),
         nb::arg("atom_cid"), nb::arg("mean_density_other_atoms"),
         get_docstring_from_xml("get_spherical_variance").c_str())
    .def("get_sum_density_for_atoms_in_residue",
         &molecules_container_t::get_sum_density_for_atoms_in_residue,
         nb::arg("imol"), nb::arg("cid"), nb::arg("atom_names"), nb::arg("imol_map"),
         get_docstring_from_xml("get_sum_density_for_atoms_in_residue").c_str())
    .def("get_svg_for_2d_ligand_environment_view",
         &molecules_container_t::get_svg_for_2d_ligand_environment_view,
         nb::arg("imol"), nb::arg("residue_cid"), nb::arg("add_key"),
         get_docstring_from_xml("get_svg_for_2d_ligand_environment_view").c_str())
    .def("get_svg_for_residue_type",
         &molecules_container_t::get_svg_for_residue_type,
         nb::arg("imol"), nb::arg("comp_id"), nb::arg("use_rdkit_svg"), nb::arg("background_type"),
         get_docstring_from_xml("get_svg_for_residue_type").c_str())
    .def("get_symmetry",
         &molecules_container_t::get_symmetry,
         nb::arg("imol"), nb::arg("symmetry_search_radius"), nb::arg("centre_x"), nb::arg("centre_y"), nb::arg("centre_z"),
         get_docstring_from_xml("get_symmetry").c_str())
    .def("get_temperature_factor_of_atom",
         &molecules_container_t::get_temperature_factor_of_atom,
         nb::arg("imol"), nb::arg("atom_cid"),
         get_docstring_from_xml("get_temperature_factor_of_atom").c_str())
    .def("get_torsion",
         &molecules_container_t::get_torsion,
         nb::arg("imol"), nb::arg("cid"), nb::arg("atom_names"),
         get_docstring_from_xml("get_torsion").c_str())
    .def("get_torsion_restraints_weight",
         &molecules_container_t::get_torsion_restraints_weight,
         get_docstring_from_xml("get_torsion_restraints_weight").c_str())
    .def("get_triangles_for_blender",
         &molecules_container_t::get_triangles_for_blender,
         nb::arg("imol"),
         get_docstring_from_xml("get_triangles_for_blender").c_str())
    .def("get_types_in_molecule",
         &molecules_container_t::get_types_in_molecule,
         nb::arg("imol"),
         get_docstring_from_xml("get_types_in_molecule").c_str())
    .def("get_use_gemmi",
         &molecules_container_t::get_use_gemmi,
         get_docstring_from_xml("get_use_gemmi").c_str())
    .def("get_use_rama_plot_restraints",
         &molecules_container_t::get_use_rama_plot_restraints,
         get_docstring_from_xml("get_use_rama_plot_restraints").c_str())
    .def("get_use_torsion_restraints",
         &molecules_container_t::get_use_torsion_restraints,
         get_docstring_from_xml("get_use_torsion_restraints").c_str())
    .def("get_validation_vs_dictionary_for_selection",
         &molecules_container_t::get_validation_vs_dictionary_for_selection,
         nb::arg("imol"), nb::arg("selection_cid"), nb::arg("include_non_bonded_contacts"),
         get_docstring_from_xml("get_validation_vs_dictionary_for_selection").c_str())
    .def("get_vertices_for_blender",
         &molecules_container_t::get_vertices_for_blender,
         nb::arg("imol"),
         get_docstring_from_xml("get_vertices_for_blender").c_str())
    .def("go_to_blob",
         &molecules_container_t::go_to_blob,
         nb::arg("x1"), nb::arg("y1"), nb::arg("z1"), nb::arg("x2"), nb::arg("y2"), nb::arg("z2"), nb::arg("contour_level"),
         get_docstring_from_xml("go_to_blob").c_str())
    .def("import_cif_dictionary",
         &molecules_container_t::import_cif_dictionary,
         nb::arg("cif_file_name"), nb::arg("imol_enc"),
         get_docstring_from_xml("import_cif_dictionary").c_str())
    .def("init_refinement_of_molecule_as_fragment_based_on_reference",
         &molecules_container_t::init_refinement_of_molecule_as_fragment_based_on_reference,
         nb::arg("imol_frag"), nb::arg("imol_ref"), nb::arg("imol_map"),
         get_docstring_from_xml("init_refinement_of_molecule_as_fragment_based_on_reference").c_str())
    .def("is_a_difference_map",
         &molecules_container_t::is_a_difference_map,
         nb::arg("imol_map"),
         get_docstring_from_xml("is_a_difference_map").c_str())
    .def("is_valid_map_molecule",
         &molecules_container_t::is_valid_map_molecule,
         nb::arg("imol_map"),
         get_docstring_from_xml("is_valid_map_molecule").c_str())
    .def("is_valid_model_molecule",
         &molecules_container_t::is_valid_model_molecule,
         nb::arg("imol"),
         get_docstring_from_xml("is_valid_model_molecule").c_str())
    .def("jed_flip",
         nb::overload_cast<int, const std::string&, bool> (&molecules_container_t::jed_flip),
         nb::arg("imol"), nb::arg("atom_cid"), nb::arg("invert_selection"),
         get_docstring_from_xml("jed_flip").c_str())
    .def("lsq_superpose",
         &molecules_container_t::lsq_superpose,
	 nb::arg("imol_ref"), nb::arg("imol_mov"),
         get_docstring_from_xml("lsq_superpose").c_str())
    .def("make_mask",
         &molecules_container_t::make_mask,
         nb::arg("imol_map_ref"), nb::arg("imol_model"), nb::arg("atom_selection_cid"), nb::arg("radius"),
         get_docstring_from_xml("make_mask").c_str())
    .def("make_mesh_for_bonds_for_blender",
         &molecules_container_t::make_mesh_for_bonds_for_blender,
         nb::arg("imol"), nb::arg("mode"), nb::arg("against_a_dark_background"), nb::arg("bond_width"),
         nb::arg("atom_radius_to_bond_width_ratio"), nb::arg("smoothness_factor"),
         get_docstring_from_xml("make_mesh_for_bonds_for_blender").c_str())
    .def("make_mesh_for_gaussian_surface_for_blender",
         &molecules_container_t::make_mesh_for_gaussian_surface_for_blender,
         nb::arg("imol"), nb::arg("sigma"), nb::arg("contour_level"),
         nb::arg("box_radius"), nb::arg("grid_scale"), nb::arg("b_factor"),
         get_docstring_from_xml("make_mesh_for_gaussian_surface_for_blender").c_str())
    .def("make_mesh_for_goodsell_style_for_blender",
         &molecules_container_t::make_mesh_for_goodsell_style_for_blender,
         nb::arg("imol"), nb::arg("colour_wheel_rotation_step"),
         nb::arg("saturation"), nb::arg("goodselliness"),
         get_docstring_from_xml("make_mesh_for_goodsell_style_for_blender").c_str())
    .def("make_mesh_for_map_contours_for_blender",
         &molecules_container_t::make_mesh_for_map_contours_for_blender,
         nb::arg("imol"), nb::arg("x"), nb::arg("y"), nb::arg("z"), nb::arg("level"), nb::arg("radius"),
         get_docstring_from_xml("make_mesh_for_map_contours_for_blender").c_str())
    .def("make_mesh_for_molecular_representation_for_blender",
         &molecules_container_t::make_mesh_for_molecular_representation_for_blender,
         nb::arg("imol"), nb::arg("cid"), nb::arg("colour_scheme"), nb::arg("style"),
         nb::arg("secondary_structure_usage_flag"),
         get_docstring_from_xml("make_mesh_for_molecular_representation_for_blender").c_str())
    .def("mask_map_by_atom_selection",
         &molecules_container_t::mask_map_by_atom_selection,
         nb::arg("imol_coords"), nb::arg("imol_map"), nb::arg("cid"),
         nb::arg("atom_radius"), nb::arg("invert_flag"),
         get_docstring_from_xml("mask_map_by_atom_selection").c_str())
    .def("make_power_scaled_map",
         &molecules_container_t::make_power_scaled_map,
         nb::arg("imol_ref"), nb::arg("imol_map_for_scaling"),
         get_docstring_from_xml("make_power_scaled_map").c_str())
    .def("merge_molecules",
         nb::overload_cast<int,const std::string &>(&molecules_container_t::merge_molecules),
         nb::arg("imol"), nb::arg("list_of_other_molecules"),
         get_docstring_from_xml("merge_molecules").c_str())
    .def("minimize_energy",
         &molecules_container_t::minimize_energy,
         nb::arg("imol"), nb::arg("atom_selection_cid"), nb::arg("n_cycles"),
         nb::arg("do_rama_plot_restraints"), nb::arg("rama_plot_weight"),
         nb::arg("do_torsion_restraints"), nb::arg("torsion_weight"),
         nb::arg("refinement_is_quiet"),
         get_docstring_from_xml("minimize_energy").c_str())
    .def("minimize",
         &molecules_container_t::minimize,
         nb::arg("imol"), nb::arg("atom_selection_cid"), nb::arg("n_cycles"),
         nb::arg("do_rama_plot_restraints"), nb::arg("rama_plot_weight"),
         nb::arg("do_torsion_restraints"), nb::arg("torsion_weight"),
         nb::arg("refinement_is_quiet"),
         get_docstring_from_xml("minimize").c_str())
    .def("mmcif_tests",
         &molecules_container_t::mmcif_tests,
         nb::arg("last_test_only"),
         get_docstring_from_xml("mmcif_tests").c_str())
    .def("mmrrcc",
         &molecules_container_t::mmrrcc,
         get_docstring_from_xml("mmrrcc").c_str())
    .def("move_molecule_to_new_centre",
         &molecules_container_t::move_molecule_to_new_centre,
         nb::arg("imol"), nb::arg("x"), nb::arg("y"), nb::arg("z"),
         get_docstring_from_xml("move_molecule_to_new_centre").c_str())
    .def("multiply_residue_temperature_factors",
         &molecules_container_t::multiply_residue_temperature_factors,
         nb::arg("imol"), nb::arg("cid"), nb::arg("factor"),
         get_docstring_from_xml("multiply_residue_temperature_factors").c_str())
    .def("mutate",
         &molecules_container_t::mutate,
         nb::arg("imol"), nb::arg("cid"), nb::arg("new_residue_type"),
         get_docstring_from_xml("mutate").c_str())
    .def("new_positions_for_atoms_in_residues",
         &molecules_container_t::new_positions_for_atoms_in_residues,
         nb::arg("imol"), nb::arg("moved_residues"),
         get_docstring_from_xml("new_positions_for_atoms_in_residues").c_str())
    .def("new_positions_for_residue_atoms",
         &molecules_container_t::new_positions_for_residue_atoms,
         nb::arg("imol"), nb::arg("residue_cid"), nb::arg("moved_atoms"),
         get_docstring_from_xml("new_positions_for_residue_atoms").c_str())
    .def("non_standard_residue_types_in_model",
         &molecules_container_t::non_standard_residue_types_in_model,
         nb::arg("imol"),
         get_docstring_from_xml("non_standard_residue_types_in_model").c_str())
    .def("package_version",
         &molecules_container_t::package_version,
         get_docstring_from_xml("package_version").c_str())
    .def("partition_map_by_chain",
         &molecules_container_t::partition_map_by_chain,
         nb::arg("imol_map"), nb::arg("imol_model"),
         get_docstring_from_xml("partition_map_by_chain").c_str())
    .def("pepflips_using_difference_map",
         &molecules_container_t::pepflips_using_difference_map,
         nb::arg("imol_coords"), nb::arg("imol_difference_map"), nb::arg("n_sigma"),
         get_docstring_from_xml("pepflips_using_difference_map").c_str())
    .def("peptide_omega_analysis",
         &molecules_container_t::peptide_omega_analysis,
         nb::arg("imol_model"),
         get_docstring_from_xml("peptide_omega_analysis").c_str())
    .def("print_secondary_structure_info",
         &molecules_container_t::print_secondary_structure_info,
         nb::arg("imol"),
         get_docstring_from_xml("print_secondary_structure_info").c_str())
    .def("rail_points_total",
         &molecules_container_t::rail_points_total,
         get_docstring_from_xml("rail_points_total").c_str())
    .def("ramachandran_analysis",
         &molecules_container_t::ramachandran_analysis,
         nb::arg("imol_model"),
         get_docstring_from_xml("ramachandran_analysis").c_str())
    .def("ramachandran_validation",
         &molecules_container_t::ramachandran_validation,
         nb::arg("imol"),
         get_docstring_from_xml("ramachandran_validation").c_str())
    .def("read_coordinates",
         &molecules_container_t::read_coordinates,
         nb::arg("file_name"),
         get_docstring_from_xml("read_coordinates").c_str())
    .def("read_ccp4_map",
         &molecules_container_t::read_ccp4_map,
         nb::arg("file_name"), nb::arg("is_a_difference_map"),
         get_docstring_from_xml("read_ccp4_map").c_str())
    .def("read_extra_restraints",
         &molecules_container_t::read_extra_restraints,
         nb::arg("imol"), nb::arg("file_name"),
         get_docstring_from_xml("read_extra_restraints").c_str())
    .def("read_mtz",
         &molecules_container_t::read_mtz,
         nb::arg("file_name"), nb::arg("f"), nb::arg("phi"), nb::arg("weight"),
         nb::arg("use_weight"), nb::arg("is_a_difference_map"),
         get_docstring_from_xml("read_mtz").c_str())
    .def("read_pdb",
         &molecules_container_t::read_pdb,
         nb::arg("file_name"),
         get_docstring_from_xml("read_pdb").c_str())
    .def("read_small_molecule_cif",
         &molecules_container_t::read_small_molecule_cif,
         nb::arg("file_name"),
         get_docstring_from_xml("read_small_molecule_cif").c_str())
    .def("redo",
         &molecules_container_t::redo,
         nb::arg("imol"),
         get_docstring_from_xml("redo").c_str())
    .def("refine",
         &molecules_container_t::refine,
         nb::arg("imol"), nb::arg("n_cycles"),
         get_docstring_from_xml("refine").c_str())
    .def("refine_residue_range",
         &molecules_container_t::refine_residue_range,
         nb::arg("imol"), nb::arg("chain_id"), nb::arg("res_no_start"),
         nb::arg("res_no_end"), nb::arg("n_cycles"),
         get_docstring_from_xml("refine_residue_range").c_str())
    .def("refine_residues",
         &molecules_container_t::refine_residues,
         nb::arg("imol"), nb::arg("chain_id"), nb::arg("res_no"), nb::arg("ins_code"),
         nb::arg("alt_conf"), nb::arg("mode"), nb::arg("n_cycles"),
         get_docstring_from_xml("refine_residues").c_str())
    .def("refine_residues_using_atom_cid",
         &molecules_container_t::refine_residues_using_atom_cid,
         nb::arg("imol"), nb::arg("cid"), nb::arg("mode"), nb::arg("n_cycles"),
         get_docstring_from_xml("refine_residues_using_atom_cid").c_str())
    .def("regen_map",
         &molecules_container_t::regen_map,
         nb::arg("imol_map"), nb::arg("imol_maps"), nb::arg("scales"),
         get_docstring_from_xml("regen_map").c_str())
    .def("replace_fragment",
         &molecules_container_t::replace_fragment,
         nb::arg("imol_base"), nb::arg("imol_reference"), nb::arg("atom_selection"),
         get_docstring_from_xml("replace_fragment").c_str())
    .def("replace_map_by_mtz_from_file",
         &molecules_container_t::replace_map_by_mtz_from_file,
         nb::arg("imol"), nb::arg("file_name"), nb::arg("f"), nb::arg("phi"),
         nb::arg("weight"), nb::arg("use_weight"),
         get_docstring_from_xml("replace_map_by_mtz_from_file").c_str())
    .def("replace_molecule_by_model_from_file",
         &molecules_container_t::replace_molecule_by_model_from_file,
         nb::arg("imol"), nb::arg("pdb_file_name"),
         get_docstring_from_xml("replace_molecule_by_model_from_file").c_str())
    .def("replace_residue",
         &molecules_container_t::replace_residue,
         nb::arg("imol"), nb::arg("residue_cid"), nb::arg("new_residue_type"), nb::arg("imol_enc"),
         get_docstring_from_xml("replace_residue").c_str())
    .def("residues_with_missing_atoms",
         &molecules_container_t::residues_with_missing_atoms,
         nb::arg("imol"),
         get_docstring_from_xml("residues_with_missing_atoms").c_str())
    .def("rigid_body_fit",
         &molecules_container_t::rigid_body_fit,
         nb::arg("imol"), nb::arg("multi_cid"), nb::arg("imol_map"),
         get_docstring_from_xml("rigid_body_fit").c_str())
    .def("rotamer_analysis",
         &molecules_container_t::rotamer_analysis,
         nb::arg("imol_model"),
         get_docstring_from_xml("rotamer_analysis").c_str())
    .def("rotate_around_bond",
         &molecules_container_t::rotate_around_bond,
         nb::arg("imol"), nb::arg("residue_cid"),
         nb::arg("atom_name_1"), nb::arg("atom_name_2"),
         nb::arg("atom_name_3"), nb::arg("atom_name_4"),
         nb::arg("torsion_angle"),
         get_docstring_from_xml("rotate_around_bond").c_str())
    .def("scale_map",
         &molecules_container_t::scale_map,
         nb::arg("imol_map"), nb::arg("scale_factor"),
         get_docstring_from_xml("scale_map").c_str())
    .def("servalcat_refine_xray",
         &molecules_container_t::servalcat_refine_xray,
         nb::arg("imol"), nb::arg("imol_map"), nb::arg("output_prefix"),
         get_docstring_from_xml("servalcat_refine_xray").c_str())
    .def("servalcat_refine_xray_with_keywords",
         &molecules_container_t::servalcat_refine_xray_with_keywords,
         get_docstring_from_xml("servalcat_refine_xray_with_keywords").c_str())
    .def("set_add_waters_sigma_cutoff",
         &molecules_container_t::set_add_waters_sigma_cutoff,
         nb::arg("d"),
         get_docstring_from_xml("set_add_waters_sigma_cutoff").c_str())
    .def("set_add_waters_variance_limit",
         &molecules_container_t::set_add_waters_variance_limit,
         nb::arg("d"),
         get_docstring_from_xml("set_add_waters_variance_limit").c_str())
    .def("set_add_waters_water_to_protein_distance_lim_min",
         &molecules_container_t::set_add_waters_water_to_protein_distance_lim_min,
         nb::arg("d"),
         get_docstring_from_xml("set_add_waters_water_to_protein_distance_lim_min").c_str())
    .def("set_add_waters_water_to_protein_distance_lim_max",
         &molecules_container_t::set_add_waters_water_to_protein_distance_lim_max,
         nb::arg("d"),
         get_docstring_from_xml("set_add_waters_water_to_protein_distance_lim_max").c_str())
    .def("set_colour_map_for_map_coloured_by_other_map",
         &molecules_container_t::set_colour_map_for_map_coloured_by_other_map,
         nb::arg("colour_table"),
         get_docstring_from_xml("set_colour_map_for_map_coloured_by_other_map").c_str())
    .def("set_colour_wheel_rotation_base",
         &molecules_container_t::set_colour_wheel_rotation_base,
         nb::arg("imol"), nb::arg("r"),
         get_docstring_from_xml("set_colour_wheel_rotation_base").c_str())
    .def("set_draw_missing_residue_loops",
         &molecules_container_t::set_draw_missing_residue_loops,
         nb::arg("state"),
         get_docstring_from_xml("set_draw_missing_residue_loops").c_str())
    .def("set_draw_missing_residue_loops",
         &molecules_container_t::set_draw_missing_residue_loops,
         nb::arg("state"),
         get_docstring_from_xml("set_draw_missing_residue_loops").c_str())
    .def("set_gltf_pbr_metalicity_factor",
         &molecules_container_t::set_gltf_pbr_metalicity_factor,
         nb::arg("imol"), nb::arg("metalicity"),
         get_docstring_from_xml("set_gltf_pbr_metalicity_factor").c_str())
    .def("set_gltf_pbr_roughness_factor",
         &molecules_container_t::set_gltf_pbr_roughness_factor,
         nb::arg("imol"), nb::arg("roughness_factor"),
         get_docstring_from_xml("set_gltf_pbr_roughness_factor").c_str())
    .def("set_imol_refinement_map",
         &molecules_container_t::set_imol_refinement_map,
         nb::arg("i"),
         get_docstring_from_xml("set_imol_refinement_map").c_str())
    .def("set_make_backups",
         &molecules_container_t::set_make_backups,
         nb::arg("state"),
         get_docstring_from_xml("set_make_backups").c_str())
    .def("set_logging_file",
         &molecules_container_t::set_logging_file,
         nb::arg("file_name"),
         get_docstring_from_xml("set_logging_level").c_str())
    .def("set_logging_level",
         &molecules_container_t::set_logging_level,
         nb::arg("level"),
         get_docstring_from_xml("set_logging_level").c_str())
    .def("set_map_sampling_rate",
         &molecules_container_t::set_map_sampling_rate,
         nb::arg("msr"),
         get_docstring_from_xml("set_map_sampling_rate").c_str())
    .def("set_map_weight",
         &molecules_container_t::set_map_weight,
         nb::arg("w"),
         get_docstring_from_xml("set_map_weight").c_str())
    .def("set_max_number_of_threads",
         &molecules_container_t::set_max_number_of_threads,
         nb::arg("n_threads"),
         get_docstring_from_xml("set_max_number_of_threads").c_str())
    .def("set_molecule_name",
         &molecules_container_t::set_molecule_name,
         nb::arg("imol"), nb::arg("new_name"),
         get_docstring_from_xml("set_molecule_name").c_str())
    .def("set_occupancy",
         &molecules_container_t::set_occupancy,
         nb::arg("imol"), nb::arg("cid"), nb::arg("occ_new"),
         get_docstring_from_xml("set_occupancy").c_str())
    .def("set_rama_plot_restraints_weight",
         &molecules_container_t::set_rama_plot_restraints_weight,
         nb::arg("f"),
         get_docstring_from_xml("set_rama_plot_restraints_weight").c_str())
    .def("set_refinement_is_verbose",
         &molecules_container_t::set_refinement_is_verbose,
         nb::arg("state"),
         get_docstring_from_xml("set_refinement_is_verbose").c_str())
    .def("set_refinement_geman_mcclure_alpha",
         &molecules_container_t::set_refinement_geman_mcclure_alpha,
         nb::arg("a"),
         get_docstring_from_xml("set_refinement_geman_mcclure_alpha").c_str())
    .def("set_residue_to_rotamer_number",
         &molecules_container_t::set_residue_to_rotamer_number,
         nb::arg("imol"), nb::arg("residue_cid"),
         nb::arg("alt_conf"), nb::arg("rotamer_number"),
         get_docstring_from_xml("set_residue_to_rotamer_number").c_str())
    .def("set_show_timings",
         &molecules_container_t::set_show_timings,
         nb::arg("s"),
         get_docstring_from_xml("set_show_timings").c_str())
    .def("set_temperature_factors_using_cid",
         &molecules_container_t::set_temperature_factors_using_cid,
         nb::arg("imol"), nb::arg("cid"), nb::arg("temp_fact"),
         get_docstring_from_xml("set_temperature_factors_using_cid").c_str())
    .def("set_torsion_restraints_weight",
         &molecules_container_t::set_torsion_restraints_weight,
         nb::arg("f"),
         get_docstring_from_xml("set_torsion_restraints_weight").c_str())
    .def("set_use_gemmi",
         &molecules_container_t::set_use_gemmi,
         nb::arg("state"),
         get_docstring_from_xml("set_use_gemmi").c_str())
    .def("set_use_rama_plot_restraints",
         &molecules_container_t::set_use_rama_plot_restraints,
         nb::arg("state"),
         get_docstring_from_xml("set_use_rama_plot_restraints").c_str())
    .def("set_use_torsion_restraints",
         &molecules_container_t::set_use_torsion_restraints,
         nb::arg("state"),
         get_docstring_from_xml("set_use_torsion_restraints").c_str())
    .def("set_user_defined_atom_colour_by_selection",
         &molecules_container_t::set_user_defined_atom_colour_by_selection,
         nb::arg("imol"), nb::arg("indexed_residues_cids"), nb::arg("colour_applies_to_non_carbon_atoms_also"),
         get_docstring_from_xml("set_user_defined_atom_colour_by_selection").c_str())
    .def("set_user_defined_bond_colours",
         &molecules_container_t::set_user_defined_bond_colours,
         nb::arg("imol"), nb::arg("colour_map"),
         get_docstring_from_xml("set_user_defined_bond_colours").c_str())
    .def("sfcalc_genmap",
         &molecules_container_t::sfcalc_genmap,
         nb::arg("imol_model"), nb::arg("imol_map_with_data_attached"), nb::arg("imol_updating_difference_map"),
         get_docstring_from_xml("sfcalc_genmap").c_str())
    .def("sfcalc_genmaps_using_bulk_solvent",
         &molecules_container_t::sfcalc_genmaps_using_bulk_solvent,
         nb::arg("imol_model"), nb::arg("imol_2fofc_map"), nb::arg("imol_updating_difference_map"), nb::arg("imol_map_with_data_attached"),
         get_docstring_from_xml("sfcalc_genmaps_using_bulk_solvent").c_str())
    .def("sharpen_blur_map",
         &molecules_container_t::sharpen_blur_map,
         nb::arg("imol_map"), nb::arg("b_factor"), nb::arg("in_place_flag"),
         get_docstring_from_xml("sharpen_blur_map").c_str())
    .def("sharpen_blur_map_with_resample",
         &molecules_container_t::sharpen_blur_map_with_resample,
         nb::arg("imol_map"), nb::arg("b_factor"), nb::arg("resample_factor"), nb::arg("in_place_flag"),
         get_docstring_from_xml("sharpen_blur_map_with_resample").c_str())
    .def("side_chain_180",
         nb::overload_cast<int, const std::string&>                         (&molecules_container_t::side_chain_180),
         nb::arg("imol"), nb::arg("atom_cid"),
         get_docstring_from_xml("side_chain_180").c_str())
    .def("split_multi_model_molecule",
         &molecules_container_t::split_multi_model_molecule,
         nb::arg("imol"),
         get_docstring_from_xml("split_multi_model_molecule").c_str())
    .def("split_residue_using_map",
         &molecules_container_t::split_residue_using_map,
         nb::arg("imol"), nb::arg("residue_cid"), nb::arg("imol_diff_map"),
         get_docstring_from_xml("split_residue_using_map").c_str())
    .def("test_function",
         &molecules_container_t::test_function,
         nb::arg("s"),
         get_docstring_from_xml("test_function").c_str())
    .def("test_origin_cube",
         &molecules_container_t::test_origin_cube,
         get_docstring_from_xml("test_origin_cube").c_str())
    .def("transform_map_using_lsq_matrix",
         &molecules_container_t::transform_map_using_lsq_matrix,
         nb::arg("imol_map"), nb::arg("lsq_matrix"), nb::arg("x"), nb::arg("y"), nb::arg("z"), nb::arg("radius"),
         get_docstring_from_xml("transform_map_using_lsq_matrix").c_str())
    .def("try_read_dictionaries_for_new_residue_types",
         &molecules_container_t::try_read_dictionaries_for_new_residue_types,
         nb::arg("imol"),
         get_docstring_from_xml("try_read_dictionaries_for_new_residue_types").c_str())
    .def("undo",
         &molecules_container_t::undo,
         nb::arg("imol"),
         get_docstring_from_xml("undo").c_str())
    .def("unmodelled_blobs",
         &molecules_container_t::unmodelled_blobs,
         nb::arg("imol_model"), nb::arg("imol_map"), nb::arg("rmsd_cut_off"),
         get_docstring_from_xml("unmodelled_blobs").c_str())
    .def("write_coordinates",
         &molecules_container_t::write_coordinates,
         nb::arg("imol"), nb::arg("file_name"),
         get_docstring_from_xml("write_coordinates").c_str())
    .def("write_map",
         &molecules_container_t::write_map,
         nb::arg("imol"), nb::arg("file_name"),
         get_docstring_from_xml("write_map").c_str())
    .def("write_png",
         &molecules_container_t::write_png,
         nb::arg("compound_id"), nb::arg("imol"), nb::arg("file_name"),
         get_docstring_from_xml("write_png").c_str())
    ;
    nb::class_<coot::chain_mutation_info_container_t>(m,"chain_mutation_info_container_t")
      .def_ro("chain_id",         &coot::chain_mutation_info_container_t::chain_id)
      .def_ro("alignedS",         &coot::chain_mutation_info_container_t::alignedS)
      .def_ro("alignedT",         &coot::chain_mutation_info_container_t::alignedT)
      .def_ro("alignedS_label",   &coot::chain_mutation_info_container_t::alignedS_label)
      .def_ro("alignedT_label",   &coot::chain_mutation_info_container_t::alignedT_label)
      .def_ro("alignment_string", &coot::chain_mutation_info_container_t::alignment_string)
      .def_ro("alignment_score",  &coot::chain_mutation_info_container_t::alignment_score)
      .def_ro("insertions",       &coot::chain_mutation_info_container_t::insertions)
      .def_ro("deletions",        &coot::chain_mutation_info_container_t::deletions)
      .def_ro("mutations",        &coot::chain_mutation_info_container_t::mutations)
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
    nb::enum_<coot::graph_data_type>(m, "graph_data_type")
       .value("Unset",          coot::graph_data_type::UNSET)
       .value("Density",        coot::graph_data_type::DENSITY)
       .value("Distortion",     coot::graph_data_type::DISTORTION)
       .value("Energy",         coot::graph_data_type::ENERGY)
       .value("Probability",    coot::graph_data_type::PROBABILITY)
       .value("Correlation",    coot::graph_data_type::CORRELATION)
       .value("LogProbability", coot::graph_data_type::LOG_PROBABILITY)
       .value("TorsionAngle",   coot::graph_data_type::TORSION_ANGLE)
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
    nb::enum_<coot::restraint_type_t>(m, "restraint_type")
       .value("Bond", coot::restraint_type_t::BOND_RESTRAINT)
       .value("Angle", coot::restraint_type_t::ANGLE_RESTRAINT)
       .value("Torsion", coot::restraint_type_t::TORSION_RESTRAINT)
       .value("Plane", coot::restraint_type_t::PLANE_RESTRAINT)
       .value("Non-Bonded-Contact", coot::restraint_type_t::NON_BONDED_CONTACT_RESTRAINT)
       .value("Chiral-Volume", coot::restraint_type_t::CHIRAL_VOLUME_RESTRAINT)
       .value("Trans-Peptide", coot::restraint_type_t::TRANS_PEPTIDE_RESTRAINT)
       .value("Geman-McClure", coot::restraint_type_t::GEMAN_MCCLURE_DISTANCE_RESTRAINT)// add start pos and target pos at some stage
       ;
    nb::class_<coot::simple_restraint>(m, "simple_restraint")
       .def_ro("restraint_type", &coot::simple_restraint::restraint_type)
       .def_ro("target_value",   &coot::simple_restraint::target_value)
    ;
    nb::class_<coot::geometry_distortion_info_t>(m, "geometry_distortion_info_t")
       .def_ro("distortion_score",  &coot::geometry_distortion_info_t::distortion_score)
       .def_ro("atom_indices",      &coot::geometry_distortion_info_t::atom_indices)
       .def_ro("atom_specs",        &coot::geometry_distortion_info_t::atom_specs)
       .def_ro("residue_spec",      &coot::geometry_distortion_info_t::residue_spec)
       .def_ro("restraint",         &coot::geometry_distortion_info_t::restraint)
       ;
    nb::class_<coot::geometry_distortion_info_container_t>(m, "geometry_distortion_info_container_t")
       .def_ro("chain_id",            &coot::geometry_distortion_info_container_t::chain_id)
       .def("distortion_sum",         &coot::geometry_distortion_info_container_t::distortion_sum)
       .def("size",                   &coot::geometry_distortion_info_container_t::size)
       .def("get_geometry_distortion_info", &coot::geometry_distortion_info_container_t::get_geometry_distortion_info)
       .def_ro("geometry_distortion", &coot::geometry_distortion_info_container_t::geometry_distortion)
       .def_ro("min_resno",           &coot::geometry_distortion_info_container_t::min_resno)
       .def_ro("max_resno",           &coot::geometry_distortion_info_container_t::max_resno)
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
    .def("format", &coot::residue_spec_t::format)
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
    .def("format", &coot::atom_spec_t::format)
    ;
    nb::class_<coot::plain_atom_overlap_t>(m,"plain_atom_overlap_t")
    .def(nb::init<>())
       .def_rw("ligand_atom_index", &coot::plain_atom_overlap_t::ligand_atom_index)
       .def_rw("atom_spec_1", &coot::plain_atom_overlap_t::atom_spec_1)
       .def_rw("atom_spec_2", &coot::plain_atom_overlap_t::atom_spec_2)
       .def_rw("overlap_volume", &coot::plain_atom_overlap_t::overlap_volume)
       .def_rw("r_1", &coot::plain_atom_overlap_t::r_1)
       .def_rw("r_2", &coot::plain_atom_overlap_t::r_2)
       .def_rw("is_h_bond", &coot::plain_atom_overlap_t::is_h_bond)
    ;
    nb::class_<positioned_atom_spec_t>(m,"positioned_atom_spec_t")
    .def(nb::init<>())
    .def_ro("atom_spec", &positioned_atom_spec_t::atom_spec)
    .def_ro("pos1", &positioned_atom_spec_t::pos1)
    .def_ro("pos2", &positioned_atom_spec_t::pos2)
    ;
    nb::class_<coot::atom_distance_t>(m,"atom_distance_t")
      .def(nb::init<>())
      .def_ro("atom_1", &coot::atom_distance_t::atom_1)
      .def_ro("atom_2", &coot::atom_distance_t::atom_2)
      .def_ro("distance", &coot::atom_distance_t::distance)
      ;
    nb::class_<coot::residue_range_t>(m,"residue_range_t")
      .def(nb::init<>())
      .def_rw("chain_id",     &coot::residue_range_t::chain_id)
      .def_rw("res_no_start", &coot::residue_range_t::res_no_start)
      .def_rw("res_no_end",   &coot::residue_range_t::res_no_end)
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
    nb::class_<coot::acedrg_types_for_bond_t>(m,"acedrg_types_for_bond_t")
       .def_ro("atom_id_1",   &coot::acedrg_types_for_bond_t::atom_id_1)
       .def_ro("atom_id_2",   &coot::acedrg_types_for_bond_t::atom_id_2)
       .def_ro("atom_type_1", &coot::acedrg_types_for_bond_t::atom_type_1)
       .def_ro("atom_type_2", &coot::acedrg_types_for_bond_t::atom_type_2)
       .def_ro("bond_length", &coot::acedrg_types_for_bond_t::bond_length)
       .def_ro("bond_is_between_atoms_in_the_same_ring", &coot::acedrg_types_for_bond_t::bond_is_between_atoms_in_the_same_ring)
    ;
    nb::class_<coot::acedrg_types_for_residue_t>(m,"acedrg_types_for_residue_t")
        .def_ro("bond_types", &coot::acedrg_types_for_residue_t::bond_types)
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
    nb::class_<moorhen::helix_t>(m,"helix_t")
       .def_ro("serNum", &moorhen::helix_t::serNum)
       .def_ro("helixID", &moorhen::helix_t::helixID)
       .def_ro("initResName", &moorhen::helix_t::initResName)
       .def_ro("initChainID", &moorhen::helix_t::initChainID)
       .def_ro("initSeqNum", &moorhen::helix_t::initSeqNum)
       .def_ro("initICode", &moorhen::helix_t::initICode)
       .def_ro("endResName", &moorhen::helix_t::endResName)
       .def_ro("endChainID", &moorhen::helix_t::endChainID)
       .def_ro("endSeqNum", &moorhen::helix_t::endSeqNum)
       .def_ro("endICode", &moorhen::helix_t::endICode)
       .def_ro("helixClass", &moorhen::helix_t::helixClass)
       .def_ro("comment", &moorhen::helix_t::comment)
       .def_ro("length", &moorhen::helix_t::length)
       ;

    nb::class_<moorhen::header_info_t>(m,"header_info_t")
       .def_ro("title", &moorhen::header_info_t::title)
       .def_ro("author_lines", &moorhen::header_info_t::author_lines)
       .def_ro("journal_lines", &moorhen::header_info_t::journal_lines)
       .def_ro("compound_lines", &moorhen::header_info_t::compound_lines)
       .def_ro("helix_info", &moorhen::header_info_t::helix_info)
       ;

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

    nb::class_<lsq_results_t>(m, "lsq_results_t")
      .def_ro("rotation_matrix", &lsq_results_t::rotation_matrix)
      .def_ro("translation",     &lsq_results_t::translation)
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
    nb::class_<coot::api::moved_atom_t>(m,"moved_atom_t")
    .def(nb::init<const std::string&, const std::string&, float, float, float>())
    .def(nb::init<const std::string&, const std::string&, float, float, float, int>())
    .def_ro("atom_name", &coot::api::moved_atom_t::atom_name)
    .def_ro("alt_conf", &coot::api::moved_atom_t::alt_conf)
    .def_ro("x", &coot::api::moved_atom_t::x)
    .def_ro("y", &coot::api::moved_atom_t::y)
    .def_ro("z", &coot::api::moved_atom_t::z)
    .def_ro("index", &coot::api::moved_atom_t::index)
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
    nb::class_<coot::api::moved_residue_t>(m,"moved_residue_t")
    .def(nb::init<const std::string&, int, const std::string&>())
    .def_ro("chain_id", &coot::api::moved_residue_t::chain_id)
    .def_ro("res_no", &coot::api::moved_residue_t::res_no)
    .def_ro("ins_code", &coot::api::moved_residue_t::ins_code)
    .def_ro("moved_atoms", &coot::api::moved_residue_t::moved_atoms)
    .def("add_atom",&coot::api::moved_residue_t::add_atom)
    ;
}
