/*
 * api/coot-molecule-moltris.cc
 * 
 * Copyright 2020 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */


#include "coot-molecule.hh"

// #include <MoleculesToTriangles/CXXClasses/RendererGL.h>
#include <MoleculesToTriangles/CXXClasses/Light.h>
#include <MoleculesToTriangles/CXXClasses/Camera.h>
#include <MoleculesToTriangles/CXXClasses/SceneSetup.h>
#include <MoleculesToTriangles/CXXClasses/ColorScheme.h>
#include <MoleculesToTriangles/CXXClasses/MyMolecule.h>
#include <MoleculesToTriangles/CXXClasses/RepresentationInstance.h>
#include <MoleculesToTriangles/CXXClasses/MolecularRepresentationInstance.h>
#include <MoleculesToTriangles/CXXClasses/VertexColorNormalPrimitive.h>
#include <MoleculesToTriangles/CXXClasses/BallsPrimitive.h>
#include <MoleculesToTriangles/CXXClasses/SurfacePrimitive.h>

#include "MoleculesToTriangles/CXXClasses/tubes.hh"

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>

#include <tuple>
#include <map>
#include <string>


//! Add a colour rule: eg. ("//A", "red")
void
coot::molecule_t::add_colour_rule(const std::string &selection, const std::string &colour_name) {

   colour_rules.push_back(std::make_pair(selection, colour_name));
}


//! delete all the colour rules
void
coot::molecule_t::delete_colour_rules() {
   colour_rules.clear();
}

void
coot::molecule_t::print_colour_rules() const {

   std::cout << "=============================" << std::endl;
   std::cout << " colour rules for molecule " << imol_no << std::endl;
   std::cout << "=============================" << std::endl;
   for (unsigned int i=0; i<colour_rules.size(); i++) {
      const auto &cr = colour_rules[i];
      std::cout << "   " << i << " " << cr.first << " " << cr.second << std::endl;
   }
   std::cout << "=============================" << std::endl;

}

void
coot::molecule_t::fill_default_colour_rules() {

   int imod = 1;
   if (! atom_sel.mol) return;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   if (model_p) {
      std::vector<std::string> chain_ids;
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int n_res = chain_p->GetNumberOfResidues();
         if (n_res > 0) {
            std::string chain_id(chain_p->GetChainID());
            chain_ids.push_back(chain_id);
         }
      }

      std::vector<std::string> colour_names = {
         "Salmon", "SandyBrown", "LightSeaGreen", "Goldenrod", "tomato",
         "limegreen", "royalblue", "gold", "aquamarine", "maroon",
         "lightcoral", "deeppink", "brown", "Burlywood", "steelblue"};

      colour_rules.clear();
      unsigned int name_index = 0;
      for (unsigned int i=0; i<chain_ids.size(); i++) {
         const auto &chain_id = chain_ids[i];
         std::string cid = std::string("//") + chain_id;
         const std::string &cn = colour_names[name_index];
         auto p = std::make_pair(cid, cn);
         colour_rules.push_back(p);

         // next round
         name_index++;
         if (name_index == colour_names.size()) name_index = 0;
      }


      // 20230204-PE try something else
      std::sort(chain_ids.begin(), chain_ids.end());
      colour_rules.clear();
      bool against_a_dark_background = true;
      for (unsigned int i=0; i<chain_ids.size(); i++) {
         const auto &chain_id = chain_ids[i];
         int col_index = i + 50; // imol_no is added in get_bond_colour_by_mol_no()
         colour_t col = get_bond_colour_by_mol_no(col_index, against_a_dark_background);
         colour_holder ch = col.to_colour_holder(); // around the houses
         std::string hex = ch.hex();
         std::string cid = std::string("//") + chain_id;
         if (false)
            std::cout << "debug:: colour_rules(); imol " << imol_no << " i " << i << " col_index " << col_index
                      << " " << hex << std::endl;
         auto p = std::make_pair(cid, hex);
         colour_rules.push_back(p);
      }
   }
}

//! get the colour rules. Preferentially return the user-defined colour rules.
//! @return If there are no user-defined colour rules, then return the stand-in rules
std::vector<std::pair<std::string, std::string> >
coot::molecule_t::get_colour_rules() const {

      return colour_rules;
}

//! Update float parameter for MoleculesToTriangles molecular mesh
void
coot::molecule_t::M2T_updateFloatParameter(const std::string &param_name, float value) {

   bool done = false;
   for (unsigned int i=0; i<M2T_float_params.size(); i++) {
      if (param_name == M2T_float_params[i].first) {
         M2T_float_params[i].second = value;
         done = true;
         break;
      }
   }
   if (! done) {
      M2T_float_params.push_back(std::make_pair(param_name, value));
   }
}

//! Update int parameter for MoleculesToTriangles molecular mesh
void
coot::molecule_t::M2T_updateIntParameter(const std::string &param_name, int value) {

   bool done = false;
   for (unsigned int i=0; i<M2T_int_params.size(); i++) {
      if (param_name == M2T_int_params[i].first) {
         M2T_int_params[i].second = value;
         done = true;
         break;
      }
   }
   if (! done) {
      M2T_int_params.push_back(std::make_pair(param_name, value));
   }
}

void
coot::molecule_t::print_M2T_FloatParameters() const {

   for (unsigned int i=0; i<M2T_float_params.size(); i++) {
      std::cout << "   " << i << " " << M2T_float_params[i].first << " " << M2T_float_params[i].second << std::endl;
   }

}


void
coot::molecule_t::print_M2T_IntParameters() const {

   for (unsigned int i=0; i<M2T_int_params.size(); i++) {
      std::cout << "   " << i << " " << M2T_int_params[i].first << " " << M2T_int_params[i].second << std::endl;
   }

}


#include "MoleculesToTriangles/CXXSurface/CXXChargeTable.h"
#include "MoleculesToTriangles/CXXSurface/CXXUtils.h"
#include "MoleculesToTriangles/CXXSurface/CXXSurface.h"
#include "MoleculesToTriangles/CXXSurface/CXXCreator.h"

#include "coot-utils/json.hpp"
using json = nlohmann::json;

//! \brief set the residue properties
//!
//! a list of propperty maps such as `{"chain-id": "A", "res-no": 34, "ins-code": "", "worm-radius": 1.2}`
//!
//! @param json_string is the properties in JSON format
//! @return true
bool
coot::molecule_t::set_residue_properties(const std::string &json_string) {

   bool status = true;

   json j = json::parse(json_string);
   for (json::iterator it=j.begin(); it!=j.end(); ++it) {
      json j_residue_properties = *it;
      json::const_iterator j_chain_id          = j_residue_properties.find("chain-id");
      json::const_iterator j_res_no            = j_residue_properties.find("res-no");
      json::const_iterator j_ins_code          = j_residue_properties.find("ins-code");
      json::const_iterator j_worm_radius       = j_residue_properties.find("worm-radius");
      json::const_iterator j_worm_radius_aniso = j_residue_properties.find("worm-radius-aniso");
      float radius = -1.1f;
      if (j_chain_id != j_residue_properties.end()) {
         if (j_res_no != j_residue_properties.end()) {
            if (j_ins_code != j_residue_properties.end()) {
               res_prop_t res_prop;

               std::string chain_id = j_chain_id.value();
               int res_no = j_res_no.value();
               std::string ins_code = j_ins_code.value();
               coot::residue_spec_t res_spec(chain_id, res_no, ins_code);

               if (j_worm_radius_aniso != j_residue_properties.end()) {
                  if (j_worm_radius_aniso->is_array()) {
                     if (j_worm_radius_aniso->size() == 3) {
                        res_prop.value_x = (*j_worm_radius_aniso)[0];
                        res_prop.value_y = (*j_worm_radius_aniso)[1];
                        res_prop.value_z = (*j_worm_radius_aniso)[2];
                     }
                  }
               }

               if (j_worm_radius != j_residue_properties.end()) {
                  float radius = j_worm_radius.value();
                  if (radius > 0.0) {
                     // std::cout << "debug:: radius for " << res_spec << " " << radius << std::endl;
                     res_prop.value = radius;
                  }
               }

               if (res_prop.value > 0 || res_prop.filled_aniso())
                  residue_properties_map[res_spec] = res_prop;

            }
         }
      }
   }
   return status;
}

void
coot::molecule_t::clear_residue_properties() {

   residue_properties_map.clear();
}


coot::simple_mesh_t
coot::molecule_t::get_molecular_representation_mesh(const std::string &atom_selection_str,
                                                    const std::string &colour_scheme,
                                                    const std::string &style,
                                                    int secondaryStructureUsageFlag) const {

   bool debug = true;

   auto get_max_resno_for_polymer = [] (mmdb::Chain *chain_p) {
      int res_no_max = -1;
      int nres = chain_p->GetNumberOfResidues();
      for (int ires=nres-1; ires>=0; ires--) {
         mmdb::Residue *residue_p = chain_p->GetResidue(ires);
         if (residue_p) {
            int seq_num = residue_p->GetSeqNum();
            if (seq_num > res_no_max) {
               if (residue_p->isAminoacid() || residue_p->isNucleotide()) {
                  res_no_max = seq_num;
               }
            }
         }
      }
      return res_no_max;
   };

   // A residue as a CID of the form /1/A/23(ALA), or /1/A/23(ALA).B with an insertion code,
   // which is what a caller will want to hand back to the rest of the API.
   //
   // The insertion code goes after the residue name, following mmdb's own
   // /mdl/chn/seq(res).ins grammar. Putting it after the number instead parses, but wrongly: a
   // reader taking the text after the dot then gets "B(ALA)" rather than "B".
   auto residue_cid = [] (mmdb::Residue *residue_p) {
      std::string cid = "/";
      cid += std::to_string(residue_p->GetModelNum());
      cid += "/";
      cid += residue_p->GetChainID();
      cid += "/";
      cid += std::to_string(residue_p->GetSeqNum());
      cid += "(";
      cid += residue_p->GetResName();
      cid += ")";
      std::string ins_code = residue_p->GetInsCode();
      if (! ins_code.empty()) {
         cid += ".";
         cid += ins_code;
      }
      return cid;
   };

   auto molecular_representation_instance_to_mesh = [&residue_cid] (std::shared_ptr<MolecularRepresentationInstance> molrepinst,
                                                        const std::vector<std::pair<std::string, float> > &M2T_float_params,
                                                        const std::vector<std::pair<std::string, int> > &M2T_int_params) {
      coot::simple_mesh_t mesh;

      std::shared_ptr<Representation> r = molrepinst->getRepresentation();

      if (! M2T_float_params.empty())
         for (const auto &par : M2T_float_params)
            r->updateFloatParameter(par.first, par.second);
      if (! M2T_int_params.empty())
         for (const auto &par : M2T_int_params)
            r->updateIntParameter(par.first, par.second);

      // testing
      // r->updateFloatParameter("ribbonStyleCoilThickness", 1.8);
      // r->updateFloatParameter("ribbonStyleHelixWidth", 2.6);

      r->redraw(); // 20221207-PE this was missing!
      std::vector<std::shared_ptr<DisplayPrimitive> > vdp = r->getDisplayPrimitives();

      auto displayPrimitiveIter = vdp.begin();
      for (displayPrimitiveIter=vdp.begin(); displayPrimitiveIter != vdp.end(); displayPrimitiveIter++) {
         DisplayPrimitive &displayPrimitive = **displayPrimitiveIter;
         if (displayPrimitive.type() == DisplayPrimitive::PrimitiveType::SurfacePrimitive    ||
             displayPrimitive.type() == DisplayPrimitive::PrimitiveType::BoxSectionPrimitive ||
             displayPrimitive.type() == DisplayPrimitive::PrimitiveType::BallsPrimitive      ||
             displayPrimitive.type() == DisplayPrimitive::PrimitiveType::CylinderPrimitive ){
            displayPrimitive.generateArrays();

            coot::simple_mesh_t submesh;
            VertexColorNormalPrimitive &surface = dynamic_cast<VertexColorNormalPrimitive &>(displayPrimitive);
            submesh.vertices.resize(surface.nVertices());

            auto vcnArray = surface.getVertexColorNormalArray();
            for (unsigned int iVertex=0; iVertex < surface.nVertices(); iVertex++){
               coot::api::vnc_vertex &gv = submesh.vertices[iVertex];
               VertexColorNormalPrimitive::VertexColorNormal &vcn = vcnArray[iVertex];
               for (int i=0; i<3; i++) {
                  gv.pos[i]    = vcn.vertex[i];
                  gv.normal[i] = vcn.normal[i];
                  gv.color[i]  = 0.0037f * vcn.color[i];
               }
               gv.color[3] = 0.00392f * vcn.color[3];
               if ((gv.color[0] + gv.color[1] + gv.color[2]) > 10.0) {
                  gv.color *= 0.00392f; // /255
                  if (gv.color[3] > 0.99)
                     gv.color[3] = 1.0;
               }
            }

            // Which residue each vertex belongs to.
            //
            // Several of these primitives already know. A surface is built from one sphere
            // patch per atom and one torus per pair, and records the generating atom against
            // every vertex; a ribbon is swept along a spline and records the alpha carbon of
            // the residue whose half-turn it is drawing. All of it survived into atomArray on
            // the base class and was dropped here, because simple_mesh_t had nowhere to put
            // it. Now it does, so a caller can tell which residue is under a triangle it hit.
            //
            // No list of primitive types: the question is whether this one recorded anything,
            // and a null array answers it. That works only because all three that allocate
            // the array now value-initialise it, so an unwritten slot is null rather than
            // heap litter - the primitives that record nothing never allocate it at all.
            {
               const mmdb::Atom **atomArray = surface.getAtomArray();
               // The other atom of a saddle and the split between the two. Only a surface has
               // these, and only a surface needs them: a ribbon vertex lies within one
               // residue's half-turn of spline, so it belongs to that residue outright.
               const mmdb::Atom **atom2Array = 0;
               const float *atomWeightArray = 0;
               if (displayPrimitive.type() == DisplayPrimitive::PrimitiveType::SurfacePrimitive) {
                  SurfacePrimitive &surfacePrimitive = dynamic_cast<SurfacePrimitive &>(displayPrimitive);
                  atom2Array = surfacePrimitive.getAtom2Array();
                  atomWeightArray = surfacePrimitive.getAtomWeightArray();
               }
               if (atomArray) {
                  submesh.vertex_owner.resize(surface.nVertices(), coot::simple_mesh_t::no_owner);
                  submesh.vertex_owner_other.resize(surface.nVertices(), coot::simple_mesh_t::no_owner);
                  submesh.vertex_owner_weight.resize(surface.nVertices(), 1.0f);
                  std::map<mmdb::Residue *, unsigned int> residue_index;

                  // Name a residue, adding it to this submesh's table the first time it is seen.
                  auto index_of_residue = [&residue_index, &submesh, &residue_cid] (mmdb::Atom *at) {
                     if (! at) return coot::simple_mesh_t::no_owner;
                     mmdb::Residue *residue_p = at->GetResidue();
                     if (! residue_p) return coot::simple_mesh_t::no_owner;
                     std::map<mmdb::Residue *, unsigned int>::const_iterator it = residue_index.find(residue_p);
                     if (it != residue_index.end()) return it->second;
                     unsigned int idx = submesh.owners.size();
                     submesh.owners.push_back(residue_cid(residue_p));
                     residue_index[residue_p] = idx;
                     return idx;
                  };

                  for (unsigned int iVertex=0; iVertex < surface.nVertices(); iVertex++) {
                     // const_cast because mmdb's accessors are not const, not because anything
                     // here modifies the atom.
                     mmdb::Atom *at = const_cast<mmdb::Atom *>(atomArray[iVertex]);
                     unsigned int owner = index_of_residue(at);
                     if (owner == coot::simple_mesh_t::no_owner) continue;

                     mmdb::Atom *at2 = atom2Array ? const_cast<mmdb::Atom *>(atom2Array[iVertex]) : 0;
                     unsigned int other = index_of_residue(at2);
                     float weight = atomWeightArray ? atomWeightArray[iVertex] : 1.0f;

                     // A saddle between two atoms of the SAME residue is not a boundary at all,
                     // and most of them are: collapse it, so the common case stays a single
                     // full-weight owner and the data only grows where a boundary really runs.
                     if (other == owner) other = coot::simple_mesh_t::no_owner;
                     if (other == coot::simple_mesh_t::no_owner) weight = 1.0f;

                     // The heavier of the two goes first, so a reader that wants one answer -
                     // which residue was clicked - can take vertex_owner without comparing.
                     if (other != coot::simple_mesh_t::no_owner && weight < 0.5f) {
                        std::swap(owner, other);
                        weight = 1.0f - weight;
                     }

                     submesh.vertex_owner[iVertex] = owner;
                     submesh.vertex_owner_other[iVertex] = other;
                     submesh.vertex_owner_weight[iVertex] = weight;
                  }
               }
            }

            auto indexArray = surface.getIndexArray();
            submesh.triangles.resize(surface.nTriangles());
            for (unsigned int iTriangle=0; iTriangle<surface.nTriangles(); iTriangle++){
               g_triangle &gt = submesh.triangles[iTriangle];
               for (int i=0; i<3; i++)
                  gt[i] = indexArray[3*iTriangle+i];
            }
            if (false)
               std::cout << "in molecular_representation_instance_to_mesh(): Here B adding "
                         << submesh.vertices.size() << " vertices and "
                         << submesh.triangles.size() << " triangles" << std::endl;
            mesh.add_submesh(submesh);
         }
      }

      return mesh;
   };

   class chain_info_t {
   public:
      mmdb::Chain *chain_p;
      int resno_min;
      int resno_max;
      chain_info_t(mmdb::Chain *c, int min, int max) : chain_p(c), resno_min(min), resno_max(max) {}
   };

   auto get_chains_in_selection = [] (std::shared_ptr<MyMolecule> my_mol,
                                const std::string &atom_selection_str) {
      std::vector<chain_info_t> ci;
      MyMolecule *mm = my_mol.get();
      mmdb::Manager *mol = my_mol.get()->mmdb;

      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            int resno_max = -999999;
            int resno_min = 999999;
            if (n_res > 0) {
               for (int ires=0; ires<n_res; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     int res_no = residue_p->GetSeqNum();
                     if (! residue_p->isSolvent()) {
                        if (res_no > resno_max) resno_max = res_no;
                        if (res_no < resno_min) resno_min = res_no;
                     }
                  }
               }
            }
            ci.push_back(chain_info_t(chain_p, resno_min, resno_max));
         }
      }

      return ci;
   };

   auto get_polymer_chains_in_selection = [] (std::shared_ptr<MyMolecule> my_mol,
                                              const std::string &atom_selection_str) {
      std::vector<chain_info_t> ci;
      MyMolecule *mm = my_mol.get();
      mmdb::Manager *mol = my_mol.get()->mmdb;

      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            int resno_max = -999999;
            int resno_min = 999999;
            if (n_res > 0) {
               for (int ires=0; ires<n_res; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     int res_no = residue_p->GetSeqNum();
                     if (residue_p->isAminoacid() || residue_p->isNucleotide()) {
                        if (res_no > resno_max) resno_max = res_no;
                        if (res_no < resno_min) resno_min = res_no;
                     }
                  }
               }
            }
            ci.push_back(chain_info_t(chain_p, resno_min, resno_max));
         }
      }

      return ci;
   };

   auto ramp_chains = [get_polymer_chains_in_selection,
                       molecular_representation_instance_to_mesh]
      (std::shared_ptr<MyMolecule> my_mol,
       const std::string &atom_selection_str,
       const std::string &style,
       const std::vector<std::pair<std::string, float> > &M2T_float_params,
       const std::vector<std::pair<std::string, int> > &M2T_int_params) {

      coot::simple_mesh_t mesh;

      std::vector<chain_info_t> ci = get_polymer_chains_in_selection(my_mol, atom_selection_str);

      for (const auto &ch : ci) {

         if (ch.resno_max < ch.resno_min) continue;
         std::string chain_sel = "//" + std::string(ch.chain_p->GetChainID());

         auto ramp_cs = std::shared_ptr<ColorScheme>(new ColorScheme());
         auto apcrr_p = std::make_shared<AtomPropertyRampColorRule>();
         int n_ramp_points = ch.resno_max - ch.resno_min;
         apcrr_p->setNumberOfRampPoints(n_ramp_points);
         apcrr_p->setStartValue(ch.resno_min);
         apcrr_p->setEndValue(ch.resno_max);
         apcrr_p->setCompoundSelection(
            std::shared_ptr<CompoundSelection>(new CompoundSelection(chain_sel + "/*.*/*:*")));
         ramp_cs->addRule(apcrr_p);

         std::shared_ptr<MolecularRepresentationInstance> molrepinst =
            MolecularRepresentationInstance::create(my_mol, ramp_cs, chain_sel, style);
         coot::simple_mesh_t submesh = molecular_representation_instance_to_mesh(molrepinst, M2T_float_params, M2T_int_params);
         mesh.add_submesh(submesh);
      }
      return mesh;
   };

   auto glm_to_clipper = [] (const glm::vec3 &gp) {
      return clipper::Coord_orth(gp[0], gp[1], gp[2]);
   };

   auto set_vertex_colour = [] (api::vnc_vertex &vertex, float potential) {
       float pot_f = min(1.0,fabs(potential)*2.);
       if (potential < 0.0){
           vertex.color = glm::vec4(1.0, 1.0-pot_f, 1.0-pot_f, 1.0);
       } else {
           vertex.color = glm::vec4(1.0-pot_f, 1.0-pot_f, 1.0, 1.0);
       }
   };

   // ----------------- main line -------------------

   coot::simple_mesh_t mesh;

   if (false)
      std::cout << "get_molecular_representation_mesh() atom_selection: " << atom_selection_str
                << " colour_scheme: " << colour_scheme << " style: " << style << std::endl;

   if (style == "Tubes") { //  bendy helices

      mmdb::Manager *mol = atom_sel.mol;
      std::string atom_selection = "//";
      std::string colour_scheme = "Helix";
      float radius_for_helices = 3.2;
      unsigned int n_slices_for_helices = 16;
      float radius_for_coil = 0.8;
      float Cn_for_coil = 2;
      int accuracy_for_coil = 12;
      unsigned int n_slices_for_coil = 12;
      // make_tubes_representation() removed until it can be compiled without using coot include/libs
      // mesh = make_tubes_representation(mol, atom_selection, colour_scheme, radius_for_coil, Cn_for_coil,
      //                                  accuracy_for_coil, n_slices_for_coil, secondaryStructureUsageFlag);

   } else {

      {

         try {

            auto my_mol = std::make_shared<MyMolecule>(atom_sel.mol, secondaryStructureUsageFlag);

            // Convert residue_properies_map to M2T format and pass to MyMolecule
            if (!residue_properties_map.empty()) {
               std::map<std::tuple<std::string, int, std::string>, std::tuple<float, float, float>> radii_map;
               for (const auto &[spec, prop] : residue_properties_map) {
                  if (prop.filled_aniso()) {
                     radii_map[std::make_tuple(spec.chain_id, spec.res_no, spec.ins_code)] =
                         std::make_tuple(prop.value_x, prop.value_y, prop.value_z);
                  } else if (prop.value > 0.0f) {
                     // Isotropic fallback: use same value for all axes
                     radii_map[std::make_tuple(spec.chain_id, spec.res_no, spec.ins_code)] =
                         std::make_tuple(prop.value, prop.value, prop.value);
                  }
               }
               if (!radii_map.empty()) {
                  my_mol->setResidueRadii(radii_map);
               }
            }

            // auto chain_cs = ColorScheme::colorChainsScheme();
            auto chain_cs = ColorScheme::colorChainsSchemeWithColourRules(colour_rules);
            if (! colour_rules.empty())
               chain_cs = ColorScheme::colorChainsSchemeWithColourRules(colour_rules);
            auto ele_cs   = ColorScheme::colorByElementScheme();
            auto ss_cs    = ColorScheme::colorBySecondaryScheme();
            auto bf_cs    = ColorScheme::colorBFactorScheme();
            auto this_cs  = chain_cs; // default
            if (colour_scheme == "Chains")    this_cs = chain_cs;
            if (colour_scheme == "Element")   this_cs = ele_cs;
            if (colour_scheme == "BFactor")   this_cs = bf_cs;
            if (colour_scheme == "Secondary") this_cs = ss_cs;
            if (colour_scheme == "RampChains" || colour_scheme == "colorRampChainsScheme") {
               mesh = ramp_chains(my_mol, atom_selection_str, style, M2T_float_params, M2T_int_params);
            } else {

               if (colour_scheme == "ByOwnPotential") {

                  std::shared_ptr<MolecularRepresentationInstance> molrepinst =
                     MolecularRepresentationInstance::create(my_mol, this_cs, atom_selection_str, style);
                  mesh = molecular_representation_instance_to_mesh(molrepinst, M2T_float_params, M2T_int_params);

                  //Instantiate an electrostatics map and cause it to calculate itself
                  CXXChargeTable theChargeTable;
                  CXXUtils::assignCharge(atom_sel.mol, atom_sel.SelectionHandle, &theChargeTable);
                  CXXCreator *theCreator = new CXXCreator(atom_sel.mol, atom_sel.SelectionHandle);
                  theCreator->calculate();
                  clipper::Cell cell;
                  clipper::NXmap<double> theClipperNXMap;
                  theClipperNXMap = theCreator->coerceToClipperMap(cell);

                  for (unsigned int i=0; i<mesh.vertices.size(); i++) {
                     clipper::Coord_orth orthogonals = glm_to_clipper(mesh.vertices[i].pos);
                     const clipper::Coord_map mapUnits(theClipperNXMap.coord_map(orthogonals));
                     float potential = theClipperNXMap.interp<clipper::Interp_cubic>( mapUnits );
                     // subSurfaceIter->setScalar(potentialHandle, i, potential);
                     // std::cout << "potential: " << potential << std::endl;  // 20240226-PE Quiet! (now that this test is run)
                     set_vertex_colour(mesh.vertices[i], potential); // change ref
                  }

                  delete theCreator;

               } else {

                  std::shared_ptr<MolecularRepresentationInstance> molrepinst =
                     MolecularRepresentationInstance::create(my_mol, this_cs, atom_selection_str, style);
                  mesh = molecular_representation_instance_to_mesh(molrepinst, M2T_float_params, M2T_int_params);

                  if (false) {
                     for (unsigned int i=0; i<mesh.vertices.size(); i++) {
                        const auto &vertex = mesh.vertices[i];
                        std::cout << i << " " << glm::to_string(vertex.pos) << " " << glm::to_string(vertex.color) << std::endl;
                     }
                  }
               }
            }
            mesh.fill_colour_map(); // for blendering
         }
         // A mesh that could not be built says so, rather than being returned empty.
         //
         // simple_mesh_t::status is documented as exactly this flag - "1 is good, 0 is bad
         // (0 is set when we get a bad_alloc)" - but nothing here was setting it, so a
         // representation that threw came back indistinguishable from one whose selection
         // legitimately matched nothing. A caller then drew no geometry and had no way to
         // tell whether that was the answer or a failure.
         //
         // Partly built is also failed: whatever was added before the throw is not the mesh
         // that was asked for, so it is cleared rather than half drawn.
         catch (const std::out_of_range &oor) {
            std::cout << "ERROR:: out of range in get_molecular_representation_mesh() " << oor.what() << std::endl;
            mesh.clear();
            mesh.status = 0;
         }
         catch (const std::runtime_error &rte) {
            std::cout << "ERROR:: runtime error in get_molecular_representation_mesh() " << rte.what() << std::endl;
            mesh.clear();
            mesh.status = 0;
         }
         catch (...) {
            std::cout << "ERROR:: unknown exception in get_molecular_representation_mesh()! " << std::endl;
            mesh.clear();
            mesh.status = 0;
         }
      }
   }
   return mesh;
}
