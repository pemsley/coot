/* src/lbg-interface.cc
 * 
 * Copyright 2012 by The University of Oxford
 * Copyright 2013 by Medical Research Council
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#if defined (USE_PYTHON)
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

// We need this high so that dcgettext() so we don't get expected unqualified-id before 'const'
// errors when we read libintl from rdkit-interface.hh
#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#include <libintl.h>
#endif // MAKE_ENHANCED_LIGAND_TOOLS

#include "compat/coot-sysdep.h"

#include <cstring>

#define ENABLE_NLS // fix dcgettext() header problems on including
		   // libintl.h (via RDKitBase.h etc (including boost
		   // stuff).

#include <iostream> // for istream?
#include <istream> // for istream?

#include "graphics-info.h"
#include "c-interface-generic-objects.h"
#include "sdf-interface.hh"

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#include "sdf-internal.hh" // has use-rdkit (because an RDKit::ROMol is in the interface)
#include <GraphMol/MolChemicalFeatures/MolChemicalFeatureFactory.h>
#endif // MAKE_ENHANCED_LIGAND_TOOLS

#include "meshed-generic-display-object.hh"
#include "colour-holder-to-glm.hh"


bool residue_to_sdf_file(int imol, const char *chain_id, int res_no, const char *ins_code, 
			 const char *sdf_file_name, bool kekulize) {

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

   bool success = true; 
   graphics_info_t g;
   if (g.is_valid_model_molecule(imol)) {
      mmdb::Residue *residue_p =
	 graphics_info_t::molecules[imol].get_residue(chain_id, res_no, ins_code);
      if (residue_p) {
	 try {
	    bool includeStereo = true;
	    int confId = 0;
	    // this can throw an exception
	    RDKit::RWMol rdkm = coot::rdkit_mol_sanitized(residue_p, imol, *g.Geom_p());
	    // maybe this can throw an exception too.
	    RDKit::MolToMolFile(rdkm, sdf_file_name, includeStereo, confId, kekulize);
	    // success = true we presume
	 }
	 catch (const std::runtime_error &coot_error) {
	    success = false;
	    std::cout << coot_error.what() << std::endl;
	    std::string m = "Residue type ";
	    m += residue_p->GetResName();
	    m += " not found in dictionary.";
	 }
	 catch (const std::exception &rdkit_error) {
	    success = false;
	    std::cout << rdkit_error.what() << std::endl;
	 }
      } else {
	 success = false;
      }
   } else {
      success = false;
   }
   return success;
#else
   std::cout << "Not compiled with MAKE_ENHANCED_LIGAND_TOOLS" << std::endl;
   return false;
#endif // MAKE_ENHANCED_LIGAND_TOOLS
}

bool residue_to_mdl_file_for_mogul(int imol, const char *chain_id,
				   int res_no, const char *ins_code, 
				   const char *mdl_file_name) {

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

   bool success = false;
   graphics_info_t g;
   if (g.is_valid_model_molecule(imol)) {
      mmdb::Residue *residue_p =
	 graphics_info_t::molecules[imol].get_residue(chain_id, res_no, ins_code);
      if (residue_p) {
	 try {
	    bool includeStereo = true;
	    int confId = 0;
	    // this can throw an exception
	    RDKit::RWMol rdkm = coot::rdkit_mol_sanitized(residue_p, imol, *g.Geom_p());

	    coot::mogulify_mol(rdkm); // convert difficult functional groups to mogul query
	                              // format (changes reference).
	    
	    // this can throw an exception too.
	    bool kekulize = false;
	    RDKit::MolToMolFile(rdkm, mdl_file_name, includeStereo, confId, kekulize);
	    success = true;
	 }
	 catch (const std::runtime_error &coot_error) {
	    std::cout << coot_error.what() << std::endl;
	    std::string m = "Residue type ";
	    m += residue_p->GetResName();
	    m += " not found in dictionary.";
	 }
	 catch (const std::exception &rdkit_error) {
	    std::cout << rdkit_error.what() << std::endl;
	 }
      }
   }
   return success;
#else
   std::cout << "Not compiled with MAKE_ENHANCED_LIGAND_TOOLS" << std::endl;
   return false;
#endif // MAKE_ENHANCED_LIGAND_TOOLS
}


// rdkit chemical features.
bool show_feats(int imol, const char *chain_id, int res_no, const char *ins_code) {

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

   bool success = false; 
   graphics_info_t g;
   if (g.is_valid_model_molecule(imol)) {
      mmdb::Residue *residue_p =
	 graphics_info_t::molecules[imol].get_residue(chain_id, res_no, ins_code);
      if (! residue_p) {
	 std::cout << "Residue not found in molecule " << imol << std::endl;
      } else { 
	 try {
	    // this can throw an exception
	    RDKit::RWMol rdkm = coot::rdkit_mol_sanitized(residue_p, imol, *g.Geom_p());
	    // create a name (used to name the  generic objects object)
	    std::string name = "Chemical Features: ";
	    name += residue_p->GetChainID();
	    name += " ";
	    name += g.int_to_string(residue_p->GetSeqNum());
	    name += " ";
	    name += residue_p->GetResName();
	    chemical_features::show(imol, rdkm, name);
	    g.graphics_draw();
	    success = true;
	 }
	 catch (const std::runtime_error &coot_error) {
	    success = false;
	    std::cout << coot_error.what() << std::endl;
	    std::string m = "Residue type ";
	    m += residue_p->GetResName();
	    m += " not found in dictionary.";
	 }
	 catch (const std::exception &rdkit_error) {
	    success = false;
	    std::cout << "RDKit molecule generation problem: "
		      << rdkit_error.what() << std::endl;
	 }
      }
   }
   return success;
#else
   std::cout << "Not compiled with MAKE_ENHANCED_LIGAND_TOOLS" << std::endl;
   return false;
#endif // MAKE_ENHANCED_LIGAND_TOOLS
}



#include "c-interface.h" // for set_display_generic_object()
#include "c-interface-widgets.hh" // for add_generic_display_object
#include "c-interface-gtk-widgets.h" // for move_molecule_to_screen_centre_internal()

//! \brief
//! import a molecule from a smiles string
//!
//! RDKit is used to interpret the SMILES string
//!
//! no dictionary is generated
//!
//! @return a molecule number or -1 on failure
int import_rdkit_mol_from_smiles(const std::string &smiles_str, const std::string &comp_id) {

   int imol = -1;

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
   try {
      RDKit::RWMol *m = RDKit::SmilesToMol(smiles_str);
      if (m) {
	 bool explicit_only = false;
	 bool add_coords = true;
	 RDKit::MolOps::addHs(*m, explicit_only, add_coords);
	 int iconf = RDKit::DGeomHelpers::EmbedMolecule(*m);
	 if(iconf < 0) {
	    std::cout << "WARNING:: RDKit::embedding failed." << std::endl;
	 } else {
	    double vdwThresh=10.0;
	    int confId = -1;
	    bool ignoreInterfragInteractions=true; // sensible?
	    ForceFields::ForceField *ff =
	       RDKit::UFF::constructForceField(*m,
					       vdwThresh, confId,
					       ignoreInterfragInteractions);
	    ff->initialize();
	    int maxIters = 500;
	    int res=ff->minimize(maxIters);
	    delete ff;

	    mmdb::Residue *residue_p = coot::make_residue(*m, iconf, comp_id);
	    if (residue_p) {
	       mmdb::Manager *mol = coot::util::create_mmdbmanager_from_residue(residue_p);
	       if (mol) {
		  graphics_info_t g;
		  imol = g.create_molecule();
		  std::string label = "Imported ";
		  label += comp_id;
		  g.molecules[imol].install_model(imol, mol, g.Geom_p(), label, 1, false, false);
		  move_molecule_to_screen_centre_internal(imol);
	       }
	       delete residue_p;
	    }
	 }
      } else {
	 std::cout << "WARNING:: BAD SMILES " << smiles_str << std::endl;
	 std::string s = "WARNING:: Bad SMILES: " + smiles_str;
	 info_dialog(s.c_str());
      }
   }
   catch (const std::exception &e) {
      std::cout << "WARNING:: exception caught " << e.what() << std::endl;
   }
   
#endif // MAKE_ENHANCED_LIGAND_TOOLS

   return imol;
}

// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------
#ifdef MAKE_ENHANCED_LIGAND_TOOLS
// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------

// ---------------------------------- internal - no public access -----------------
// 
void chemical_features::show(int imol, const RDKit::ROMol &rdkm, std::string name) {

   graphics_info_t g;

   RDKit::MolChemicalFeatureFactory *factory = get_feature_factory();
   if (! factory) {
      std::cout << "WARNING:: no factory" << std::endl;
      return;
   }

   RDKit::FeatSPtrList features = factory->getFeaturesForMol(rdkm);

   int new_obj_idx = new_generic_object_number(name);
   meshed_generic_display_object &features_obj = g.generic_display_objects[new_obj_idx];

   RDKit::Conformer conf = rdkm.getConformer(0); // iconf?

   std::list<RDKit::FeatSPtr>::const_iterator it;
   for (it=features.begin(); it!=features.end(); it++) {
      RDKit::FeatSPtr feat_ptr = *it;
      boost::shared_ptr<RDKit::MolChemicalFeature> sp = *it;
      RDGeom::Point3D pos = sp.get()->getPos();
      clipper::Coord_orth centre(pos.x, pos.y, pos.z);
      meshed_generic_display_object::sphere_t sphere(centre, 0.4); // was 0.5 20230616-PE
      std::string family = sp.get()->getFamily();
      coot::colour_holder col;
      if (family == "Hydrophobe")
	 col = coot::colour_holder(0.4, 0.6, 0.4);
      if (family == "LumpedHydrophobe")
	 col = coot::colour_holder(0.4, 0.4, 0.5);
      if (family == "Aromatic")
	 col = coot::colour_holder(0.6, 0.6, 0.3);
      if (family == "Donor")
	 col = coot::colour_holder(0.2, 0.6, 0.7);
      if (family == "Acceptor")
	 col = coot::colour_holder(0.7, 0.2, 0.7);
      if (family == "PosIonizable")
	 col = coot::colour_holder(0.2, 0.2, 0.7);
      if (family == "NegIonizable")
	 col = coot::colour_holder(0.7, 0.2, 0.2);

      sphere.col = colour_holder_to_glm(col);

      // make the lumped sphere be smaller for aesthetic reasons (more
      // easily distinguished)
      // 
      if (family == "LumpedHydrophobe")
	 sphere.radius = 0.32; // was 0.38;

      // don't show a sphere for an aromatic - but do for everything else.
      // 
      if (family != "Aromatic")
         features_obj.add_sphere(sphere);

      if (family == "Donor") {
	 std::pair<bool, clipper::Coord_orth> normal = get_normal_info(sp.get(), rdkm, conf);
	 if (normal.first) {
	    clipper::Coord_orth p1(centre + 1.3 * normal.second);
 	    meshed_generic_display_object::arrow_t arrow(centre, p1);
	    arrow.col = col;
 	    features_obj.add_arrow(arrow);
	 }
      }

      // "pawn" shape
      if (family == "Acceptor") {
	 std::pair<bool, clipper::Coord_orth> normal = get_normal_info(sp.get(), rdkm, conf);

         if (normal.first) {
            float base_radius = 0.9f;
            float  top_radius = 0.2f;
            unsigned int n_slices = 24;
            col = coot::colour_holder(0.6, 0.1, 0.7);

            clipper::Coord_orth p1(centre + 1.3 * normal.second);
            std::pair<glm::vec3, glm::vec3> start_end(coord_orth_to_glm(centre),
                                                      coord_orth_to_glm(p1));
            features_obj.add_cone(start_end, col, base_radius, top_radius, n_slices, false, true,
                                  meshed_generic_display_object::FLAT_CAP,
                                  meshed_generic_display_object::FLAT_CAP);
         }
      }

      if (family == "Aromatic") {
	 std::pair<bool, clipper::Coord_orth> normal = get_normal_info(sp.get(), rdkm, conf);
	 if (normal.first) {
	    clipper::Coord_orth p1(centre + 1.3 * normal.second);
	    clipper::Coord_orth p2(centre - 1.3 * normal.second);
	    
            float r_1 =   1.0f;
            if (sp.get()->getNumAtoms() == 5) r_1 = 0.8f;
 	    meshed_generic_display_object::torus_t torus_1(p1, normal.second, r_1, 0.2);
 	    meshed_generic_display_object::torus_t torus_2(p2, normal.second, r_1, 0.2);
            torus_1.height_scale = 0.66;
            torus_2.height_scale = 0.66;
 	    torus_1.col = col;
 	    torus_2.col = col;
	    if (sp.get()->getNumAtoms() == 5) {
	       torus_1.n_ring_atoms = 5;
	       torus_2.n_ring_atoms = 5;
	    }
	    features_obj.add_torus(torus_1);
 	    features_obj.add_torus(torus_2);
	 }
      }
   }

   attach_generic_object_to_molecule(new_obj_idx, imol);

   // more or less boilerplate when using meshed_generic_display_object
   graphics_info_t::attach_buffers();
   Material material;
   material.do_specularity = true;
   material.specular_strength = 0.9;
   features_obj.mesh.setup(material);
   set_display_generic_object(new_obj_idx, 1);
   graphics_draw();
   
   set_display_generic_object(new_obj_idx, 1);
   graphics_draw();
}



std::pair<bool, clipper::Coord_orth>
chemical_features::get_normal_info(RDKit::MolChemicalFeature *feat,
				   const RDKit::ROMol &mol,
				   const RDKit::Conformer &conf) {
   
   bool r = false;
   clipper::Coord_orth v;
   if (feat->getFamily() == "Aromatic") {
      return get_normal_info_aromatic(feat, conf);
   } 
   if (feat->getFamily() == "Donor") {
      return get_normal_info_donor(feat, mol, conf);
   } 
   if (feat->getFamily() == "Acceptor") {
      return get_normal_info_donor(feat, mol, conf);
   } 
   return std::pair<bool, clipper::Coord_orth>(false, v); // fail
}

std::pair<bool, clipper::Coord_orth>
chemical_features::get_normal_info_aromatic(RDKit::MolChemicalFeature *feat, const RDKit::Conformer &conf) {


   bool r = false;
   clipper::Coord_orth v(0,0,0);
   if (feat->getNumAtoms() > 1) {
      RDGeom::Point3D centre_pos = feat->getPos();
      clipper::Coord_orth centre(centre_pos.x, centre_pos.y, centre_pos.z);
      const std::vector<const RDKit::Atom *> &feat_atoms = feat->getAtoms();
      int idx0 = feat_atoms[0]->getIdx();
      int idx1 = feat_atoms[1]->getIdx();
      const RDGeom::Point3D &r_pos_0 = conf.getAtomPos(idx0);
      const RDGeom::Point3D &r_pos_1 = conf.getAtomPos(idx1);
      clipper::Coord_orth pt_0(r_pos_0.x, r_pos_0.y, r_pos_0.z);
      clipper::Coord_orth pt_1(r_pos_1.x, r_pos_1.y, r_pos_1.z);
      clipper::Coord_orth v0 = pt_0 - centre;
      clipper::Coord_orth v1 = pt_1 - centre;
      clipper::Coord_orth cp(clipper::Coord_orth::cross(v0,v1).unit());
      r = true;
      v = cp;
   }
   return std::pair<bool, clipper::Coord_orth>(r, v);
 }


std::pair<bool, clipper::Coord_orth>
chemical_features::get_normal_info_donor(RDKit::MolChemicalFeature *feat,
					 const RDKit::ROMol &mol,
					 const RDKit::Conformer &conf) {

   bool r = false;
   clipper::Coord_orth v(0,0,0);
   if (feat->getNumAtoms() == 1) {
      RDGeom::Point3D centre_pos = feat->getPos();
      clipper::Coord_orth centre(centre_pos.x, centre_pos.y, centre_pos.z);
      const std::vector<const RDKit::Atom *> &feat_atoms = feat->getAtoms();
      int atom_index = feat_atoms[0]->getIdx();
      // what we do now depends on how many non-hydrogen neighbours this atom has.
      RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
      boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(feat_atoms[0]);
      std::vector<clipper::Coord_orth> neighbour_positions;
      while(nbrIdx != endNbrs){
	 const RDKit::Atom *at = mol[*nbrIdx];
	 if (at->getAtomicNum() != 1) {
	    RDGeom::Point3D r_pos = conf.getAtomPos(*nbrIdx);
	    neighbour_positions.push_back(clipper::Coord_orth(r_pos.x, r_pos.y, r_pos.z));
	 } 
	 ++nbrIdx;
      }

      if (neighbour_positions.size()) {
	 // normalize all the vectors from the neighbours to the
	 // centre atom, add them, normalize the sum - and that is the
	 // direction of the feature.
	 clipper::Coord_orth sum_vec(0,0,0);
	 for (unsigned int i=0; i<neighbour_positions.size(); i++) { 
	    sum_vec += (centre-neighbour_positions[i]);
	 }
	 if (sum_vec.lengthsq() > 0.0001) { 
	    v = clipper::Coord_orth(sum_vec.unit());
	    r = true;
	 }
      } 
   }
   return std::pair<bool, clipper::Coord_orth>(r, v);
}


// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------
#endif // MAKE_ENHANCED_LIGAND_TOOLS
// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------

