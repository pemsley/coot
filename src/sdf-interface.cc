/* src/lbg-interface.cc
 * 
 * Copyright 2012 by The University of Oxford
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

#include <cstring>
#define ENABLE_NLS // fix dcgettext() header problems on including
		   // libintl.h (via RDKitBase.h etc (including boost
		   // stuff).

#include <iostream> // for istream?
#include <istream> // for istream?

#include "graphics-info.h"
#include "sdf-interface.hh"
#include "sdf-internal.hh" // has use-rdkit (because an RDKit::ROMol is in the interface)

#include <GraphMol/MolChemicalFeatures/MolChemicalFeatureFactory.h>

bool residue_to_sdf_file(int imol, const char *chain_id, int res_no, const char *ins_code, 
			 const char *sdf_file_name, bool kekulize) {

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

   bool success = true; 
   graphics_info_t g;
   if (g.is_valid_model_molecule(imol)) {
      CResidue *residue_p =
	 graphics_info_t::molecules[imol].get_residue(chain_id, res_no, ins_code);
      if (residue_p) {
	 try {
	    bool includeStereo = true;
	    int confId = 0;
	    // this can throw an exception
	    RDKit::RWMol rdkm = coot::rdkit_mol_sanitized(residue_p, *g.Geom_p());
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
      CResidue *residue_p =
	 graphics_info_t::molecules[imol].get_residue(chain_id, res_no, ins_code);
      if (residue_p) {
	 try {
	    bool includeStereo = true;
	    int confId = 0;
	    // this can throw an exception
	    RDKit::RWMol rdkm = coot::rdkit_mol_sanitized(residue_p, *g.Geom_p());

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
      CResidue *residue_p =
	 graphics_info_t::molecules[imol].get_residue(chain_id, res_no, ins_code);
      if (! residue_p) {
	 std::cout << "Residue not found in molecule " << imol << std::endl;
      } else { 
	 try {
	    // this can throw an exception
	    RDKit::RWMol rdkm = coot::rdkit_mol_sanitized(residue_p, *g.Geom_p());
	    // create a name (used to name the  generic objects object)
	    std::string name = "Chemical Features: ";
	    name += residue_p->GetChainID();
	    name += " ";
	    name += g.int_to_string(residue_p->GetSeqNum());
	    name += " ";
	    name += residue_p->GetResName();
	    chemical_features::show(rdkm, name);
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



// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------
#ifdef MAKE_ENHANCED_LIGAND_TOOLS
// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------


// internal - no public access
// 
void chemical_features::show(const RDKit::ROMol &rdkm, std::string name) {

   graphics_info_t g;

   RDKit::MolChemicalFeatureFactory *factory = get_feature_factory();
   if (! factory) {
      std::cout << "WARNING:: no factory" << std::endl;
      return;
   }

   RDKit::FeatSPtrList features = factory->getFeaturesForMol(rdkm);
   coot::generic_display_object_t features_obj(name);
   RDKit::Conformer conf = rdkm.getConformer(0); // iconf?

   std::list<RDKit::FeatSPtr>::const_iterator it;
   for (it=features.begin(); it!=features.end(); it++) {
      RDKit::FeatSPtr feat_ptr = *it;
      boost::shared_ptr<RDKit::MolChemicalFeature> sp = *it;
      RDGeom::Point3D pos = sp.get()->getPos();
      // std::cout << "feature pos: " << pos << std::endl;
      clipper::Coord_orth centre(pos.x, pos.y, pos.z);
      coot::generic_display_object_t::sphere_t sphere(centre, 0.5);
      std::string family = sp.get()->getFamily();
      // std::cout << "family: " << sp.get()->getFamily() << std::endl;
      // std::cout << "   type: " << sp.get()->getType() << std::endl;
      coot::colour_t col;
      if (family == "Hydrophobe")
	 col = coot::colour_t(0.4, 0.6, 0.4);
      if (family == "LumpedHydrophobe")
	 col = coot::colour_t(0.4, 0.4, 0.5);
      if (family == "Aromatic")
	 col = coot::colour_t(0.6, 0.6, 0.3);
      if (family == "Donor")
	 col = coot::colour_t(0.2, 0.6, 0.7);
      if (family == "Acceptor")
	 col = coot::colour_t(0.7, 0.2, 0.7);
      if (family == "PosIonizable")
	 col = coot::colour_t(0.2, 0.2, 0.7);
      if (family == "NegIonizable")
	 col = coot::colour_t(0.7, 0.2, 0.2);

      sphere.col = col;

      // make the lumped sphere be smaller for aesthetic reasons (more
      // easily distinguished)
      // 
      if (family == "LumpedHydrophobe")
	 sphere.radius = 0.38;

      // don't show a sphere for an aromatic - but do for everything else.
      // 
      if (family != "Aromatic")
	 features_obj.spheres.push_back(sphere);

      if (family == "Donor" || family == "Acceptor") {
	 std::pair<bool, clipper::Coord_orth> normal = get_normal_info(sp.get(), rdkm, conf);
	 if (normal.first) {
	    clipper::Coord_orth p1(centre + 1.3 * normal.second);
 	    coot::generic_display_object_t::arrow_t arrow(centre, p1);
	    arrow.col = col;
 	    features_obj.arrows.push_back(arrow);
	 }
      }
      
      if (family == "Aromatic") {
	 std::pair<bool, clipper::Coord_orth> normal = get_normal_info(sp.get(), rdkm, conf);
	 if (normal.first) {
	    clipper::Coord_orth p1(centre + 1.3 * normal.second);
	    clipper::Coord_orth p2(centre - 1.3 * normal.second);
	    
 	    coot::generic_display_object_t::torus_t torus_1(centre, p1, 0.18, 1.1);
 	    coot::generic_display_object_t::torus_t torus_2(centre, p2, 0.18, 1.1);
 	    torus_1.col = col;
 	    torus_2.col = col;
	    if (sp.get()->getNumAtoms() == 5) {
	       torus_1.n_ring_atoms = 5;
	       torus_2.n_ring_atoms = 5;
	    }
	    features_obj.tori.push_back(torus_1);
 	    features_obj.tori.push_back(torus_2);
	 }
      }
   }
   g.generic_objects_p->push_back(features_obj);
   // this is hacky.
   int iobj = g.generic_objects_p->size() -1;
   g.generic_objects_p->back().is_displayed_flag = true;
}

RDKit::MolChemicalFeatureFactory *
chemical_features::get_feature_factory() {

   RDKit::MolChemicalFeatureFactory *factory = NULL;
   
   std::string features_file_name("BaseFeatures.fdef");
   std::string data_dir = "Data";
   std::string dir =
      coot::util::append_dir_dir(coot::rdkit_package_data_dir(), data_dir);
   std::string full_name =
      coot::util::append_dir_file(dir, features_file_name);

   // override that with an environment variable
   const char *e = getenv("COOT_CHEMICAL_FEATURES_DEF");
   if (e) full_name = e;
   
   if (coot::file_exists(full_name)) {
      std::ifstream inStream(full_name.c_str());
      std::istream &instrm = static_cast<std::istream &>(inStream);
      factory = RDKit::buildFeatureFactory(instrm);

   } else { 
      std::cout << "WARNING:: " << full_name << " does not exist. "
		<< "Stoping now." << std::endl;
   }
   return factory;
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
	 const RDKit::ATOM_SPTR at = mol[*nbrIdx];
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

