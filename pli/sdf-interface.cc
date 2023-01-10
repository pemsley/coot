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
#define ENABLE_NLS // 20220126-PE Charles says this is needed to fix dcgettext() problems
                   // when including libintl.h - hmm!

// #include "graphics-info.h"
#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#include <libintl.h>
#endif // MAKE_ENHANCED_LIGAND_TOOLS

#include "compat/coot-sysdep.h"

#include <cstring>
#include <iostream> // for istream?
#include <istream>  // for istream?

// #include "c-interface-generic-objects.h" // no longer in src

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

// rename these ideally.
#include "sdf-interface.hh" // internal
#include "sdf-interface-for-export.hh"

#include "coot-utils/shape-types.hh"
#include "coot-utils/shapes.hh"

#include "sdf-internal.hh" // has use-rdkit (because an RDKit::ROMol is in the interface)
#include <GraphMol/MolChemicalFeatures/MolChemicalFeatureFactory.h>
#endif // MAKE_ENHANCED_LIGAND_TOOLS

// kekulize is an optional argument, default true.
//
bool residue_to_sdf_file(int imol, mmdb::Residue *residue_p, const char *sdf_file_name,
                         const coot::protein_geometry &geom, bool kekulize) {

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

   bool success = true; 
   // graphics_info_t g;
   if (true) {
      // mmdb::Residue *residue_p = graphics_info_t::molecules[imol].get_residue(chain_id, res_no, ins_code);
      if (residue_p) {
	 try {
	    bool includeStereo = true;
	    int confId = 0;
	    // this can throw an exception
	    RDKit::RWMol rdkm = coot::rdkit_mol_sanitized(residue_p, imol, geom);
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

bool residue_to_mdl_file_for_mogul(int imol, mmdb::Residue *residue_p,
				   const std::string &mdl_file_name,
                                   const coot::protein_geometry &geom) {

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

   bool success = false;
   // graphics_info_t g;
   if (true) {
      // mmdb::Residue *residue_p =
      // graphics_info_t::molecules[imol].get_residue(chain_id, res_no, ins_code);
      if (residue_p) {
	 try {
	    bool includeStereo = true;
	    int confId = 0;
	    // this can throw an exception
	    RDKit::RWMol rdkm = coot::rdkit_mol_sanitized(residue_p, imol, geom);

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

#include "utils/coot-utils.hh"

// rdkit chemical features.
std::vector<coot::simple_mesh_t>
chemical_features::generate_meshes(int imol, mmdb::Residue *residue_p, const coot::protein_geometry &geom) {


   std::vector<coot::simple_mesh_t> meshes;

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

   // graphics_info_t g;
   // if (g.is_valid_model_molecule(imol)) {
   if (true) {
      // mmdb::Residue *residue_p = graphics_info_t::molecules[imol].get_residue(chain_id, res_no, ins_code);
      if (! residue_p) {
	 std::cout << "Residue not found in molecule " << imol << std::endl;
      } else { 
	 try {
	    // this can throw an exception
            int iconf = 0;
	    RDKit::RWMol rdkm = coot::rdkit_mol_sanitized(residue_p, imol, geom);
	    // create a name (used to name the  generic objects object)
	    std::string name = "Chemical Features: ";
	    name += residue_p->GetChainID();
	    name += " ";
	    name += coot::util::int_to_string(residue_p->GetSeqNum());
	    name += " ";
	    name += residue_p->GetResName();
            meshes = generate_meshes(imol, rdkm, iconf, name);
	 }
	 catch (const std::runtime_error &coot_error) {
	    std::cout << coot_error.what() << std::endl;
	    std::string m = "Residue type ";
	    m += residue_p->GetResName();
	    m += " not found in dictionary.";
	 }
	 catch (const std::exception &rdkit_error) {
	    std::cout << "RDKit molecule generation problem: " << rdkit_error.what() << std::endl;
	 }
      }
   }
   return meshes;
#else
   std::cout << "Not compiled with MAKE_ENHANCED_LIGAND_TOOLS" << std::endl;
   return meshes;
#endif // MAKE_ENHANCED_LIGAND_TOOLS
}

// std::vector<coot::simple_mesh_t>
// chemical_features::generate_meshes(int imol, const RDKit::ROMol &rdkm, const std::string &name) {

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

std::vector<coot::simple_mesh_t>
chemical_features::generate_meshes(int imol, const RDKit::ROMol &rdkm, int iconf, const std::string &name) {

   std::vector<coot::simple_mesh_t> meshes;

   RDKit::MolChemicalFeatureFactory *factory = get_feature_factory();
   if (! factory) {
      std::cout << "WARNING:: no factory" << std::endl;
      return meshes;
   }

   RDKit::FeatSPtrList features = factory->getFeaturesForMol(rdkm);
   //coot::generic_display_object_t features_obj(name);

   RDKit::Conformer conf = rdkm.getConformer(iconf); // typically 0

   auto make_h_bond_arrow = [] (boost::shared_ptr<RDKit::MolChemicalFeature> &sp,
                                const RDKit::ROMol &rdkm, const RDKit::Conformer &conf,
                                const clipper::Coord_orth &centre) {
      coot::simple_mesh_t m;
      std::pair<bool, clipper::Coord_orth> normal = get_normal_info(sp.get(), rdkm, conf);
      clipper::Coord_orth p1(centre + 1.3 * normal.second);
      // meshed_generic_display_object::arrow_t arrow(centre, p1); // 20230109-PE code in src need moving to coot-utils
      return m;
   };

   auto make_aromatic_rings = [] (boost::shared_ptr<RDKit::MolChemicalFeature> &sp,
                                const RDKit::ROMol &rdkm, const RDKit::Conformer &conf,
                                const clipper::Coord_orth &centre) {
      coot::simple_mesh_t m;
      std::pair<bool, clipper::Coord_orth> normal = get_normal_info(sp.get(), rdkm, conf);
      if (normal.first) {
         clipper::Coord_orth p1(centre + 1.3 * normal.second);
         clipper::Coord_orth p2(centre - 1.3 * normal.second);
         float r_1 =   1.0f;
         if (sp.get()->getNumAtoms() == 5) r_1 = 0.8f;
         shapes::torus_t t1(p1, normal.second, r_1, 0.2f);
         shapes::torus_t t2(p2, normal.second, r_1, 0.2f);
         t1.height_scale = 0.66;
         t2.height_scale = 0.66;
         coot::simple_mesh_t m1 = coot::torus_mesh(t1);
         coot::simple_mesh_t m2 = coot::torus_mesh(t2);
         m.add_submesh(m1);
         m.add_submesh(m2);
      } else {
         std::cout << "make_aromatic_rings(): no normal " << std::endl;
      }

      glm::vec4 col(0.7, 0.7, 0.3, 1.0);
      m.change_colour(col);
      return m;
   };

   std::list<RDKit::FeatSPtr>::const_iterator it;
   for (it=features.begin(); it!=features.end(); ++it) {
      RDKit::FeatSPtr feat_ptr = *it;
      boost::shared_ptr<RDKit::MolChemicalFeature> sp = *it;
      RDGeom::Point3D pos = sp.get()->getPos();
      clipper::Coord_orth centre(pos.x, pos.y, pos.z);
      glm::vec3 centre_glm(pos.x, pos.y, pos.z);
      coot::simple_mesh_t mesh(name);

      // make a sphere for everything (but don't add it if family is "Aromatic"

      coot::simple_mesh_t sphere = coot::simple_mesh_t::make_sphere();
      glm::vec4 col(0.5, 0.5, 0.5, 1.0);

      std::string family = sp.get()->getFamily();

      if (family == "Hydrophobe")       col = glm::vec4(0.4, 0.6, 0.4, 1.0);
      if (family == "LumpedHydrophobe") col = glm::vec4(0.4, 0.4, 0.5, 1.0);
      if (family == "Aromatic")         col = glm::vec4(0.6, 0.6, 0.3, 1.0);
      if (family == "Donor")            col = glm::vec4(0.2, 0.6, 0.7, 1.0);
      if (family == "Acceptor")         col = glm::vec4(0.7, 0.2, 0.7, 1.0);
      if (family == "PosIonizable")     col = glm::vec4(0.2, 0.2, 0.7, 1.0);
      if (family == "NegIonizable")     col = glm::vec4(0.7, 0.2, 0.2, 1.0);
      // add more, like halogens
      sphere.change_colour(col);

      if (family == "LumpedHydrophobe") 
         sphere.scale(0.38f);
      else
         sphere.scale(0.5f);

      sphere.translate(centre_glm);

      if (family != "Aromatic")
         meshes.push_back(sphere);

      if (family == "Donor" || family == "Acceptor") {
         coot::simple_mesh_t arrow = make_h_bond_arrow(sp, rdkm, conf, centre);
         meshes.push_back(arrow);
      }

      if (family == "Aromatic") {
         coot::simple_mesh_t rings = make_aromatic_rings(sp, rdkm, conf, centre);
         meshes.push_back(rings);
      }
   }

   return meshes;
}

#endif




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

