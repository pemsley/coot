/* lidia-core/chemical-feature-clusters.cc
 * 
 * Copyright 2012 by the University of Oxford
 * Copyright 2016 by Medical Research Council
 * Author: Paul Emsley
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

// header-here

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#include <iostream>
#include <fstream>
#include <cstring>

#include "utils/coot-utils.hh"
// #include "coot-utils/coot-coord-utils.hh" out of order
#include "rdkit-interface.hh"
#include "chemical-feature-clusters.hh"
#include <GraphMol/MolChemicalFeatures/MolChemicalFeature.h>

#include "get-residue.hh"


// do_alignment is an optional arg - default false
// 
coot::chem_feat_clust::chem_feat_clust(const std::vector<residue_spec_t> &protein_residues,
				       const std::vector<chem_feat_solvated_ligand_spec> &ligands_in,
				       double water_dist_cutoff_in,
				       const protein_geometry *geometry_p_in,
				       bool do_alignment) {

   water_dist_cutoff = water_dist_cutoff_in;
   geometry_p = geometry_p_in;
   setup_success = fill_ligands(ligands_in);
   if (setup_success) {
      setup_success = check_dictionaries();

      if (setup_success) { // setup_success
	 if (do_alignment)
	    align();
	 
	 fill_waters();
      }
   }
}



void
coot::chem_feat_solvated_ligand::init_residue() {

   residue = lidia_utils::get_residue(ligand_spec, mol);

   if (!residue) {
      std::cout << "WARNING:: null residue from spec " << ligand_spec << std::endl;
   }
}


std::vector<coot::simple_chemical_feature_attributes>
coot::chem_feat_clust::get_chemical_features(int imol,
					     residue_spec_t lig_spec,
					     mmdb::Manager *mol) {

   std::vector<coot::simple_chemical_feature_attributes> v;
   bool success = false;

   if (! setup_success)
      return v;   
	 
   mmdb::Residue *residue_p = lidia_utils::get_residue(lig_spec, mol);

   if (! residue_p) {
      std::cout << "WARNING:: failed to get ligand for molecule " << imol
		<< " at " << lig_spec << std::endl;
   } else {

      try {
	 // this can throw an exception
	 RDKit::RWMol rdkm = coot::rdkit_mol_sanitized(residue_p, imol, *geometry_p);
	 
	 RDKit::MolChemicalFeatureFactory *factory = chemical_features::get_feature_factory();
	 if (! factory) {
	    std::cout << "WARNING:: no factory" << std::endl;
	    return v;
	 }

	 RDKit::FeatSPtrList features = factory->getFeaturesForMol(rdkm);

	 RDKit::Conformer conf = rdkm.getConformer(0); // iconf?

	 std::list<RDKit::FeatSPtr>::const_iterator it;
	 for (it=features.begin(); it!=features.end(); it++) {
	    RDKit::FeatSPtr feat_ptr = *it;
	    boost::shared_ptr<RDKit::MolChemicalFeature> sp = *it;
	    RDGeom::Point3D pos = sp.get()->getPos();
	    clipper::Coord_orth centre(pos.x, pos.y, pos.z);
	    std::string family = sp.get()->getFamily();
	    simple_chemical_feature_attributes cf(family, centre, imol, lig_spec);
	    v.push_back(cf);
	 }
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "WARNING:: " << rte.what() << std::endl;
      }
   }

   return v;
}


// return all chemical features
// 
std::vector<coot::simple_chemical_feature_attributes>
coot::chem_feat_clust::get_chemical_features() {

   std::vector<coot::simple_chemical_feature_attributes> v;

   for (unsigned int i=0; i<ligands.size(); i++) {
      std::vector<coot::simple_chemical_feature_attributes> cfs = 
	 get_chemical_features(ligands[i].imol, ligands[i].ligand_spec, ligands[i].mol);
      for (unsigned int j=0; j<cfs.size(); j++)
	 v.push_back(cfs[j]);
   }
   return v;
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

bool
coot::chem_feat_clust::fill_ligands(const std::vector<chem_feat_solvated_ligand_spec> &ligands_in) {

   bool success = true;
   for (unsigned int i=0; i<ligands_in.size(); i++) { 
      chem_feat_solvated_ligand lig(ligands_in[i]);
      if (! lig.residue) {
	 success = false;
	 break;
      }
      ligands.push_back(lig);
   }
   return success;
}

bool
coot::chem_feat_clust::check_dictionaries() {

   // I think that imol should be part of the class chem_feat_clust
   // for now I will fake up an imol
   int imol_fake = 0;
   bool success = true;
   for (unsigned int ilig=0; ilig<ligands.size(); ilig++) {
      mmdb::Residue *res = ligands[ilig].residue;
      if (res) {
	 // it should be set by now
	 std::string res_name = res->GetResName();

	 bool have = geometry_p->have_at_least_minimal_dictionary_for_residue_type(res_name, imol_fake);

	 if (! have) {
	    success = false;
	    break;
	 } 
      } 
   }
   return success;
}


void
coot::chem_feat_clust::fill_waters() {

   // This needs to be changed so that we only add waters that are
   // some cut-off (4.2?) from the position of any of the atoms of any
   // of the ligands
   //
   // we need to make a big vector of the positions of the ligand atoms
   // and check that we are close enough
   //
   // so first get a big vector of ligand coords:
   //
   std::vector<clipper::Coord_orth> lig_coords = get_ligands_coords();
   // 
   // double dist_cutoff = 4.2; // a class member data

   for (unsigned int ilig=0; ilig<ligands.size(); ilig++) {
      for (unsigned int iw=0; iw<ligands[ilig].waters.size(); iw++) {
	 mmdb::Residue *res = lidia_utils::get_residue(ligands[ilig].waters[iw], ligands[ilig].mol);
	 if (res) {
	    std::string res_name = res->GetResName();
	    if (res_name == "HOH") {
	       mmdb::Atom *at = res->GetAtom("O");
	       if (! at) {
		  std::cout << "Missing O at HOH in " << ligands[ilig].waters[iw]
			    << std::endl;
	       } else {
		  clipper::Coord_orth pt = lidia_utils::co(at);

		  if (is_close_to_a_ligand_atom(pt, lig_coords)) {

		     water_attribs wa(ilig, iw, at, pt);

		     if (false)
			std::cout << "fill_waters(): adding " << ilig << " "
				  << iw << " " << atom_spec_t(at) << std::endl;

		     water_positions.push_back(wa);
		  }
	       }
	    }
	 }
      }
   }
}

std::vector<clipper::Coord_orth>
coot::chem_feat_clust::get_ligands_coords() const {

   std::vector<clipper::Coord_orth> v;
   
   for (unsigned int i=0; i<ligands.size(); i++) {
      mmdb::Residue *res = ligands[i].residue;
      if (res) {
	 mmdb::PAtom *residue_atoms = 0;
	 int n_residue_atoms;
	 res->GetAtomTable(residue_atoms, n_residue_atoms);
	 for (int iat=0; iat<n_residue_atoms; iat++) { 
	    clipper::Coord_orth pt = lidia_utils::co(residue_atoms[iat]);
	    v.push_back(pt);
	 }
      }
   }

   return v;
}

bool
coot::chem_feat_clust::is_close_to_a_ligand_atom(const clipper::Coord_orth &pt_test,
						 const std::vector<clipper::Coord_orth> &ligand_atom_positions) const {

   bool s = false;
   double dc_sqrd = water_dist_cutoff * water_dist_cutoff;
   for (unsigned int i=0; i<ligand_atom_positions.size(); i++) { 
      if (clipper::Coord_orth(ligand_atom_positions[i] - pt_test).lengthsq() < dc_sqrd) {
	 s = true;
	 break;
      }
   }
   return s;
}




void
coot::chem_feat_clust::align() {

   std::cout << "missing alignment fuction " << std::endl;

}
   

#endif // MAKE_ENHANCED_LIGAND_TOOLS
