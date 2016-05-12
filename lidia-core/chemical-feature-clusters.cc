
// header-here

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#include <iostream>
#include <fstream>
#include <cstring>

#include "utils/coot-utils.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "rdkit-interface.hh"
#include "chemical-feature-clusters.hh"
#include <GraphMol/MolChemicalFeatures/MolChemicalFeature.h>


// do_alignment is an optional arg - default false
// 
coot::chem_feat_clust::chem_feat_clust(const std::vector<residue_spec_t> &protein_residues,
				       const std::vector<chem_feat_solvated_ligand_spec> &ligands_in,
				       const protein_geometry *geometry_p_in,
				       bool do_alignment) {

   geometry_p = geometry_p_in;
   setup_success = fill_ligands(ligands_in);
   if (setup_success) {
      setup_success = check_dictionaries();

      if (true) { // setup_success
	 if (do_alignment)
	    align();
	 
	 fill_waters();
      }
   }
}

void
coot::chem_feat_solvated_ligand::init_residue() {

   residue = util::get_residue(ligand_spec, mol);

   if (!residue) {
      std::cout << "WARNING:: null residue from spec " << ligand_spec << std::endl;
   }
}


bool
coot::chem_feat_clust::get_chemical_features(unsigned int idx,
					     residue_spec_t lig_spec,
					     mmdb::Manager *mol) {

   bool success = false;

   if (! setup_success)
      return success;   
	 
   mmdb::Residue *residue_p = util::get_residue(lig_spec, mol);

   if (! residue_p) {
      std::cout << "WARNING:: failed to get ligand " << idx << " at " << lig_spec << std::endl;
   } else {

      try {
	 // this can throw an exception
	 RDKit::RWMol rdkm = coot::rdkit_mol_sanitized(residue_p, *geometry_p);
	 
	 RDKit::MolChemicalFeatureFactory *factory = chemical_features::get_feature_factory();
	 if (! factory) {
	    std::cout << "WARNING:: no factory" << std::endl;
	    return success;
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
      
	 }
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "WARNING:: " << rte.what() << std::endl;
      }
   }
   return success;
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

   bool success = true;
   for (unsigned int ilig=0; ilig<ligands.size(); ilig++) {
      mmdb::Residue *res = ligands[ilig].residue;
      if (res) {
	 // it should be set by now
	 std::string res_name = res->GetResName();

	 bool have = geometry_p->have_at_least_minimal_dictionary_for_residue_type(res_name);

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

   for (unsigned int ilig=0; ilig<ligands.size(); ilig++) {
      for (unsigned int iw=0; iw<ligands[ilig].waters.size(); iw++) {
	 mmdb::Residue *res = util::get_residue(ligands[ilig].waters[iw],
						ligands[ilig].mol);
	 if (res) {
	    std::string res_name = res->GetResName();
	    if (res_name == "HOH") {
	       mmdb::Atom *at = res->GetAtom("O");
	       if (! at) {
		  std::cout << "Missing O at HOH in " << ligands[ilig].waters[iw]
			    << std::endl;
	       } else {
		  clipper::Coord_orth pt = co(at);
		  water_attribs wa(ilig, iw, at, pt);
		  water_positions.push_back(wa);
	       }
	    }
	 }
      }
   }
}

void
coot::chem_feat_clust::cluster_waters() {

   std::cout << "INFO:: clustering " << water_positions.size()
	     << " water positions" << std::endl;


}

void
coot::chem_feat_clust::cluster_ligand_chemical_features() {

   for (unsigned int i=0; i<ligands.size(); i++) { 
      get_chemical_features(i, ligands[i].ligand_spec, ligands[i].mol);
   }

}

void
coot::chem_feat_clust::align() {

   std::cout << "missing alignment fuction " << std::endl;

}
   

#endif // MAKE_ENHANCED_LIGAND_TOOLS
