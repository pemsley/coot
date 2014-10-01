
#include <GraphMol/GraphMol.h>

#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <ForceField/ForceField.h>
#include <ForceField/MMFF/BondStretch.h>
#include <ForceField/MMFF/AngleBend.h>

#include <boost/python.hpp>
using namespace boost::python;

#define HAVE_GSL
#include <ideal/simple-restraint.hh>
#include <coot-utils/coot-coord-utils.hh>
#include <lidia-core/rdkit-interface.hh>

#include "py-restraints.hh"
#include "restraints-private.hh" // for bond-order conversion

namespace coot {

   // this uses atom indices
   class mmff_bond_restraint_info_t {
   public:
      mmff_bond_restraint_info_t() { sigma = -1;}
      mmff_bond_restraint_info_t(unsigned int idx_1_in,
				 unsigned int idx_2_in,
				 const std::string &type_in,
				 const double bl,
				 const double sigma_in) {
	 idx_1 = idx_1_in;
	 idx_2 = idx_2_in;
	 type = type_in;
	 resting_bond_length = bl;
	 sigma = sigma_in;
      }
      unsigned int idx_1;
      unsigned int idx_2;
      std::string type;
      double resting_bond_length;
      double sigma; // pseudo sigma based on K_{bond}
      unsigned int get_idx_1() const { return idx_1; } 
      unsigned int get_idx_2() const { return idx_2; }
      std::string get_type() const { return type; }
      double get_resting_bond_length() const { return resting_bond_length; }
      double get_sigma() const { return sigma; } 
   };

   // this uses atom indices
   class mmff_angle_restraint_info_t {
   public:
      mmff_angle_restraint_info_t() { sigma = -1;}
      mmff_angle_restraint_info_t(unsigned int idx_1_in,
				  unsigned int idx_2_in,
				  unsigned int idx_3_in,
				  const double angle,
				  const double sigma_in) {
	 idx_1 = idx_1_in;
	 idx_2 = idx_2_in;
	 idx_3 = idx_3_in;
	 resting_angle = angle;
	 sigma = sigma_in;
      }
      unsigned int idx_1;
      unsigned int idx_2;
      unsigned int idx_3;
      double resting_angle;
      double sigma; // pseudo sigma based on K_{angle}
      unsigned int get_idx_1() const { return idx_1; } 
      unsigned int get_idx_2() const { return idx_2; }
      unsigned int get_idx_3() const { return idx_3; }
      double get_resting_angle() const { return resting_angle; }
      double get_sigma() const { return sigma; } 
   };


   class mmff_b_a_restraints_container_t {
   public:
      std::vector<mmff_bond_restraint_info_t>  bonds;
      std::vector<mmff_angle_restraint_info_t> angles;
      mmff_b_a_restraints_container_t() { }
      unsigned int bonds_size() const { return bonds.size(); }
      unsigned int angles_size() const { return angles.size(); }

      // these will crash if you feed them an out-of-bounds index
      mmff_bond_restraint_info_t get_bond(const unsigned int i) {
	 return bonds[i];
      } 
      mmff_angle_restraint_info_t get_angle(const unsigned int i) {
	 return angles[i];
      }
      
   };

   RDKit::ROMol *regularize(RDKit::ROMol &r);
   RDKit::ROMol *regularize_with_dict(RDKit::ROMol &r,
				      PyObject *py_restraints,
				      const std::string &comp_id);
   // This tries to get a residue from the dictionary using the
   // model_Cartn of the dict_atoms.  If that is not available, return
   // a molecule with atoms that do not have positions.
   RDKit::ROMol *rdkit_mol_chem_comp_pdbx(const std::string &chem_comp_dict_file_name,
					  const std::string &comp_id);
   RDKit::ROMol *hydrogen_transformations(const RDKit::ROMol &r);
   RDKit::ROMol *mogulify(const RDKit::ROMol &r);
   mmff_b_a_restraints_container_t *mmff_bonds_and_angles(RDKit::ROMol &mol_in);

}


BOOST_PYTHON_MODULE(libpyrogen_boost) {
   def("regularize",               coot::regularize,               return_value_policy<manage_new_object>());
   def("regularize_with_dict",     coot::regularize_with_dict,     return_value_policy<manage_new_object>());
   def("rdkit_mol_chem_comp_pdbx", coot::rdkit_mol_chem_comp_pdbx, return_value_policy<manage_new_object>());
   def("hydrogen_transformations", coot::hydrogen_transformations, return_value_policy<manage_new_object>());
   def("mogulify",                 coot::mogulify,                 return_value_policy<manage_new_object>());
   def("mmff_bonds_and_angles",    coot::mmff_bonds_and_angles,    return_value_policy<manage_new_object>());

   class_<coot::mmff_bond_restraint_info_t>("mmff_bond_restraint_info_t")
      .def("get_idx_1",         &coot::mmff_bond_restraint_info_t::get_idx_1)
      .def("get_idx_2",         &coot::mmff_bond_restraint_info_t::get_idx_2)
      .def("get_type",          &coot::mmff_bond_restraint_info_t::get_type)
      .def("get_resting_bond_length", &coot::mmff_bond_restraint_info_t::get_resting_bond_length)
      .def("get_sigma",               &coot::mmff_bond_restraint_info_t::get_sigma)
      ;

   class_<coot::mmff_angle_restraint_info_t>("mmff_angle_restraint_info_t")
      .def("get_idx_1",    &coot::mmff_angle_restraint_info_t::get_idx_1)
      .def("get_idx_2",    &coot::mmff_angle_restraint_info_t::get_idx_2)
      .def("get_idx_3",    &coot::mmff_angle_restraint_info_t::get_idx_3)
      .def("get_resting_angle",  &coot::mmff_angle_restraint_info_t::get_resting_angle)
      .def("get_sigma",          &coot::mmff_angle_restraint_info_t::get_sigma)
      ;


   // established coot class, works with atom names though - so not useful ATM.
   class_<coot::dict_bond_restraint_t>("dict_bond_restraint_t")
      .def("atom_id_1",  &coot::dict_bond_restraint_t::atom_id_1)
      .def("atom_id_2",  &coot::dict_bond_restraint_t::atom_id_2)
      .def("type",       &coot::dict_bond_restraint_t::type)
      .def("value_dist", &coot::dict_bond_restraint_t::value_dist)
      .def("value_esd",  &coot::dict_bond_restraint_t::value_esd)
      ;
   
   
   class_<coot::mmff_b_a_restraints_container_t>("mmff_b_a_restraints_container_t")
      .def("bonds_size",  &coot::mmff_b_a_restraints_container_t::bonds_size)
      .def("angles_size", &coot::mmff_b_a_restraints_container_t::angles_size)
      .def("get_bond",    &coot::mmff_b_a_restraints_container_t::get_bond)
      .def("get_angle",   &coot::mmff_b_a_restraints_container_t::get_angle)
      ;


}


RDKit::ROMol*
coot::mogulify(const RDKit::ROMol &mol) {
   
   RDKit::RWMol rw(mol);
   coot::mogulify_mol(rw);
   RDKit::ROMol *ro = new RDKit::ROMol(rw);
   return ro;
}

coot::mmff_b_a_restraints_container_t *
coot::mmff_bonds_and_angles(RDKit::ROMol &mol) {

   mmff_b_a_restraints_container_t *r = new mmff_b_a_restraints_container_t;

   RDKit::MMFF::MMFFMolProperties *mmffMolProperties = new RDKit::MMFF::MMFFMolProperties(mol);
   if (! mmffMolProperties->isValid()) {
      std::cout << "invalid properties " << std::endl;
   } else {
      // happy path

      // iterate over bonds - simple
      // 
      ForceFields::MMFF::MMFFBondCollection *mmff_bonds =
	 ForceFields::MMFF::MMFFBondCollection::getMMFFBond();
      RDKit::ROMol::BondIterator bondIt;
      for (bondIt=mol.beginBonds(); bondIt!=mol.endBonds(); bondIt++) {
	 unsigned int idx_1 = (*bondIt)->getBeginAtomIdx();
	 unsigned int idx_2 = (*bondIt)->getEndAtomIdx();
	 unsigned int iAtomType_1 = mmffMolProperties->getMMFFAtomType(idx_1);
	 unsigned int iAtomType_2 = mmffMolProperties->getMMFFAtomType(idx_2);
	 unsigned int bondType  = mmffMolProperties->getMMFFBondType(*bondIt);
	 const ForceFields::MMFF::MMFFBond *mmffBondParams =
	    (*mmff_bonds)(bondType, iAtomType_1, iAtomType_2);
	 if (mmffBondParams) { 
	    double r0 = ForceFields::MMFF::Utils::calcBondRestLength(mmffBondParams);
	    double kb = ForceFields::MMFF::Utils::calcBondForceConstant(mmffBondParams);
	    double sigma = 0.04/sqrt(kb);
	    std::string order = convert_to_energy_lib_bond_type((*bondIt)->getBondType());
	    mmff_bond_restraint_info_t br(idx_1, idx_2, order, r0, sigma);
	    r->bonds.push_back(br);
	 }
      }

      
      // iterate over angles
      // 
      ForceFields::MMFF::MMFFAngleCollection *mmff_angles =
	 ForceFields::MMFF::MMFFAngleCollection::getMMFFAngle();
      unsigned int n_atoms = mol.getNumAtoms();
      std::map<unsigned long long, bool> done_angle;
      for (unsigned int iat_1=0; iat_1<n_atoms; iat_1++) { 
	 RDKit::ATOM_SPTR at_1 = mol[iat_1];
	 RDKit::ROMol::ADJ_ITER nbr_idx_1, end_nbrs_1;
	 boost::tie(nbr_idx_1, end_nbrs_1) = mol.getAtomNeighbors(at_1);
	 while(nbr_idx_1 != end_nbrs_1){
	    const RDKit::ATOM_SPTR at_2 = mol[*nbr_idx_1];

	    RDKit::ROMol::ADJ_ITER nbr_idx_2, end_nbrs_2;
	    boost::tie(nbr_idx_2, end_nbrs_2) = mol.getAtomNeighbors(at_2);
	    while(nbr_idx_2 != end_nbrs_2){
	       const RDKit::ATOM_SPTR at_3 = mol[*nbr_idx_2];
	       if (at_3 != at_1) {

		  unsigned int idx_1 = at_1->getIdx();
		  unsigned int idx_2 = at_2->getIdx();
		  unsigned int idx_3 = at_3->getIdx();

		  unsigned int m = 10000;
		  unsigned long long angle_key_1 = idx_1 * m * m + idx_2 * m + idx_3;
		  unsigned long long angle_key_2 = idx_3 * m * m + idx_2 * m + idx_1;

		  if (done_angle.find(angle_key_1) == done_angle.end() &&
		      done_angle.find(angle_key_2) == done_angle.end()) {

		     done_angle[m] = true;

		     unsigned int iAtomType_1 = mmffMolProperties->getMMFFAtomType(idx_1);
		     unsigned int iAtomType_2 = mmffMolProperties->getMMFFAtomType(idx_2);
		     unsigned int iAtomType_3 = mmffMolProperties->getMMFFAtomType(idx_3);

		     unsigned int angle_type =
			mmffMolProperties->getMMFFAngleType(mol, idx_1, idx_2, idx_3);

 		     const ForceFields::MMFF::MMFFAngle *mmffAngleParams =
 			(*mmff_angles)(angle_type, iAtomType_1, iAtomType_2, iAtomType_3);
		     
		     if (mmffAngleParams) {
			double a = ForceFields::MMFF::Utils::calcAngleRestValue(mmffAngleParams);
			double k = ForceFields::MMFF::Utils::calcAngleForceConstant(mmffAngleParams);
			double esd = 3.0/sqrt(k);
			if (0)
			   std::cout << idx_1 << " " << idx_2 << " " << idx_3 << "    "
				     << a << " " << k << std::endl;
			mmff_angle_restraint_info_t angle(idx_1, idx_2, idx_3, a, esd);
			r->angles.push_back(angle);
		     }
		  }
	       }
	       nbr_idx_2++;
	    }
	    nbr_idx_1++;
	 }
      }
   }
   return r;
}

RDKit::ROMol *
coot::regularize(RDKit::ROMol &mol_in) {

   RDKit::ROMol *m = new RDKit::ROMol(mol_in);
   return m;
}

RDKit::ROMol *
coot::regularize_with_dict(RDKit::ROMol &mol_in, PyObject *restraints_py, const std::string &comp_id) {

   coot::dictionary_residue_restraints_t dict_restraints = 
      monomer_restraints_from_python(restraints_py);
   RDKit::RWMol *m = new RDKit::RWMol(mol_in);
   mmdb::Residue *residue_p = make_residue(mol_in, 0, comp_id);
   if (! residue_p) {
      std::cout << "WARNING:: bad residue " << std::endl;
   } else {
      // deep copy residue and add to new molecule.
      mmdb::Manager *cmmdbmanager = util::create_mmdbmanager_from_residue(residue_p);
      mmdb::Residue *new_residue_p = coot::util::get_first_residue(cmmdbmanager);
      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      new_residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      std::cout << "------------------ simple_refine() called from "
		<< "restraints-boost.cc:regularize_with_dict()" << std::endl;
      simple_refine(new_residue_p, cmmdbmanager, dict_restraints);
      std::cout << "------------------ simple_refine() finished" << std::endl;
      int iconf = 0;
      update_coords(m, iconf, new_residue_p);
   }
   return static_cast<RDKit::ROMol *>(m);
}


RDKit::ROMol *
coot::rdkit_mol_chem_comp_pdbx(const std::string &chem_comp_dict_file_name,
			       const std::string &comp_id) {

   RDKit::ROMol *mol = new RDKit::ROMol;
   coot::protein_geometry geom;
   geom.set_verbose(false);
   int read_number = 0;
   geom.init_refmac_mon_lib(chem_comp_dict_file_name, read_number);
   bool idealized = false;
   bool try_autoload_if_needed = false;
   
   mmdb::Residue *r = geom.get_residue(comp_id, idealized, try_autoload_if_needed);

   std::pair<bool, dictionary_residue_restraints_t> rest = geom.get_monomer_restraints(comp_id);
   if (rest.first) {

      if (r) {
	 // makes a 3d conformer

	 // 20140618: This was false - now try true so that it processes 110 (i.e. no more
	 //           valence error).
	 bool sanitize = true;
	 RDKit::RWMol mol_rw = coot::rdkit_mol(r, rest.second, "", sanitize);

	 RDKit::MolOps::assignStereochemistry(mol_rw);
	 RDKit::ROMol *m = new RDKit::ROMol(mol_rw);

	 // Here test that the propety mmcif_chiral_volume_sign has
	 // been set on atoms of mol_rw and m
	 //
	 if (0) { 
	    for (unsigned int iat=0; iat<m->getNumAtoms(); iat++) { 
	       RDKit::ATOM_SPTR at_p = (*m)[iat];
	       try {
		  std::string name;
		  std::string cv;
		  at_p->getProp("name", name);
		  at_p->getProp("mmcif_chiral_volume_sign", cv);
		  std::cout << "m: name " << name << " " << cv << std::endl;
	       }
	       catch (const KeyErrorException &err) {
		  std::cout << "m: no name or cv " << std::endl;
	       }
	    }
	 }

	 // debug.  OK, so the bond orders are undelocalized here.
	 // debug_rdkit_molecule(&mol_rw);
      
	 return m;
      } else {

	 // Are you here unexpectedly?  That's because the input cif dictionary doesn't have
	 // _chem_comp_atom.x{yz} or model_Cartn_x{yz} or pdbx_model_Cartn_x{yz}_ideal for
	 // the given comp_id.
	 // 
	 std::cout << "INFO:: No 3d coords in dictionary : using 2d from dictionary bonds and angles."
		   << std::endl;
	 // makes a 2d conformer

	 // Bleugh.  We deal with mol vs m badly here.  
	 // 
	 RDKit::RWMol mol_rw = coot::rdkit_mol(rest.second);
	 RDKit::MolOps::sanitizeMol(mol_rw);
	 RDKit::ROMol *m = new RDKit::ROMol(mol_rw);
	 bool canon_orient = false;
	 bool clear_confs  = false;
	 int iconf = RDDepict::compute2DCoords(*m, NULL, canon_orient, clear_confs, 10, 20);
	 delete mol;
	 return m;
      }
   } 
   return mol;
}
     
RDKit::ROMol *
coot::hydrogen_transformations(const RDKit::ROMol &mol) {

   RDKit::RWMol *r = new RDKit::RWMol(mol);
   
   RDKit::ROMol *query_cooh = RDKit::SmartsToMol("[C^2](=O)O[H]");
   RDKit::ROMol *query_n = RDKit::SmartsToMol("[N^3;H2]");
   std::vector<RDKit::MatchVectType>  matches_cooh;
   std::vector<RDKit::MatchVectType>  matches_n;
   bool recursionPossible = true;
   bool useChirality = true;
   bool uniquify = true;

   // Note that we get the atom indexing using mol, but we change the
   // molecule r (I hope that the atoms are equivalently indexed).
   
   int matched_cooh = RDKit::SubstructMatch(mol,*query_cooh,matches_cooh,uniquify,recursionPossible, useChirality);
   int matched_n    = RDKit::SubstructMatch(mol,*query_n,   matches_n,   uniquify,recursionPossible, useChirality);

   if (1)
      std::cout << "COOH-NH2 hydrogen_transformations(): matches:  num COOH matches: "
		<< matches_cooh.size() << "  num N matches: " << matches_n.size() << std::endl;

   for (unsigned int imatch_cooh=0; imatch_cooh<matches_cooh.size(); imatch_cooh++) {
      std::cout << "INFO:: hydrogen exchanges matches_cooh: ";
      for (unsigned int i=0; i<matches_cooh[imatch_cooh].size(); i++) { 
	 std::cout << " [" << matches_cooh[imatch_cooh][i].first 
		   << ": "  << matches_cooh[imatch_cooh][i].second
		   << "]";
      }
      std::cout << std::endl;
   }
	 
   for (unsigned int imatch_cooh=0; imatch_cooh<matches_cooh.size(); imatch_cooh++) {
      
      RDKit::ATOM_SPTR at_c  = (*r)[matches_cooh[imatch_cooh][0].second];
      RDKit::ATOM_SPTR at_o1 = (*r)[matches_cooh[imatch_cooh][1].second];
      RDKit::ATOM_SPTR at_o2 = (*r)[matches_cooh[imatch_cooh][2].second];
      RDKit::ATOM_SPTR at_h  = (*r)[matches_cooh[imatch_cooh][3].second];

      RDKit::Bond *bond_1 = r->getBondBetweenAtoms(at_c.get()->getIdx(), at_o1.get()->getIdx());
      RDKit::Bond *bond_2 = r->getBondBetweenAtoms(at_c.get()->getIdx(), at_o2.get()->getIdx());
      RDKit::Bond *bond_3 = r->getBondBetweenAtoms(at_h.get()->getIdx(), at_o2.get()->getIdx());

      if (bond_1) {
	 if (bond_2) {
	 // bond_1->setBondType(RDKit::Bond::ONEANDAHALF);
	 // bond_2->setBondType(RDKit::Bond::ONEANDAHALF);
	 }
      }

      at_o2->setFormalCharge(-1);
      if (bond_3)
	 r->removeBond(at_o2.get()->getIdx(), at_h.get()->getIdx());
      else
	 std::cout << "no bond_3 found!" << std::endl;
      r->removeAtom(at_h.get());
   }

   // do we neet to sanitize? Yes, we do because we go on to minimize this molecule
   RDKit::MolOps::sanitizeMol(*r);
   
   RDKit::ROMol *ro_mol = new RDKit::ROMol(*r);
   if (0)
      std::cout << "hydrogen_transformations returns mol: " << RDKit::MolToSmiles(*ro_mol) << std::endl;
   delete r;
   return ro_mol;

}
