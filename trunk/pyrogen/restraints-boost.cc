
#include <GraphMol/GraphMol.h>

#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <ForceField/ForceField.h>

#include <boost/python.hpp>
using namespace boost::python;

#define HAVE_GSL
#include <ideal/simple-restraint.hh>
#include <coot-utils/coot-coord-utils.hh>
#include <lidia-core/rdkit-interface.hh>

#include "py-restraints.hh"


namespace coot {

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
   void mmff_stuff(RDKit::ROMol &mol_in);

}


BOOST_PYTHON_MODULE(libpyrogen_boost) {
   def("regularize",               coot::regularize,               return_value_policy<manage_new_object>());
   def("regularize_with_dict",     coot::regularize_with_dict,     return_value_policy<manage_new_object>());
   def("rdkit_mol_chem_comp_pdbx", coot::rdkit_mol_chem_comp_pdbx, return_value_policy<manage_new_object>());
   def("hydrogen_transformations", coot::hydrogen_transformations, return_value_policy<manage_new_object>());
   def("mogulify",                 coot::mogulify,                 return_value_policy<manage_new_object>());
   def("mmff_stuff",               coot::mmff_stuff);
}


RDKit::ROMol*
coot::mogulify(const RDKit::ROMol &mol) {
   
   RDKit::RWMol rw(mol);
   coot::mogulify_mol(rw);
   RDKit::ROMol *ro = new RDKit::ROMol(rw);
   return ro;
}

void coot::mmff_stuff(RDKit::ROMol &mol_in) {

   RDKit::MMFF::MMFFMolProperties *mmffMolProperties = new RDKit::MMFF::MMFFMolProperties(mol_in);
   ForceFields::ForceField *field = new ForceFields::ForceField();
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

      std::cout << "debug at_h " << at_h << std::endl;
      std::cout << "debug at_h.get() " << at_h.get() << std::endl;

      if (bond_1) {
	 if (bond_2) {
	 std::cout << "found bond 1 and bond 2  " << bond_1 << std::endl;
	 // bond_1->setBondType(RDKit::Bond::ONEANDAHALF);
	 // bond_2->setBondType(RDKit::Bond::ONEANDAHALF);
	 }
      }

      std::cout << "done bond type setting " << std::endl;
      
      at_o2->setFormalCharge(-1);
      std::cout << "done formal charge setting " << std::endl;
      if (bond_3)
	 r->removeBond(at_o2.get()->getIdx(), at_h.get()->getIdx());
      else
	 std::cout << "no bond_3 found!" << std::endl;
      std::cout << "done remove bond_3" << std::endl;
      
      r->removeAtom(at_h.get());
      std::cout << "done remove atom " << std::endl;

   }

   // do we neet to sanitize? Yes, we do because we go on to minimize this molecule
   std::cout << "sanitizing r... " << std::endl;
   RDKit::MolOps::sanitizeMol(*r);
   std::cout << "done sanitizing r... " << std::endl;
   
   RDKit::ROMol *ro_mol = new RDKit::ROMol(*r);
   if (0)
      std::cout << "hydrogen_transformations returns mol: " << RDKit::MolToSmiles(*ro_mol) << std::endl;
   delete r;
   return ro_mol;

}
