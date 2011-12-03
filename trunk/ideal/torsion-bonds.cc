
#include <algorithm>
#include "simple-restraint.hh"
#include "coot-coord-extras.hh"

// this can throw an exception
// 
std::vector<std::pair<CAtom *, CAtom *> >
coot::torsionable_bonds(CMMDBManager *mol, PPCAtom atom_selection,
			int n_selected_atoms,
			coot::protein_geometry *geom_p) { 

   std::vector<std::pair<CAtom *, CAtom *> > v;
   bool include_pyranose_ring_torsions_flag = false;

   std::vector<CResidue *> residues;
   std::map<CResidue *, std::vector<int> > atoms_in_residue;
   // fill residues and atoms_in_residue
   for (unsigned int i=0; i<n_selected_atoms; i++) {
      CResidue *r = atom_selection[i]->residue;
      if (std::find(residues.begin(), residues.end(), r) == residues.end())
	 residues.push_back(r);
      atoms_in_residue[r].push_back(i);
   }
   std::map<CResidue *, coot::dictionary_residue_restraints_t> res_restraints;
   for (unsigned int ires=0; ires<residues.size(); ires++) { 
      std::string rn = residues[ires]->GetResName();
      std::pair<bool, coot::dictionary_residue_restraints_t> rest =
	 geom_p->get_monomer_restraints(rn);
      if (! rest.first) {
	 std::string m = "Restraints not found for type ";
	 m += rn;
	 throw std::runtime_error(m);
      }
      res_restraints[residues[ires]] = rest.second;
   }

   for (unsigned int ires=0; ires<residues.size(); ires++) {
      // a coot-coord-extras function
      std::vector<std::pair<CAtom *, CAtom *> > v_inner =
	 coot::torsionable_bonds_monomer_internal(residues[ires], atom_selection, n_selected_atoms,
						  include_pyranose_ring_torsions_flag, geom_p);
      // std::cout << "found " << v_inner.size() << " monomer internal torsions for "
      // << residues[ires]->GetResName() << std::endl;
      for (unsigned int ip=0; ip<v_inner.size(); ip++)
	 v.push_back(v_inner[ip]);
   }

   std::vector<std::pair<CAtom *, CAtom *> > v_link =    
      coot::torsionable_link_bonds(residues, mol, geom_p);
   for (unsigned int il=0; il<v_link.size(); il++)
      v.push_back(v_link[il]);
   
   for (unsigned int ipair=0; ipair<v.size(); ipair++) { 
      std::cout << "   torsionable bond: "
		<< coot::atom_spec_t(v[ipair].first) << "  "
		<< coot::atom_spec_t(v[ipair].second)
		<< std::endl;
   }
   return v;
}

std::vector<std::pair<CAtom *, CAtom *> >
coot::torsionable_link_bonds(std::vector<CResidue *> residues_in,
			     CMMDBManager *mol, coot::protein_geometry *geom_p) {

   std::vector<std::pair<CAtom *, CAtom *> > v;
   
   std::vector<std::pair<bool, CResidue *> > residues(residues_in.size());
   for (unsigned int i=0; i<residues_in.size(); i++)
      residues[i] = std::pair<bool, CResidue *> (0, residues_in[i]);

   std::vector<coot::atom_spec_t> dummy_fixed_atom_specs;
   coot::restraints_container_t restraints(residues, *geom_p, mol, dummy_fixed_atom_specs);
   coot::bonded_pair_container_t bpc = restraints.bonded_residues_from_res_vec(*geom_p);

   // add in the torsion: CB-CG-ND2-C1 (Psi-N)
   // add in the torsion: CG-ND2-C1-O5 (Phi-N)
   // t.geom.link_add_torsion("NAG-ASN", 1, 2, 2, 2, " C1 ", "ND2 ", " CG ", " CB ", 180, 40, 3, "Psi-N");
   // t.geom.link_add_torsion("NAG-ASN", 1, 1, 2, 2, " O5 ", " C1 ", "ND2 ", " CG ", 180, 40, 3, "Psi-N");
      
   // std::cout << "found LINKR linked pairs:\n   " <<  bpc;

   for (unsigned int i=0; i<bpc.bonded_residues.size(); i++) { 
      coot::dictionary_residue_link_restraints_t link = geom_p->link(bpc[i].link_type);
      if (link.link_id != "") {
	 if (0) 
	    std::cout << "   dictionary link found " << link.link_id << " with "
		      << link.link_bond_restraint.size() << " bond restraints and "
		      << link.link_torsion_restraint.size() << " link torsions " << std::endl;
	 for (unsigned int ib=0; ib<link.link_bond_restraint.size(); ib++) { 
	    // we need to get the atoms and add them to "pairs".

	    if (0) 
	       std::cout << "   "
			 << link.link_bond_restraint[ib].atom_id_1_4c() << " "
			 << link.link_bond_restraint[ib].atom_1_comp_id << " to "
			 << link.link_bond_restraint[ib].atom_id_2_4c() << " " 
			 << link.link_bond_restraint[ib].atom_2_comp_id << " "
			 << " of "
			 << bpc[i].res_1->GetResName() << " to "
			 << bpc[i].res_2->GetResName()
			 << std::endl;
	    CAtom *link_atom_1 = bpc[i].res_1->GetAtom(link.link_bond_restraint[ib].atom_id_1_4c().c_str());
	    CAtom *link_atom_2 = bpc[i].res_2->GetAtom(link.link_bond_restraint[ib].atom_id_2_4c().c_str());
	    if (link_atom_1 && link_atom_2) { 
	       std::pair<CAtom *, CAtom *> pair(link_atom_1, link_atom_2);
	       v.push_back(pair);
	    }
	 }
      }
   }
   return v;
} 


// this can throw an exception.
void
coot::multi_residue_torsion_fit_map(CMMDBManager *mol,
				    const clipper::Xmap<float> &xmap,
				    coot::protein_geometry *geom_p) {

   try { 
      PPCAtom atom_selection = 0;
      int n_selected_atoms;
      int selhnd = mol->NewSelection();
      mol->SelectAtoms(selhnd, 0, "*",
		       ANY_RES, "*",
		       ANY_RES, "*",
		       "*", // residue name
		       "*",
		       "*", 
		       "*"); // alt-loc
      mol->GetSelIndex(selhnd, atom_selection, n_selected_atoms);
      std::vector<std::pair<CAtom *, CAtom *> > pairs = 
 	 coot::torsionable_bonds(mol, atom_selection, n_selected_atoms, geom_p);
      mol->DeleteSelection(selhnd);
   } 
   catch (std::runtime_error rte) {
      std::cout << "WARNING:: " << rte.what() << std::endl;
   } 
   // coot::atom_tree_t tree();
   
} 
