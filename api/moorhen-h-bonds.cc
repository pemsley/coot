
#include "molecules_container.hh"

#include "coot-utils/coot-h-bonds.hh"


std::vector<moorhen::h_bond>
molecules_container_t::get_h_bonds(int imol, const std::string &cid_str) const {

   mmdb::realtype max_dist = 3.8; // pass this

   std::vector<moorhen::h_bond> m_hbs;

   if (! is_valid_model_molecule(imol)) return m_hbs;

   mmdb::Manager *mol = get_mol(imol);
   coot::h_bonds hb;
   int SelHnd_all = mol->NewSelection();
   int SelHnd_lig = mol->NewSelection();
   mol->SelectAtoms(SelHnd_all, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", "*", "*");
   mol->Select(SelHnd_lig, mmdb::STYPE_ATOM, cid_str.c_str(), mmdb::SKEY_NEW);

   if (true) { // debug
      int nSelAtoms;
      mmdb::PPAtom local_SelAtom;
      mol->GetSelIndex(SelHnd_lig, local_SelAtom, nSelAtoms);
      for (int i=0; i<nSelAtoms; i++) {
         std::cout << i << "   " << local_SelAtom[i] << " " << coot::atom_spec_t(local_SelAtom[i]) << std::endl;
      }
   }

   // copy the attributes of atom_in into m_at
   auto copy_to_mh_atom = [] (mmdb::Atom *atom_in, moorhen::h_bond_atom &m_at) {

      if (atom_in) { // can be null (strange).

         m_at.serial = atom_in->serNum;
         m_at.x = atom_in->x;
         m_at.y = atom_in->y;
         m_at.z = atom_in->z;
         m_at.charge = atom_in->charge;
         m_at.occ = atom_in->occupancy;
         m_at.b_iso = atom_in->tempFactor;
         m_at.element = std::string(atom_in->element);
         m_at.name = std::string(atom_in->name);
         m_at.model = atom_in->GetModelNum();
         m_at.chain = std::string(atom_in->GetChainID());
         m_at.res_no = atom_in->GetSeqNum();
         m_at.residue_name = std::string(atom_in->GetResidue()->name);
         m_at.altLoc = std::string(atom_in->altLoc);
         m_at.charge = atom_in->charge;
         m_at.occ = atom_in->occupancy;
         m_at.b_iso = atom_in->tempFactor;
         m_at.element = std::string(atom_in->element);
      }
   };

   std::vector<coot::h_bond> hbonds = hb.get_mcdonald_and_thornton(SelHnd_lig, SelHnd_all, mol, geom, max_dist);

   for(unsigned ib=0;ib<hbonds.size();ib++) {

      moorhen::h_bond mhb;

      copy_to_mh_atom(hbonds[ib].hb_hydrogen,    mhb.hb_hydrogen);
      copy_to_mh_atom(hbonds[ib].donor,          mhb.donor);
      copy_to_mh_atom(hbonds[ib].donor_neigh,    mhb.donor_neigh);
      copy_to_mh_atom(hbonds[ib].acceptor,       mhb.acceptor);
      copy_to_mh_atom(hbonds[ib].acceptor_neigh, mhb.acceptor_neigh);

      mhb.angle_1 = hbonds[ib].angle_1;
      mhb.angle_2 = hbonds[ib].angle_2;
      mhb.angle_3 = hbonds[ib].angle_3;
      mhb.dist    = hbonds[ib].dist;
      mhb.ligand_atom_is_donor    = hbonds[ib].ligand_atom_is_donor;
      mhb.hydrogen_is_ligand_atom = hbonds[ib].hydrogen_is_ligand_atom;
      mhb.bond_has_hydrogen_flag  = hbonds[ib].bond_has_hydrogen_flag;
      m_hbs.push_back(mhb);
   }

   mol->DeleteSelection(SelHnd_lig);
   mol->DeleteSelection(SelHnd_all);

   return m_hbs;

}

