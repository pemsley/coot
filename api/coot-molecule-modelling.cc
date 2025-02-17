
#include "coot-utils/atom-tree.hh"
#include "coot-utils/coot-coord-extras.hh"
#include "coot-molecule.hh"

int
coot::molecule_t::replace_residue(const std::string &residue_cid, const std::string &new_residue_type,
                                  int imol_enc, const coot::protein_geometry &geom) {

   int status = 0;

   mmdb::Residue *residue_p = cid_to_residue(residue_cid);
   if (residue_p) {

      std::pair<bool, dictionary_residue_restraints_t> rp = geom.get_monomer_restraints(new_residue_type, imol_enc);
      if (rp.first) {
      const auto &restraints_new_type = rp.second;

         std::string current_residue_type = residue_p->GetResName();
         std::pair<bool, dictionary_residue_restraints_t> rp_current = geom.get_monomer_restraints(current_residue_type, imol_enc);
         if (rp_current.first) {
            const auto &restraints_current_type = rp_current.second;
            status = util::mutate_by_overlap(residue_p, atom_sel.mol, restraints_current_type, restraints_new_type);
         }
      }
   }
   return status;
}


int
coot::molecule_t::rotate_around_bond(const std::string &residue_cid,
                                     const std::string &alt_conf,
                                     coot::atom_name_quad quad,
                                     double torsion_angle,
                                     coot::protein_geometry &geom) {

   int status = 0;
   double r = -999.9;

   mmdb::Residue *residue_p = cid_to_residue(residue_cid);
   if (residue_p) {
      std::string res_name(residue_p->GetResName());
      std::pair<bool, coot::dictionary_residue_restraints_t> restraints_info =
         geom.get_monomer_restraints(res_name, imol_no);
      if (! restraints_info.first) {
         std::cout << "WARNING:: set_torsion: No restraints for " << res_name << std::endl;
      } else {
         coot::atom_tree_t tree(restraints_info.second, residue_p, alt_conf);

         try {
            r = tree.set_dihedral(quad.atom_name(0),
                                  quad.atom_name(1),
                                  quad.atom_name(2),
                                  quad.atom_name(3),
                                  torsion_angle);
            atom_sel.mol->FinishStructEdit();
         }
         catch(const std::runtime_error &rte) {
            std::cout << "in set_torsion:: set_dihedral() error: " << rte.what() << std::endl;
         }

      }
   } else {
      std::cout << "failed to find residue " << residue_cid << std::endl;
   }
   return status;
}


//! copy chain using NCS matrix
bool
coot::molecule_t::copy_ncs_chain(const std::string &from_chain_id, const std::string &to_chain_id) {

   bool status = false;

   return status;

}


void
coot::molecule_t::execute_simple_nucleotide_addition(const std::string &cid,
                                                     mmdb::Manager *standard_residues_asc_mol) {

   mmdb::Residue *residue_p = cid_to_residue(cid);
   if (residue_p) {
      std::string chain_id = residue_p->GetChainID();
   } else {
      std::cout << "WARNING:: failed to find residue " << cid << std::endl;
   }

}

void
coot::molecule_t::execute_simple_nucleotide_addition(mmdb::Residue *residue_p,
                                                     mmdb::Manager *standard_residues_asc_mol) {

   if (! residue_p) {
      std::cout << "WARNING:: " << __FUNCTION__ << "() missing-residue " << std::endl;
   } else {
      int res_no = residue_p->GetSeqNum();
      std::string term_type = "";
      std::string chain_id = residue_p->GetChainID();
      coot::residue_spec_t r_p_spec(chain_id, res_no - 1, "");
      coot::residue_spec_t r_n_spec(chain_id, res_no + 1, "");
      mmdb::Residue *r_p = get_residue(r_p_spec);
      mmdb::Residue *r_n = get_residue(r_n_spec);
      if (r_p  && ! r_n) term_type = "C";
      if (r_n  && ! r_p) term_type = "N";
      if (!r_n && ! r_p) term_type = "MC";
      execute_simple_nucleotide_addition(term_type, residue_p, chain_id, standard_residues_asc_mol);
   }
}

#include "ligand/ideal-rna.hh"

void
coot::molecule_t::execute_simple_nucleotide_addition(const std::string &term_type,
                                                     mmdb::Residue *res_p, const std::string &chain_id,
                                                     mmdb::Manager *standard_residues_asc_mol) {
   // If it's RNA beam it in in ideal A form,
   // If it's DNA beam it in in ideal B form

   // What's the plan?
   //
   // OK the plan is to generate a 2 residue molecule of
   // single-stranded RNA (or DNA).
   //
   // Depending on if this is N or C terminal type, we define the
   // sequence, adding a "base" residue (that we'll use to match to
   // res_p);

   if (term_type == "not-terminal-residue") {

      // maybe return this
      std::cout << "WARNING:: That was not a terminal residue (check for neighbour solvent residues maybe) "
                << coot::residue_spec_t(res_p) << std::endl;
   } else {

      std::string seq = "aa";
      std::string RNA_or_DNA_str = "RNA";
      std::string form_str = "A";
      bool single_stranded_flag = 1;

      if (coot::util::nucleotide_is_DNA(res_p)) {
	 RNA_or_DNA_str = "DNA";
	 form_str = "B";
      }

      coot::ideal_rna ir(RNA_or_DNA_str, form_str, single_stranded_flag,
			 seq, standard_residues_asc_mol);
      ir.use_v3_names();
      mmdb::Manager *mol = ir.make_molecule();

      int match_resno;
      int interesting_resno;
      if (term_type == "C" || term_type == "MC") {
	 match_resno = 1;
	 interesting_resno = 2;
      } else {
	 interesting_resno = 1;
	 match_resno = 2;
      }

      mmdb::Residue *moving_residue_p = NULL;
      mmdb::Residue *interesting_residue_p = NULL;
      int imod = 1;
      // now set moving_residue_p and interesting_residue_p:
      mmdb::Model *model_p = mol->GetModel(imod);
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 mmdb::PResidue residue_p;
	 for (int ires=0; ires<nres; ires++) {
	    residue_p = chain_p->GetResidue(ires);
	    // 	 std::cout << "testing vs resno " << residue_p->GetSeqNum()
	    // 		   << std::endl;
	    if (residue_p->GetSeqNum() == match_resno) {
	       moving_residue_p = residue_p;
	    }
	    if (residue_p->GetSeqNum() == interesting_resno) {
	       interesting_residue_p = residue_p;
	    }
	    if (moving_residue_p && interesting_residue_p)
	       break;
	 }
	 if (moving_residue_p && interesting_residue_p)
	    break;
      }

      if (interesting_residue_p) {
	 if (moving_residue_p) {
	    bool use_old_names = false;
	    std::pair<bool, clipper::RTop_orth> rtop_pair =
	       coot::util::nucleotide_to_nucleotide(res_p, moving_residue_p,
						    use_old_names);

	    // now apply rtop to mol:
	    if (rtop_pair.first) {
	       // fix up the residue number and chain id to match the clicked atom
	       int new_resno = res_p->GetSeqNum() + interesting_resno - match_resno;
	       interesting_residue_p->seqNum = new_resno;

               // we always want to remove OP3 from the residue to which a new residue
               // is added when we add to the "N-terminus"
               //
               if (term_type == "N" || term_type == "MN") {
                  mmdb::Atom **residue_atoms = 0;
                  int n_residue_atoms = 0;
                  bool deleted = false;
                  res_p->GetAtomTable(residue_atoms, n_residue_atoms);
                  for (int iat=0; iat<n_residue_atoms; iat++) {
                     mmdb::Atom *at = residue_atoms[iat];
                     if (at) {
                        std::string at_name(at->name);
                        if (at_name == " OP3") {  // PDBv3 FIXME
                           delete at;
                           at = NULL;
                           deleted = true;
                           break;
                        }
                     }
                  }
                  if (deleted)
                     res_p->TrimAtomTable();
               }

	       coot::util::transform_mol(mol, rtop_pair.second);
	       // byte gz = GZM_NONE;
	       // mol->WritePDBASCII("overlapped.pdb", gz);
	       mmdb::Manager *residue_mol =
		  coot::util::create_mmdbmanager_from_residue(interesting_residue_p);

	       atom_selection_container_t asc = make_asc(residue_mol);
	       // set the chain id of the chain that contains interesting_residue_p:
	       model_p = residue_mol->GetModel(imod);
	       // run over chains of the existing mol
	       nchains = model_p->GetNumberOfChains();
	       for (int ichain=0; ichain<nchains; ichain++) {
                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
		  int nres = chain_p->GetNumberOfResidues();
		  for (int ires=0; ires<nres; ires++) {
                     mmdb::PResidue residue_p = chain_p->GetResidue(ires);
		     if (residue_p->GetSeqNum() == interesting_residue_p->GetSeqNum()) {
			chain_p->SetChainID(chain_id.c_str());
		     }
		  }
	       }

               // now include asc into this atom selection
               insert_coords_internal(asc);
	    }
	 }
      } else {
	 std::cout << "Failed to find interesting residue (with resno " << interesting_resno
		   << ")" << std::endl;
      }
      delete mol;
   }

}


// insert coords - c.f function in molecule-class-info_t
void
coot::molecule_t::insert_coords_internal(const atom_selection_container_t &asc) {

   // run over each chain, residue of the asc (if terminal residue
   // fit only one chain, one residue, of course).

   bool inserted = false; // not inserted yet
   int imod = 1;
   mmdb::Model *asc_model_p = asc.mol->GetModel(imod);
   int asc_n_chains = asc_model_p->GetNumberOfChains();
   for (int i_asc_chain=0; i_asc_chain<asc_n_chains; i_asc_chain++) {
      mmdb::Chain *asc_chain = asc.mol->GetChain(1,i_asc_chain);
      int nres_asc = asc_chain->GetNumberOfResidues();
      int udd_atom_index = asc.UDDAtomIndexHandle;
      for (int ires_asc=0; ires_asc<nres_asc; ires_asc++) {
         mmdb::Residue *asc_residue = asc_chain->GetResidue(ires_asc);

         // Now find the corresponding chain in our atom_sel.mol:
         int imodel = 1;
         int n_chains = atom_sel.mol->GetNumberOfChains(imodel);
         for (int i_chain=0; i_chain<n_chains; i_chain++) {
            mmdb::Chain *chain = atom_sel.mol->GetChain(1,i_chain);

            // test chains
            std::string asc_chain_str(asc_chain->GetChainID());
            std::string mol_chain_str(    chain->GetChainID());
            if (asc_chain_str == mol_chain_str) {

               bool embed_in_chain_flag = false;
               // use the coot-utils function, not this one
               mmdb::PResidue res = coot::deep_copy_this_residue_old_style(asc_residue, "", 1,
                                                                           udd_atom_index, embed_in_chain_flag);
               std::pair<int, mmdb::Residue *> serial_number =
                  find_serial_number_for_insert(asc_residue->GetSeqNum(),
                                                asc_residue->GetInsCode(),
                                                mol_chain_str);

               if (res) {

                  if (serial_number.first != -1) {
                     // insert at this position (other residues are
                     // shifted up).
                     chain->InsResidue(res, serial_number.first);
                     coot::copy_segid(serial_number.second, res);
                     inserted = true;
                  } else {

                     // std::cout << "DEBUG:: insert_coords_internal() add residue\n";
                     mmdb::Residue *last_residue = util::get_last_residue_in_chain(chain);
                     if (last_residue) {
                        chain->AddResidue(res);
                        coot::copy_segid(last_residue, res);
                        inserted = 1;
                     }
                  }
               }
            }
         }


         if (! inserted) {
            // OK, there was no chain in the current mol that matches
            // the chain of the asc.
            // Let's copy the asc chain and add it to atom_sel.mol
            mmdb::Chain *new_chain = new mmdb::Chain;
            int imodel = 1;
            mmdb::Model *this_model = atom_sel.mol->GetModel(imodel);
            this_model->AddChain(new_chain);
            new_chain->SetChainID(asc_chain->GetChainID());

            std::cout << "DEBUG:: Creating a new chain " << asc_chain->GetChainID() << std::endl;
            bool embed_in_chain_flag = false;
            // use the coot-utils function, not this one
            mmdb::Residue *res = coot::deep_copy_this_residue_old_style(asc_residue, "", 1,
                                                                        udd_atom_index, embed_in_chain_flag);
            if (res) {
               new_chain->AddResidue(res);
               atom_sel.mol->FinishStructEdit(); // so that we don't keep adding a
                                                 // new Chain to atom_sel.mol
            }

         }
      }
   }
   atom_sel.mol->FinishStructEdit();
}
