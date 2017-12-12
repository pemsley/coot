/* ligand/multi-peptide.cc
 * 
 * Copyright 2010 by the University of Oxford
 * Copyright 2015 by Medical Research Council
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

#include "geometry/residue-and-atom-specs.hh"
#include "mini-mol/mini-mol-utils.hh"
#include "ideal/simple-restraint.hh"
#include "multi-peptide.hh"
#include "residue_by_phi_psi.hh"
#include "coot-utils/coot-map-utils.hh"
#include "analysis/stats.hh"


coot::minimol::fragment
coot::multi_build_N_terminal_ALA(mmdb::Residue *res_p,
				 const std::string &chain_id,
				 float b_factor_in,
				 int n_trials,
				 const coot::protein_geometry &geom,
				 const clipper::Xmap<float> &xmap,
				 std::pair<float, float> mv) {

   return multi_build_terminal_ALA(-1, res_p, chain_id, b_factor_in, n_trials, geom, xmap, mv);

}

coot::minimol::fragment
coot::multi_build_C_terminal_ALA(mmdb::Residue *res_p,
				 const std::string &chain_id,
				 float b_factor_in,
				 int n_trials,
				 const coot::protein_geometry &geom,
				 const clipper::Xmap<float> &xmap,
				 std::pair<float, float> mv) {

   return multi_build_terminal_ALA(1, res_p, chain_id, b_factor_in, n_trials, geom, xmap, mv);

}


coot::minimol::fragment
coot::multi_build_terminal_ALA(int offset, // direction
			       mmdb::Residue *res_p,
			       const std::string &chain_id,
			       float b_factor_in,
			       int n_trials,
			       const coot::protein_geometry &geom,
			       const clipper::Xmap<float> &xmap,
			       std::pair<float, float> mv) {


   
   // res_p might not be a member of hierarchy, so pass the chain id too.

   minimol::fragment many_residues;
   // std::string chain_id = res_p->GetChainID();

   bool prev_happy_fit = true;
   bool happy_fit = true; // starting
   int icount = 0;

   while (happy_fit) {

      std::string dir = "N";
      if (offset == 1)
	 dir = "C";

      std::string residue_type = "ALA";
      residue_type = "GLY"; // does good?
      
      residue_by_phi_psi rpp(dir, res_p, chain_id, residue_type, 20);
      rpp.import_map_from(xmap);

      // res from f is what we really want
      minimol::residue res;
      minimol::fragment f = rpp.best_fit_phi_psi(n_trials, offset);
      if (offset == -1) { 
	 for (int ires=f.max_residue_number(); ires>=f.min_res_no(); ires--) {
	    if (f[ires].atoms.size() > 0) {
	       res = f[ires];
	       break;
	    }
	 }
      } else {
	 // for C-terminal builds, get the first residue.
	 for (int ires=f.min_res_no(); ires<=f.max_residue_number(); ires++) {
	    if (f[ires].atoms.size() > 0) {
	       res = f[ires];
	       break;
	    }
	 }
      } 

      if (res.is_empty()) {
	 happy_fit = false;
      } else {

	 std::cout << "   debug::multi_build_terminal_ALA() with offset "
		   << offset << " was passed residue seqnum " << res_p->GetSeqNum()
		   << " and found res with seqnum " << res.seqnum << std::endl;
	 
	 std::pair<bool, clipper::Coord_orth> cbeta_info = cbeta_position(res);
	 if (cbeta_info.first) { 
	    res.addatom(" CB ", " C", cbeta_info.second, "", 1.0, 20.0);
	    res.name = "ALA";
	 }

	 if (offset == -1)
	    res.seqnum = res_p->GetSeqNum() -1;
	 if (offset ==  1)
	    res.seqnum = res_p->GetSeqNum() + 1;
	    

	 bool this_happy_fit = does_residue_fit(res, xmap, mv);

 	 if (! this_happy_fit && ! prev_happy_fit) {
 	    happy_fit = false;
	    std::cout << "   info:: 2 unhappy fits - stopping now " << std::endl;
	 } 

	 if (happy_fit) { 
	    icount++;
	    // add res to mol and refine
	    many_residues.addresidue(res, 20);

	    if (false) {
	       // debug:  Let's write a pdb file for this fragment
	       coot::minimol::molecule m_tmp;
	       m_tmp.fragments.push_back(many_residues);
	       std::string tmp_filename = "phi-psi-";
	       tmp_filename += util::int_to_string(icount); 
	       tmp_filename += ".pdb";
	       m_tmp.write_file(tmp_filename, 10);
	    }

	    // refine res and its neighbour, update the atoms of many_residues
	    int seqnum=res.seqnum;
	    refine_end(&many_residues, seqnum, offset, geom, xmap);

	    // ------------- For next round ----------------
	    // 
	    // set res_p to be the newly constructed residue
	    res_p = res.make_residue();
	    prev_happy_fit = this_happy_fit;
	 }
      }

      // synthetic stop
      // 
      // if (icount == 40)
      // happy_fit = false;
      
   }
   return many_residues;
}

// move the atoms (the end residue and its neighbour) of many_residue
// 
// offset tells us which direction the new residue was added (-1 means
// new N terminus)
void
coot::refine_end(coot::minimol::fragment *many_residues,
		 int seqnum, int offset,
		 const coot::protein_geometry &geom,
		 const clipper::Xmap<float> &xmap) {


   mmdb::Manager *mol = new mmdb::Manager;
   mmdb::Model *model_p = new mmdb::Model;
   int ires_new = -1;

   if (offset == -1)
      ires_new = many_residues->first_residue(); // with atoms

   mmdb::Chain *chain_p = many_residues->make_chain();

   model_p->AddChain(chain_p);
   mol->AddModel(model_p);
 
   mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
   mol->FinishStructEdit();

   std::vector<std::pair<bool, mmdb::Residue *> > residues;
   std::vector<mmdb::Residue *> moving_residues; 

   int nres = chain_p->GetNumberOfResidues();

   mmdb::Residue *residue_p;
   for (int ires=0; ires<nres; ires++) {
      residue_p = chain_p->GetResidue(ires);
      std::pair<bool, mmdb::Residue *> p(false, residue_p);
      int res_seq_num = residue_p->GetSeqNum();

      if (res_seq_num == ires_new) { 
	 moving_residues.push_back(residue_p);
	 residues.push_back(p);
      } 
      if (offset == -1) { 
	 if (res_seq_num == (ires_new + 1)) { // the next residue after the N-terminal residue
	    moving_residues.push_back(residue_p);
	    residues.push_back(p);
	 }
	 if (res_seq_num == (ires_new + 2)) {
	    std::pair<bool, mmdb::Residue *> fixed_p(false, residue_p);
	    residues.push_back(fixed_p);
	 }
      }
      if (offset == 1) { 
	 if (res_seq_num == (ires_new - 1)) { // the residue before the C-terminal residue
	    moving_residues.push_back(residue_p);
	    residues.push_back(p);
	 }
	 if (res_seq_num == (ires_new - 2)) {
	    std::pair<bool, mmdb::Residue *> fixed_p(false, residue_p);
	    residues.push_back(fixed_p);
	 }
      }
   }

   // std::cout << "   chain_p has " << nres << " residues" << std::endl;
   
   std::vector<mmdb::Link> links;
   std::vector<atom_spec_t> fixed_atom_specs;
   restraint_usage_Flags flags = TYPICAL_RESTRAINTS;

   restraints_container_t restraints(residues, links, geom, mol, fixed_atom_specs, xmap);
   pseudo_restraint_bond_type pseudos = NO_PSEUDO_BONDS;
   bool do_internal_torsions = false;
   float weight = 60;

   restraints.set_quiet_reporting();
   restraints.add_map(weight);
   bool do_trans_peptide_restraints = true;
   int imol = 0;
   restraints.make_restraints(imol, geom, flags, do_internal_torsions, do_internal_torsions, 0, 0, pseudos);
   restraints.minimize(flags);

   for (unsigned int ii=0; ii<moving_residues.size(); ii++) {
      int seqnum_l = moving_residues[ii]->GetSeqNum();
      (*many_residues)[seqnum_l].update_positions_from(moving_residues[ii]);
   }

   delete mol;
} 

bool
coot::does_residue_fit(const coot::minimol::residue &res, const clipper::Xmap<float> &xmap,
		       std::pair<float, float> mv) {

   bool r = true;
   double z_crit = 1.3;
   
   double rmsd = sqrt(mv.second);

   std::vector<double> rho(res.atoms.size());
   for (unsigned int iat=0; iat<res.atoms.size(); iat++) { 
      float d = util::density_at_point(xmap, res.atoms[iat].pos);
      rho[iat] = d;
   }

   stats::single st(rho);
   for (unsigned int iat=0; iat<res.atoms.size(); iat++) {
      if (rho[iat] < (mv.first + rmsd * z_crit)) {

	 if (res.atoms[iat].name != " CB ") {  // PDBv3 FIXME
	    std::cout << "   low density for atom residue: " << res.seqnum
		      << " atom: " << res.atoms[iat].name
		      << rho[iat] << " vs " << mv.first <<  " + " << rmsd << " * " << z_crit << " at "
		      << res.atoms[iat].pos.format()
		      << std::endl;
	    r = false;
	    break;
	 }
      } 
   }

   return r;
} 
