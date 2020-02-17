
// ----------------------------------------
// How to use the thread pool for threading
// ----------------------------------------

// convert this:

      for(int i=0; i<n_atoms_max; i++) {
	 mmdb::Atom *at = asc.atom_selection[i];
	 clipper::Coord_orth pt = coot::co(at);
	 clipper::Coord_map cm_try_2(rtop_og * pt);
	 float dn = coot::util::density_at_point_by_cubic_interp(nxmap, cm_try_2);
	 sum += dn;
      }
      auto tp_2 = std::chrono::high_resolution_clock::now();

// to this:

#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
      unsigned int n_threads = 3; // some number less than max proc
      ctpl::thread_pool thread_pool(n_threads);
      std::atomic<unsigned int> done_count_for_threads(0);
      std::vector<float> dv(n_threads, 0.0);
      std::vector<std::pair<unsigned int, unsigned int> > ranges =
	 coot::atom_index_ranges(n_atoms_max, n_threads);
      for (std::size_t i=0; i<ranges.size(); i++) {
	 thread_pool.push(density_for_atoms_multithread,
			  std::cref(asc),
			  std::cref(rtop_og),
			  std::cref(ranges[i]),
			  std::cref(nxmap),
			  &dv[i],
			  std::ref(done_count_for_threads));
      }
      while (done_count_for_threads < ranges.size())
	 std::this_thread::sleep_for(std::chrono::microseconds(1));
      for (std::size_t i=0; i<ranges.size(); i++) sum += dv[i];
#endif // HAVE_BOOST_BASED_THREAD_POOL_LIBRARY

// now define the function that is used for threading:

#define HAVE_BOOST_BASED_THREAD_POOL_LIBRARY

#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
#include "utils/split-indices.hh"
#include "utils/ctpl.h"

void density_for_atoms_multithread(int thread_index,
				   const atom_selection_container_t &asc,
				   const clipper::RTop<> &rtop_og,
				   const std::pair<unsigned int, unsigned int> &atom_index_range,
				   const clipper::NXmap<float> &nxmap,
				   float *dv,
				   std::atomic<unsigned int> &done_count_for_threads) {

   for (unsigned int i=atom_index_range.first; i<atom_index_range.second; i++) {
      mmdb::Atom *at = asc.atom_selection[i];
      clipper::Coord_orth pt = coot::co(at);
      clipper::Coord_map cm_try_2(rtop_og * pt);
      float dn = coot::util::density_at_point_by_cubic_interp(nxmap, cm_try_2);
      *dv += dn;
   }

   done_count_for_threads++;
}



// --------------------------------------------------------------------------------------

   int n_residue_atoms;
   mmdb::PPAtom residue_atoms;

   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      mmdb::Atom *at = residue_atoms[iat];
   }


         clipper::HKL_info::HKL_reference_index hri;
	 for (hri=fphidata.first(); !hri.last(); hri.next()) {
	    std::cout << " MTZ fphi: " << hri.hkl().h() << " "
		      << hri.hkl().k() << " " << hri.hkl().l() << " "
		      << fphidata[hri].f() << " "
		      << clipper::Util::rad2d(fphidata[hri].phi()) << std::endl;
	 }

// template code for history addition.
   std::string cmd = "";
   std::vector<coot::command_arg_t> args;
   args.push_back();
   add_to_history_typed(cmd, args);

   add_to_history_simple("");

//    if (atom_sel.n_selected_atoms > 0) {

      int n_models = atom_sel.mol->GetNumberOfModels();
      for (int imod=1; imod<=n_models; imod++) {

	 mmdb::Model *model_p = mol->GetModel(imod);
	 mmdb::Chain *chain_p;
	 // run over chains of the existing mol
	 int nchains = model_p->GetNumberOfChains();
	 if (nchains <= 0) {
	    std::cout << "bad nchains in molecule " << nchains
		      << std::endl;
	 } else {
	    for (int ichain=0; ichain<nchains; ichain++) {
	       chain_p = model_p->GetChain(ichain);
	       if (chain_p == NULL) {
		  // This should not be necessary. It seem to be a
		  // result of mmdb corruption elsewhere - possibly
		  // DeleteChain in update_molecule_to().
		  std::cout << "NULL chain in ... " << std::endl;
	       } else {
		  int nres = chain_p->GetNumberOfResidues();
		  mmdb::PResidue residue_p;
		  mmdb::Atom *at;
		  for (int ires=0; ires<nres; ires++) {
		     residue_p = chain_p->GetResidue(ires);
		     int n_atoms = residue_p->GetNumberOfAtoms();

		     for (int iat=0; iat<n_atoms; iat++) {
			at = residue_p->GetAtom(iat);



 // ----
void check_chiral_volumes(int imol) {
   graphics_info_t g;
   if (imol < graphics_info_t::n_molecules) {
      if (graphics_info_t::molecule[imol].has_model()) {
	 // my function here
      } else {
	 std::cout << "WARNING:: molecule " << imol
		   <<  " does not have coordinates\n";
      }
   } else {
      std::cout << "WARNING:: no such molecule " << imol << std::endl;
   }
}


// debug a mol

       {
	  int imod = 1;
	  mmdb::Model *model_p = flat_mol->GetModel(imod);
	  mmdb::Chain *chain_p;
	  int nchains = model_p->GetNumberOfChains();
	  for (int ichain=0; ichain<nchains; ichain++) {
	     chain_p = model_p->GetChain(ichain);
	     std::cout << "%%%%%%%%%% DEBUG chain :" << chain_p->GetChainID() << ":" << std::endl;
	     int nres = chain_p->GetNumberOfResidues();
	     mmdb::Residue *residue_p;
	     mmdb::Atom *atom_p;
	     for (int ires=0; ires<nres; ires++) {
		residue_p = chain_p->GetResidue(ires);
		std::cout << "%%%%%%%%%% DEBUG    residue number " << residue_p->GetSeqNum()
			  << std::endl;
		int n_atoms = residue_p->GetNumberOfAtoms();
		for (int iat=0; iat<n_atoms; iat++) {
		   atom_p = residue_p->GetAtom(iat);
		   std::cout << "%%%%%%%%%% DEBUG       atom :" << atom_p->GetAtomName() << ":"
			     << std::endl;
		}
	     }
	  }
       }


      // debug
      std::cout << "------------ molecule from residue selection ---- " << std::endl;
      int imod = 1;
      mmdb::Model *model_p = x->GetModel(imod);
      mmdb::Chain *chain_p;
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 mmdb::Residue *residue_p;
	 for (int ires=0; ires<nres; ires++) {
	    residue_p = chain_p->GetResidue(ires);
	    int n_atoms = residue_p->GetNumberOfAtoms();
	    std::cout << "   :" << chain_p->GetChainID() << ": " << residue_p->GetSeqNum()
		      << " :" << residue_p->GetInsCode() << ":    " << n_atoms << " atoms"
		      << std::endl;
	 }
      }


   // ---- simple version with protection

   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (! model_p) {
      std::cout << "Null model" << std::endl;
   } else {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         if (! chain_p) {
            std::cout << "Null chain" << std::endl;
         } else {
            int nres = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<nres; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (! residue_p) {
                  std::cout << "Null residue" << std::endl;
               } else {
                  int n_atoms = residue_p->GetNumberOfAtoms();
                  std::cout << "residue has " << n_atoms << " atoms " << std::endl;
                  for (int iat=0; iat<n_atoms; iat++) {
                     mmdb::Atom *at = residue_p->GetAtom(iat);
                     if (at)
                     std::cout << "   " << iat << " " << coot::atom_spec_t(at) << std::endl;
                  }
               }
            }
         }
      }
   }

   // ---- print mol

   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         std::cout << "   Chain " << chain_p->GetChainID() << std::endl;
         int nres = chain_p->GetNumberOfResidues();
         for (int ires=0; ires<nres; ires++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(ires);
            std::cout << "      " << residue_spec_t(residue_p) << std::endl;
            int n_atoms = residue_p->GetNumberOfAtoms();
            for (int iat=0; iat<n_atoms; iat++) {
               mmdb::Atom *at = residue_p->GetAtom(iat);
               std::cout << "       " << atom_spec_t(at) << std::endl;
            }
         }
      }
   }


   // ---- simple version


   // for(int imod = 1; imod<=asc.mol->GetNumberOfModels(); imod++) {
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int nres = chain_p->GetNumberOfResidues();
         for (int ires=0; ires<nres; ires++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(ires);
            int n_atoms = residue_p->GetNumberOfAtoms();
            for (int iat=0; iat<n_atoms; iat++) {
               mmdb::Atom *at = residue_p->GetAtom(iat);
