
#include <thread>

#include "utils/coot-utils.hh"
#include "utils/split-indices.hh"
#include "mmdb-extras.h"
#include "mmdb-crystal.h"
#include "Bond_lines.h"

#include "ligand/rotamer.hh" // do we have the directory hierachy correct?

std::vector<rotamer_markup_container_t>
Bond_lines_container::get_rotamer_dodecs(const atom_selection_container_t &asc) const {

   std::vector<rotamer_markup_container_t> dodecs;

#ifdef HAVE_CXX_THREAD
   unsigned int n_threads = coot::get_max_number_of_threads();

   // std::vector<coot::generic_display_object_t::dodec_t> dodecs;

   if (rotamer_probability_tables_p) {

      if (asc.mol) {
	 std::vector<std::pair<mmdb::Residue *, mmdb::Atom *> > residues;
	 int imod = 1;

	 mmdb::Model *model_p = asc.mol->GetModel(imod);
	 if (model_p) {
	    mmdb::Chain *chain_p;
	    int n_chains = model_p->GetNumberOfChains();
	    for (int ichain=0; ichain<n_chains; ichain++) {
	       chain_p = model_p->GetChain(ichain);
	       int nres = chain_p->GetNumberOfResidues();
	       for (int ires=0; ires<nres; ires++) {
		  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
		  if (residue_p) {
		     mmdb::Atom *atom_p = coot::util::intelligent_this_residue_mmdb_atom(residue_p);
		     if (atom_p) {
			std::pair<mmdb::Residue *, mmdb::Atom *> p(residue_p, atom_p);
			residues.push_back(p);
		     }
		  }
	       }
	    }
	 }

	 std::vector<std::vector<unsigned int> > splits;
	 coot::split_indices(&splits, residues.size(), n_threads);
	 // now each thread has its own splits vector
	 dodecs.resize(residues.size());

	 std::vector<std::thread> threads;
	 for (std::size_t i=0; i<n_threads; i++) {
	    // each thread partially fills dodecs
	    threads.push_back(std::thread(add_rotamer_markups, std::cref(splits[i]), std::cref(residues),
					  rotamer_probability_tables_p, &dodecs));
	 }
	 for (std::size_t i=0; i<n_threads; i++)
	    threads[i].join();
      }
   }
#endif // HAVE_CXX_THREAD
   return dodecs;
}

// partially fill dodecs
// static
void
Bond_lines_container::add_rotamer_markups(const std::vector<unsigned int> &indices,
					  const std::vector<std::pair<mmdb::Residue *, mmdb::Atom *> > &residues,
					  coot::rotamer_probability_tables *rpt,
					  std::vector<rotamer_markup_container_t> *dodecs) {
   for (std::size_t i=0; i<indices.size(); i++) {
      const unsigned int &idx = indices[i];
      rotamer_markup_container_t rmc = get_rotamer_probability(residues[idx], rpt);
      dodecs->at(idx) = rmc;
   }
}

// static
rotamer_markup_container_t
Bond_lines_container::get_rotamer_probability(const std::pair<mmdb::Residue *, mmdb::Atom *> &ra,
					      coot::rotamer_probability_tables *rpt) {

   rotamer_markup_container_t rmc;
   coot::residue_spec_t res_spec(ra.first);
   rmc.spec = res_spec;

   // old: integer probabilities
   // coot::rotamer rot(residue_p);
   // coot::rotamer_probability_info_t pr = rot.probability_of_this_rotamer();

   mmdb::Residue *residue_p = ra.first;
   std::string res_name(residue_p->GetResName());
   if (coot::util::is_standard_amino_acid_name(res_name)) {

      try {

	 std::vector<coot::rotamer_probability_info_t> pr_v = rpt->probability_this_rotamer(residue_p);

	 if (pr_v.size() > 0) {
	    const coot::rotamer_probability_info_t &pr = pr_v[0]; // hack

	    if (pr.state != coot::rotamer_probability_info_t::RESIDUE_IS_GLY_OR_ALA) {
	       // OK or MISSING_ATOMS or ROTAMER_NOT_FOUND
	       clipper::Coord_orth pos = coot::co(ra.second);
	       double z = 0;
               bool use_deuteranomaly_mode = false;
	       coot::colour_holder col(z, 0.0, 1.0, use_deuteranomaly_mode, std::string(""));
	       if (pr.state == coot::rotamer_probability_info_t::OK) {

                  if (false)
                     std::cout << "in get_rotamer_probability() OK "
                               << res_spec << " " << pr.probability << std::endl;

                  // pr should be between 0 and 100.
                  //
		  // pr is high, z low, -> green
		  // pr is ~0, z is ~1 -> red
                  //
		  z = 1.0 - sqrt(pr.probability*0.01);

		  // args fraction, min, max, dummy-not-colour-triple-flag
		  col = coot::colour_holder(z, 0.0, 1.0, use_deuteranomaly_mode, std::string(""));
                  col.brighten(0.1);
	       }
	       if (pr.state == coot::rotamer_probability_info_t::MISSING_ATOMS)
		  col = coot::colour_holder("#bb22bb"); // purple
	       if (pr.state == coot::rotamer_probability_info_t::ROTAMER_NOT_FOUND)
		  col = coot::colour_holder("#999999"); // grey

	       if (false)
		  std::cout << "debug:: atom ra.second " << coot::atom_spec_t(ra.second)
			    << " hal pr " << pr.probability
			    << " has col " << col << std::endl;

	       rmc = rotamer_markup_container_t(res_spec, pos, col, pr);
	    }
	 }
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "exception caught in get_rotamer_dodecs() " << std::endl;
	 std::cout << "    " << rte.what() << std::endl;
      }
   }

   return rmc;
}
