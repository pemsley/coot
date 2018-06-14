
#include "mmdb-extras.h"
#include "mmdb-crystal.h"
#include "Bond_lines.h"

#include "ligand/rotamer.hh" // do we have the directory hierachy correct?

std::vector<rotamer_markup_container_t>
Bond_lines_container::get_rotamer_dodecs(const atom_selection_container_t &asc) const {

   // std::vector<coot::generic_display_object_t::dodec_t> dodecs;
   std::vector<rotamer_markup_container_t> dodecs;

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

	 for (unsigned int i=0; i<residues.size(); i++) {

	    mmdb::Residue *residue_p = residues[i].first;

	    // old: integer probabilities
	    // coot::rotamer rot(residue_p);
	    // coot::rotamer_probability_info_t pr = rot.probability_of_this_rotamer();

	    std::string res_name(residue_p->GetResName());

	    if (is_standard_amino_acid_name(res_name)) {

	       try {

		  std::vector<coot::rotamer_probability_info_t> pr_v =
		     rotamer_probability_tables_p->probability_this_rotamer(residue_p);

		  if (pr_v.size() > 0) {
		     const coot::rotamer_probability_info_t &pr = pr_v[0]; // hack

		     std::cout << "   " << coot::residue_spec_t(residue_p) << " " << pr << std::endl;
		     if (pr.state != coot::rotamer_probability_info_t::RESIDUE_IS_GLY_OR_ALA) {
			// OK or MISSING_ATOMS or ROTAMER_NOT_FOUND
			double size = 0.5;
			clipper::Coord_orth pos = coot::co(residues[i].second);
			double z = 0;
			coot::colour_holder col(z, 0.0, 1.0, std::string(""));
			if (pr.state == coot::rotamer_probability_info_t::OK) {
			   if (false)
			      std::cout << coot::residue_spec_t(residues[i].first) << " " << pr << "\n";
			   z = sqrt(sqrt(pr.probability*0.01));
			   col = coot::colour_holder(z, 0.0, 1.0, std::string(""));
			}
			if (pr.state == coot::rotamer_probability_info_t::MISSING_ATOMS)
			   col = coot::colour_holder("#bb22bb"); // purple
			if (pr.state == coot::rotamer_probability_info_t::ROTAMER_NOT_FOUND)
			   col = coot::colour_holder("#22eeee"); // cyan

			rotamer_markup_container_t rc(pos, col);
			dodecs.push_back(rc);
		     }
		  }
	       }
	       catch (const std::runtime_error &rte) {
		  std::cout << "exception caught in get_rotamer_dodecs() " << std::endl;
		  std::cout << "    " << rte.what() << std::endl;
	       }
	    }
	 }
      }
   }
   return dodecs;
}
