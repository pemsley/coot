/* ligand/wligand.hh
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006 by Paul Emsley, The University of York
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

#include "monomer-utils.hh"
#include "ligand.hh"


namespace coot {

   class wligand : public ligand, public monomer_utils {

      // monomer_type, residue_type, comp_id are synonymous.
      // 


       std::vector<coot::atom_name_pair> get_torsion_bonds_atom_pairs(const std::string &monomer_type,
 								     const std::vector <coot::dict_torsion_restraint_t> &monomer_torsions) const;

      
      std::vector<coot::atom_index_pair>
      get_atom_index_pairs(std::vector<coot::atom_name_pair>atom_name_pairs,
			   const coot::minimol::molecule &ligand) const;

      std::vector<coot::atom_index_pair> get_atom_index_pairs(std::vector<coot::atom_name_pair> atom_name_pairs,
							      const coot::minimol::molecule &ligand,
							      int imono) const;

      std::vector<std::vector<int> > getcontacts(const coot::minimol::molecule &mol) const;
      std::vector<std::vector<int> > getcontacts(const coot::minimol::molecule &mol,
						 const coot::dictionary_residue_restraints_t &dict) const;


      std::string get_monomer_type_from_mol(const coot::minimol::molecule &mol) const;

      std::vector <float>
      get_torsions_by_random(const std::vector <coot::dict_torsion_restraint_t> &m_torsions) const;

      float probability_of_torsions(const std::vector <coot::dict_torsion_restraint_t> &m_torsions,
				    const std::vector <float> &r) const;

      std::vector<float> cumulative;
      float cumulative_step;
      bool debug_wiggly_ligands; 
      void setup_normal_cumulative_table() {

	 // set up cumulative 
	 float cumulative_step=0.01;
	 float sigma_lim = 4.0; // from -4.0 to +4.0 sigma.
	 float sum = 0.0;
	 for (float v= -sigma_lim; v<sigma_lim; v+=cumulative_step) {
	    double f = (1.0/sqrt(2.0*M_PI)) * exp(-0.5*v*v);
	    sum += f;
	    cumulative.push_back(sum);
	 }

	 // This looks fine.
	 //for (unsigned int i=1; i<cumulative.size(); i++) {
	 // std::cout << " cumulative: " << i << " " << cumulative[i] << std::endl;
	 //  std::cout << " cumulative: " << i << " "
	 //	      << cumulative[i] - cumulative[i-1] << std::endl;
	 //}
      }
      
   public:

      wligand() {
	 cumulative_step=0.01;
	 setup_normal_cumulative_table();
	 debug_wiggly_ligands = 0;
      }

      // This should be a util function, I think.
      double get_random_normal_value() const;

      // There are 2 types of installation of wiggly ligands.
      // 
      // This is because I had a re-think about how to implement
      // wiggly-ligands.  At the first attempt to write wiggly ligands
      // with linked residues, I realised it was very complex
      // (first-go).
      //
      // I realised later that we have a situation similar to the
      // refinement, i.e. we have torsions between arbitary atom
      // indices.  So when I re-write linked-wiggly-ligands, I will
      // use the restraints_container_t and methods as a model.  Maybe
      // I will even use mmdb format, rather than mini-mol.  Perhaps I
      // will derive from restraints_container_t or whatever is
      // appropriate.
      //
      // For now, (20030707) I will write
      // install_simple_wiggly_ligands using a minimol.  The ligand in
      // this case is only one residue, and the function fails
      // (returns 0) if it is not.
      // 

      // update pg with try_dynamic_add if needed.
      // 
      
      // Throw an exception if there is failure to install the ligands.
      // [used to return std::pair<short int, std::string>]
      //
      std::vector<minimol::molecule>
      install_simple_wiggly_ligands(coot::protein_geometry *pg,
				    const coot::minimol::molecule &ligand,
				    int n_samples,
				    bool optimize_geometry_flag,
				    bool fill_returned_molecules_vector_flag);

      short int install_linked_wiggly_ligands(const coot::protein_geometry &pg,
					      const coot::minimol::molecule &ligand,
					      int n_samples);

      void set_debug_wiggly_ligands() { debug_wiggly_ligands = 1; } 
   }; 

}
