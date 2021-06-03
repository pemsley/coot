/* ligand/wligand.hh
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006 by The University of York
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

#ifndef WLIGAND_HH
#define WLIGAND_HH

#include "monomer-utils.hh"
#include "ligand.hh"


namespace coot {

   class torsioned_atoms_info_t {
   public:
      torsioned_atoms_info_t(dict_torsion_restraint_t rest, double tors_in) {
	 quad = atom_name_quad(rest.atom_id_1_4c(),
			       rest.atom_id_2_4c(),
			       rest.atom_id_3_4c(),
			       rest.atom_id_4_4c());
	 torsion_from_restraint = rest.angle();
	 esd_from_restraint = rest.esd();
	 period_from_restraint = rest.periodicity();
	 torsion = tors_in;
      } 
      atom_name_quad quad;
      double torsion_from_restraint;
      double esd_from_restraint;
      int period_from_restraint;
      double torsion; // torsion fed into mgtree class
   };
      
      
   class installed_wiggly_ligand_info_t {
   public:
      minimol::molecule mol;
      std::vector<torsioned_atoms_info_t> torsioned_atoms;
      void add_torsion(const dict_torsion_restraint_t &rest, double torsion);
      void add_torsions(const std::vector<dict_torsion_restraint_t> &rests,
			const std::vector<float> &torsions);
      unsigned int n_torsions() const;
      // throw an exception if not possible
      std::pair<float, float> get_set_and_ideal_torsions(int i) const;
      // throw an exception if not possible
      std::pair<float, float> get_set_and_real_torsions(int i) const;
   };

   class wligand : public ligand, public monomer_utils {

      // monomer_type, residue_type, comp_id are synonymous.
      // 


       std::vector<atom_name_pair> get_torsion_bonds_atom_pairs(const std::string &monomer_type,
 								     const std::vector <dict_torsion_restraint_t> &monomer_torsions) const;
       std::vector<atom_name_quad> get_torsion_bonds_atom_quads(const std::string &monomer_type,
 								     const std::vector <dict_torsion_restraint_t> &monomer_torsions) const;

      
      std::vector<atom_index_pair>
      get_atom_index_pairs(std::vector<atom_name_pair>atom_name_pairs,
			   const minimol::molecule &ligand) const;

      std::vector<atom_index_quad>
      get_atom_index_quads(const std::vector<atom_name_quad> &atom_name_quads,
			   const minimol::molecule &ligand) const;

      std::vector<atom_index_pair>
      get_atom_index_pairs(std::vector<atom_name_pair> atom_name_pairs,
			   const minimol::molecule &ligand,
			   int imono) const;

      std::vector<std::vector<int> > getcontacts(const minimol::molecule &mol) const;
      std::vector<std::vector<int> > getcontacts(const minimol::molecule &mol,
						 const dictionary_residue_restraints_t &dict) const;


      std::string get_monomer_type_from_mol(const minimol::molecule &mol) const;

      std::vector <float>
      get_torsions_by_random(const std::vector <dict_torsion_restraint_t> &m_torsions) const;

      float probability_of_torsions(const std::vector <dict_torsion_restraint_t> &m_torsions,
				    const std::vector <float> &r) const;

      std::vector<float> cumulative;
      float cumulative_step;
      bool debug_wiggly_ligands; 
      void setup_normal_cumulative_table() {

	 // set up cumulative 
	 cumulative_step=0.01;
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

      bool is_unique_conformer(const coot::minimol::molecule &mol) const;

      static installed_wiggly_ligand_info_t
      optimize(const coot::minimol::residue &wiggled_ligand_residue,
               const coot::protein_geometry &pg,
               const std::vector <dict_torsion_restraint_t> &non_const_torsions,
               const std::vector<float> &torsion_set,
               const std::string &ligand_chain_id,
               int isample);

      // can throw a std::runtime_error
      // (if the dictionary/atoms are not found)
      //
      installed_wiggly_ligand_info_t
      optimize_and_install_if_unique(const minimol::residue &wiggled_ligand_residue,
				     const coot::protein_geometry &pg,
				     const std::vector <dict_torsion_restraint_t> &non_const_torsions, 
				     const std::vector<float> &torsion_set,
				     const std::string &chain_id,
				     int isample,
				     bool optimize_geometry_flag,
				     bool fill_returned_molecules_vector_flag);
      
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
      std::vector<installed_wiggly_ligand_info_t>
      install_simple_wiggly_ligands(protein_geometry *pg,
				    const minimol::molecule &ligand,
				    int imol_ligand,
				    int n_samples,
				    bool optimize_geometry_flag,
				    bool fill_returned_molecules_vector_flag);

      // install one by one for dialog updating
      installed_wiggly_ligand_info_t
      install_simple_wiggly_ligand(protein_geometry *pg,
				   const minimol::molecule &ligand,
				   int imol_ligand,
				   int isample,
				   bool optimize_geometry_flag);

      short int install_linked_wiggly_ligands(const protein_geometry &pg,
					      const minimol::molecule &ligand,
					      int n_samples);

      void set_debug_wiggly_ligands() { debug_wiggly_ligands = 1; }
      std::vector<minimol::molecule> get_conformers() const {
	 return initial_ligand;
      }
   }; 

}

#endif // WLIGAND_HH
