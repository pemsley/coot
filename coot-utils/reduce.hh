
#ifndef REDUCE_HH
#define REDUCE_HH

#include "clipper/core/coords.h"
#include "mmdb2/mmdb_manager.h"

#include "geometry/protein-geometry.hh" // needed to get H-bond types when checking HIS protonatino
                                        // (and spin searching OHs).

namespace coot {

   class reduce {

      class torsion_info_t {
      public:
	 torsion_info_t(const std::string &at_name_1_in, // CA 
			const std::string &at_name_2_in, // CB
			const std::string &at_name_3_in, // CG1 (i.e. other way round)
			double bond_length_in,
			double angle_deg_in,
			double torsion_deg_in) {
	    bond_length = bond_length_in;
	    at_name_1 = at_name_1_in;
	    at_name_2 = at_name_2_in;
	    at_name_3 = at_name_3_in;
	    angle_deg = angle_deg_in;
	    torsion_deg = torsion_deg_in;
	 }
	 // can throw a std::runtime_error
	 torsion_info_t(const std::string &first_neighb,
			const std::vector<std::string> &second_neighb_vec,
			const std::map<std::string, std::vector<std::string> > &third_neighbour_map,
			double bl,
			double angle_in,
			double torsion_in) {
	    if (second_neighb_vec.size() > 0) {
	       std::string second = second_neighb_vec[0];
	       std::map<std::string, std::vector<std::string> >::const_iterator it;
	       it = third_neighbour_map.find(second);
	       if (it != third_neighbour_map.end()) {
		  std::vector<std::string> thirds = it->second;
		  if (thirds.size() > 0) {
		     const std::string &third = thirds[0];
		     at_name_1 = third;
		     at_name_2 = second;
		     at_name_3 = first_neighb;
		     angle_deg = angle_in;
		     torsion_deg = torsion_in;
		     bond_length = bl;
		  }
	       }
	    }
	 }
	 std::string at_name_1; // CA (not sure if this atom ordering is sensible)
	 std::string at_name_2; // CB
	 std::string at_name_3; // CG1
	 double bond_length;
	 double angle_deg;
	 double torsion_deg;
      };

      class atom_with_attached_Hs {
      public:
	 enum hydrogen_t { METHYL, HYDROXYL, SULFHYDRYL };

	 atom_with_attached_Hs(mmdb::Atom *at_in, hydrogen_t type_in,
			       const std::vector<mmdb::Atom *> &hydrogens_in) :
	    at(at_in), type(type_in), hydrogens(hydrogens_in) {}
	 mmdb::Atom* at;
	 hydrogen_t type;
	 std::vector<mmdb::Atom *> hydrogens;
      };

      class atoms_with_spinnable_Hs {
      public:

	 std::map<std::string, std::vector<atom_with_attached_Hs> > typed_atoms;
	 std::vector<std::vector<atom_with_attached_Hs> > cliques;

      public:
	 atoms_with_spinnable_Hs() { typed_atoms[""].reserve(200); }

	 void add(mmdb::Atom* at, atom_with_attached_Hs::hydrogen_t type,
		  const std::vector<mmdb::Atom *> &attahed_hydrogen_atoms);
	 void add(mmdb::Atom* at, atom_with_attached_Hs::hydrogen_t type,
		  mmdb::Atom *attahed_hydrogen_atom);
	 void cliquize();
	 std::vector<std::vector<atom_with_attached_Hs> > get_cliques() const {
	    return cliques;
	 }
	 void resolve_clashing_cliques() {
	    for (std::size_t i=0; i<cliques.size(); i++)
	       resolve_clashing_clique(cliques[i]);
	 }
	 // Move the H atoms of the clique if there are classes.
	 // (clashes may still remain after having done so, of course).
	 void resolve_clashing_clique(const std::vector<atom_with_attached_Hs> &clique);
      };

      clipper::Coord_orth position_by_bond_length_angle_torsion(mmdb::Atom *at_1,  // CA
								mmdb::Atom *at_2,  // CB
								mmdb::Atom *at_3,  // CG
								double bl,
								double angle_rad,
								double torsion_rad) const;
      clipper::Coord_orth position_by_bisection(mmdb::Atom *at_1, // for Hs on PHE etc
						mmdb::Atom *at_2,
						mmdb::Atom *at_3,
						double bl) const;
      std::pair<clipper::Coord_orth, clipper::Coord_orth>
      position_pair_by_bisection(mmdb::Atom *at_1,  // CA
				 mmdb::Atom *at_2,  // CB
				 mmdb::Atom *at_3,  // CG
				 double bl,
				 double angle_separation // how far apart are the two H atoms?
				 ) const;
      clipper::Coord_orth position_by_tetrahedron(mmdb::Atom *at_central,
						  mmdb::Atom *at_n_1,
						  mmdb::Atom *at_n_2,
						  mmdb::Atom *at_n_3,
						  double bl) const;
      void add_methyl_Hs(const std::string &at_name_1,
			 const std::string &at_name_2,
			 const std::string &at_name_3,
			 torsion_info_t torsion_1,
			 mmdb::Residue *residue_p);
      void add_methyl_Hs(const std::string &at_name_1,
			 const std::string &at_name_2,
			 const std::string &at_name_3,
			 torsion_info_t torsion_1, torsion_info_t torsion_2,
			 mmdb::Residue *residue_p);
      void add_2_sp3_hydrogens(const std::string &H_at_name_1,
			       const std::string &H_at_name_2,
			       const std::string &at_name_1,
			       const std::string &at_name_2,
			       const std::string &at_name_3,
			       double bond_length,
			       double angle_between_Hs, // in degrees
			       mmdb::Residue *residue_p,
			       bool choose_only_farthest_position=false);
      void add_2_sp3_hydrogens(const std::string &H_at_name_1,
			       const std::string &H_at_name_2,
			       const std::string &first_neighb,
			       const std::vector<std::string> &second_neighb_vec,
			       double bond_length,
			       double angle_between_Hs, // in degrees
			       mmdb::Residue *residue_p,
			       bool choose_only_farthest_position=false);
      void add_tetrahedral_hydrogen(const std::string &H_at_name,
				    const std::string &central_name,
				    const std::string &neighb_at_name_1,
				    const std::string &neighb_at_name_2,
				    const std::string &neighb_at_name_3,
				    double bond_length,
				    mmdb::Residue *residue_p);
      void add_tetrahedral_hydrogen(const std::string &H_at_name,
				    const std::string &first_neighb,
				    const std::vector<std::string> &second_neighb_vec,
				    double bl, mmdb::Residue *residue_p);
      
      void add_aromatic_hydrogen(const std::string &H_at_name,
				 const std::string &neighb_at_name_1,
				 const std::string &neighb_at_name_2, // add to this
				 const std::string &neighb_at_name_3,
				 double bl, mmdb::Residue *residue_p);
      void add_amino_hydrogens(const std::string &H_at_name_1,
			       const std::string &H_at_name_2,
			       const std::string &at_name_1,
			       const std::string &at_name_2,
			       const std::string &at_name_3,
			       double bl_amino, // angle is 120, torsions are 180 and 0
			       mmdb::Residue *residue_p);
      void add_amino_hydrogens(const std::string &H_at_name_1,
			       const std::string &H_at_name_2,
			       const std::string &first_neighb,
			       const std::vector<std::string> &second_neighb_vec,
			       const std::map<std::string, std::vector<std::string> > &third_neighb_vec,
			       double bl_amino,
			       mmdb::Residue *residue_p);

      void add_guanidinium_hydrogens(mmdb::Residue *residue_p);
      void add_trp_indole_hydrogens(mmdb::Residue *residue_p);
      void add_trp_indole_hydrogen(const std::string &H_name,
				   const std::string &at_name_1,
				   const std::string &at_name_2,
				   const std::string &at_name_3,
				   double bl,
				   mmdb::Residue *residue_p);

      atoms_with_spinnable_Hs spinables;

      // this will need a spin-search
      void add_OH_H(const std::string &H_name,
		    const std::string &at_name_1,
		    const std::string &at_name_2,
		    const std::string &at_name_3,
		    double bl,
		    double angle,      // deg
		    double tor_inital, // deg
		    mmdb::Residue *residue_p);
      // this will need a spin-search
      std::vector<mmdb::Atom *> add_SH_H(const std::string &H_name,
					 const std::string &at_name_1,
					 const std::string &at_name_2,
					 const std::string &at_name_3,
					 double bl,
					 double angle,      // deg
					 double tor_inital, // deg
					 mmdb::Residue *residue_p);
      // both of the above wrap this:
      std::vector<mmdb::Atom *> add_xH_H(const std::string &H_name,
					 const std::string &at_name_1,
					 const std::string &at_name_2,
					 const std::string &at_name_3,
					 double bl,
					 double angle,      // deg
					 double tor_inital, // deg
					 mmdb::Residue *residue_p);

      void add_OH_H(const std::string &H_at_name,
		    const std::string &first_neighb,
		    const std::vector<std::string> &second_neighb_vec,
		    const std::map<std::string, std::vector<std::string> > &third_neighb_vec,
		    double bond_length,
		    double ang_deg,
		    double torsion_deg,
		    mmdb::Residue *residue_p);
      
      
      void add_his_ring_C_Hs(mmdb::Residue *residue_p);
      std::vector<mmdb::Atom *> add_his_ring_H(const std::string &H_name,
					       const std::string &at_name_1,
					       const std::string &at_name_2,
					       const std::string &at_name_3,
					       double bl,
					       mmdb::Residue *residue_p);

      void add_his_ring_H(const std::string &H_at_name,
			  const std::string &first_neigh,
			  const std::vector<std::string> second_neighb_vec,
			  double bl,
			  mmdb::Residue *residue_p);

      void add_aromatic_hydrogen(const std::string &H_at_name,
				 const std::string &first_neigh,
				 const std::vector<std::string> second_neighb_vec,
				 double bl,
				 mmdb::Residue *residue_p);

      mmdb::Manager *mol;
      int imol; // for dictionary lookups.
      protein_geometry *geom_p;
      void add_riding_hydrogens(); // non-spin-search
      bool add_riding_hydrogens(mmdb::Residue *residue_p, mmdb::Residue *residue_prev_p);
      void add_main_chain_hydrogens(mmdb::Residue *residue_p, mmdb::Residue *residue_prev_p,
				    bool is_gly=false);
      void add_main_chain_HA(mmdb::Residue *residue_p);
      void add_main_chain_H(mmdb::Residue *residue_p, mmdb::Residue *residue_prev_p);
      // like above but for ligands (both second neighbours come from ligand residue)
      void add_amino_single_H(const std::string &H_at_name,
			      const std::string &first_neighb,
			      const std::vector<std::string> &second_neighb_vec,
			      double bl,
			      mmdb::Residue *residue_p);
      // add H to second atom by bisection
      void add_amino_single_H(const std::string H_at_name,
			      const std::string &at_name_1,
			      const std::string &at_name_2,
			      const std::string &at_name_3,
			      double bl,
			      mmdb::Residue *residue_p);

      mmdb::Atom *add_hydrogen_atom(std::string atom_name, clipper::Coord_orth &pos,
				    mmdb::realtype bf,
				    const std::string &altconf,
				    mmdb::Residue *residue_p);
      // score hypotheses and convert to the best scoring one.
      void find_best_his_protonation_orientation(mmdb::Residue *residue_p);
      void delete_atom_by_name(const std::string &at_name, mmdb::Residue *residue_p);
      void hydrogen_placement_by_dictionary(mmdb::Residue *residue_p);
      void hydrogen_placement_by_dictionary(const dictionary_residue_restraints_t &rest,
					    mmdb::Residue *residue_p);
      // return a list of names of placed atoms
      std::vector<std::string>
      place_hydrogen_by_connected_atom_energy_type(unsigned int iat,
						   unsigned int iat_neighb,
						   const dictionary_residue_restraints_t &rest,
						   mmdb::Residue *residue_p);
      // which calls
      std::vector<std::string>
      place_hydrogen_by_connected_atom_energy_type(const std::string &energy_type,
						   unsigned int iat,
						   unsigned int iat_neighb,
						   const dictionary_residue_restraints_t &rest,
						   mmdb::Residue *residue_p);
      void place_hydrogen_by_connected_2nd_neighbours(unsigned int iat,
						      unsigned int iat_neighb,
						      const dictionary_residue_restraints_t &rest,
						      mmdb::Residue *residue_p);

      // atoms that are connected to the second neighbour that aren't the first neighbour.
      //
      std::map<std::string, std::vector<std::string> > 
      third_neighbour_map(const std::string &first_neighb,
			  const std::vector<std::string> second_neighb_vec,
			  const dictionary_residue_restraints_t &rest) const;

      std::string get_other_H_name(const std::string &first_neighb,
				   const std::string &H_at_name,
				   const dictionary_residue_restraints_t &dict) const;
      std::vector<std::string> get_other_H_names(const std::string &first_neighb,
						 const std::string &H_at_name,
						 const dictionary_residue_restraints_t &dict) const;

      bool is_ss_bonded(mmdb::Residue *residue_p) const;
      // No SSBOND interogation yet, so...
      bool hack_ss_bond_test(mmdb::Residue *residue_p, mmdb::Model *model_p) const;

      // as a rule, we don't want to add hydrogen to atoms that are linked
      //
      bool is_linked(const std::string &atom_name, mmdb::Residue *residue_p) const;

      bool verbose_output;

   public:
      reduce(mmdb::Manager *mol_in, int imol_in) {
	 mol = mol_in;
	 imol = imol_in;
         verbose_output = true;
      }
      void add_hydrogen_atoms(); // changes mol
      void delete_hydrogen_atoms();
      void add_geometry(protein_geometry *geom_p_in) { geom_p = geom_p_in; }
      // change HE2 to HD1 and vice versa
      void switch_his_protonation(mmdb::Residue *residue_p, mmdb::Atom *current_H_atom);
      void set_verbose_output(bool flag) { verbose_output = flag; }
      
   };

}

#endif // REDUCE_HH
