
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
	 torsion_info_t(const std::string &at_name_1_in,
			const std::string &at_name_2_in,
			const std::string &at_name_3_in,
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
	 std::string at_name_1;
	 std::string at_name_2;
	 std::string at_name_3;
	 double bond_length;
	 double angle_deg;
	 double torsion_deg;
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
			       mmdb::Residue *residue_p);
      void add_tetrahedral_hydrogen(const std::string &H_at_name,
				    const std::string &central_name,
				    const std::string &neighb_at_name_1,
				    const std::string &neighb_at_name_2,
				    const std::string &neighb_at_name_3,
				    double bond_length,
				    mmdb::Residue *residue_p);
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
      void add_guanidinium_hydrogens(mmdb::Residue *residue_p);
      void add_trp_indole_hydrogens(mmdb::Residue *residue_p);
      void add_trp_indole_hydrogen(const std::string &H_name,
				   const std::string &at_name_1,
				   const std::string &at_name_2,
				   const std::string &at_name_3,
				   double bl,
				   mmdb::Residue *residue_p);
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
      void add_SH_H(const std::string &H_name,
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
      void add_main_chain_H(mmdb::Residue *residue_p, mmdb::Residue *residue_prev_p);
      void add_main_chain_HA(mmdb::Residue *residue_p);
      mmdb::Atom *add_hydrogen_atom(std::string atom_name, clipper::Coord_orth &pos,
				    mmdb::realtype bf, mmdb::Residue *residue_p);
      // score hypotheses and convert to the best scoring one.
      void find_best_his_protonation_orientation(mmdb::Residue *residue_p);
      void delete_atom_by_name(const std::string &at_name, mmdb::Residue *residue_p);
      void hydrogen_placement_by_dictionary(mmdb::Residue *residue_p);
      void hydrogen_placement_by_dictionary(const dictionary_residue_restraints_t &rest,
					    mmdb::Residue *residue_p);
      void place_hydrogen_by_connected_atom_energy_type(unsigned int iat,
							unsigned int iat_neighb,
							const dictionary_residue_restraints_t &rest,
							mmdb::Residue *residue_p);
      void place_hydrogen_by_connected_2nd_neighbours(unsigned int iat,
						      unsigned int iat_neighb,
						      const dictionary_residue_restraints_t &rest,
						      mmdb::Residue *residue_p);
      
   public:
      reduce(mmdb::Manager *mol_in, int imol_in) {
	 mol = mol_in;
	 imol = imol_in;
      }
      void add_hydrogen_atoms(); // changes mol
      void delete_hydrogen_atoms();
      void add_geometry(protein_geometry *geom_p_in) { geom_p = geom_p_in; }
      // change HE2 to HD1 and vice versa
      void switch_his_protonation(mmdb::Residue *residue_p, mmdb::Atom *current_H_atom);
      
   };

}

#endif // REDUCE_HH
