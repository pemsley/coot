
#include <map>

#include "mmdb_manager.h"
#include "clipper/core/coords.h"
#include "coot-coord-utils.hh"
#include "protein-geometry.hh"
#include "lbg-shared.hh"


namespace coot { 

   class fle_residues_helper_t {
   public:
      bool is_set;
      clipper::Coord_orth centre;
      residue_spec_t spec;
      std::string residue_name;
      fle_residues_helper_t() { is_set = 0; }
      fle_residues_helper_t(const clipper::Coord_orth &pt,
			    const residue_spec_t &spec_in,
			    const std::string &res_name_in) {
	 centre = pt;
	 spec = spec_in;
	 residue_name = res_name_in;
	 is_set = 1;
      }
   };

   class pi_stacking_instance_t {
   public:

      // CATION_PI_STACKING sets ligand_cationic_atom_name, not the
      // ligand_ring_atom_names vector.
      //
      enum { PI_PI_STACKING,
	     PI_CATION_STACKING, // for cations on the protein residues (ligand pi)
	     CATION_PI_STACKING, // for cations on the ligand (protein TRY, PRO, TRP)
      };
      CResidue *res;
      int type; // pi-pi or pi-cation
      std::vector<std::string> ligand_ring_atom_names;
      float overlap_score; 
      std::string ligand_cationic_atom_name; // for cations on the ligand
      pi_stacking_instance_t(CResidue *res_in, int type_in,
			     const std::vector<std::string> &ring_atoms) {
	 res = res_in;
	 type = type_in;
	 ligand_ring_atom_names = ring_atoms;
      }
      
      // and the constructor for CATION_PI_STACKING
      // 
      pi_stacking_instance_t(CResidue *residue_in,
			     const std::string &ligand_atom_name_in) {
	 type = CATION_PI_STACKING;
	 res = residue_in;
	 ligand_cationic_atom_name = ligand_atom_name_in;
	 overlap_score = 0;
      }
      friend std::ostream& operator<< (std::ostream& s, const pi_stacking_instance_t &spec);
   };
   std::ostream& operator<< (std::ostream& s, const pi_stacking_instance_t &spec);

   class pi_stacking_container_t {
   private: 
      // can throw an exception
      std::pair<float, int> get_pi_overlap_to_ligand_ring(CResidue *res,
							  const clipper::Coord_orth &pt) const;

      float get_pi_overlap_to_ligand_cation(CResidue *res, const clipper::Coord_orth &pt) const;
      
      std::pair<clipper::Coord_orth, clipper::Coord_orth>
      get_ring_pi_centre_points(const std::vector<std::string> &ring_atom_names,
				CResidue *res_ref) const;
      
      // can throw an exception if not enough points found in pts.
      std::pair<clipper::Coord_orth, clipper::Coord_orth>
      ring_centre_and_normal(const std::vector<clipper::Coord_orth> &pts) const;

      // TRP has 2 rings, so we have to return a vector
      // 
      std::vector<std::vector<std::string> >
      ring_atom_names(const std::string &residue_name) const;
      
      float overlap_of_pi_spheres(const clipper::Coord_orth &pt1,
				  const clipper::Coord_orth &pt2,
				  const double &m1_pt_1, const double &m2_pt_1,
				  const double &m1_pt_2, const double &m2_pt_2) const;

      float
      overlap_of_cation_pi(const clipper::Coord_orth &ligand_pi_point,
			   const clipper::Coord_orth &cation_atom_point) const;

      std::vector<clipper::Coord_orth> get_cation_atom_positions(CResidue *res) const;
      // by search through res_ref
      std::vector<std::pair<std::string, clipper::Coord_orth> > get_ligand_cations(CResidue *res, const coot::dictionary_residue_restraints_t &monomer_restraints) const;

      std::vector<std::vector<std::string> >
      get_ligand_aromatic_ring_list(const dictionary_residue_restraints_t &monomer_restraints) const;
      
      
   public:
      // a vector of residues and types
      std::vector<pi_stacking_instance_t> stackings;
      pi_stacking_container_t (const dictionary_residue_restraints_t &monomer_restraints,
			       const std::vector<CResidue *> &filtered_residues,
			       CResidue *res_ref);
   };


   void write_solvent_accessibilities(const std::vector<std::pair<coot::atom_spec_t, float> > &sav,
				      CResidue *reference_residue);

   bool is_a_metal(CResidue *res);

   // The bonds from the protein to the ligand which contain
   // ligand-atom-name residue-spec and bond type (acceptor/donor).
   // These (ligand atom names) will have to be mapped to x y position
   // of the flat ligand.
   // 
   class fle_ligand_bond_t {
   public:
      enum { H_BOND_DONOR_MAINCHAIN,
	     H_BOND_DONOR_SIDECHAIN,
	     H_BOND_ACCEPTOR_MAINCHAIN, 
	     H_BOND_ACCEPTOR_SIDECHAIN,
	     METAL_CONTACT_BOND,
	     BOND_COVALENT,
	     BOND_OTHER };  // must sync this to lbg.hh
      std::string ligand_atom_name;
      int bond_type; // acceptor/donor
      residue_spec_t res_spec;
      double bond_length;  // from residue atom to ligand atom
      double water_protein_length; // if residue is a water, this is the closest
                                   // distance to protein (100 if very far).
      fle_ligand_bond_t(const std::string &name_in,
			int bond_type_in,
			double bl_in,
			residue_spec_t spec_in) {
	 ligand_atom_name = name_in;
	 bond_type = bond_type_in;
	 res_spec = spec_in;
	 bond_length = bl_in;
      }
      static int get_bond_type(CAtom *at_donor, CAtom *at_acceptor, bool ligand_atom_is_donor_flag) {
	 int r_bond_type = BOND_OTHER;

	 CAtom *ligand_atom = at_donor;
	 CAtom *residue_atom = at_acceptor;

	 if (! ligand_atom_is_donor_flag)
	    std::swap(ligand_atom, residue_atom);

	 if (is_a_metal(residue_atom->residue)) {
	    r_bond_type = METAL_CONTACT_BOND;
	 } else { 

	    if (ligand_atom_is_donor_flag) { 
	       if (coot::is_main_chain_p(residue_atom))
		  r_bond_type = H_BOND_ACCEPTOR_MAINCHAIN;
	       else
		  r_bond_type = H_BOND_ACCEPTOR_SIDECHAIN;
	       
	    } else {
	       if (coot::is_main_chain_p(residue_atom))
		  r_bond_type = H_BOND_DONOR_MAINCHAIN;
	       else
		  r_bond_type = H_BOND_DONOR_SIDECHAIN;
	    }
	 }
	 return r_bond_type;
      } 
   };
   
   std::map<std::string, std::string> make_flat_ligand_name_map(CResidue *flat_res);

   // return 100 if no other contact found (strange!)
   // 
   double find_water_protein_length(CResidue *ligand_residue, CMMDBManager *mol);


   std::vector<fle_ligand_bond_t> get_covalent_bonds(CMMDBManager *mol,
						     int SelHnd_lig,
						     int SelHnd_all,
						     const residue_spec_t &ligand_spec,
						     const protein_geometry &geom);

   std::vector<fle_ligand_bond_t> get_metal_bonds(CResidue *ligand_res,
						  const std::vector<CResidue *> &residues);


   std::vector<fle_ligand_bond_t> get_fle_ligand_bonds(CResidue *res_ref,
						       const std::vector<CResidue *> &residues,
						       const std::map<std::string, std::string> &name_map);
   // uses the coot::h_bond class (which uses the dictionary).
   // 
   std::vector<fle_ligand_bond_t> get_fle_ligand_bonds(CResidue *res_ref,
						       const std::vector<CResidue *> &residues,
						       CMMDBManager *mol,
						       const std::map<std::string, std::string> &name_map,
						       const protein_geometry &geom);
   
   void write_fle_centres(const std::vector<fle_residues_helper_t> &v,
			  const std::vector<coot::fle_ligand_bond_t> &bonds_to_ligand,
			  const std::vector<coot::solvent_exposure_difference_helper_t> &sed,
			  const pi_stacking_container_t &stacking,
			  CResidue *res_flat);

   bool standard_residue_name_p(const std::string &rn);

   enum { SP_HYBRIDIZATION, SP2_HYBRIDIZATION, SP3_HYBRIDIZATION };

   std::vector<std::pair<std::string, int> >
   get_prodrg_hybridizations(const std::string &prodrg_log_file_name);

   std::vector<std::pair<CAtom *, std::vector<clipper::Coord_orth> > >
   get_cannonball_vectors(CResidue *ligand_res_3d,
			  const coot::dictionary_residue_restraints_t &monomer_restraints);

   enum { H_IS_RIDING, H_IS_ROTATABLE }; // shared between named_torsion_t and flev_attached_hydrogens_t.

   
   // we need to map (the hydrogens torsions) between ideal prodrg
   // ligand atoms and the atoms in the residue/ligand of interest
   // (the reference ligand).
   class named_torsion_t {
   public:
      double torsion;
      double angle;
      double dist;
      std::string base_atom_name;
      std::string atom_name_2;
      std::string atom_name_bonded_to_H;
      int hydrogen_type;
      named_torsion_t(const std::string &base_name_in,
		      const std::string &a2,
		      const std::string &anbtoH,
		      double dist_in,
		      double angle_in,
		      double torsion_in,
		      int hydrogen_type_in) {
	 torsion = torsion_in;
	 angle = angle_in;
	 dist = dist_in;
	 base_atom_name = base_name_in;
	 atom_name_2 = a2;
	 atom_name_bonded_to_H = anbtoH;
	 hydrogen_type = hydrogen_type_in; // H_IS_ROTATABLE or H_IS_RIDING.
      } 
   };

   class flev_attached_hydrogens_t {
      std::vector<std::string> atoms_with_riding_hydrogens;
      std::vector<std::string> atoms_with_rotating_hydrogens;
      bool add_named_torsion(CAtom *h_at, CAtom *at,
			     const dictionary_residue_restraints_t &restraints,
			     CMMDBManager *mol,
			     int hydrogen_type); // fill named_torsions
      std::vector<std::pair<CAtom *, std::vector<clipper::Coord_orth> > >
      named_hydrogens_to_reference_ligand(CResidue *ligand_residue_3d,
					  const dictionary_residue_restraints_t &restraints) const;

      // Can throw an exception
      // 
      // Return the position of the H-ligand atom (the atom to which
      // the H is attached) and the hydrogen position - in that order.
      // 
      std::pair<clipper::Coord_orth, clipper::Coord_orth>
      hydrogen_pos(const coot::named_torsion_t &named_tor, CResidue *res) const;
      
      std::vector<CAtom *> close_atoms(const clipper::Coord_orth &pt,
				       const std::vector<CResidue *> &env_residues) const;

      bash_distance_t find_bash_distance(const clipper::Coord_orth &ligand_atom_pos,
					 const clipper::Coord_orth &hydrogen_pos,
					 const std::vector<CAtom *> &close_residue_atoms) const;
      double get_radius(const std::string &ele) const;
      
   public:
      flev_attached_hydrogens_t(const dictionary_residue_restraints_t &restraints);

      std::vector<named_torsion_t> named_torsions;
      
      // fill the named_torsions vector
      void cannonballs(CResidue *ligand_residue_3d,
		       const std::string &prodrg_3d_ligand_file_name,
		       const coot::dictionary_residue_restraints_t &restraints);

      // apply those cannonball direction onto the real reference ligand:
      void distances_to_protein(CResidue *residue_reference,
				CMMDBManager *mol_reference);

      std::map<std::string, std::vector<coot::bash_distance_t> > atom_bashes;
   
   };

   // the reference_residue is the flat residue
   void write_ligand_atom_accessibilities(const std::vector<std::pair<coot::atom_spec_t, float> > &sav,
					  const coot::flev_attached_hydrogens_t &attached_hydrogens,
					  CResidue *reference_residue);


}
