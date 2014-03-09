
#ifndef FLEV_ANNOTATIONS_HH
#define FLEV_ANNOTATIONS_HH

#include "geometry/protein-geometry.hh"
#include "lidia-core/lbg-shared.hh" // bash_distance_t


namespace coot {
      
   // we need to map (the hydrogens torsions) between ideal prodrg
   // ligand atoms and the atoms in the residue/ligand of interest
   // (the reference ligand).
   // 
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
      // the "base" (heavy) atom name in first and the H name in second.
      std::vector<std::pair<std::string, std::string> > atoms_with_riding_hydrogens;
      std::vector<std::pair<std::string, std::string> > atoms_with_rotating_hydrogens;
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

      // find an atom (the atom, perhaps) bonded to lig_at that is not H_at.
      // Return its position. Can throw a std::runtime_error if not found.
      // 
      clipper::Coord_orth get_atom_pos_bonded_to_atom(CAtom *lig_at, CAtom *H_at, // not H_at
						      CResidue *ligand_residue,
						      const protein_geometry &geom) const;
      
      
   public:
      flev_attached_hydrogens_t(const dictionary_residue_restraints_t &restraints);

      std::vector<named_torsion_t> named_torsions;
      
      // fill the named_torsions vector, a trivial wrapper to the below function
      void cannonballs(CResidue *ligand_residue_3d,
		       const std::string &prodrg_3d_ligand_file_name,
		       const coot::dictionary_residue_restraints_t &restraints);

      // fill the named_torsions vector
      void cannonballs(CResidue *ligand_residue_3d,
		       CMMDBManager *mol,
		       const coot::dictionary_residue_restraints_t &restraints);
      
      // apply those cannonball direction onto the real reference ligand:
      void distances_to_protein(CResidue *residue_reference,
				CMMDBManager *mol_reference);
      void distances_to_protein_using_correct_Hs(CResidue *residue_reference,
						 CMMDBManager *mol_reference,
						 const protein_geometry &geom);

      std::map<std::string, std::vector<coot::bash_distance_t> > atom_bashes;
   
   };

   
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
   std::ostream& operator<<(std::ostream &s, fle_residues_helper_t fler);


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
      atom_spec_t ligand_atom_spec;
      int bond_type; // acceptor/donor

      residue_spec_t res_spec;
      atom_spec_t interacting_residue_atom_spec; // contains res_spec obviously.
      
      double bond_length;  // from residue atom to ligand atom
      double water_protein_length; // if residue is a water, this is the closest
                                   // distance to protein (100 if very far).
      fle_ligand_bond_t(const atom_spec_t &ligand_atom_spec_in,
			const atom_spec_t &interacting_residue_atom_spec_in,
			int bond_type_in,
			double bl_in) {
	 ligand_atom_spec = ligand_atom_spec_in;
	 interacting_residue_atom_spec = interacting_residue_atom_spec_in;
	 res_spec = residue_spec_t(interacting_residue_atom_spec_in);
	 bond_type = bond_type_in;
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
      bool operator==(const fle_ligand_bond_t &in) const {
	 bool status = false;
	 if (in.bond_type == bond_type) {
	    if (in.ligand_atom_spec == ligand_atom_spec) {
	       if (in.interacting_residue_atom_spec == interacting_residue_atom_spec) {
		  status = true;
	       }
	    }
	 }
	 return status;
      }
   };
   std::ostream& operator<<(std::ostream &s, fle_ligand_bond_t flb);

}


#endif // FLEV_ANNOTATIONS_HH
