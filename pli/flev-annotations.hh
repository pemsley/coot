/* pli/flev-annotations.hh
 * 
 * Copyright 2010 by The University of Oxford
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

#ifndef FLEV_ANNOTATIONS_HH
#define FLEV_ANNOTATIONS_HH

#include "geometry/protein-geometry.hh"
#include "geometry/main-chain.hh"
#include "geometry/residue-and-atom-specs.hh"
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
      bool add_named_torsion(mmdb::Atom *h_at, mmdb::Atom *at,
			     const dictionary_residue_restraints_t &restraints,
			     mmdb::Manager *mol,
			     int hydrogen_type); // fill named_torsions
      std::vector<std::pair<mmdb::Atom *, std::vector<clipper::Coord_orth> > >
      named_hydrogens_to_reference_ligand(mmdb::Residue *ligand_residue_3d,
					  const dictionary_residue_restraints_t &restraints) const;

      // Can throw an exception
      // 
      // Return the position of the H-ligand atom (the atom to which
      // the H is attached) and the hydrogen position - in that order.
      // 
      std::pair<clipper::Coord_orth, clipper::Coord_orth>
      hydrogen_pos(const coot::named_torsion_t &named_tor, mmdb::Residue *res) const;
      
      std::vector<mmdb::Atom *> close_atoms(const clipper::Coord_orth &pt,
				       const std::vector<mmdb::Residue *> &env_residues) const;

      bash_distance_t find_bash_distance(const clipper::Coord_orth &ligand_atom_pos,
					 const clipper::Coord_orth &hydrogen_pos,
					 const std::vector<mmdb::Atom *> &close_residue_atoms) const;
      double get_radius(const std::string &ele) const;

      // find an atom (the atom, perhaps) bonded to lig_at that is not H_at.
      // Return its position. Can throw a std::runtime_error if not found.
      // 
      clipper::Coord_orth get_atom_pos_bonded_to_atom(mmdb::Atom *lig_at, mmdb::Atom *H_at, // not H_at
						      mmdb::Residue *ligand_residue,
						      const protein_geometry &geom) const;
      
      
   public:
      flev_attached_hydrogens_t(const dictionary_residue_restraints_t &restraints);

      std::vector<named_torsion_t> named_torsions;
      
      // fill the named_torsions vector, a trivial wrapper to the below function
      void cannonballs(mmdb::Residue *ligand_residue_3d,
		       const std::string &prodrg_3d_ligand_file_name,
		       const coot::dictionary_residue_restraints_t &restraints);

      // fill the named_torsions vector
      void cannonballs(mmdb::Residue *ligand_residue_3d,
		       mmdb::Manager *mol,
		       const coot::dictionary_residue_restraints_t &restraints);
      
      // apply those cannonball direction onto the real reference ligand:
      void distances_to_protein(mmdb::Residue *residue_reference,
				mmdb::Manager *mol_reference);
      void distances_to_protein_using_correct_Hs(mmdb::Residue *residue_reference,
						 mmdb::Manager *mol_reference,
						 const protein_geometry &geom);

      std::map<std::string, std::vector<coot::bash_distance_t> > atom_bashes;
   
   };

   
   class fle_residues_helper_t {
   public:
      bool is_set;
      clipper::Coord_orth transformed_relative_centre;
      clipper::Coord_orth interaction_position; // when the user
						// clicks on the
						// interaction of this
						// residue, this is
						// where we go.
      residue_spec_t spec;
      std::string residue_name;
      fle_residues_helper_t() { is_set = 0; }
      fle_residues_helper_t(const clipper::Coord_orth &t_r_pt,
			    const residue_spec_t &spec_in,
			    const std::string &res_name_in) {
	 transformed_relative_centre = t_r_pt;
	 spec = spec_in;
	 residue_name = res_name_in;
	 is_set = 1;
	 interaction_position = clipper::Coord_orth(0,0,0);
      }
      void set_interaction_position(const clipper::Coord_orth &p) {
	 interaction_position = p;
      } 
   };
   std::ostream& operator<<(std::ostream &s, fle_residues_helper_t fler);


   bool is_a_metal(mmdb::Residue *res);

   // The bonds from the protein to the ligand which contain
   // ligand-atom-name residue-spec and bond type (acceptor/donor).
   // These (ligand atom names) will have to be mapped to x y position
   // of the flat ligand.
   // 
   class fle_ligand_bond_t {
   public:
      enum ligand_bond_t {
         H_BOND_DONOR_MAINCHAIN,
         H_BOND_DONOR_SIDECHAIN,
         H_BOND_ACCEPTOR_MAINCHAIN,
         H_BOND_ACCEPTOR_SIDECHAIN,
         METAL_CONTACT_BOND,
         BOND_COVALENT,
         BOND_OTHER };  // must sync this to lbg.hh (why not extract it? (you can do it now))
      atom_spec_t ligand_atom_spec;
      int bond_type; // acceptor/donor

      residue_spec_t res_spec;
      atom_spec_t interacting_residue_atom_spec; // contains res_spec obviously.

      bool is_H_bond_to_water;
      double bond_length;  // from residue atom to ligand atom
      double water_protein_length; // if residue is a water, this is the closest
                                   // distance to protein (100 if very far).
      fle_ligand_bond_t(const atom_spec_t &ligand_atom_spec_in,
                        const atom_spec_t &interacting_residue_atom_spec_in,
                        int bond_type_in,
                        double bl_in,
                        bool is_water) {
         ligand_atom_spec = ligand_atom_spec_in;
         interacting_residue_atom_spec = interacting_residue_atom_spec_in;
         res_spec = residue_spec_t(interacting_residue_atom_spec_in);
         bond_type = bond_type_in;
         bond_length = bl_in;
         is_H_bond_to_water = is_water;
      }
      static int get_bond_type(mmdb::Atom *at_donor, mmdb::Atom *at_acceptor, bool ligand_atom_is_donor_flag) {
         int r_bond_type = BOND_OTHER;

         mmdb::Atom *ligand_atom = at_donor;
         mmdb::Atom *residue_atom = at_acceptor;

         if (at_donor) {
            if (at_acceptor) {

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
