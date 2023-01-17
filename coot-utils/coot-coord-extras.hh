/* coot-utils/coot-coord-extras.cc
 * 
 * Copyright 2004, 2005, 2006, 2007 by The University of York
 * Copyright 2009 by the University of Oxford
 * Copyright 2011 by the University of Oxford
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

#ifndef HAVE_COOT_COORD_EXTRAS_HH
#define HAVE_COOT_COORD_EXTRAS_HH

#include <fstream>
#include <map>
#include "tree.hh" // Kasper Peeters tree

#include "coot-coord-utils.hh"
#include "geometry/protein-geometry.hh"
#include "mini-mol/atom-quads.hh"
#include "mini-mol/mini-mol.hh"
#include "bonded-pairs.hh"

// functions and classes here (in extras) can use protein_geometry.


namespace coot {

   namespace util { 
      // geom_p gets updated to include the residue restraints if necessary
      // 
      std::pair<int, std::vector<std::string> >
      check_dictionary_for_residues(mmdb::PResidue *SelResidues, int nSelResidues,
				    protein_geometry *geom_p, int read_number);

      // We also now pass regular_residue_flag so that the indexing of the
      // contacts is inverted in the case of not regular residue.  I don't
      // know why this is necessary, but I have stared at it for hours, this
      // is a quick (ugly hack) fix that works.  I suspect that there is
      // some atom order dependency in mgtree that I don't understand.
      // Please fix (remove the necessity of depending on
      // regular_residue_flag) if you know how.
      // 
      std::vector<std::vector<int> >
      get_contact_indices_from_restraints(mmdb::Residue *residue,
					  protein_geometry *geom_p,
					  bool regular_residue_flag,
					  bool add_reverse_contacts);

      std::vector<std::vector<int> >
      get_contact_indices_from_restraints(mmdb::Residue *residue,
					  const dictionary_residue_restraints_t &restraints,
					  bool regular_residue_flag,
					  bool add_reverse_contacts);


      std::vector<std::vector<int> >
      get_contact_indices_for_PRO_residue(mmdb::PPAtom residue_atom,
					  int nResidueAtoms, 
					  protein_geometry *geom_p);

      // class definition is here but functionality is in molecule-class-info-other
      class missing_atom_info {
      public:
	 std::vector<std::string> residues_with_no_dictionary;
	 std::vector<mmdb::Residue *> residues_with_missing_atoms;
         std::map<mmdb::Residue *, std::vector<std::string> > residue_missing_atom_names_map;
	 std::vector<std::pair<mmdb::Residue *, std::vector<mmdb::Atom *> > > atoms_in_coords_but_not_in_dict;
	 missing_atom_info() {}
	 missing_atom_info(const std::vector<std::string> &residues_with_no_dictionary_in,
			   const std::vector<mmdb::Residue *>  &residues_with_missing_atoms_in,
			   const std::vector<std::pair<mmdb::Residue *, std::vector<mmdb::Atom *> > > &atoms_in_coords_but_not_in_dict_in) :
            residues_with_no_dictionary(residues_with_no_dictionary_in),
            residues_with_missing_atoms(residues_with_missing_atoms_in),
            atoms_in_coords_but_not_in_dict(atoms_in_coords_but_not_in_dict_in) {}
      };

      missing_atom_info missing_atoms(mmdb::Manager *mol,
                                      bool do_missing_hydrogen_atoms_flag,
                                      protein_geometry *geom_p);

      class dict_atom_info_t {
      public:
	 std::string name;
	 bool is_Hydrogen_flag;
	 dict_atom_info_t(const std::string &name_in, bool is_Hydrogen_flag_in) : name(name_in) {
	    is_Hydrogen_flag = is_Hydrogen_flag_in;
	 }
      };

      // a trivial helper class
      class dict_residue_atom_info_t {
      public:
	 std::string residue_name;
	 std::vector<dict_atom_info_t> atom_info;
	 dict_residue_atom_info_t(const std::string &residue_name_in,
				  const std::vector<dict_atom_info_t> &atom_info_in) :
            residue_name(residue_name_in), atom_info(atom_info_in)  {
	 }
	 // Here is the clever stuff, get the restraints info for a
	 // particular residue and from that set the atom_info.
	 dict_residue_atom_info_t(const std::string &residue_name,
				  protein_geometry *geom_p);
	 bool is_empty_p() const {
	    return (atom_info.size() == 0);
	 }
      };

      // do we need to pass read number to this function too?
      bool is_nucleotide_by_dict_dynamic_add(mmdb::Residue *residue_p, protein_geometry *geom_p);
      bool is_nucleotide_by_dict(mmdb::Residue *residue_p, const protein_geometry &geom);
   }

   mmdb::Residue *GetResidue(const minimol::residue &r); // For use with wiggly
					   // ligands, constructed from
					   // a minimol residue, the
					   // get_contact_indices_from_restraints()
					   // needs a mmdb::Residue *.  Caller disposes.
      


   class recursive_forwards_container_t {
   public:
      bool done;
      std::vector<int> forwards;
      recursive_forwards_container_t(bool done_in, const std::vector<int> &forwards_in) : forwards(forwards_in) {
	 done = done_in;
      }
      recursive_forwards_container_t() {
	 done = 0;
      } 
   };

   // Move the atoms in res_moving.
   // 
   // Return the number of rotated torsions.
   //
   class match_torsions {
      enum { REFERENCE_TORSION, MOVING_TORSION};
      mmdb::Residue *res_moving;
      mmdb::Residue *res_ref;
      dictionary_residue_restraints_t moving_residue_restraints;
      std::pair<bool, double> apply_torsion(const atom_name_quad &quad_moving,
					    const atom_name_quad &reference,
					    const std::string &alt_conf);
      std::pair<bool, double> apply_torsion_by_contacts(const atom_name_quad &quad_moving,
							const atom_name_quad &reference,
							const std::string &alt_conf);
      // return in radians
      std::pair<bool, double> get_torsion(int torsion_type, const atom_name_quad &quad) const;
      // return in radians
      std::pair<bool, double> get_torsion(mmdb::Residue *res, const atom_name_quad &quad) const;
      
   public: 
      match_torsions(mmdb::Residue *res_moving, mmdb::Residue *res_ref,
		     const dictionary_residue_restraints_t &moving_residue_restraints_in);
      int match (const std::vector<dict_torsion_restraint_t>  &tr_moving,
		 const std::vector<dict_torsion_restraint_t>  &tr_ref);

   };


   // which uses:
   //
   // (just the torsionable bonds (middle atom pairs)) by looking at
   // the monomer restraints
   // 
   // Don't return any hydrogen torsions - perhaps we should make that a
   // passed parameter.
   //
   // These functions should be part of protein_geometry, should they not?
   // 
   std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> >
   torsionable_bonds_monomer_internal(mmdb::Residue *residue_p,
				      mmdb::PPAtom atom_selection, int n_selected_atoms,
				      bool include_pyranose_ring_torsions_flag,
				      protein_geometry *geom_p);

   // The quad version of this (for actually setting torsions)
   // 
   std::vector<torsion_atom_quad>
   torsionable_bonds_monomer_internal_quads(mmdb::Residue *residue_p,
				      mmdb::PPAtom atom_selection, int n_selected_atoms,
				      bool include_pyranose_ring_torsions_flag,
				      protein_geometry *geom_p);

   // only uses the LINKR records in the first model.
   //
   bonded_pair_container_t    
   linkrs_in_atom_selection(mmdb::Manager *mol, mmdb::PPAtom atom_selection, int n_selected_atoms,
			    protein_geometry *geom_p);

   

   // Note: this is a simple-minded hack.  The right way of doing this
   // is to define a bonding tree that includes atoms from both
   // residues.  Then we don't need reference structures - the
   // "moving" residue atoms will get placed by internal coordinates.
   //
   // This can be viewed as starting (or test material) for the Proper
   // Way code.
   // 
   // Scenario Simple Beam-in:
   //    User has an ASN onto which they want to beam in a NAG.
   // 
   //    setup:
   //    Make mmdb::Residue *s and/or molecule for the N-linked NAG reference residues.
   // 
   ///   Get the mmdb::Residue * for the user residue ASN 
   //    Get the mmdb::Atoms *s for the OD1, ND2, GC and CB in user residue [1]
   //    Get the mmdb::Atoms *s for the OD1, ND2, GC and CB in N-linked ASN molecule [2]
   //
   //    LSQ fit the NAG residue from the reference ASN+NAG pair using
   //    matrix that rotates [2] onto [1].  (We don't need the ASN from
   //    that pair).  Now we can add that rotated NAG mmdb::Residue * to user
   //    molecule.  we have N-linked-NAG template - consider renaming to
   //    ASN-NAG-via-NAG-ASN i.e. the general case
   //    {ResType1}-{ResType2}-via-{LinkName} where ResType1 and ResType2
   //    are comp-ids, of course.  Actually, NAG-ASN is a pyranose-ASN
   //    link (group to comp_id). Hmm...
   //
   // This can throw a std::runtime_error if we can't find the group
   // of the input residues (for example).
   // 
   class beam_in_linked_residue {
      mmdb::Residue *residue_ref; // in user-defined molecule
      mmdb::Residue *template_res_ref;
      mmdb::Residue *template_res_mov;
      std::string comp_id_ref;
      std::string comp_id_new;
      protein_geometry *geom_p;
      std::string link_type;
      
      bool have_template;

      // return success status (0 = fail).
      std::vector<std::string> make_reference_atom_names(const std::string &comp_id) const;
      // get the given residue from the template coordinates
      mmdb::Residue *get_residue(const std::string &comp_id, mmdb::Manager*mol) const;
      std::vector<mmdb::Atom *> get_atoms(mmdb::Residue *residue_p,
				     const std::vector<std::string> &names) const;
      bool setup_by_comp_id(const std::string &comp_id_ref,
			    const std::string &new_res_type);
      bool setup_by_comp_id_group(const std::string &comp_id_ref,
				  const std::string &group_new);
      bool setup_by_group_group(const std::string &group_ref,
				const std::string &group_new);
      // move the residues of mov_res, don't change object variables.
      bool lsq_fit(mmdb::Residue *ref_res,
		   mmdb::Residue *matcher_res,
		   mmdb::Residue *mov_res,
		   const std::vector<std::string> &lsq_atom_names_ref,
		   const std::vector<std::string> &lsq_atom_names_match) const;
      // apply the chem mod (specifically, the CHEM_MOD_FUNCTION_DELETE
      // (delete all atoms with the given name)
      void delete_atom(mmdb::Residue *res, const std::string &atom_name) const;
      std::string atom_id_mmdb_expand(const std::string &atom_id,
				      const std::string &res_name,
				      int imol) const; 

      // If the link is a BETA1-6 or an ALPHA1-6 then the linked
      // residue (and the O6 of the residue to which we are adding)
      // can rotate about the C5-C6 bond (the template is just one of
      // the many options).
      clipper::Coord_orth get_O6_position_from_template() const;
      //
      // simply get the attached residue, don't handle the positioning
      // of the O6 on the residue to which we are adding.
      mmdb::Residue *get_residue_raw() const; 
      
   public:
      beam_in_linked_residue(mmdb::Residue *residue_ref,
			     const std::string &link_type_in,
			     const std::string &new_residue_type,
			     protein_geometry *geom_p);
      // This can return NULL if we were unable to make the residue to be attached.
      mmdb::Residue *get_residue() const;
   };


   // the class that is the template for the glyco_tree
   class linked_residue_t {
   public:
      mmdb::Residue * residue;
      std::string residue_name;
      std::string link_type; // to parent (root node has this as "")
      bool order_switch; // should be false most of the time
      linked_residue_t(mmdb::Residue *residue_in, const std::string &link_in) {
	 residue = residue_in;
	 if (residue)
	    residue_name = residue->GetResName();
	 link_type = link_in;
	 order_switch = false;
      }
      linked_residue_t() {
	 residue = NULL;
	 order_switch = false;
      }
      linked_residue_t(const std::string &res_name_in, const std::string &link_in) {
	 residue = NULL;
	 residue_name = res_name_in;
	 link_type = link_in;
	 order_switch = false;
      }
      std::string res_name() const {
	 if (residue)
	    return std::string(residue->GetResName());
	 else
	    return residue_name;
      } 
      bool operator==(const linked_residue_t &test_lr) const {

	 // should we test order switch here too?
	 
	 if (test_lr.link_type == link_type)
	    if (test_lr.res_name() == res_name())
	       return true;
	    else
	       return false; 
	 else
	    return false;
      }
      // needs testing
      bool operator<(const linked_residue_t &test_lr) const {
	 return (residue_spec_t(residue) < residue_spec_t(test_lr.residue));
      }
      friend std::ostream& operator<<(std::ostream &o, const linked_residue_t &lr);
   };
   std::ostream& operator<<(std::ostream &o, const linked_residue_t &lr);

   // use residues-near-residue to find linked residues
   std::vector<mmdb::Residue *> simple_residue_tree(mmdb::Residue *, mmdb::Manager *mol, float dist_max);

   class glyco_tree_t {
   public:
      class residue_id_t {
      public:
	 enum prime_arm_flag_t { UNSET, PRIME, NON_PRIME };
	 std::string res_type; // this is tested as empty to see if this object is filled
	 std::string link_type;
	 std::string parent_res_type;
	 residue_spec_t parent_res_spec;
	 unsigned int level;
	 prime_arm_flag_t prime_arm_flag; // are we in the (4') arm?
	 residue_id_t() {}
	 residue_id_t(int level_in,
		      prime_arm_flag_t prime_flag_in,
		      const std::string &res_type_in,
		      const std::string &link_type_in,
		      const std::string &parent_res_type_in,
		      const residue_spec_t &parent_res_spec_in) : res_type(res_type_in),
								  link_type(link_type_in),
								  parent_res_type(parent_res_type_in),
								  parent_res_spec(parent_res_spec_in),
								  level(level_in),
								  prime_arm_flag(prime_flag_in) {}
      };
   private:
      protein_geometry *geom_p;
      std::vector<mmdb::Residue *> linked_residues;
      tree<linked_residue_t> glyco_tree; // 20170503 the constructor now stores its own tree

      bool is_pyranose(mmdb::Residue *r) const; 
      tree<linked_residue_t> find_rooted_tree(mmdb::Residue *residue_root_p,
					      const std::vector<mmdb::Residue *> &residues) const;
      tree<linked_residue_t> find_ASN_rooted_tree(mmdb::Residue *residue_p,
						  const std::vector<mmdb::Residue *> &residues) const;
      tree<linked_residue_t> find_stand_alone_tree(const std::vector<mmdb::Residue *> &residues) const;
      void compare_vs_allowed_trees(const tree<linked_residue_t> &tr) const;
      bool compare_trees(const tree<linked_residue_t> &tree_for_testing,
			 const tree<linked_residue_t> &tree_reference) const;
      tree<linked_residue_t> oligomannose_tree() const;
      tree<linked_residue_t>      complex_tree() const;
      tree<linked_residue_t>       hybrid_tree() const;
      static bool residue_comparitor(mmdb::Residue *res1, mmdb::Residue *res2) {
	 return (residue_spec_t(res1) < residue_spec_t(res2));
      }
      void print(const tree<linked_residue_t> &glyco_tree) const;
      std::vector<mmdb::Residue *> residues(const tree<linked_residue_t> &glyco_tree) const;
      void output_internal_distances(mmdb::Residue *residue_p,
				     mmdb::Residue *parent_residue_p,
				     double dist_lim,
				     std::ofstream &f) const;
      // old, all-residue
      void output_internal_distances(mmdb::Residue *residue_p,
				     std::vector<mmdb::Residue *> residues,
				     double dist_lim,
				     std::ofstream &f) const;
      residue_id_t::prime_arm_flag_t get_prime(mmdb::Residue *residue_p) const;
      int get_level(mmdb::Residue *residue_p) const;

   public:
      glyco_tree_t(mmdb::Residue *residue_p, mmdb::Manager *mol, protein_geometry *geom_p_in);
      std::vector<mmdb::Residue *> residues(const coot::residue_spec_t &containing_res_spec) const;
      void internal_distances(double dist_lim, const std::string &file_name) const;
      residue_id_t get_id(mmdb::Residue *residue_p) const;
      // for tree comparison
      tree<linked_residue_t> get_glyco_tree() const { return glyco_tree; }
      bool compare_trees(const tree<linked_residue_t> &tree_in) const;
      std::vector<std::pair<coot::residue_spec_t, coot::residue_spec_t> > matched_pairs(const tree<linked_residue_t> &t_in) const;
   };
}

#endif // HAVE_COOT_COORD_EXTRAS_HH
