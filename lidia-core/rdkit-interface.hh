/* lidia-core/rdkit-interface.hh
 * 
 * Copyright 2010, 2011 by The University of Oxford
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

#ifndef RDKIT_INTERFACE_HH
#define RDKIT_INTERFACE_HH

#include "compat/coot-sysdep.h"
#include "use-rdkit.hh"

#include <mmdb2/mmdb_manager.h>
#include "geometry/protein-geometry.hh"
#include "lbg-molfile.hh"

namespace coot { 

   // can throw an runtime_error exception (e.g. residue not in dictionary)
   //
   RDKit::RWMol rdkit_mol_sanitized(mmdb::Residue *residue_p, int imol_enc,
				    const protein_geometry &geom);

   // sets properties _Name, ResName ResName ChainID alt_id
   //
   RDKit::RWMol rdkit_mol(mmdb::Residue *residue_p, int imol_enc,
			  const protein_geometry &geom);

   // For this function, should we do this at end?
   // 
   // for (unsigned int iat=0; iat<m.getNumAtoms(); iat++)
   //     m[iat]->calcImplicitValence(true);
   // 
   RDKit::RWMol rdkit_mol(const dictionary_residue_restraints_t &restraints); // fill the coord from
                                                                              // the dictionary if
                                                                              // you can.

   // Fill a conformer with the 2D depiction if you can. The first is the conformer id, -1 means failure
   std::pair<int, RDKit::RWMol> rdkit_mol_with_2d_depiction(const dictionary_residue_restraints_t &restraints);

   // this can return a molecule with 0 atoms (e.g. atom is a ZN and there is 1 atom in the dictionary -
   // so no bonded atoms - hence empty molecule is returned).
   //
   // If non-empty, the `alt_conf` means "only select conformers with alt-confs that match `alt_conf`."
   //
    RDKit::RWMol rdkit_mol(mmdb::Residue *residue_p,
			  const dictionary_residue_restraints_t &restraints,
			  const std::string &alt_conf="",
			  bool undelocalise=true);
   // return a kekulized molecule
   RDKit::RWMol remove_Hs_and_clean(const RDKit::ROMol &m, bool set_aromaticity=false);

   void set_atom_chirality(RDKit::Atom *rdkit_at,
			   mmdb::Atom *atom_name,
			   mmdb::Residue *residue_p,
			   const coot::dictionary_residue_restraints_t &restraints);
   void set_atom_chirality(RDKit::Atom *rdkit_at, const dict_atom &dict_atom);

   // tinker with mol
   void set_3d_conformer_state(RDKit::RWMol *mol); // hack the setting of 3D state, seems not to
                                                   // be done for mdl files when zs are 0.   
   bool has_zero_coords(RDKit::RWMol *mol, unsigned int iconf); // e.g. reading from a MolFile,
                                            // all coords are 0.0.
                                            // in such a case we need to do a Compute2DCoords()
                                            // before showing the molecule in lidia.
   void rdkit_mol_sanitize(RDKit::RWMol &mol);
   // tinker with mol
   void mogulify_mol(RDKit::RWMol &mol);
   void charge_guanidinos(RDKit::RWMol *rdkm);
   void mogulify_nitro_groups(RDKit::RWMol *rdkm);
   bool chiral_check_order_swap(RDKit::Atom *at_1, RDKit::Atom *at_2,
				const std::vector<dict_chiral_restraint_t>  &chiral_restraints);
   bool chiral_check_order_swap(RDKit::Atom *at_1, RDKit::Atom *at_2);
   bool chiral_check_order_swap_singleton(RDKit::Atom *at_1, RDKit::Atom *at_2,
					  const dictionary_residue_restraints_t &restraints);


   // tweaking function used by above (change mol maybe).
   // @return the added hydrogen name - or "" if nothing was added.
   //
   // When adding and atom, try to find the name of the Hydrogen from
   // the bond restraints.  If not found, add an atom called "-".
   // 
   std::string add_H_to_ring_N_as_needed(RDKit::RWMol *mol,
				  int idx, const std::string &atom_name,
				  const dictionary_residue_restraints_t &restraints); 

   // can throw an RDKit::ConformerException (std::exception)
   // can return -1 if current conformer is 3D.
   int add_2d_conformer(RDKit::ROMol *rdkmol_in, double weight_for_3d_distances); // tweak rdkmol_in
   RDKit::Bond::BondType convert_bond_type(const std::string &t);

   lig_build::molfile_molecule_t make_molfile_molecule(const RDKit::ROMol &rdkm, int iconf);

   // Returns NULL on fail.
   //
   // Caller deletes.
   //
   // resulting residues has chain id of "" and residue number 1.
   //
   mmdb::Residue *make_residue(const RDKit::ROMol &rdkm, int iconf, const std::string &res_name);
   dictionary_residue_restraints_t make_dictionary(const RDKit::ROMol &rdkm, int iconf, const std::string &res_name);

   
   lig_build::bond_t::bond_type_t convert_bond_type(const RDKit::Bond::BondType &type);
   
   // This can throw a std::exception
   void remove_non_polar_Hs(RDKit::RWMol *rdkm); // fiddle with rdkm


   // a wrapper for the above, matching hydrogens names to the
   // dictionary.  Add atoms to residue_p, return success status and error
   // message pair
   // 
   std::pair<bool, std::string>
   add_hydrogens_with_rdkit(mmdb::Residue *residue_p,
			    const dictionary_residue_restraints_t &restraints);

   std::string infer_H_name(int iat,
			    RDKit::Atom *atom_p,
			    const RDKit::ROMol *mol,
			    const dictionary_residue_restraints_t &restraints,
			    const std::vector<std::string> &H_names_already_added);

   // print out the atom names RDKit bond types
   void debug_rdkit_molecule(const RDKit::ROMol *rdkm);
   //
   void undelocalise(RDKit::RWMol *rdkm); // fiddle with rdkm
   void undelocalise_phosphates(RDKit::ROMol *rdkm);
   void undelocalise_sulphates(RDKit::ROMol *rdkm);
   void undelocalise_aminos(RDKit::RWMol *rdkm);
   void undelocalise_nitros(RDKit::RWMol *rdkm);
   void undelocalise_carboxylates(RDKit::RWMol *rdkm);
   void undelocalise_methyl_carboxylates(RDKit::RWMol *rdkm);
   void charge_metals(RDKit::RWMol *rdkm);
   void charge_sp3_borons(RDKit::RWMol *rdkm);
   void charge_undelocalized_guanidinos(RDKit::RWMol *rdkm); // A + on the  C.

   // which calls (not for public use)
   // fiddle with the bonds in rdkm as needed.
   void deloc_O_check_inner(RDKit::RWMol *rdkm, RDKit::Atom *central_C,
			    RDKit::Atom *O1, RDKit::Atom *O2,
			    RDKit::Bond *b1, RDKit::Bond *b2);


   // try to add (for instance) a +1 to Ns with 4 bonds.
   //
   enum { FORMAL_CHARGE_UNKNOWN = -1 };
   
   // account for Ns with "too many" hydrogens by assigning a formal
   // charge to the N.
   void assign_formal_charges(RDKit::RWMol *rdkm);

   // when using a refmac cif dictionary, we construct a molecule with
   // deloc bonds (e.g. on a phosphate)
   // valence on P: (1 1/2) * 3 + 1 -> 6 => problem
   void charge_phosphates(RDKit::RWMol *rdkm);

   // These work on the valence model of sulphate and phospates: i.e. those
   // that have single bonds and a double bond to their oxygens from the P/S.
   // 
   int remove_phosphate_hydrogens(RDKit::RWMol *m, bool deloc_bonds); 
   int remove_sulphate_hydrogens (RDKit::RWMol *m, bool deloc_bonds); 
   // which are wrappers for the non-user function
   int remove_PO4_SO4_hydrogens(RDKit::RWMol *m, 
				unsigned int atomic_num, 
				bool deloc_bonds);

   int remove_carboxylate_hydrogens(RDKit::RWMol *m, bool deloc_bonds); 

   // account for Ns with "too many" hydrogens by assigning deleting a
   // hydrogen.
   // 
   void delete_excessive_hydrogens(RDKit::RWMol *rdkm);

   // used in the rdkit_mol() "constructor".
   // 
   RDKit::Atom::ChiralType get_chiral_tag(mmdb::Residue *residue_p,
					  const dictionary_residue_restraints_t &restraints,
					  mmdb::Atom *atom_p);
   
   RDKit::Atom::ChiralType get_chiral_tag_v2(mmdb::Residue *residue_p,
					     const dictionary_residue_restraints_t &restraints,
					     mmdb::Atom *atom_p);

   bool cip_rank_sorter(const std::pair<const RDKit::Atom *, unsigned int> &at_1,
			const std::pair<const RDKit::Atom *, unsigned int> &at_2);

   // are all the bonds between the atoms (in the vector) all aromatic?
   //
   bool is_aromatic_ring(const std::vector<int> &ring_atom_indices,
			 RDKit::ROMol &rdkm);


   // add "energy_type" properties to the atoms
   void set_energy_lib_atom_types(RDKit::ROMol *mol);

   // dicitionaries from CCDs don't have energy atom types. We need them for ligand
   // environment analysis (flev).  It is presumed that mol has energy_type atom types
   // (e.g. set from the above function). mol is not modified.
   //
   void set_dictionary_atom_types_from_mol(dictionary_residue_restraints_t *dictionary,
					   const RDKit::ROMol *mol);

   // make an rdkit molecule from dictionary and use the above function to the energy types
   // why is this here and not in protein-geometry? (e.g. after having parsed a cif file
   // check for the existance of atom types and if they are not there, add them)
   // Currently protein-geometry doesn't know about the RDKit.  This is the lowest
   // level of Coot that know about the RDKit.
   void set_dictionary_atom_types(dictionary_residue_restraints_t *dictionary);

   // update the atom positions of the rdkit_molecule from residue_p
   // 
   void update_coords(RDKit::RWMol *mol, int iconf, mmdb::Residue *residue_p);

   // return a copy of the input molecule having deleted the atoms of the R-group
   //
   // bond_index is the index of the bond to be replace.  atom_index
   // is the index of the atom (i.e. one of the atoms in bond of
   // bond_index), the R-group of which the atoms are marked for
   // deletion and the atom_index atom is changed to "*" for later
   // modification with various R-groups.
   // 
   RDKit::ROMol *split_molecule(const RDKit::ROMol &mol_in, int bond_index, int atom_index); 

   // Caller deletes the molecules of the vector.
   // 
   std::vector<RDKit::ROMol *> join_molecules(const RDKit::ROMol &mol, int atom_index,
					      const RDKit::ROMol &trial_fragment);

   // test calling with bad conf_id
   mmdb::Residue *residue_from_rdkit_mol(const RDKit::ROMol &mol_in, int conf_id, const std::string &new_comp_id);

   // like the above but make the dictionary for the above
   // All atoms must have atom names or this will return std::nullopt
   std::optional<dictionary_residue_restraints_t> dictionary_from_rdkit_mol(const RDKit::ROMol &mol_in, int conf_id, const std::string &new_comp_id);
}

#endif // RDKIT_INTERFACE_HH
