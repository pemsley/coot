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

#include "use-rdkit.hh"

#include <mmdb/mmdb_manager.h>
#include "protein-geometry.hh"
#include "lbg-molfile.hh"

namespace coot { 

   // can throw an runtime_error exception (e.g. residue not in dictionary)
   //
   RDKit::RWMol rdkit_mol_sanitized(CResidue *residue_p, const protein_geometry &geom);
   RDKit::RWMol rdkit_mol(CResidue *residue_p, const protein_geometry &geom);
   RDKit::RWMol rdkit_mol(CResidue *residue_p, const coot::dictionary_residue_restraints_t &restraints);
   // tinker with mol
   void rdkit_mol_sanitize(RDKit::RWMol &mol);

   // tweaking function used by above (change mol maybe).
   // @return the added hydrogen name - or "" if nothing was added.
   //
   // When adding and atom, try to find the name of the Hydrogen from
   // the bond restraints.  If not found, add an atom called "-".
   // 
   std::string add_H_to_ring_N_as_needed(RDKit::RWMol *mol,
				  int idx, const std::string &atom_name,
				  const coot::dictionary_residue_restraints_t &restraints); 
   
   int add_2d_conformer(RDKit::ROMol *rdkmol_in, double weight_for_3d_distances); // tweak rdkmol_in
   RDKit::Bond::BondType convert_bond_type(const std::string &t);

   lig_build::molfile_molecule_t make_molfile_molecule(const RDKit::ROMol &rdkm, int iconf);

   // Returns NULL on fail.
   //
   // Caller deletes.
   //
   // resulting residues has chain id of "" and residue number 1.
   //
   CResidue *make_residue(const RDKit::ROMol &rdkm, int iconf, const std::string &res_name);
   
   lig_build::bond_t::bond_type_t convert_bond_type(const RDKit::Bond::BondType &type);
   
   // This can throw a std::exception
   void remove_non_polar_Hs(RDKit::RWMol *rdkm); // fiddle with rdkm


   // a wrapper for the above, matching hydrogens names to the
   // dictionary.  Add atoms to residue_p, return success status and error
   // message pair
   // 
   std::pair<bool, std::string>
   add_hydrogens_with_rdkit(CResidue *residue_p,
			    const dictionary_residue_restraints_t &restraints);

   std::string infer_H_name(int iat,
			    RDKit::ATOM_SPTR atom_p,
			    const RDKit::ROMol *mol,
			    const dictionary_residue_restraints_t &restraints,
			    const std::vector<std::string> &H_names_already_added);

   //
   void undelocalise(RDKit::RWMol *rdkm); // fiddle with rdkm
   void undelocalise_phosphates(RDKit::ROMol *rdkm);
   void undelocalise_aminos(RDKit::RWMol *rdkm);
   void undelocalise_methyl_carboxylates(RDKit::RWMol *rdkm);
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

   // account for Ns with "too many" hydrogens by assigning deleting a
   // hydrogen.
   // 
   void delete_excessive_hydrogens(RDKit::RWMol *rdkm);

   // used in the rdkit_mol() "constructor".
   // 
   RDKit::Atom::ChiralType get_chiral_tag(CResidue *residue_p,
					  const dictionary_residue_restraints_t &restraints,
					  CAtom *atom_p);

   // update the atom positions of the rdkit_molecule from residue_p
   // 
   void update_coords(RDKit::RWMol *mol, int iconf, CResidue *residue_p);

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
}

#endif // RDKIT_INTERFACE_HH
