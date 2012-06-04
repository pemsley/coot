
namespace coot {

   // private (no SWIG interface)
   // 
   // the engine for the above calls
   std::pair<CMMDBManager *, CResidue *>
   regularize_inner(PyObject *rdkit_mol,
		    PyObject *restraints_py,
		    const std::string &res_name);

   std::pair<CMMDBManager *, CResidue *>
   regularize_inner(RDKit::ROMol &mol,
		    PyObject *restraints_py,
		    const std::string &res_name);
   

   bool is_const_torsion(const RDKit::ROMol &mol,
			 const RDKit::Atom *at_2,
			 const RDKit::Atom *at_3);
   // also private (no interface)
   // 
   // now update the atom positions of the rdkit_molecule from residue_p
   // (perhaps this should be in rdkit-interface.hh?)
   // 
   void update_coords(RDKit::RWMol *mol, int iconf, CResidue *residue_p);

      // alter restraints
   int assign_chirals(const RDKit::ROMol &mol, coot::dictionary_residue_restraints_t *restraints);
   // alter restraints
   void add_chem_comp_atoms(const RDKit::ROMol &mol, coot::dictionary_residue_restraints_t *restraints);
   // alter restraints
   void add_chem_comp_planes(const RDKit::ROMol &mol, coot::dictionary_residue_restraints_t *restraints);
   // alter restraints
   void add_chem_comp_aromatic_planes(const RDKit::ROMol &mol, coot::dictionary_residue_restraints_t *restraints);
   // alter restraints
   void add_chem_comp_deloc_planes(const RDKit::ROMol &mol, coot::dictionary_residue_restraints_t *restraints);


   // alter restraints
   void fill_with_energy_lib_bonds(const RDKit::ROMol &mol,
				   const coot::energy_lib_t &energy_lib,
				   coot::dictionary_residue_restraints_t *restraints);
   // alter restraints
   void fill_with_energy_lib_angles(const RDKit::ROMol &mol,
				    const coot::energy_lib_t &energy_lib,
				    coot::dictionary_residue_restraints_t *restraints);

   // alter restraints
   void fill_with_energy_lib_torsions(const RDKit::ROMol &mol,
				      const coot::energy_lib_t &energy_lib,
				      coot::dictionary_residue_restraints_t *restraints);

   std::string convert_to_energy_lib_bond_type(RDKit::Bond::BondType bt);


}

