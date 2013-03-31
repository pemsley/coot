#include "testing.hh"
#include "testing-data.hh"


int test_dreiding_torsion_energy() {

   int r = 0;
   std::string filename = greg_test("tutorial-modern.pdb");
   atom_selection_container_t atom_sel = get_atom_selection(filename, 1);
   bool ifound = 0;
   testing_data t;

   int imod = 1;
   CResidue *residue_p = test_get_residue(atom_sel.mol, "B", 1);
   if (residue_p) {
      PPCAtom residue_atoms = 0;
      int n_residue_atoms;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      std::string comp_id = residue_p->GetResName();
      std::vector<coot::torsion_atom_quad> quads =
	 torsionable_bonds_monomer_internal_quads(residue_p, residue_atoms, n_residue_atoms,
						  false, &t.geom);
      std::cout << "# torsionable quads: " << quads.size() << std::endl;
      for (unsigned int i=0; i<quads.size(); i++) {
	 double tors = quads[i].torsion();
// 	 double d = t.geom.dreiding_torsion_energy(comp_id, quads[i]);
// 	 std::cout << "   " << i << " " << quads[i] << " " << tors << " " << d << std::endl;
      }
   }
   return r;
}
