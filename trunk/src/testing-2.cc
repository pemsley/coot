
#include "clipper/ccp4/ccp4_map_io.h"

#include "testing.hh"
#include "testing-data.hh"
#include "ideal/simple-restraint.hh"
#include "coot-utils/coot-map-utils.hh"
#include "coot-utils/emma.hh"
#include "coords/mmdb-crystal.h"

int test_map_tools() {

   int r = 0;

   clipper::Xmap<float> m;

   std::string mtz_file_name("test.mtz");
   std::string f_col = "FWT";
   std::string phi_col = "PHWT";
   coot::util::map_fill_from_mtz(&m, mtz_file_name, f_col, phi_col, "", 0, 0);
   clipper::Coord_orth c(-15,-4,21);
   coot::util::map_fragment_info_t n = coot::util::map_from_map_fragment(m, c, 20);

   clipper::CCP4MAPfile mapout;
   mapout.open_write("fragment.map");
   mapout.export_xmap(n.xmap);
   mapout.close_write();

   std::string pdb_file_name = "4c62-fragment.pdb";
   pdb_file_name = "4cgf.ent-coot-0.pdb";
   atom_selection_container_t asc = get_atom_selection(pdb_file_name, false);

   if (asc.mol) {
      coot::util::emma sphd(asc.mol, 5); // 5 is border
      sphd.integrate(n.xmap);
   } 
   return r;
}

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


int test_parallel_plane_restraints() {

   int status = 0;
   testing_data t;
   CMMDBManager *mol = new CMMDBManager;
   int ierr = mol->ReadPDBASCII("3tu4-test-37,38.pdb");
   std::cout << "ReadPDBASCII() returned " << ierr << std::endl;

   short int have_flanking_residue_at_start = 1;
   short int have_flanking_residue_at_end = 1;
   short int have_disulfide_residues = 0;
   std::string alt_conf = "";
   std::string chain_id("I");
   std::vector<coot::atom_spec_t> fixed_atom_specs;
   
   coot::restraints_container_t restraints(37,
					   38,
					   have_flanking_residue_at_start,
					   have_flanking_residue_at_end,
					   have_disulfide_residues,
					   alt_conf,
					   chain_id.c_str(),
					   mol,
					   fixed_atom_specs);

   short int do_rama_restraints = 0;
   short int do_residue_internal_torsions = 1;
   short int do_link_torsions = 0;
   float rama_plot_restraint_weight = 1.0;

   coot::restraint_usage_Flags flags =
      coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_PARALLEL_PLANES;
   coot::pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
   int nrestraints = 
      restraints.make_restraints(t.geom, flags,
				 do_residue_internal_torsions,
				 rama_plot_restraint_weight,
				 do_rama_restraints,
				 pseudos);

   std::string extra_restraints_file_name("test-base-pairing-extras-I-chain.txt");
   coot::extra_restraints_t er;
   er.read_refmac_extra_restraints(extra_restraints_file_name);
   restraints.add_extra_restraints(er, t.geom); // we need a geom to look up expansions for restraint atom names.
   restraints.minimize(flags);
   mol->WritePDBASCII("3tu4-test-37,38-par-planes-out.pdb");

   delete mol;
   return status;
} 
