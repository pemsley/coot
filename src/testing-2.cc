
#include "clipper/ccp4/ccp4_map_io.h"

#include "ideal/simple-restraint.hh"
#include "coot-utils/coot-map-utils.hh"
#include "coords/mmdb-crystal.h"

#include "testing.hh"
#include "testing-data.hh"


int test_map_tools() {

   int r = 0;
   return r;
}

int test_dreiding_torsion_energy() {

   int r = 0;
   std::string filename = greg_test("tutorial-modern.pdb");
   atom_selection_container_t atom_sel = get_atom_selection(filename, true, true);
   bool ifound = 0;
   testing_data t;

   int imod = 1;
   mmdb::Residue *residue_p = test_get_residue(atom_sel.mol, "B", 1);
   if (residue_p) {
      mmdb::PPAtom residue_atoms = 0;
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
   mmdb::Manager *mol = new mmdb::Manager;
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
   restraints.add_extra_restraints(er, t.geom); // we need a geom to look up expansions
                                                // for restraint atom names.
   restraints.minimize(flags);
   mol->WritePDBASCII("3tu4-test-37,38-par-planes-out.pdb");

   delete mol;
   return status;
} 

int test_string_splitting() {

   int status = 1;

   std::string s_1 = "HEADER    DNA-RNA HYBRID                          05-DEC-94   100D              ";
   std::string s_2 = "/xx/pemsley/ligand-analysis/output/5c/coot-ligand-analysis.log:metrics-for-ligand: \"/net/nfs5/gmssd/share/databases/pdb_data/pdb/5c/35c8.pdb\" \"L\" 212 \"\" \"NOX\" corr: 0.842709898948669 mogul: 3.51116991043091 bumps: 3 2 35.0000000  67 118 diff-map-stats: -0.0315599167321964 0.257134500774758 0.0189810063062419 1834.0 693.190333617851 23.9969999967143 0.0478254141277309 0.0379304815326336 -4.4906369112141e-6 0.0130845147201278 0.119391269981861 0.137771572925048 b-factor-metrics: 1.2905129105334 14.2150001525879 11.0150003433228 0.428217821782178 ";

   std::vector<std::string> bits_1 = coot::util::split_string_no_blanks(s_1);
   std::vector<std::string> bits_2 = coot::util::split_string_no_blanks(s_2);

   for (unsigned int i=0; i<bits_1.size(); i++) { 
      if (bits_1[i] == "" || bits_1[i] == " ") {
	 status = 0;
	 std::cout << "Fail on split_string_no_blanks() test-1" << std::endl;
	 break;
      }
   }
   for (unsigned int i=0; i<bits_2.size(); i++) { 
      if (bits_2[i] == "" || bits_2[i] == " ") {
	 status = 0;
	 std::cout << "Fail on split_string_no_blanks() test-2" << std::endl;
	 break;
      }
   }
   return status;
}
