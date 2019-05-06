/* src/testing-2.cc
 * 
 * Copyright 2015 by Medical Research Council
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
   clipper::Xmap<float> dummy_xmap;

   coot::restraints_container_t restraints(37,
					   38,
					   have_flanking_residue_at_start,
					   have_flanking_residue_at_end,
					   have_disulfide_residues,
					   alt_conf,
					   chain_id.c_str(),
					   mol,
					   fixed_atom_specs,
					   &dummy_xmap);

   short int do_rama_restraints = 0;
   short int do_residue_internal_torsions = 1;
   short int do_link_torsions = 0;
   float rama_plot_restraint_weight = 1.0;

   coot::restraint_usage_Flags flags =
      coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_PARALLEL_PLANES;
   coot::pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
   bool do_trans_peptide_restraints = false;
   int imol = 0; // dummy
   int nrestraints = 
      restraints.make_restraints(imol, t.geom, flags,
				 do_residue_internal_torsions,
				 do_trans_peptide_restraints,
				 rama_plot_restraint_weight,
				 do_rama_restraints, false, false,
				 pseudos);

   std::string extra_restraints_file_name("test-base-pairing-extras-I-chain.txt");
   coot::extra_restraints_t er;
   er.read_refmac_extra_restraints(extra_restraints_file_name);
   restraints.add_extra_restraints(imol, er, t.geom); // we need a geom to look up expansions
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

#include "utils/split-indices.hh"

int test_index_splitting() {

   int status = 1;

   unsigned int n_atoms = 303;
   unsigned int n_threads = 20;

   std::vector<std::pair<unsigned int, unsigned int> > air = coot::atom_index_ranges(n_atoms, n_threads);

   std::cout << "DEBUG:: test_index_splitting:: n_atoms " << n_atoms << " nt " << n_threads << std::endl;

   air = coot::atom_index_ranges(20, 20);
   for (unsigned int i=0; i<air.size(); i++)
      std::cout << "test-2:  atom_index_ranges  " << i << " : "
		<< air[i].first << " " << air[i].second << std::endl;
   air = coot::atom_index_ranges(20, 2);
   for (unsigned int i=0; i<air.size(); i++)
      std::cout << "test-2b: atom_index_ranges  " << i << " : "
		<< air[i].first << " " << air[i].second << std::endl;
   air = coot::atom_index_ranges(24, 6);
   for (unsigned int i=0; i<air.size(); i++)
      std::cout << "test-3:  atom_index_ranges  " << i << " : "
		<< air[i].first << " " << air[i].second << std::endl;
   air = coot::atom_index_ranges(n_atoms, n_threads);
   for (unsigned int i=0; i<air.size(); i++)
      std::cout << "test-1:  atom_index_ranges  " << i << " : "
		<< air[i].first << " " << air[i].second << std::endl;
   air = coot::atom_index_ranges(25, 6);
   for (unsigned int i=0; i<air.size(); i++)
      std::cout << "test-4:  atom_index_ranges  " << i << " : "
		<< air[i].first << " " << air[i].second << std::endl;
   air = coot::atom_index_ranges(5, 3);
   for (unsigned int i=0; i<air.size(); i++)
      std::cout << "test-5:  atom_index_ranges  " << i << " : "
		<< air[i].first << " " << air[i].second << std::endl;
   air = coot::atom_index_ranges(5, 2);
   for (unsigned int i=0; i<air.size(); i++)
      std::cout << "test-6:  atom_index_ranges  " << i << " : "
		<< air[i].first << " " << air[i].second << std::endl;

   // if (air[19].first == 304) status = 0;

   return status;
}
