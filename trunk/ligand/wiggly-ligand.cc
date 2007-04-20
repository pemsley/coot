/* ligand/wiggly-ligand.cc 
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006 by Paul Emsley, The University of York
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

#include "coot-utils.hh"  // needed, it seems for the case where GSL
			  // is not defined.
#include "mgtree.h"

#include "wligand.hh"
#include "regularize-minimol.hh"

// Ligand installation from coordinates (an minimol), number of trials.
//
// In the constructor we check all atoms of the coordinates are
// in the restraints of that residue type.
//
// We need to modify restraints_container_t to read the atom list too
// (see notes).
// 
//  done  (it was half done anyway).
//

std::pair<short int, std::string>
coot::wligand::install_simple_wiggly_ligands(coot::protein_geometry *pg,
					     const coot::minimol::molecule &ligand_in,
					     int n_samples) {

   short int istat = 0;
   std::string m = ""; 
   
   // m_torsions: a set of torsions for this monomer
  // 
   std::string monomer_type = get_monomer_type_from_mol(ligand_in);
//    std::cout << "DEBUG:: in install_simple_wiggly_ligands: "
// 	     << "get_monomer_type_from_mol returns :"
// 	     << monomer_type << ":" << std::endl;
   short int do_hydrogen_torsions_flag = 0;
   std::vector <coot::dict_torsion_restraint_t> m_torsions =
      pg->get_monomer_torsions_from_geometry(monomer_type, do_hydrogen_torsions_flag);


   if (m_torsions.size() == 0) {
      std::cout << "Requested flexible molecule for ligand " << monomer_type
		<< " but no non-Hydrogen rotatable bonds found.\n"
		<< " Did you forget to read the dictionary?\n";
      
      m = "Requested flexible molecule for ligand ";
      m += monomer_type;
      m += "\n"; 
      m += " but no non-Hydrogen rotatable bonds found.\n";
      m += " Did you forget to read the dictionary?\n";
      
      istat = 0;
      return std::pair<short int, std::string> (istat, m);

   } else {
      istat = 1; // OK.... so far.
   }

   // This should be inside the else (then we can remove the above
   // return) , it's just a mess to do.
   // 
   for (int isample=0; isample<n_samples; isample++) {

      coot::minimol::molecule ligand = ligand_in; // local changable copy
      
      // Let's make the coordinates:
      // 
      std::vector< ::Cartesian> coords;
      for(int ifrag=0; ifrag<ligand.fragments.size(); ifrag++) {
	 for (int ires=ligand[ifrag].min_res_no(); ires<=ligand[ifrag].max_residue_number(); ires++) {
	    for (unsigned int iat=0; iat<ligand[ifrag][ires].atoms.size(); iat++) {
	       // coords.push_back(coord_orth_to_cartesian(ligand[ifrag][ires][iat].pos));
	       coords.push_back( ::Cartesian(ligand[ifrag][ires][iat].pos.x(),
					     ligand[ifrag][ires][iat].pos.y(),
					     ligand[ifrag][ires][iat].pos.z()));
	    }
	 }
      }
	 
      //  Rotatable bonds -> coord indices
      // 
      std::vector<coot::atom_name_pair> atom_name_pairs =
	 get_torsion_bonds_atom_pairs(monomer_type, m_torsions);

      std::vector<coot::atom_index_pair> atom_index_pairs =
	 get_atom_index_pairs(atom_name_pairs, ligand);

      if (atom_name_pairs.size() != atom_index_pairs.size()) { 
	 std::cout << "DISASTER:: atom_name_pairs.size() != atom_index_pairs.size()\n";
	 std::cout << "atom_name_pairs.size()   " << atom_name_pairs.size()  << std::endl;
	 std::cout << "atom_index_pairs.size()  " << atom_index_pairs.size() << std::endl;

      }

      std::vector<std::vector<int> > contacts = getcontacts(ligand);
//       for (int i=0; i<contacts.size(); i++) {
// 	 for (int j=0; j<contacts[i].size(); j++) {
// 	    std::cout << "contacts " << i << " to " << contacts[i][j] << std::endl;
// 	 }
//       }

      Tree tree;
      // std::cout << "Filling tree with " << coords.size() << " coordinates \n";

      tree.SetCoords(coords, 0, contacts); // compile XXX

      std::vector<float> torsion_set = get_torsions_by_random(m_torsions);

      // Here: check that torsion_set, atom_index_pairs and
      // atom_name_pairs all have the same size.
      
//       std::cout << " There are " << atom_index_pairs.size() 
// 	        << " atom_index_pairs" << std::endl;
      for (unsigned int ibond=0; ibond< atom_index_pairs.size(); ibond++) {
         if (m_torsions[ibond].periodicity() > 0)  {  
// 	    std::cout << "rotating " << torsion_set[ibond] << " round " 
// 		      << atom_index_pairs[ibond].index1 << " "
// 		      << atom_index_pairs[ibond].index2 << "\n";
	    tree.SetDihedralAngle(atom_index_pairs[ibond].index1,
				  atom_index_pairs[ibond].index2,
				  clipper::Util::d2rad(torsion_set[ibond]));
         }
      }

      int ncoord = 0;
      std::vector< ::Cartesian > rotated_coords = tree.GetAllCartesians();
      for(int ifrag=0; ifrag<ligand.fragments.size(); ifrag++) {
	 for (int ires=ligand[ifrag].min_res_no(); ires<=ligand[ifrag].max_residue_number(); ires++) {
	    for (unsigned int iat=0; iat<ligand[ifrag][ires].atoms.size(); iat++) {
	       ligand[ifrag][ires][iat].pos =
		  clipper::Coord_orth(rotated_coords[ncoord].get_x(),
				      rotated_coords[ncoord].get_y(),
				      rotated_coords[ncoord].get_z());
// 	       std::cout << "ipdating position to " << ligand[ifrag][ires][iat].pos.format()
// 		         << std::endl;
	       ncoord++;
	    }
	 }
      }


      if (1) {  // debugging wiggly ligands
	 std::string filename = "wligand-";
	 filename += int_to_string(isample);
	 filename += ".pdb";
	 ligand.write_file(filename, default_b_factor);
      }
      
#ifdef HAVE_GSL      
      coot::minimol::molecule reg_ligand = coot::regularize_minimol_molecule(ligand, *pg);
      // coot::minimol::molecule reg_ligand;
      install_ligand(reg_ligand);
#else
      install_ligand(ligand);
#endif
   }
   return std::pair<short int, std::string> (istat, m);
}

std::string
coot::wligand::get_monomer_type_from_mol(const coot::minimol::molecule &mol) const {

   std::string r;
   short int ifound = 0;
   
   for(int ifrag=0; ifrag<mol.fragments.size(); ifrag++) {
      for (int ires=mol[ifrag].min_res_no(); ires<=mol[ifrag].max_residue_number(); ires++) {
	 if (mol[ifrag][ires].n_atoms() > 0) {
	    r = mol[ifrag][ires].name;
	    ifound = 1;
	    break;
	 }
      }
      if (ifound)
	 break;
   }

   std::cout << "get_monomer type: " << r << std::endl;
   return r;
}


std::vector<std::vector<int> >
coot::wligand::getcontacts(const coot::minimol::molecule &mol) {

   std::vector<coot::minimol::atom *> atoms = mol.select_atoms_serial();
   std::vector<std::vector<int> > contacts;

   for (unsigned int i=0; i<atoms.size(); i++) {
      std::vector<int> v;
      contacts.push_back(v);
      for (unsigned int j=0; j<atoms.size(); j++) {
	 if (j != i) {
	    if (clipper::Coord_orth::length(atoms[i]->pos, atoms[j]->pos) < 1.85) {
	       contacts[i].push_back(j);
	    }
	 }
      }
   }
   return contacts;
}

std::vector <float>
coot::wligand::get_torsions_by_random(const std::vector <coot::dict_torsion_restraint_t> &m_torsions) const { 
   
   std::vector <float> sample_tors(m_torsions.size());
   float prob;
   float r;

   float max_prob = 1.0;
   float non_rotating_torsion_cut_off = 11.0;

   if (0) { 
      for(int itor=0; itor<m_torsions.size(); itor++) {
	 std::cout << "DEBUG input torsion: " << itor << " "
		   << m_torsions[itor].atom_id_2_4c() << " "
		   << m_torsions[itor].atom_id_3_4c() << " "
		   << m_torsions[itor].angle() << " " << m_torsions[itor].esd()
		   << std::endl;
      }
   }
   
   for(int itor=0; itor<m_torsions.size(); itor++) {
      if (m_torsions[itor].periodicity() == 1) {
	 if (m_torsions[itor].esd() < non_rotating_torsion_cut_off) {
	    sample_tors[itor] = m_torsions[itor].angle(); // don't sample it.
	 } else {
	    sample_tors[itor] = m_torsions[itor].angle();
	    double v = get_random_normal_value();
	    sample_tors[itor] += v*m_torsions[itor].esd();
	 } 
      } else {
	 // Get a random periodicity sample, e.g. 1 from {0,1,2}
	 float frac = float (coot::util::random())/float (RAND_MAX);
	 float p_f = float(m_torsions[itor].periodicity());
	 // do floor()?
	 int irandom_periodicity_incidence = int(p_f * frac);
	 float random_periodicity_incidence = float (irandom_periodicity_incidence);
	 // e.g. random_periodicity_incidence is 1.0 from {0.0, 1.0, 2.0}
	 float periodicity_angle = 360.0 * (random_periodicity_incidence/p_f);
	 sample_tors[itor] = m_torsions[itor].angle();
	 sample_tors[itor] += periodicity_angle;

	 if (m_torsions[itor].esd() < non_rotating_torsion_cut_off) {
	    // add nothing
	 } else {
	    // add random selection of distribution using angle esd
	    double v = get_random_normal_value();
	    sample_tors[itor] += v*m_torsions[itor].esd();
	 } 
      }
   }

   if (0) { // debugging torsions
      for(int itor=0; itor<m_torsions.size(); itor++)
	 std::cout << "DEBUG:: sampled torsion " << itor << "  "
		   <<  sample_tors[itor] << std::endl;
   }

   return sample_tors; 
}

// should this be a coot util function?
//
// get random selection from a normal distribution, mean 0, std: 1
// 
double
coot::wligand::get_random_normal_value() const {

   // now index a random number into that cumlative vector
   //
   float x = -16.0; // some nonsense initial value
   float sum = cumulative.back();
   float r = sum*coot::util::random()/float(RAND_MAX);

   for (unsigned int i=0; i<cumulative.size(); i++) {
      if (r < cumulative[i]) {
	 x = (cumulative_step*float(i))-4.0;
	 if (i>0) { 
	    float cum_i_1 = cumulative[i-1];
	    float cum_i   = cumulative[i];
	    float i_interp = (float(i)-1) + (r-cum_i_1)/(cum_i - cum_i_1);
	    x = (cumulative_step*i_interp)-4.0;
	 }
	 break;
      }
   }
   return x;
}


#ifdef COMPILE_OLD_STUFF

std::vector <float>
coot::wligand::get_torsions_by_random_old(const std::vector <coot::dict_torsion_restraint_t> &m_torsions) const { 
   
   std::vector <float> sample_tors(m_torsions.size());
   float prob;
   float r;

   float max_prob = 1.0;
   for (int itor=0; itor<m_torsions.size(); itor++)
      if (m_torsions[itor].periodicity() > 0) 
         max_prob *= 1/(m_torsions[itor].esd() * sqrt(2.0 * M_PI));  // 1/[2s sqrt(2pi)], what is that?
   std::cout << "max_prob: " << max_prob << std::endl;
   
   for (;;) {
      // Only sample torsions with eds below
      // non_rotating_torsion_cut_off at exact positions,
      // corresponding to repeats due to periodicity
      // 
      float non_rotating_torsion_cut_off = 10.0; 
      for(int itor=0; itor<m_torsions.size(); itor++) {
         if (m_torsions[itor].periodicity() > 0) {
	    if (m_torsions[itor].periodicity() == 1) {
	       sample_tors[itor] = m_torsions[itor].angle(); // don't sample it.
	    } else { 
	       if (m_torsions[itor].esd() > non_rotating_torsion_cut_off) { 
		  sample_tors[itor] = m_torsions[itor].angle() + 
		     360.0 * float (coot::util::random())/float (RAND_MAX);
	       } else {
		  sample_tors[itor] = m_torsions[itor].angle();
		  float frac = float (coot::util::random())/float (RAND_MAX);
		  float p_f = float(m_torsions[itor].periodicity());
		  int n_random_periodicity_incidence = int(p_f * frac);
		  float fnr = float(n_random_periodicity_incidence);
		  sample_tors[itor] += fnr * 360/p_f;
	       }
	       if (sample_tors[itor] > 360.0) {
		  sample_tors[itor] -= 360.0;
	       }
	    }
         }
      }
      r = max_prob * coot::util::random()/ float (RAND_MAX);
      prob = probability_of_torsions(m_torsions, sample_tors);

      std::cout << "max_prob: " << max_prob << std::endl;
      std::cout << "comparing " << prob << " and " << r << std::endl;
      if (prob > r) break;
   }
   std::cout << "get_torsions by random returns : " << std::endl;
   for(int itor=0; itor<m_torsions.size(); itor++)
      std::cout << "    " << itor << "  " <<  sample_tors[itor] << std::endl;

    return sample_tors; 
}

#endif // 0000, don't compile

// Possibly not specific to wligand:
float
coot::wligand::probability_of_torsions(const std::vector <coot::dict_torsion_restraint_t> &m_torsions,
				       const std::vector <float> &r) const {
   double pr = 1.0;
   if (m_torsions.size() != r.size()) {
      std::cout << "ERROR: this should never happen in wligand::probability_of_torsions"
		<< std::endl;
      return -999.0;
   } else {

      //for(int i=0; i<m_torsions.size(); i++) {
         //std::cout << "torsion " << i << " " 
		 //<< m_torsions[i].atom_id_1_4c() << " " 
		 //<< m_torsions[i].atom_id_2_4c() << " " 
		 //<< m_torsions[i].atom_id_3_4c() << " " 
		 //<< m_torsions[i].atom_id_4_4c() << " " 
		 //<< m_torsions[i].angle() << " " 
		 //<< m_torsions[i].esd() << " " 
		 //<< m_torsions[i].periodicity() << " " 
		 //<<std::endl;
      //}

      double z;
      double s;
      for(int i=0; i<r.size(); i++) {
	 int per = m_torsions[i].periodicity();

	 if (per > 0 ) { 
		 // i.e. actually is a torsion - not some kludged up plane restraint...

	    double diff = 99999.9; 
	    double tdiff;
	    double trial_target;
	    for(int iper=0; iper<per; iper++) { 
	       trial_target =  m_torsions[i].angle()+ double(iper)*360.0/double(per);
	       if (trial_target > 360.0) trial_target -= 360.0; 
	       tdiff = r[i] - trial_target; 
	       if (abs(tdiff) < abs(diff)) { 
	          diff = tdiff;
	       }
	    }
	    if (diff == 99999.9) { 
	       std::cout << "Error in periodicity (" << per << ") check" << std::endl;
	       std::cout << "target_value: " << m_torsions[i].angle()
		         << ", theta: " << r[i] << std::endl;
	    }
   
	    z = diff/m_torsions[i].esd();
	    s = 1/(m_torsions[i].esd() * sqrt(2.0 * M_PI));
 	    std::cout << "DEBUG:: torsion " << i << " " << diff << "/" << m_torsions[i].esd()
 		      << " multiplying " << pr << " by " << s
 		      << " and " << exp( -0.5 * z *z ) << " z is " << z << "\n";
	    pr *= s * exp( -0.5 * z *z );
	 }
      }
   }
   return pr;
} 
				       
 
std::vector<coot::atom_index_pair>
coot::wligand::get_atom_index_pairs(std::vector<coot::atom_name_pair>atom_name_pairs,
				    const coot::minimol::molecule &ligand) const {
   
   int nres = 0;
   int i_store_index; 
   std::vector<coot::atom_index_pair> index_pairs;


   for(int ifrag=0; ifrag<ligand.fragments.size(); ifrag++) {
      for (int ires=ligand[ifrag].min_res_no(); ires<=ligand[ifrag].max_residue_number(); ires++) {
	 for (unsigned int ipair=0; ipair<atom_name_pairs.size(); ipair++) {
	    i_store_index = -1;
	    for (unsigned int iat=0; iat<ligand[ifrag][ires].atoms.size(); iat++) {
	       //
	       if (ligand[ifrag][ires][iat].name == atom_name_pairs[ipair].atom1) {
		  i_store_index = iat;
	       }
	    }
	    for (unsigned int iat=0; iat<ligand[ifrag][ires].atoms.size(); iat++) {
	       if (ligand[ifrag][ires][iat].name == atom_name_pairs[ipair].atom2) {
		  if (i_store_index > -1) { 
		     index_pairs.push_back(coot::atom_index_pair(i_store_index, iat));
		  }
	       }
	    }
	 }
      }
   }
   
   return index_pairs; 
}


std::vector<coot::atom_name_pair>
coot::wligand::get_torsion_bonds_atom_pairs(const std::string &monomer_type,
					    const std::vector <coot::dict_torsion_restraint_t> &monomer_torsions) const {

   std::vector<coot::atom_name_pair> atom_pairs;
   for(int i=0; i<monomer_torsions.size(); i++) {
      coot::atom_name_pair pair(monomer_torsions[i].atom_id_2_4c(),
				monomer_torsions[i].atom_id_3_4c());
      atom_pairs.push_back(pair);
   }
   return atom_pairs;
} 
