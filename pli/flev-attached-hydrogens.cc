
#include "coot-utils/coot-coord-utils.hh"
#include "flev-annotations.hh" // break out flev_attached_hydrogens_t from here?

// RDKit version (well, the version that is used when the hydrogens
// are correctly named (according to the dictionary) and placed on the
// ligand of interest.
//
// This fills the atom_bashes vector.
//
void
pli::flev_attached_hydrogens_t::distances_to_protein_using_correct_Hs(mmdb::Residue *ligand_residue,
                                                                      mmdb::Manager *mol,
                                                                      const coot::protein_geometry &geom) {

   // the constructor (called just before this) should fill
   // atoms_with_rotating_hydrogens and atoms_with_riding_hydrogens
   // vectors (using the restraints).

   float radius = 6.0;
   std::vector<mmdb::Residue *> env_residues =
      coot::residues_near_residue(ligand_residue, mol, radius);

   // -------------------------------------------------------------------
   //                    riding hydrogens
   // -------------------------------------------------------------------
   //
   mmdb::PPAtom residue_atoms = 0;
   int n_ligand_atoms;
   ligand_residue->GetAtomTable(residue_atoms, n_ligand_atoms);
   for (unsigned int irh=0; irh<atoms_with_riding_hydrogens.size(); irh++) {
      mmdb::Atom *lig_at = NULL;
      mmdb::Atom *H_at = NULL;
      for (int iat=0; iat<n_ligand_atoms; iat++) {
         std::string atom_name(residue_atoms[iat]->name);
         if (atom_name == atoms_with_riding_hydrogens[irh].first)
            lig_at = residue_atoms[iat];
         if (atom_name == atoms_with_riding_hydrogens[irh].second)
            H_at = residue_atoms[iat];
         if (lig_at && H_at)
            break;
      }
      if (lig_at && H_at) {
         clipper::Coord_orth H_pt(H_at->x, H_at->y, H_at->z);
         clipper::Coord_orth lig_atom_pt(lig_at->x, lig_at->y, lig_at->z);

         std::vector<mmdb::Atom *> atoms = close_atoms(H_pt, env_residues);
         coot::bash_distance_t bash = find_bash_distance(lig_atom_pt, H_pt, atoms);
         atom_bashes[atoms_with_riding_hydrogens[irh].first].push_back(bash);
         if (true)
            std::cout << " adding bash distance " << bash << " to atom "
                      << atoms_with_riding_hydrogens[irh].first << std::endl;
      }
   }

   // -------------------------------------------------------------------
   //                 rotatable hydrogens (more complex)
   // -------------------------------------------------------------------
   //

   for (unsigned int irh=0; irh<atoms_with_rotating_hydrogens.size(); irh++) {
      mmdb::Atom *lig_at = NULL;
      mmdb::Atom *H_at = NULL;
      for (int iat=0; iat<n_ligand_atoms; iat++) {
         std::string atom_name(residue_atoms[iat]->name);
         if (atom_name == atoms_with_rotating_hydrogens[irh].first)
            lig_at = residue_atoms[iat];
         if (atom_name == atoms_with_rotating_hydrogens[irh].second)
            H_at = residue_atoms[iat];
         if (lig_at && H_at)
            break;
      }
      if (lig_at && H_at) {
         clipper::Coord_orth H_pt(H_at->x, H_at->y, H_at->z);
         clipper::Coord_orth lig_atom_pt(lig_at->x, lig_at->y, lig_at->z);

         std::vector<mmdb::Atom *> atoms = close_atoms(H_pt, env_residues);

         try {
            clipper::Coord_orth vector_pt = get_atom_pos_bonded_to_atom(lig_at, H_at, // not H_at
                                                                        ligand_residue, geom);
            clipper::Coord_orth base_ref_pt(0,0,0);
            double tors = clipper::Coord_orth::torsion(base_ref_pt, vector_pt, lig_atom_pt, H_pt);
            double dist = sqrt((lig_atom_pt - H_pt).lengthsq());
            double angle = clipper::Coord_orth::angle(vector_pt, lig_atom_pt, H_pt);

            int n_torsion_samples = 8;
            for (int itor=0; itor<n_torsion_samples; itor++) {

               double tmp_tor_d =  double(itor) * 360.0/double(n_torsion_samples);
               double tmp_tor = clipper::Util::d2rad(tmp_tor_d);
               tmp_tor += tors;
               clipper::Coord_orth new_pt =
                  clipper::Coord_orth(base_ref_pt, vector_pt, lig_atom_pt, dist, angle, tmp_tor);
               coot::bash_distance_t bash = find_bash_distance(lig_atom_pt, new_pt, atoms);
               atom_bashes[atoms_with_rotating_hydrogens[irh].first].push_back(bash);
            }
         }
         catch (const std::runtime_error &rte) {
            std::cout << rte.what() << std::endl;
         }
      }
   }
}



coot::bash_distance_t
pli::flev_attached_hydrogens_t::find_bash_distance(const clipper::Coord_orth &ligand_atom_pos,
                                                    const clipper::Coord_orth &hydrogen_pos,
                                                    const std::vector<mmdb::Atom *> &close_residue_atoms) const {

   // find the residue from the close residue atoms and cache the
   // dictionaries here so that we can then call
   // dictionary_map[residue_p].type_energy(atom_name) and use that
   // type energy (checking for not "") to find if the atom is a
   // hydrogen bond accetor (or both) to use
   // energy_lib_t::some_new_accessor_function_hb_type(type_energy).
   // If hb_type is acceptor, then decrease bash distance,
   // atom_radius_plus_cbr by 0.8A or so.
   //
   std::map<mmdb::Residue *, coot::dictionary_residue_restraints_t> dictionary_map;

   double cannonball_radius = 0.8; // radius of the cannonball, c.f. at least a hydrogen.

   double max_dist = 4.05; // if we have travelled 4A without bashing
                           // into a protein atom then this has
                           // practically unlimited substitution distance.

   double max_dist_squared = max_dist * max_dist;
   clipper::Coord_orth h_vector((hydrogen_pos - ligand_atom_pos).unit());

   if (0)
      std::cout << "h_vector: " << h_vector.format() << " from hydrogen pos: "
                << hydrogen_pos.format() << "and ligand atom pos: " << ligand_atom_pos.format()
                << std::endl;

   // set the atomic radii:
   //
   std::vector<double> radius(close_residue_atoms.size());
   for (unsigned int iat=0; iat<close_residue_atoms.size(); iat++) {
      std::string ele(close_residue_atoms[iat]->element);
      radius[iat] = get_radius(ele);
   }

   coot::bash_distance_t dd;

   std::vector<clipper::Coord_orth> atom_positions(close_residue_atoms.size());
   // likewise set the atom positions so that we don't have to keep doing it.
   for (unsigned int i=0; i<close_residue_atoms.size(); i++)
      atom_positions[i] = clipper::Coord_orth(close_residue_atoms[i]->x,
                                              close_residue_atoms[i]->y,
                                              close_residue_atoms[i]->z);

   for (double slide=0; slide<=max_dist; slide+=0.04) {
      clipper::Coord_orth test_pt = ligand_atom_pos + slide * h_vector;
      if (true)
         std::cout << "   bash distance for ligand atom at " << ligand_atom_pos.format() << " "
                   << "determined from " << atom_positions.size() << " atom positions"
                   << std::endl;
      for (unsigned int iat=0; iat<atom_positions.size(); iat++) {
         double atom_radius_plus_cbr = radius[iat] + cannonball_radius;
         double d_squared = (test_pt - atom_positions[iat]).lengthsq();
         if (true)
            std::cout << "   atom " << iat << " "
                      << close_residue_atoms[iat]->GetChainID() << " "
                      << close_residue_atoms[iat]->GetSeqNum() << " "
                      << close_residue_atoms[iat]->GetAtomName() << " "
                      << " slide: " << slide
                      << " comparing " << sqrt(d_squared) << "^2  and "
                    << atom_radius_plus_cbr << "^2" << std::endl;
         if (d_squared < atom_radius_plus_cbr*atom_radius_plus_cbr) {
            dd = coot::bash_distance_t(slide);
            break;
         }
      }
      if (dd.limited) {
         break;
      }
   }
   return dd;
}



// What are the atoms that are close (distance < 6A) to pt?
//
// waters are not counted as close atoms.
//
std::vector<mmdb::Atom *>
pli::flev_attached_hydrogens_t::close_atoms(const clipper::Coord_orth &pt,
                                            const std::vector<mmdb::Residue *> &env_residues) const {

   std::vector<mmdb::Atom *> v;
   double dist_crit = 6.0;
   double dist_crit_squared = dist_crit * dist_crit;

   for (unsigned int i=0; i<env_residues.size(); i++) {
      mmdb::Residue *residue_p = env_residues[i];
      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         clipper::Coord_orth atom_pos(residue_atoms[iat]->x, residue_atoms[iat]->y, residue_atoms[iat]->z);
         double d_squared = (pt - atom_pos).lengthsq();
         if (d_squared < dist_crit_squared) {
            std::string rn(residue_atoms[iat]->GetResName());
            if (rn != "HOH")
               v.push_back(residue_atoms[iat]);
         }
      }
   }
   return v;
}

double
pli::flev_attached_hydrogens_t::get_radius(const std::string &ele) const {

   double radius = 1.70;
   if (ele == " H")
      radius = 1.20;
   if (ele == " N")
      radius = 1.55;
   if (ele == " O")
      radius = 1.52;
   if (ele == " S")
      radius = 1.8;
   return radius;
}
