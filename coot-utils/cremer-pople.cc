
#include <vector>
#include <string>
#include <map>
#include <algorithm>
// why don't I need iostream here?
#include <iomanip>
#include "coot-coord-utils.hh"
#include "cremer-pople.hh"

coot::cremer_pople_t::cremer_pople_t(mmdb::Residue *residue_p) {

   auto is_member = [] (const std::string &an, const std::vector<std::string> &atom_names) {
      return std::find(atom_names.begin(), atom_names.end(), an) != atom_names.end();
   };


   bool debug = false;

   if (debug)
      std::cout << "::: Residue " << coot::residue_spec_t(residue_p) << " "
                << residue_p->GetResName() << std::endl;

   std::vector<std::string> atom_names = {" C1 ", " C2 ", " C3 ", " C4 ", " C5 ", " O5 "};

   // we will deal wih alt conf residues another time
   std::map<std::string, mmdb::Atom *> atoms_map;

   mmdb::Atom **residue_atoms = 0;
   int n_residue_atoms = 0;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      mmdb::Atom *at = residue_atoms[iat];
      if (! at->isTer()) {
         std::string atom_name(at->GetAtomName());
         if (is_member(atom_name, atom_names)) {
            atoms_map[atom_name] = at;
         }
      }
   }
   if (atoms_map.size() == 6) {
      std::vector<clipper::Coord_orth> positions;
      for (auto atom : atoms_map) {
         clipper::Coord_orth c(atom.second->x, atom.second->y, atom.second->z);
         positions.push_back(c);
      }
      lsq_plane_info_t lsq_plane_info(positions);

      // the majority of the carbons should have a positive displacement
      // so we test here if we need a sign flip:
      unsigned int n_positive = 0;
      for (auto atom : atoms_map) {
         std::string ele = atom.second->element;
         if (ele == " C") {
            clipper::Coord_orth c(atom.second->x, atom.second->y, atom.second->z);
            double d = lsq_plane_info.plane_deviation(c);
            if (d> 0.0)
               n_positive++;
         }
      }
      double sign_flip = 1.0;
      if (n_positive < 3)
         sign_flip = -1.0;

      // We need indexing now, so convert the atoms_map to a vector
      unsigned int hash = 0; // sort of a hash
      std::vector<std::pair<std::string, mmdb::Atom *>> atoms_vector(6);
      for (unsigned int j=0;j<6; j++) {
         const std::string &atom_name = atom_names[j];
         auto it = atoms_map.find(atom_name);
         if (it != atoms_map.end()) {
            atoms_vector[j] = *it;
            hash += (j + 1);
         } else {
            std::cout << "ERROR:: missing atom " << atom_name << " in residue "
                      << coot::residue_spec_t(residue_p) << std::endl;
         }
      }

      if (hash == 21) { // 1 + 2 + 3 + 4 + 5 + 6

         const double k1 = - std::sqrt(2.0/6.0);
         const double k2 =   std::sqrt(1.0/6.0);
         double c2_sum = 0.0;
         double s2_sum = 0.0;
         double q3_sum = 0.0;
         for (unsigned int jj=0;jj<6; jj++) {
            mmdb:: Atom *atom = atoms_vector[jj].second;
            unsigned int j = jj + 1; // j is 1 to 6
            clipper::Coord_orth c(atom->x, atom->y, atom->z);
            double z_j = sign_flip * lsq_plane_info.plane_deviation(c);
            if (debug)
               std::cout << "   for atom " << coot::atom_spec_t(atom) << " z_j: " << z_j << std::endl;
            double angle = 2.0 * M_PI * (static_cast<double>(j-1)) / 3.0;
            double c_part = z_j * cos(angle);
            double s_part = z_j * sin(angle);
            c2_sum += c_part;
            s2_sum += s_part;
            double q3_part = z_j * std::pow(-1, (j-1));
            q3_sum += q3_part;
         }
         double c2 = k1 * c2_sum;
         double s2 = k1 * s2_sum;
         double q3 = k2 * q3_sum;
         double q2 = std::sqrt(c2*c2 + s2*s2);
         double Q  = std::sqrt(q2*q2 + q3*q3);
         double theta = atan2(q2, q3);
         double phi   = atan2(-s2, c2);
         double theta_d = clipper::Util::rad2d(theta);
         double   phi_d = clipper::Util::rad2d(phi);
         phi_d -= 120.0; // to match privateer
         if (theta_d < 0.0) theta_d += 360.0;
         if (phi_d   < 0.0)   phi_d += 360.0;
         std::string space = " ";
         std::string conf = "-";
         if (fabs(theta_d) < 20.0)         conf = "4c1";
         if (fabs(theta_d - 360.0) < 20.0) conf = "4c1";
         if (fabs(theta_d - 180.0) < 20.0) conf = "1c4";
         if (fabs(theta_d + 180.0) < 20.0) conf = "1c4";
         if (residue_p->GetSeqNum() > 9) space = "";
         std::cout << "Residue " << coot::residue_spec_t(residue_p) << " " << space
                   << residue_p->GetResName()
                   << "   Q: "    << std::setw(8) << Q
                   << "  theta: " << std::setw(8) << theta_d
                   << "  phi: "   << std::setw(8) << phi_d
                   << "  conf: "  << std::setw(8) << conf << std::endl;
      }
   }
}
