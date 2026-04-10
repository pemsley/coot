/* coot-utils/coot-coord-utils-gemmi.cc
 *
 * Copyright 2026 by Paul Emsley
 * Copyright 2026 by Medical Research Council
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

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <set>
#include <numeric>
#include "coot-coord-utils.hh"
#include "coot-coord-utils-gemmi.hh"

// ==================== trim_atom_names ====================

static std::string trim_spaces(const std::string &s) {
   auto start = s.find_first_not_of(' ');
   if (start == std::string::npos) return "";
   auto end = s.find_last_not_of(' ');
   return s.substr(start, end - start + 1);
}

void
coot::trim_atom_names(gemmi::Structure &st) {
   for (auto &model : st.models)
      for (auto &chain : model.chains)
         for (auto &res : chain.residues)
            for (auto &at : res.atoms)
               at.name = trim_spaces(at.name);
}

// ==================== atom-level ====================

double
coot::distance(const gemmi::Atom &at_1, const gemmi::Atom &at_2) {
   double dx = at_1.pos.x - at_2.pos.x;
   double dy = at_1.pos.y - at_2.pos.y;
   double dz = at_1.pos.z - at_2.pos.z;
   return std::sqrt(dx*dx + dy*dy + dz*dz);
}

double
coot::angle(const gemmi::Atom &at_1, const gemmi::Atom &at_2, const gemmi::Atom &at_3) {
   double ax = at_1.pos.x - at_2.pos.x;
   double ay = at_1.pos.y - at_2.pos.y;
   double az = at_1.pos.z - at_2.pos.z;
   double bx = at_3.pos.x - at_2.pos.x;
   double by = at_3.pos.y - at_2.pos.y;
   double bz = at_3.pos.z - at_2.pos.z;
   double dot = ax*bx + ay*by + az*bz;
   double la = std::sqrt(ax*ax + ay*ay + az*az);
   double lb = std::sqrt(bx*bx + by*by + bz*bz);
   if (la < 1e-12 || lb < 1e-12)
      throw std::runtime_error("angle: degenerate atom positions");
   double cos_theta = dot / (la * lb);
   cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
   return std::acos(cos_theta) * 180.0 / M_PI;
}

clipper::Coord_orth
coot::co(const gemmi::Atom &at) {
   return clipper::Coord_orth(at.pos.x, at.pos.y, at.pos.z);
}

bool
coot::is_hydrogen_atom(const gemmi::Atom &at) {
   return at.is_hydrogen();
}

// ==================== residue-level ====================

std::pair<bool, clipper::Coord_orth>
coot::util::get_residue_centre(const gemmi::Residue &res) {
   if (res.atoms.empty())
      return std::make_pair(false, clipper::Coord_orth(0,0,0));
   double sx = 0, sy = 0, sz = 0;
   int n = 0;
   for (const auto &at : res.atoms) {
      sx += at.pos.x;
      sy += at.pos.y;
      sz += at.pos.z;
      n++;
   }
   if (n == 0)
      return std::make_pair(false, clipper::Coord_orth(0,0,0));
   return std::make_pair(true, clipper::Coord_orth(sx/n, sy/n, sz/n));
}

std::pair<bool, clipper::Coord_orth>
coot::util::get_CA_position_in_residue(const gemmi::Residue &res) {
   for (const auto &at : res.atoms) {
      if (at.name == "CA")
         return std::make_pair(true, clipper::Coord_orth(at.pos.x, at.pos.y, at.pos.z));
   }
   return std::make_pair(false, clipper::Coord_orth(0,0,0));
}

std::pair<bool, clipper::Coord_orth>
coot::util::get_CB_position_in_residue(const gemmi::Residue &res) {
   for (const auto &at : res.atoms) {
      if (at.name == "CB")
         return std::make_pair(true, clipper::Coord_orth(at.pos.x, at.pos.y, at.pos.z));
   }
   return std::make_pair(false, clipper::Coord_orth(0,0,0));
}

short int
coot::util::is_nucleotide(const gemmi::Residue &r) {
   // match the mmdb version logic: check by residue name
   std::string rn = r.name;
   if (rn == "A"  || rn == "C"  || rn == "G"  || rn == "U"  || rn == "T"  ||
       rn == "DA" || rn == "DC" || rn == "DG" || rn == "DT" ||
       rn == "Ar" || rn == "Cr" || rn == "Gr" || rn == "Ur" || rn == "Tr" ||
       rn == "Ad" || rn == "Cd" || rn == "Gd" || rn == "Td")
      return 1;
   return 0;
}

bool
coot::util::residue_has_hydrogens_p(const gemmi::Residue &res) {
   for (const auto &at : res.atoms)
      if (at.is_hydrogen())
         return true;
   return false;
}

int
coot::util::residue_has_hetatms(const gemmi::Residue &res) {
   if (res.atoms.empty())
      return -1;
   if (res.het_flag == 'H')
      return 1;
   return 0;
}

// ==================== chain-level ====================

std::pair<int, int>
coot::util::min_and_max_residues(const gemmi::Chain &chain) {
   int min_res = 9999;
   int max_res = -9999;
   for (const auto &res : chain.residues) {
      int seqnum = res.seqid.num.value;
      if (seqnum < min_res) min_res = seqnum;
      if (seqnum > max_res) max_res = seqnum;
   }
   return std::make_pair(min_res, max_res);
}

std::pair<bool, int>
coot::util::min_resno_in_chain(const gemmi::Chain &chain) {
   if (chain.residues.empty())
      return std::make_pair(false, 0);
   int min_res = chain.residues[0].seqid.num.value;
   for (const auto &res : chain.residues) {
      int seqnum = res.seqid.num.value;
      if (seqnum < min_res) min_res = seqnum;
   }
   return std::make_pair(true, min_res);
}

std::pair<bool, int>
coot::util::max_resno_in_chain(const gemmi::Chain &chain) {
   if (chain.residues.empty())
      return std::make_pair(false, 0);
   int max_res = chain.residues[0].seqid.num.value;
   for (const auto &res : chain.residues) {
      int seqnum = res.seqid.num.value;
      if (seqnum > max_res) max_res = seqnum;
   }
   return std::make_pair(true, max_res);
}

std::vector<std::string>
coot::util::residue_types_in_chain(const gemmi::Chain &chain) {
   std::set<std::string> types_set;
   for (const auto &res : chain.residues)
      types_set.insert(res.name);
   return std::vector<std::string>(types_set.begin(), types_set.end());
}

std::pair<unsigned int, unsigned int>
coot::util::get_number_of_protein_or_nucleotides(const gemmi::Chain &chain) {
   unsigned int n_prot = 0;
   unsigned int n_nuc = 0;
   for (const auto &res : chain.residues) {
      if (is_nucleotide(res))
         n_nuc++;
      else if (res.entity_type == gemmi::EntityType::Polymer)
         n_prot++;
   }
   return std::make_pair(n_prot, n_nuc);
}

// ==================== structure-level ====================

std::pair<bool, clipper::Coord_orth>
coot::centre_of_molecule(const gemmi::Structure &st) {
   double sx = 0, sy = 0, sz = 0;
   int n = 0;
   for (const auto &model : st.models) {
      for (const auto &chain : model.chains) {
         for (const auto &res : chain.residues) {
            for (const auto &at : res.atoms) {
               sx += at.pos.x;
               sy += at.pos.y;
               sz += at.pos.z;
               n++;
            }
         }
      }
      break; // first model only, matching mmdb version
   }
   if (n == 0)
      return std::make_pair(false, clipper::Coord_orth(0,0,0));
   return std::make_pair(true, clipper::Coord_orth(sx/n, sy/n, sz/n));
}

std::pair<bool, double>
coot::radius_of_gyration(const gemmi::Structure &st) {
   auto centre = coot::centre_of_molecule(st);
   if (!centre.first)
      return std::make_pair(false, 0.0);
   double sum_sq = 0;
   int n = 0;
   for (const auto &model : st.models) {
      for (const auto &chain : model.chains) {
         for (const auto &res : chain.residues) {
            for (const auto &at : res.atoms) {
               double dx = at.pos.x - centre.second.x();
               double dy = at.pos.y - centre.second.y();
               double dz = at.pos.z - centre.second.z();
               sum_sq += dx*dx + dy*dy + dz*dz;
               n++;
            }
         }
      }
      break;
   }
   if (n == 0)
      return std::make_pair(false, 0.0);
   return std::make_pair(true, std::sqrt(sum_sq / n));
}

bool
coot::mol_has_symmetry(const gemmi::Structure &st) {
   // match mmdb version: check if symmetry operations are available,
   // not just whether spacegroup_hm is set
   const gemmi::SpaceGroup *sg = st.find_spacegroup();
   if (sg)
      return sg->operations().order() > 1; // more than just identity
   return false;
}

bool
coot::mol_is_anisotropic(const gemmi::Structure &st) {
   for (const auto &model : st.models)
      for (const auto &chain : model.chains)
         for (const auto &res : chain.residues)
            for (const auto &at : res.atoms)
               if (at.aniso.nonzero())
                  return true;
   return false;
}

float
coot::get_position_hash(const gemmi::Structure &st) {
   // match mmdb version: sum of (x - x_prev) differences, skipping TER atoms
   float hash = 0.0;
   for (const auto &model : st.models) {
      for (const auto &chain : model.chains) {
         float x_prev = 0.0;
         unsigned int atom_count = 0;
         for (const auto &res : chain.residues) {
            for (const auto &at : res.atoms) {
               if (atom_count > 0) {
                  hash += static_cast<float>(at.pos.x) - x_prev;
               }
               atom_count++;
               x_prev = static_cast<float>(at.pos.x);
            }
         }
      }
      break;
   }
   return hash;
}

std::vector<std::string>
coot::util::residue_types_in_molecule(const gemmi::Structure &st) {
   std::set<std::string> types_set;
   for (const auto &model : st.models) {
      for (const auto &chain : model.chains)
         for (const auto &res : chain.residues)
            types_set.insert(res.name);
      break;
   }
   return std::vector<std::string>(types_set.begin(), types_set.end());
}

std::vector<std::string>
coot::util::non_standard_residue_types_in_molecule(const gemmi::Structure &st) {
   std::vector<std::string> standard = standard_residue_types();
   std::set<std::string> standard_set(standard.begin(), standard.end());
   std::set<std::string> non_standard_set;
   for (const auto &model : st.models) {
      for (const auto &chain : model.chains)
         for (const auto &res : chain.residues)
            if (standard_set.find(res.name) == standard_set.end())
               non_standard_set.insert(res.name);
      break;
   }
   return std::vector<std::string>(non_standard_set.begin(), non_standard_set.end());
}

std::vector<std::string>
coot::util::chains_in_molecule(const gemmi::Structure &st) {
   std::vector<std::string> r;
   for (const auto &model : st.models) {
      for (const auto &chain : model.chains)
         r.push_back(chain.name);
      break;
   }
   return r;
}

int
coot::util::number_of_residues_in_molecule(const gemmi::Structure &st) {
   int n = 0;
   for (const auto &model : st.models) {
      for (const auto &chain : model.chains)
         n += chain.residues.size();
      break;
   }
   return n;
}

int
coot::util::max_number_of_residues_in_chain(const gemmi::Structure &st) {
   int max_n = -1;
   for (const auto &model : st.models) {
      for (const auto &chain : model.chains) {
         int n = chain.residues.size();
         if (n > max_n) max_n = n;
      }
      break;
   }
   return max_n;
}

int
coot::util::number_of_chains(const gemmi::Structure &st) {
   if (st.models.empty()) return -1;
   return st.models[0].chains.size();
}

std::pair<bool, int>
coot::util::max_resno_in_molecule(const gemmi::Structure &st) {
   bool found = false;
   int max_res = -9999;
   for (const auto &model : st.models) {
      for (const auto &chain : model.chains) {
         for (const auto &res : chain.residues) {
            int seqnum = res.seqid.num.value;
            if (seqnum > max_res) {
               max_res = seqnum;
               found = true;
            }
         }
      }
      break;
   }
   return std::make_pair(found, max_res);
}

int
coot::util::max_min_max_residue_range(const gemmi::Structure &st) {
   int max_range = -1;
   for (const auto &model : st.models) {
      for (const auto &chain : model.chains) {
         if (chain.residues.empty()) continue;
         auto mm = min_and_max_residues(chain);
         int range = mm.second - mm.first + 1;
         if (range > max_range) max_range = range;
      }
      break;
   }
   return max_range;
}

std::vector<std::string>
coot::util::alt_confs_in_molecule(const gemmi::Structure &st) {
   // match mmdb version: include empty string "" for non-alt atoms,
   // and single-char strings for alt confs
   std::set<std::string> alts_set;
   for (const auto &model : st.models) {
      for (const auto &chain : model.chains)
         for (const auto &res : chain.residues)
            for (const auto &at : res.atoms) {
               if (at.altloc == '\0' || at.altloc == ' ')
                  alts_set.insert(std::string(""));
               else
                  alts_set.insert(std::string(1, at.altloc));
            }
      break;
   }
   return std::vector<std::string>(alts_set.begin(), alts_set.end());
}

std::pair<clipper::Coord_orth, clipper::Coord_orth>
coot::util::extents(const gemmi::Structure &st) {
   double xmin = 1e30, ymin = 1e30, zmin = 1e30;
   double xmax = -1e30, ymax = -1e30, zmax = -1e30;
   for (const auto &model : st.models) {
      for (const auto &chain : model.chains) {
         for (const auto &res : chain.residues) {
            for (const auto &at : res.atoms) {
               if (at.pos.x < xmin) xmin = at.pos.x;
               if (at.pos.y < ymin) ymin = at.pos.y;
               if (at.pos.z < zmin) zmin = at.pos.z;
               if (at.pos.x > xmax) xmax = at.pos.x;
               if (at.pos.y > ymax) ymax = at.pos.y;
               if (at.pos.z > zmax) zmax = at.pos.z;
            }
         }
      }
      break;
   }
   return std::make_pair(clipper::Coord_orth(xmin, ymin, zmin),
                         clipper::Coord_orth(xmax, ymax, zmax));
}

clipper::Coord_orth
coot::util::median_position(const gemmi::Structure &st) {
   std::vector<double> xs, ys, zs;
   for (const auto &model : st.models) {
      for (const auto &chain : model.chains) {
         for (const auto &res : chain.residues) {
            for (const auto &at : res.atoms) {
               xs.push_back(at.pos.x);
               ys.push_back(at.pos.y);
               zs.push_back(at.pos.z);
            }
         }
      }
      break;
   }
   if (xs.empty())
      throw std::runtime_error("median_position: no atoms");
   std::sort(xs.begin(), xs.end());
   std::sort(ys.begin(), ys.end());
   std::sort(zs.begin(), zs.end());
   size_t n = xs.size();
   size_t mid = n / 2;
   return clipper::Coord_orth(xs[mid], ys[mid], zs[mid]);
}

int
coot::util::count_cis_peptides(const gemmi::Structure &st) {
   int n_cis = 0;
   for (const auto &model : st.models) {
      for (const auto &chain : model.chains) {
         for (size_t i=0; i+1<chain.residues.size(); i++) {
            const auto &res1 = chain.residues[i];
            const auto &res2 = chain.residues[i+1];
            const gemmi::Atom *ca1 = nullptr, *c1 = nullptr;
            const gemmi::Atom *n2 = nullptr, *ca2 = nullptr;
            for (const auto &at : res1.atoms) {
               if (at.name == "CA") ca1 = &at;
               if (at.name == "C")  c1 = &at;
            }
            for (const auto &at : res2.atoms) {
               if (at.name == "N")  n2 = &at;
               if (at.name == "CA") ca2 = &at;
            }
            if (ca1 && c1 && n2 && ca2) {
               // match mmdb version: check C-N distance < 3.0 A
               clipper::Coord_orth pc(c1->pos.x, c1->pos.y, c1->pos.z);
               clipper::Coord_orth pn(n2->pos.x, n2->pos.y, n2->pos.z);
               double d = std::sqrt((pc - pn).lengthsq());
               if (d < 3.0) {
                  clipper::Coord_orth p1(ca1->pos.x, ca1->pos.y, ca1->pos.z);
                  clipper::Coord_orth p4(ca2->pos.x, ca2->pos.y, ca2->pos.z);
                  double tors = clipper::Coord_orth::torsion(p1, pc, pn, p4);
                  double torsion = clipper::Util::rad2d(tors);
                  double pos_torsion = (torsion > 0.0) ? torsion : 360.0 + torsion;
                  double distortion = std::fabs(180.0 - pos_torsion);
                  if (distortion > 90.0)
                     n_cis++;
               }
            }
         }
      }
      break;
   }
   return n_cis;
}
