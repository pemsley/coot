// Author Fei Long
// Author Paul Emsley
// 3 July, 2014

// This is a port of the Acedrg code to coot.

#include "cod-atom-types.hh"

// static
//
// Replicates acedrg's desSortMapKey2 (kernel/utility.cpp:374) for ordering an
// atom's codClass neighbour groups: string length descending, then neighbour
// degree (nNB) descending. The final fall-back on the string keeps the order
// deterministic when length and nNB are equal (acedrg leaves those "equal").
bool
cod::atom_types_t::neighbour_pair_sorter(const std::pair<std::string, int> &a,
                                        const std::pair<std::string, int> &b) {
   if (a.first.length() != b.first.length())
      return a.first.length() > b.first.length();
   if (a.second != b.second)
      return a.second > b.second;
   return (a.first < b.first);
}


// Port of acedrg's strTransSP() (chemPropSet.cpp:1352). The hybridisation
// string is a pure function of the bondingIdx, so it is derived from the
// "acedrg_bondingIdx" property computed by set_acedrg_bonding_indices().
//
// static
std::string
cod::atom_types_t::strTransSP(int tSP) {

   if (tSP == 1) return "SP1";
   if (tSP == 2) return "SP2";
   if (tSP == 3) return "SP3";
   if (tSP > 3)  return "SPD" + std::to_string(tSP);   // SPD5, SPD6, SPD24 ...
   return "SP-NON";                                    // 0 / unset
}


// Port of acedrg's setAtomsBondingAndChiralCenter() bondingIdx logic
// (acedrg/src/kernel/chemPropSet.cpp:125) plus the only reachable bondingIdx
// change of modAtomsBondingAndChiralCenter() with tMode==1 (the non-ring,
// 2-coordinate N geometry rule). Result is stored per atom as the int property
// "acedrg_bondingIdx".
//
// Notes on fidelity:
//  - acedrg's t_len counts non-metal neighbours and t_m_len the metal ones.
//    Metals are out of scope here, so every neighbour is treated as non-metal
//    (t_m_len == 0) and t_len == getDegree().
//  - acedrg's parCharge (partial charge) is 0 in this path, so the O sp2
//    detection condition (parCharge == 0) always holds.
//  - H atoms and terminal halogens have no element branch, so they keep the
//    default bondingIdx 0 (matching libmol, e.g. the "_0_0_0" tails).
//
// static
void
cod::atom_types_t::set_acedrg_bonding_indices(RDKit::ROMol &rdkm) {

   const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
   unsigned int n_atoms = rdkm.getNumAtoms();

   auto element_of = [tbl](const RDKit::Atom *a) -> std::string {
      return tbl->getElementSymbol(a->getAtomicNum());
   };

   auto n_oxy_connect = [&element_of](const RDKit::Atom *a, const RDKit::ROMol &m) -> int {
      int nO = 0;
      RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
      boost::tie(nbrIdx, endNbrs) = m.getAtomNeighbors(a);
      while (nbrIdx != endNbrs) {
         if (element_of(m[*nbrIdx]) == "O") nO++;
         ++nbrIdx;
      }
      return nO;
   };

   std::vector<int> bonding(n_atoms, 0);

   // ---- Pass 1: per-element connectivity rules (bondingIdx only) ----
   for (unsigned int idx=0; idx<n_atoms; idx++) {
      const RDKit::Atom *at = rdkm[idx];
      std::string ele = element_of(at);
      int t_len = at->getDegree();        // t_m_len == 0 (no metals modelled)
      int charge = at->getFormalCharge();
      int &b = bonding[idx];

      if (t_len > 4) {
         b = t_len;
      } else if (ele == "C" || ele == "Si" || ele == "Ge") {
         if (t_len == 4) {
            b = 3;
         } else if (t_len == 3) {
            b = 2;                         // all charge sub-cases collapse to 2
         } else if (t_len == 2) {
            if      (n_oxy_connect(at, rdkm) == 1) b = 1;
            else if (charge == -1)                 b = 2;
            else                                   b = 1;
         }
      } else if (ele == "N" || ele == "As") {
         if      (t_len == 4 || t_len == 3) b = 3;
         else if (t_len == 2)               b = (charge == 1) ? 1 : 2;
         else if (t_len == 1)               b = 1;
      } else if (ele == "B") {
         if      (t_len == 4) b = 3;
         else if (t_len == 3) b = 2;
         else if (t_len == 2) b = (charge == 1) ? 1 : 2;
         else if (t_len == 1) b = 1;
      } else if (ele == "O") {
         unsigned int nconn = at->getDegree();   // acedrg uses raw connAtoms.size()
         if (nconn == 2) {
            b = 3;
         } else if (nconn == 1) {
            const RDKit::Atom *nb = 0;
            RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
            boost::tie(nbrIdx, endNbrs) = rdkm.getAtomNeighbors(at);
            if (nbrIdx != endNbrs) nb = rdkm[*nbrIdx];
            if (nb && nb->getDegree() != 1) b = 2;
            else                            b = 3;
         } else {
            b = 3;
         }
      } else if (ele == "P") {
         if (t_len == 2 || t_len == 3 || t_len == 4 || t_len == 5) b = 3;
      } else if (ele == "S") {
         if      (t_len == 2 || t_len == 3 || t_len == 4) b = 3;
         else if (t_len == 6)                             b = 5;
         else if (t_len == 1)                             b = 3;
      } else if (ele == "Se") {
         if      (t_len == 4 || t_len == 3 || t_len == 2) b = (t_len == 3) ? 2 : 3;
         else if (t_len == 6)                             b = 5;
         else if (t_len == 1)                             b = 3;
      } else if (ele == "Br") {
         if (t_len == 3) b = 3;
      }
      // H, terminal halogens, etc. keep b == 0
   }

   // ---- Pass 2: oxygen sp2 detection (parCharge == 0 always holds here) ----
   for (unsigned int idx=0; idx<n_atoms; idx++) {
      const RDKit::Atom *at = rdkm[idx];
      if (element_of(at) != "O") continue;
      if (at->getDegree() != 2)  continue;
      bool l_sp2 = false;
      RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
      boost::tie(nbrIdx, endNbrs) = rdkm.getAtomNeighbors(at);
      while (nbrIdx != endNbrs) {
         if (bonding[*nbrIdx] == 2) { l_sp2 = true; break; }
         ++nbrIdx;
      }
      if (l_sp2) bonding[idx] = 2;
   }

   // ---- Pass 3: N / S / C sp2 propagation. N and S read a pre-pass snapshot;
   //      C reads the live values (mirroring the acedrg source exactly). ----
   std::vector<int> pre = bonding;
   for (unsigned int idx=0; idx<n_atoms; idx++) {
      const RDKit::Atom *at = rdkm[idx];
      std::string ele = element_of(at);
      int t_len = at->getDegree();
      int charge = at->getFormalCharge();

      if (ele == "N" || ele == "As") {
         if (t_len == 3) {
            if (charge == 0) {
               bool l_sp2 = false;
               RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
               boost::tie(nbrIdx, endNbrs) = rdkm.getAtomNeighbors(at);
               while (nbrIdx != endNbrs) {
                  if (pre[*nbrIdx] == 2 && element_of(rdkm[*nbrIdx]) != "O") { l_sp2 = true; break; }
                  ++nbrIdx;
               }
               bonding[idx] = l_sp2 ? 2 : 3;   // t_m_len == 0 -> sp2 when l_sp2
            } else if (charge == 1) {
               bonding[idx] = 2;
            }
         }
      } else if (ele == "S") {
         if (t_len == 2 && charge == 0) {
            bool l_sp2 = false;
            RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
            boost::tie(nbrIdx, endNbrs) = rdkm.getAtomNeighbors(at);
            while (nbrIdx != endNbrs) {
               if (pre[*nbrIdx] == 2 && element_of(rdkm[*nbrIdx]) != "O") { l_sp2 = true; break; }
               ++nbrIdx;
            }
            if (l_sp2) bonding[idx] = 2;
         }
      } else if (ele == "C") {
         if (t_len == 3 && charge == -1) {
            int sp2count = 0;
            RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
            boost::tie(nbrIdx, endNbrs) = rdkm.getAtomNeighbors(at);
            while (nbrIdx != endNbrs) {
               if (bonding[*nbrIdx] == 2) sp2count++;   // live values, per acedrg
               ++nbrIdx;
            }
            bonding[idx] = (sp2count == 2) ? 2 : 3;
         }
      }
   }

   // ---- modAtomsBondingAndChiralCenter(tMode==1): the only reachable
   //      bondingIdx change is the non-ring, 2-coordinate N geometry rule
   //      (bond angle > 150 deg -> sp(1), else sp2(2)). Needs a conformer. ----
   if (rdkm.getNumConformers() > 0) {
      const RDKit::Conformer &conf = rdkm.getConformer();
      RDKit::RingInfo *ri = rdkm.getRingInfo();
      for (unsigned int idx=0; idx<n_atoms; idx++) {
         const RDKit::Atom *at = rdkm[idx];
         if (element_of(at) != "N") continue;
         if (at->getDegree() != 2)  continue;
         if (ri->numAtomRings(idx) != 0) continue;
         std::vector<unsigned int> nbr;
         RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
         boost::tie(nbrIdx, endNbrs) = rdkm.getAtomNeighbors(at);
         while (nbrIdx != endNbrs) { nbr.push_back(*nbrIdx); ++nbrIdx; }
         if (nbr.size() != 2) continue;
         RDGeom::Point3D pN = conf.getAtomPos(idx);
         RDGeom::Point3D v1 = conf.getAtomPos(nbr[0]) - pN;
         RDGeom::Point3D v2 = conf.getAtomPos(nbr[1]) - pN;
         double denom = v1.length() * v2.length();
         double ang_deg = 0.0;
         if (denom > 1e-9) {
            double cosang = v1.dotProduct(v2) / denom;
            if (cosang >  1.0) cosang =  1.0;
            if (cosang < -1.0) cosang = -1.0;
            ang_deg = std::acos(cosang) * 180.0 / M_PI;
         }
         bonding[idx] = (ang_deg > 150.0) ? 1 : 2;
      }
   }

   for (unsigned int idx=0; idx<n_atoms; idx++)
      rdkm[idx]->setProp("acedrg_bondingIdx", bonding[idx]);
}

