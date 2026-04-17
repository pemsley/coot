/* geometry/residue-and-atom-specs-gemmi.cc
 *
 * Copyright 2026 by Paul Emsley
 * Copyright 2026 by Medical Research Council
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

#include "residue-and-atom-specs-gemmi.hh"

namespace {

   std::string icode_to_string(char icode) {
      if (icode != ' ' && icode != '\0')
         return std::string(1, icode);
      return std::string();
   }

   std::string altloc_to_string(char altloc) {
      if (altloc != '\0' && altloc != ' ')
         return std::string(1, altloc);
      return std::string();
   }

} // anon namespace

// ==================== atom_spec_t from gemmi ====================

coot::atom_spec_t
coot::atom_spec_from_gemmi(const gemmi::Chain &chain,
                           const gemmi::Residue &res,
                           const gemmi::Atom &at,
                           int model_number) {

   atom_spec_t spec(chain.name, res.seqid.num.value,
                    icode_to_string(res.seqid.icode),
                    at.name, altloc_to_string(at.altloc));
   spec.model_number = model_number;
   return spec;
}

bool
coot::matches_spec(const atom_spec_t &spec,
                   const gemmi::Chain &chain,
                   const gemmi::Residue &res,
                   const gemmi::Atom &at) {

   if (spec.chain_id != chain.name) return false;
   if (spec.res_no != res.seqid.num.value) return false;
   if (spec.ins_code != icode_to_string(res.seqid.icode)) return false;
   if (spec.atom_name != at.name) return false;
   if (spec.alt_conf != altloc_to_string(at.altloc)) return false;
   return true;
}

const gemmi::Atom *
coot::get_atom(const atom_spec_t &spec, const gemmi::Structure &st) {

   if (st.models.empty()) return nullptr;
   for (const auto &chain : st.models[0].chains) {
      if (chain.name != spec.chain_id) continue;
      for (const auto &res : chain.residues) {
         if (res.seqid.num.value != spec.res_no) continue;
         if (icode_to_string(res.seqid.icode) != spec.ins_code) continue;
         for (const auto &at : res.atoms) {
            if (at.name == spec.atom_name && altloc_to_string(at.altloc) == spec.alt_conf)
               return &at;
         }
      }
   }
   return nullptr;
}

// ==================== residue_spec_t from gemmi ====================

coot::residue_spec_t
coot::residue_spec_from_gemmi(const gemmi::Chain &chain,
                              const gemmi::Residue &res,
                              int model_number) {

   return residue_spec_t(model_number, chain.name, res.seqid.num.value,
                         icode_to_string(res.seqid.icode));
}

const gemmi::Residue *
coot::get_residue(const residue_spec_t &spec, const gemmi::Structure &st) {

   if (st.models.empty()) return nullptr;
   for (const auto &chain : st.models[0].chains) {
      if (chain.name != spec.chain_id) continue;
      for (const auto &res : chain.residues) {
         if (res.seqid.num.value != spec.res_no) continue;
         if (icode_to_string(res.seqid.icode) == spec.ins_code)
            return &res;
      }
   }
   return nullptr;
}

// ==================== link_atoms from gemmi ====================

static coot::atom_spec_t
spec_from_atom_address(const gemmi::AtomAddress &addr, int model_number) {

   coot::atom_spec_t spec(addr.chain_name,
                          addr.res_id.seqid.num.value,
                          icode_to_string(addr.res_id.seqid.icode),
                          addr.atom_name,
                          altloc_to_string(addr.altloc));
   spec.model_number = model_number;
   return spec;
}

std::pair<coot::atom_spec_t, coot::atom_spec_t>
coot::link_atoms_from_gemmi(const gemmi::Connection &con, int model_number) {

   return std::make_pair(spec_from_atom_address(con.partner1, model_number),
                         spec_from_atom_address(con.partner2, model_number));
}
