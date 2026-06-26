/*
 * MoleculesToTriangles/CXXClasses/gemmi-sse.cc
 *
 * gemmi secondary-structure assignment. Stamps Residue::flag with 'H'/'E'/'L'.
 *
 * NOTE: gemmi's DSSP (dssp.hpp) is, in the current gemmi, a work in progress -
 * DsspCalculator::calculate_secondary_structure() has its body commented out and
 * returns "" (only the H-bond pass is implemented). So we try DSSP first (this
 * will start working automatically when gemmi's DSSP is completed, or when a port
 * of mmdb's CalcSecStructure is dropped in here), and fall back to the HELIX/SHEET
 * header records (the equivalent of the old secondary_structure_header_to_residue_sse).
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#include "gemmi-sse.hh"
#include <gemmi/dssp.hpp>
#include <gemmi/neighbor.hpp>
#include <gemmi/topo.hpp>
#include <gemmi/polyheur.hpp>  // setup_entities, is_polypeptide

namespace {
   inline bool seqid_num(const gemmi::SeqId &s, int &out) {
      if (!s.num) return false;
      out = *s.num;
      return true;
   }
   // stamp residues of chain `chain_name` with seqnum in [lo,hi] with `code`.
   void stamp_range(gemmi::Model &model, const std::string &chain_name, int lo, int hi, char code) {
      if (lo > hi) std::swap(lo, hi);
      for (gemmi::Chain &chain : model.chains) {
         if (chain.name != chain_name) continue;
         for (gemmi::Residue &res : chain.residues) {
            int n;
            if (!seqid_num(res.seqid, n)) continue;
            if (n >= lo && n <= hi) res.flag = code;
         }
      }
   }
}

void coot::m2t::assign_secondary_structure(gemmi::Structure &st) {

   if (st.models.empty()) return;
   gemmi::setup_entities(st);
   gemmi::Model &model = st.models[0];

   gemmi::NeighborSearch ns(model, st.cell, 5.0);
   ns.populate(false);
   gemmi::DsspOptions opts;

   bool any_dssp = false;
   for (gemmi::Chain &chain : model.chains) {
      gemmi::ResidueSpan polymer = chain.get_polymer();
      if (polymer.size() == 0) continue;
      const gemmi::Entity *ent = st.get_entity_of(polymer);
      gemmi::Topo::ChainInfo ci(polymer, chain, ent);
      if (!ci.polymer) continue;
      if (!gemmi::is_polypeptide(ci.polymer_type)) continue;

      // peptide residues default to loop
      for (gemmi::Topo::ResInfo &ri : ci.res_infos)
         if (ri.res) ri.res->flag = 'L';

      std::string ss = gemmi::calculate_dssp(ns, ci, opts); // currently returns "" (gemmi WIP)
      if (!ss.empty()) {
         for (size_t i = 0; i < ci.res_infos.size() && i < ss.size(); i++) {
            gemmi::Residue *res = ci.res_infos[i].res;
            if (!res) continue;
            char c = ss[i];
            if (c == 'H' || c == 'G' || c == 'I')   res->flag = 'H';
            else if (c == 'E' || c == 'B')           res->flag = 'E';
            else                                     res->flag = 'L';
         }
         any_dssp = true;
      }
   }

   // Fallback: HELIX / SHEET header records (works today; DSSP above is a gemmi stub).
   if (!any_dssp) {
      for (const gemmi::Helix &h : st.helices) {
         int lo, hi;
         if (seqid_num(h.start.res_id.seqid, lo) && seqid_num(h.end.res_id.seqid, hi))
            stamp_range(model, h.start.chain_name, lo, hi, 'H');
      }
      for (const gemmi::Sheet &sheet : st.sheets) {
         for (const gemmi::Sheet::Strand &strand : sheet.strands) {
            int lo, hi;
            if (seqid_num(strand.start.res_id.seqid, lo) && seqid_num(strand.end.res_id.seqid, hi))
               stamp_range(model, strand.start.chain_name, lo, hi, 'E');
         }
      }
   }
}
