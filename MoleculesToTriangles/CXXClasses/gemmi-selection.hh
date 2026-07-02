/*
 * MoleculesToTriangles/CXXClasses/gemmi-selection.hh
 *
 * gemmi-native replacement for the mmdb-based CompoundSelection.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef MOLECULESTOTRIANGLES_GEMMI_SELECTION_HH
#define MOLECULESTOTRIANGLES_GEMMI_SELECTION_HH

#include <string>
#include <vector>
#include <memory>

#include <gemmi/model.hpp>   // Model, Chain, Residue, Atom, CRA
#include <gemmi/select.hpp>  // gemmi::Selection

namespace coot {
   namespace m2t {

      // How a clause combines with the running result, mirroring mmdb's SKEY_*.
      enum class combine_t { NEW, AND, OR };

      // One clause of a selection. Each leaf applies its own `invert` itself, so
      // containers can combine raw boolean results uniformly.
      class selection_leaf_t {
      public:
         bool invert = false;
         virtual ~selection_leaf_t() = default;
         virtual bool matches(const gemmi::CRA &cra) const = 0;
         virtual std::string format() const = 0;
      };

      // A CID / named-subset clause: thin wrapper over a gemmi::Selection.
      // gemmi::Selection parses the mmdb-style CID strings M2T uses verbatim
      // (verified: see MoleculesToTriangles/gemmi-migration-plan.md §4A).
      class cid_leaf_t : public selection_leaf_t {
         gemmi::Selection sel;
         std::string cid;
      public:
         explicit cid_leaf_t(const std::string &cid_in);
         bool matches(const gemmi::CRA &cra) const override;
         std::string format() const override { return cid; }
      };

      // A secondary-structure clause (SSE_Helix / SSE_Strand / SSE_None).
      //
      // TODO GEMMI-SSE: secondary structure is not yet computed gemmi-side
      // (mmdb's model->CalcSecStructure / residue->SSE was removed by decision).
      // Until a replacement (gemmi dssp.hpp or a CA-geometry calc) is wired in,
      // this leaf matches nothing. drawRibbon has the corresponding stub.
      class sse_leaf_t : public selection_leaf_t {
         int sse_type; // mirrors mmdb SSE codes (see secondary_types())
         std::string name;
      public:
         explicit sse_leaf_t(const std::string &name_in);
         bool matches(const gemmi::CRA &cra) const override; // returns false (stub)
         std::string format() const override { return name; }
      };

      // A compound selection: an ordered list of (combine-rule, leaf) pairs, parsed
      // from the same '& | ! { }' mini-language as the old CompoundSelection.
      class compound_selection_t : public selection_leaf_t {
         std::vector<std::pair<combine_t, std::shared_ptr<selection_leaf_t> > > pairs;
         std::string selection_text;
         void parse(const std::string &text);
      public:
         explicit compound_selection_t(const std::string &selection_text_in);
         bool matches(const gemmi::CRA &cra) const override;
         std::string format() const override;

         // Membership over a whole model, in iteration order. Also stamps each
         // matched/unmatched atom's serial with its iteration index so callers
         // can map an atom back to its bitset slot.
         std::vector<bool> matching_atoms(gemmi::Model &model) const;
         int count_matching_atoms(gemmi::Model &model) const;
      };

      // The named subset / secondary-structure tables (MAIN, SIDE, WATER, ...).
      // Exposed for testing and reuse.
      const std::vector<std::pair<std::string, std::string> > &subset_types();
      const std::vector<std::pair<std::string, int> > &secondary_types();
   }
}

#endif // MOLECULESTOTRIANGLES_GEMMI_SELECTION_HH
