/*
 * coot-utils/pdbqt.hh
 *
 * Write PDBQT (AutoDock/Vina input) files. PDBQT is PDB with two extra pieces of
 * per-atom data - a partial charge (columns 71-76) and an AutoDock atom type
 * (columns 78-79) - and, for flexible ligands, a torsion tree.
 *
 * This file holds the shared, RDKit-free machinery: the rigid receptor writer and
 * the per-atom formatting/typing helpers. The flexible-ligand writer (which needs
 * RDKit) lives in lidia-core/pdbqt-ligand.hh.
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option) any
 * later version.
 */

#ifndef COOT_UTILS_PDBQT_HH
#define COOT_UTILS_PDBQT_HH

#include <string>
#include <vector>
#include <mmdb2/mmdb_manager.h>
#include "geometry/protein-geometry.hh"
#include "geometry/residue-and-atom-specs.hh"   // residue_spec_t
#include "lidia-core/pdbqt-common.hh"   // writable_atom_t, atom_line, ad_type_from_element

namespace coot {

   namespace pdbqt {

      //! AutoDock Vina scoring information for one docked pose (one model),
      //! recovered from the per-model UDData written by read().
      class pose_score_t {
      public:
         int model_no;      // 1-based model number
         float affinity;    // predicted binding affinity (kcal/mol)
         float rmsd_lb;     // RMSD lower bound to the best mode
         float rmsd_ub;     // RMSD upper bound to the best mode
         float inter;       // intermolecular energy term
         float intra;       // intramolecular energy term
         float unbound;     // unbound-state energy term
         pose_score_t() : model_no(0), affinity(0.0f), rmsd_lb(0.0f), rmsd_ub(0.0f),
                          inter(0.0f), intra(0.0f), unbound(0.0f) {}
      };

      //! Recover the per-model AutoDock Vina scores stored as UDData by read().
      //! @return one entry per model; empty if the molecule carries no Vina scores.
      std::vector<pose_score_t> get_scores(mmdb::Manager *mol);

      //! Write a rigid, whole-molecule (receptor) PDBQT file. Native (no RDKit):
      //! AutoDock types come from the CCP4 energy-library type and hydrogen
      //! positions; charges from the reference table (standard residues) or a
      //! per-energy-type fallback. Non-polar hydrogens are merged; waters dropped.
      //!
      //! @return the number of atoms written (0 on failure)
      int write_receptor(mmdb::Manager *mol, protein_geometry &geom, int imol_enc,
                         const std::string &file_name);

      //! Write a flexible receptor as the AutoDock/Vina pair of files. The rigid
      //! file holds everything except the movable side-chain atoms of the chosen
      //! residues (their backbone N/C/O and amide H stay rigid); the flex file
      //! holds each chosen side chain as a CA-rooted torsion tree. Rotatable bonds
      //! are the dictionary (non-const, non-ring) side-chain torsions - no RDKit.
      //! Waters are dropped. Give the two files to Vina as --receptor and --flex.
      //!
      //! @return the number of flexible side chains written to the flex file
      //!         (0 on failure or if none of flex_residues was usable)
      int write_flexible_receptor(mmdb::Manager *mol, protein_geometry &geom, int imol_enc,
                                  const std::vector<residue_spec_t> &flex_residues,
                                  const std::string &rigid_file_name,
                                  const std::string &flex_file_name);

      //! Convenience selector: the polymer residues that have a rotatable side
      //! chain and lie within `radius` (any atom) of the point (x,y,z). Skips
      //! GLY/ALA/PRO (no useful rotatable side chain) and waters. Feed the result
      //! to write_flexible_receptor().
      std::vector<residue_spec_t>
      flexible_side_chain_residues_near(mmdb::Manager *mol, double x, double y, double z,
                                        double radius);

      //! Read a PDBQT file (e.g. an AutoDock/Vina docking result) into a new
      //! mmdb::Manager. This is a lightweight reader, not a full PDB parser: it
      //! reads ATOM/HETATM records and skips the torsion-tree records. Each MODEL
      //! becomes a separate model. AutoDock atom types (columns 78-79) are mapped
      //! to chemical elements. The "REMARK VINA RESULT" affinity and RMSD values,
      //! and the INTER/INTRA/UNBOUND energy terms, are stored as per-model UDData
      //! reals ("vina_affinity", "vina_rmsd_lb", "vina_rmsd_ub", "vina_inter",
      //! "vina_intra", "vina_unbound").
      //!
      //! @return a new mmdb::Manager (caller owns it), or nullptr on failure
      mmdb::Manager *read(const std::string &file_name);

   }
}

#endif // COOT_UTILS_PDBQT_HH
