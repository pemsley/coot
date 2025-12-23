/*
 * src/rsr-functions.hh
 *
 * Copyright 2022 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
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
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */
#ifndef RSR_FUNCTIONS_HH
#define RSR_FUNCTIONS_HH

/**
 * @brief Regularize the currently active residue.
 *
 * @details
 * Finds the currently active atom (via graphics_info_t::active_atom_spec()),
 * obtains its residue and active alternate-conformation identifier (altLoc),
 * builds a one-element residue vector and calls the internal
 * regularize_residues_vec operation to regularize the residue coordinates.
 * The function clears the flag that indicates a user-picked residue-range so
 * that the regularization operation is treated as program-initiated.
 *
 * @note No parameters. Uses the application-wide active atom selection.
 *
 * @return void
 */
void regularize_residue();

/**
 * @brief Regularize a contiguous block of up to 7 residues centred on the active residue.
 *
 * @details
 * Finds the active atom and its residue, then collects up to three preceding
 * and three following residues (if present) together with the active residue,
 * preserving the active atom's alternate-conformation identifier (altLoc).
 * Calls regularize_residues_vec on that residue list to regularize the fragment.
 *
 * @note No parameters. Uses the application-wide active atom selection.
 *
 * @return void
 */
void regularize_tandem_3();

/**
 * @brief Regularize residues in a spherical neighbourhood around the active residue.
 *
 * @details
 * Uses the active atom's residue as the centre and finds residues within a
 * fixed radius (6.6 Angstroms). The active residue is included explicitly.
 * The set of residues is converted to mmdb::Residue* and passed to
 * regularize_residues_vec together with the current altLoc.
 *
 * @note No parameters. Uses the application-wide active atom selection.
 *
 * @return void
 */
void regularize_sphere();

/**
 * @brief Regularize a small connected fragment around the active atom using a distance-based tree.
 *
 * @details
 * Obtains the active atom's residue and constructs a simple residue tree
 * (coot::simple_residue_tree) around that residue using a close-distance
 * cutoff (2.0 Angstroms). The returned residues are regularized via
 * regularize_residues_vec. The active atom's altLoc is used to restrict the
 * operation to the selected alternate-conformation where applicable.
 *
 * @note No parameters. Uses the application-wide active atom selection.
 *
 * @return void
 */
void regularize_fragment_active_atom();

/**
 * @brief Regularize all residues in the chain containing the active residue.
 *
 * @details
 * Locates the chain for the active atom's residue and collects all residues
 * in that chain (coot::util::residues_in_chain). The collected residues are
 * passed to regularize_residues_vec for regularization. The active atom's
 * alternate-conformation identifier (altLoc) is forwarded so that only the
 * selected conformation is affected if appropriate.
 *
 * @note No parameters. Uses the application-wide active atom selection.
 *
 * @return void
 */
void regularize_chain();

/**
 * @brief Perform real-space refinement (RSR) on the currently active residue.
 *
 * @details
 * Builds a one-element vector containing the residue specification for the
 * active atom and invokes refine_residues_with_alt_conf() with the active
 * atom's altLoc to perform real-space refinement on that residue only.
 *
 * @note No parameters. Uses the application-wide active atom selection.
 *
 * @return void
 */
void rsr_refine_residue();

/**
 * @brief Perform real-space refinement (RSR) on the chain containing the active residue.
 *
 * @details
 * Collects all residues in the chain of the active residue and converts each
 * to a coot::residue_spec_t. Passes the vector of residue specifications and
 * the active atom's altLoc to refine_residues_with_alt_conf() to refine the
 * entire chain (for the selected alt-conformation if applicable).
 *
 * @note No parameters. Uses the application-wide active atom selection.
 *
 * @return void
 */
void rsr_refine_chain();

/**
 * @brief Perform real-space refinement (RSR) on all residues in the current molecule.
 *
 * @details
 * Gathers all residues in the molecule (coot::util::residues_in_molecule),
 * converts them to coot::residue_spec_t and calls
 * refine_residues_with_alt_conf() with the active atom's altLoc. Intended to
 * refine the whole molecule (the function contains a comment noting this
 * implementation may need reworking, particularly with respect to chain fixes).
 *
 * @note No parameters. Uses the application-wide active atom selection.
 *
 * @return void
 */
void rsr_refine_all_atoms();

/**
 * @brief Perform real-space refinement (RSR) on a block of up to 11 residues centred on the active residue.
 *
 * @details
 * Locates up to five preceding and five following residues around the active
 * residue (if present), constructs a vector of coot::residue_spec_t for them
 * and calls refine_residues_with_alt_conf() with the active atom's altLoc.
 * This is useful for refining a larger local stretch of residues together.
 *
 * @note No parameters. Uses the application-wide active atom selection.
 *
 * @return void
 */
void rsr_refine_tandem_5();

/**
 * @brief Perform real-space refinement (RSR) on a block of up to 7 residues centred on the active residue.
 *
 * @details
 * Similar to rsr_refine_tandem_5 but collects up to three preceding and three
 * following residues around the active residue, converts them to
 * coot::residue_spec_t and invokes refine_residues_with_alt_conf() using the
 * active atom's altLoc.
 *
 * @note No parameters. Uses the application-wide active atom selection.
 *
 * @return void
 */
void rsr_refine_tandem_3();

/**
 * @brief Perform real-space refinement (RSR) on residues within a larger spherical neighbourhood (radius 6.6 Å).
 *
 * @details
 * Finds residues near the active residue using a radius of 6.6 Å (via
 * graphics_info_t::molecules[imol].residues_near_residue), includes the
 * active residue, and calls refine_residues_with_alt_conf() with the active
 * atom's altLoc to refine that neighbourhood.
 *
 * @note No parameters. Uses the application-wide active atom selection.
 *
 * @return void
 */
void rsr_sphere_refine_plus();

/**
 * @brief Perform real-space refinement (RSR) on residues within a spherical neighbourhood (radius 4.3 Å).
 *
 * @details
 * As rsr_sphere_refine_plus but uses a smaller radius of 4.3 Å to collect a
 * tighter neighbourhood around the active residue before calling
 * refine_residues_with_alt_conf() with the active atom's altLoc.
 *
 * @note No parameters. Uses the application-wide active atom selection.
 *
 * @return void
 */
void rsr_sphere_refine();

/**
 * @brief Perform real-space refinement (RSR) on a distance-based fragment around the active residue.
 *
 * @details
 * Builds a fragment by calling coot::simple_residue_tree() with a close
 * distance cutoff of 2.0 Å around the active residue, converts the returned
 * mmdb::Residue* list to coot::residue_spec_t and invokes
 * refine_residues_with_alt_conf() with the active atom's altLoc.
 *
 * @note No parameters. Uses the application-wide active atom selection.
 *
 * @return void
 */
void rsr_refine_fragment_active_residue();

#endif // RSR_FUNCTIONS_HH
