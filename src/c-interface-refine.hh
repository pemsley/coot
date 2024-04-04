/*
 * src/c-interface-refine.hh
 *
 * Copyright 2013 by Medical Research Council
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

#include <string>
#include <vector>
#include "geometry/residue-and-atom-specs.hh"

/*! \file
  \brief Coot Scripting Interface - Refinement utilities
*/

void    add_initial_position_restraints(int imol, const std::vector<coot::residue_spec_t> &residue_specs, double weight);
// this removed all initial position restraints, not just those listed.
void remove_initial_position_restraints(int imol, const std::vector<coot::residue_spec_t> &residue_specs);

//! \name  More Refinement Functions
//! \{

//! \brief use unimodal ring torsion restraints (e.g. for carbohydrate pyranose)
//         uses built-in list of of torsions for specific residue types
void use_unimodal_ring_torsion_restraints(const std::string &res_name);

#ifdef USE_PYTHON
//! \brief use unimodal ring torsion restraints (e.g. for carbohydrate pyranose)
//
//         allow user definition of torsions for given residue
//  @var{torsions_info_list} is a list of item that are of the form
//  @var{[atom_name_1, atom_name_2, atom_name_3, atom_name_4, double torsion_1234]}
void use_unimodal_ring_torsion_restraints_for_residue(const std::string &res_name, PyObject *torsions_info_list);
#endif

//! \brief set the Geman-McClure distance alpha value (weight)
void set_refinement_geman_mcclure_alpha(float alpha);

//! \brief get the Geman-McClure distance alpha value (weight)
float get_refinement_geman_mcclure_alpha();

//! \brief set the Lennard Jones epsilon parameter
void set_refinement_lennard_jones_epsilon(float epsilon);

//! \brief set the log cosh scale factor for target position restraints
void set_log_cosh_target_distance_scale_factor(float sf);

#ifdef USE_GUILE
//! \brief Apply crankshaft peptide rotation optimization to the specified residue
void crankshaft_peptide_rotation_optimization_scm(int imol, SCM residue_spec_smc);
#endif

#ifdef USE_PYTHON
//! \brief Apply crankshaft peptide rotation optimization to the specified residue
void crankshaft_peptide_rotation_optimization_py(int imol, PyObject *residue_spec_py);
#endif

void convert_dictionary_planes_to_improper_dihedrals();


//! \}
