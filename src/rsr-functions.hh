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

void regularize_residue();
void regularize_tandem_3();
void regularize_sphere();
void regularize_fragment_active_atom();
void regularize_chain();
void rsr_refine_residue();
void rsr_refine_chain();
void rsr_refine_all_atoms();
void rsr_refine_tandem_5();
void rsr_refine_tandem_3();
void rsr_sphere_refine_plus();
void rsr_sphere_refine();
void rsr_refine_fragment_active_residue();

#endif // RSR_FUNCTIONS_HH
