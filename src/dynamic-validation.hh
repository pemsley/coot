/*
 * src/dynamic-validation.hh
 *
 * Copyright 2023 by Medical Research Council
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
#ifndef DYNAMIC_VALIDATION_DIALOG_HH
#define DYNAMIC_VALIDATION_DIALOG_HH

#include "coords/phenix-geo.hh"

void dynamic_validation_internal(int imol, int imol_map);

void overlaps_peptides_cbeta_ramas_and_rotas_internal(int imol);

void update_dynamic_validation();

void update_dynamic_validation_for_molecule(int imol);

void phenix_geo_validation_buttons(int imol,
                                   const coot::phenix_geo::phenix_geometry &pg,
                                   double residual_cutoff);

#endif // DYNAMIC_VALIDATION_DIALOG_HH
