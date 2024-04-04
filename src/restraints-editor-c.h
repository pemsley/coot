/*
 * src/restraints-editor-c.h
 *
 * Copyright 2008 by University of York
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

#ifndef RESTRAINTS_EDITOR_C_H
#define RESTRAINTS_EDITOR_C_H

#ifndef BEGIN_C_DECLS
#ifdef __cplusplus
#define BEGIN_C_DECLS extern "C" {
#define END_C_DECLS }
#else
#define BEGIN_C_DECLS
#define END_C_DECLS
#endif
#endif

BEGIN_C_DECLS

void apply_restraint_by_widget(GtkWidget *w);
void restraints_editor_save_restraint_by_widget(GtkWidget *w);
void restraints_editor_delete_restraint_by_widget(GtkWidget *w);
void restraints_editor_add_restraint_by_widget(GtkWidget *w);
void save_monomer_restraints_by_widget(GtkDialog *chooser);

END_C_DECLS

#endif /* RESTRAINTS_EDITOR_C_H */
