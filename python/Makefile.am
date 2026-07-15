# Makefile.am
# 
# Copyright 2005, 2006, 2007 The University of York
# Copyright 2005, 2006, 2007 Bernhard Lohkamp
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA

dist_pkgpython_PYTHON = \
                        coot_utils.py        get_ebi.py      mutate.py        \
   background_demo.py   hello.py             raster3d.py                      \
   brute_lsqman.py                                           refmac.py        \
   coot_gui.py          filter.py            libcheck.py                      \
   coot_load_modules.py                      local_code.py                    \
   gap.py               remote_control.py    povray.py       coot_lsq.py      \
   shelx.py             group_settings.py    coot_fitting.py americanisms.py  \
   tips.py              tips_gui.py          generic_objects.py  cns2coot.py  \
   extensions.py        coot_ncs.py          create_server.py                 \
   python_coot_docs.py  clear_backup.py      bones_reader.py                  \
   shelx_extensions.py  parse_pisa_xml.py    		                      \
   check_for_updates.py coot_generator_3d_import.py     gui_hole.py     gui_prosmart.py  \
   get_recent_pdbe.py   jligand.py	     	 jligand_gui.py	              \
   coot_redefine_functions.py user_define_restraints.py template_key_bindings.py   \
   gui_add_linked_cho.py coot_ligand_check.py run_mogul.py                     \
   coot_contact_score_isolated_ligand.py          gui_contact_score_isolated_ligand.py \
   dictionary_generators.py cho_restraints_from_models.py interactive_nudge_residues.py \
   pdbe_validation_data.py  pathology_plots.py ligand_validation_sliders.py \
   enhanced_ligand.py	gui_ligand_sliders.py coot_sharpen_blur.py    \
   coot_acedrg_link.py coot_toolbuttons.py rcrane_loader.py add_linked_cho.py coot_find_baddies.py \
   coot_chemical_feature_clusters.py coot_res_spec_utils.py coot_output_to_window.py \
   populate_python_menus.py coot_hole.py refine.py what_check.py coot_servalcat_outliers.py \
   coot_testing_utils.py point_group_symmetry.py gphl_screen.py \
   coot_command_interface.py \
   $(TEXI_FILE)

# Natural-language command interface package. The nobase_ prefix preserves the
# coot_commands/ directory layout under site-packages so it stays importable.
# Add each new command module here (e.g. coot_commands/commands/maps.py).
nobase_dist_pkgpython_PYTHON = \
   coot_commands/__init__.py                  \
   coot_commands/registry.py                  \
   coot_commands/completion.py                \
   coot_commands/speech.py                    \
   coot_commands/types.py                     \
   coot_commands/gui.py                       \
   coot_commands/docs.py                      \
   coot_commands/tools.py                     \
   coot_commands/context_tools.py             \
   coot_commands/agent.py                     \
   coot_commands/retrieval.py                 \
   coot_commands/socket_client.py             \
   coot_commands/agent_serve.py               \
   coot_commands/new.py                       \
   coot_commands/try.py                       \
   coot_commands/commands/__init__.py         \
   coot_commands/commands/delete.py           \
   coot_commands/commands/display.py          \
   coot_commands/commands/fetch.py            \
   coot_commands/commands/view.py             \
   coot_commands/commands/maps.py             \
   coot_commands/commands/navigation.py       \
   coot_commands/commands/refine.py           \
   coot_commands/commands/representation.py   \
   coot_commands/commands/labels.py           \
   coot_commands/commands/model_edit.py       \
   coot_commands/commands/build.py            \
   coot_commands/commands/ligand.py           \
   coot_commands/commands/validation.py       \
   coot_commands/commands/session.py          \
   coot_commands/commands/state.py            \
   coot_commands/commands/models.py           \
   coot_commands/commands/settings.py         \
   coot_commands/commands/help.py

# VTE REPL helper script — installed to pkgdatadir/python/ (not site-packages)
# because it is executed as a subprocess, not imported as a module.
vtepythondir = $(pkgdatadir)/python
dist_vtepython_DATA = coot_vte_repl.py

