/* layla/similar_ligands.cpp
 *
 * Copyright 2026 by Medical Research Council
 * Author: Paul Emsley
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include <fstream>
#include <thread>
#include <mutex>
#include <vector>
#include <cctype>
#include <filesystem>
#include <sys/stat.h>

#include "similar_ligands.hpp"
#include "network_utils.hpp"

#include "utils/coot-utils.hh"
#include "geometry/protein-geometry.hh"
#include "lidia-core/cod-types-similarity.hh"
#include "lidia-core/rdkit-interface.hh"   // rdkit_mol, rdkit_mol_with_2d_depiction

#include <GraphMol/MolOps.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/FMCS/FMCS.h>
#include <boost/range/iterator_range.hpp>

#if RDKIT_HAS_CAIRO_SUPPORT
#include <GraphMol/MolDraw2D/MolDraw2DCairo.h>
#include <GraphMol/Depictor/RDDepictor.h>  // compute2DCoords fallback
#endif

namespace {

   namespace fs = std::filesystem;

   // ---- cache directories ------------------------------------------------

   std::string home_dir() {
      const char *h = getenv("HOME");
      return h ? std::string(h) : std::string(".");
   }

   std::string monomer_images_dir() {
      std::string d = home_dir() + "/.cache/Coot/monomer-images";
      std::error_code ec;
      fs::create_directories(d, ec);
      return d;
   }

   std::string monomer_cif_cache_dir() {
      std::string d = home_dir() + "/.cache/Coot/monomer-images/cif";
      std::error_code ec;
      fs::create_directories(d, ec);
      return d;
   }

   bool file_exists_and_nonempty(const std::string &path, long min_bytes = 1) {
      struct stat st;
      if (stat(path.c_str(), &st) != 0) return false;
      return st.st_size >= min_bytes;
   }

   // ---- ensure a monomer dictionary is available in geom -----------------
   //
   // Order: local monomer library (try_dynamic_add) -> GitHub monomers -> PDBe CCD.
   // Returns true if geom ends up with restraints for comp_id.
   //
   // The whole fetch/draw path is serialized by the caller's mutex, so the
   // (non-thread-safe) protein_geometry is only touched by one worker at a time.
   bool ensure_monomer_in_geom(const std::string &comp_id, coot::protein_geometry *geom_p) {

      static int read_number = 50; // arbitrary, kept away from the app's own numbers

      auto have_it = [&]() {
         return geom_p->get_monomer_restraints(comp_id, coot::protein_geometry::IMOL_ENC_ANY).first;
      };

      if (have_it()) return true;

      // local monomer library
      geom_p->try_dynamic_add(comp_id, read_number++);
      if (have_it()) return true;

      std::string lc(1, static_cast<char>(std::tolower(static_cast<unsigned char>(comp_id[0]))));
      std::string cif_path = monomer_cif_cache_dir() + "/" + comp_id + ".cif";

      // GitHub monomers (raw)
      if (! file_exists_and_nonempty(cif_path, 200)) {
         std::string url = "https://raw.githubusercontent.com/MonomerLibrary/monomers/master/"
                           + lc + "/" + comp_id + ".cif";
         coot::layla::get_url(url, cif_path, 200);
      }
      if (file_exists_and_nonempty(cif_path, 200)) {
         geom_p->init_refmac_mon_lib(cif_path, read_number++, coot::protein_geometry::IMOL_ENC_ANY);
         if (have_it()) return true;
      }

      // PDBe CCD
      {
         std::string url = "https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/" + comp_id + ".cif";
         if (coot::layla::get_url(url, cif_path, 200) == 0) {
            geom_p->init_refmac_mon_lib(cif_path, read_number++, coot::protein_geometry::IMOL_ENC_ANY);
            if (have_it()) return true;
         }
      }

      return false;
   }

   // Draw the monomer to a PNG in the cache. Returns the png path, or "" on failure.
   std::string draw_monomer_png(const std::string &comp_id, coot::protein_geometry *geom_p) {

      std::string png_path = monomer_images_dir() + "/" + comp_id + ".png";

#if RDKIT_HAS_CAIRO_SUPPORT
      std::pair<bool, coot::dictionary_residue_restraints_t> r =
         geom_p->get_monomer_restraints(comp_id, coot::protein_geometry::IMOL_ENC_ANY);
      if (! r.first)
         return std::string();

      try {
         std::pair<int, RDKit::RWMol> mol_pair = coot::rdkit_mol_with_2d_depiction(r.second);
         RDKit::RWMol mol = mol_pair.second;
         int conf_id = mol_pair.first;
         if (conf_id < 0) {
            // No stored 2D depiction in the dictionary - compute one (as
            // dict-to-image-grid.py does with Compute2DCoords). The mol from
            // rdkit_mol_with_2d_depiction() still has hydrogens in this branch.
            try { RDKit::MolOps::removeHs(mol); } catch (...) { }
            conf_id = static_cast<int>(RDDepict::compute2DCoords(mol));
         }
         if (conf_id < 0)
            return std::string();
         RDKit::MolDraw2DCairo drawer(300, 300);
         drawer.drawMolecule(mol, nullptr, nullptr, nullptr, conf_id);
         drawer.finishDrawing();
         std::string dt = drawer.getDrawingText();
         std::ofstream f(png_path.c_str(), std::ios::binary);
         f << dt;
         f.close();
         if (file_exists_and_nonempty(png_path, 100))
            return png_path;
      }
      catch (const std::exception &e) {
         std::cout << "WARNING:: draw_monomer_png() " << comp_id << " " << e.what() << std::endl;
      }
#endif
      return std::string();
   }

   // fetch (if needed) + draw, returning a usable png path or "".
   // Serialized across worker threads to keep protein_geometry single-threaded.
   std::string ensure_and_draw(const std::string &comp_id, coot::protein_geometry *geom_p) {

      std::string png_path = monomer_images_dir() + "/" + comp_id + ".png";
      if (file_exists_and_nonempty(png_path, 100))
         return png_path; // cache hit - no lock needed

      static std::mutex geom_mutex;
      std::lock_guard<std::mutex> lock(geom_mutex);

      // another thread may have produced it while we waited
      if (file_exists_and_nonempty(png_path, 100))
         return png_path;

      if (! ensure_monomer_in_geom(comp_id, geom_p))
         return std::string();
      return draw_monomer_png(comp_id, geom_p);
   }

   // ---- results window + lazy thumbnail tiles ----------------------------

   struct tile_task_t {
      std::string comp_id;
      coot::protein_geometry *geom_p;
      GtkWidget *child;     // the GtkFlowBoxChild
      GtkPicture *picture;
      bool started;
   };

   struct results_window_t {
      GtkWidget *window;                 // the results window (dialog parent)
      GtkWidget *scrolled;
      CootLigandEditorCanvas *canvas;    // the sketch's canvas
      unsigned int mol_idx;              // index of the sketched molecule
      coot::protein_geometry *geom_p;
      std::vector<tile_task_t *> tiles;
   };

   void results_window_free(gpointer data) {
      results_window_t *rw = static_cast<results_window_t *>(data);
      for (tile_task_t *t : rw->tiles) delete t;
      delete rw;
   }

   struct draw_done_t {
      GtkPicture *picture; // holds a ref
      std::string png_path;
   };

   gboolean set_picture_idle(gpointer data) {
      draw_done_t *d = static_cast<draw_done_t *>(data);
      if (! d->png_path.empty())
         gtk_picture_set_filename(d->picture, d->png_path.c_str());
      g_object_unref(d->picture);
      delete d;
      return G_SOURCE_REMOVE;
   }

   void start_tile(tile_task_t *tile) {
      if (tile->started) return;
      tile->started = true;

      GtkPicture *pic = tile->picture;
      g_object_ref(pic); // keep alive across the worker thread
      std::string comp_id = tile->comp_id;
      coot::protein_geometry *geom_p = tile->geom_p;

      std::thread([pic, comp_id, geom_p]() {
         std::string png = ensure_and_draw(comp_id, geom_p);
         draw_done_t *d = new draw_done_t{pic, png};
         g_idle_add(set_picture_idle, d);
      }).detach();
   }

   // Start loading any tile whose allocation intersects the visible viewport
   // (plus a prefetch margin). This is what makes the grid lazy.
   void load_visible_tiles(results_window_t *rw) {
      double vh = gtk_widget_get_height(rw->scrolled);
      if (vh <= 0) vh = 600.0; // not yet allocated - be generous
      const double margin = vh; // prefetch one viewport ahead

      for (tile_task_t *tile : rw->tiles) {
         if (tile->started) continue;
         graphene_rect_t bounds;
         if (gtk_widget_compute_bounds(tile->child, rw->scrolled, &bounds)) {
            double top = bounds.origin.y;
            double bot = bounds.origin.y + bounds.size.height;
            if (bot >= -margin && top <= vh + margin)
               start_tile(tile);
         } else {
            // could not compute bounds yet - load it to be safe
            start_tile(tile);
         }
      }
   }

   gboolean load_visible_idle(gpointer data) {
      load_visible_tiles(static_cast<results_window_t *>(data));
      return G_SOURCE_REMOVE;
   }

   void on_vadjustment_value_changed(GtkAdjustment * /*adj*/, gpointer data) {
      load_visible_tiles(static_cast<results_window_t *>(data));
   }

   void info_dialog(GtkWindow *parent, const std::string &message) {
      GtkWidget *d = gtk_message_dialog_new(parent,
                                            GTK_DIALOG_MODAL,
                                            GTK_MESSAGE_INFO, GTK_BUTTONS_CLOSE,
                                            "%s", message.c_str());
      g_signal_connect(d, "response", G_CALLBACK(+[](GtkDialog *dd, int, gpointer){
         gtk_window_destroy(GTK_WINDOW(dd));
      }), nullptr);
      gtk_window_present(GTK_WINDOW(d));
   }

   std::string strip_spaces(const std::string &s) {
      std::size_t a = s.find_first_not_of(' ');
      std::size_t b = s.find_last_not_of(' ');
      return (a == std::string::npos) ? std::string() : s.substr(a, b - a + 1);
   }

   // Match the chosen reference monomer to the sketched molecule and write the
   // reference's dictionary atom names onto the sketch, then switch the canvas to
   // the "Atom Names" display mode so the user sees them.
   void assign_atom_names_from_monomer(GtkWindow *parent,
                                       CootLigandEditorCanvas *canvas,
                                       unsigned int mol_idx,
                                       coot::protein_geometry *geom_p,
                                       const std::string &comp_id) {

      if (! ensure_monomer_in_geom(comp_id, geom_p)) {
         info_dialog(parent, "Could not fetch the dictionary for " + comp_id + ".");
         return;
      }
      std::pair<bool, coot::dictionary_residue_restraints_t> r =
         geom_p->get_monomer_restraints(comp_id, coot::protein_geometry::IMOL_ENC_ANY);
      if (! r.first) {
         info_dialog(parent, "No restraints found for " + comp_id + ".");
         return;
      }

      // Reference: the dictionary molecule, WITH hydrogens, carrying the
      // dictionary atom names (heavy and H). We keep the hydrogens so we can
      // read the names of each heavy atom's bonded hydrogens.
      RDKit::RWMol ref;
      try {
         ref = coot::rdkit_mol(r.second);   // sets the "name" property per atom
         coot::rdkit_mol_sanitize(ref);
      }
      catch (const std::exception &e) {
         info_dialog(parent, std::string("Could not build a molecule for ") + comp_id + ": " + e.what());
         return;
      }

      // Query: the sketched molecule. Keep its atom indexing so names can be
      // written back by index; sketches carry implicit hydrogens.
      const RDKit::ROMol *q_ro = coot_ligand_editor_canvas_get_rdkit_molecule(canvas, mol_idx);
      if (! q_ro) {
         info_dialog(parent, "There is no molecule to name.");
         return;
      }
      RDKit::RWMol query(*q_ro);
      try { coot::rdkit_mol_sanitize(query); } catch (...) { /* matching may still work */ }

      // Map query atom index -> reference atom index.
      std::map<int, int> q_to_r;

      RDKit::MatchVectType match; // pairs of (queryAtomIdx, molAtomIdx)
      if (RDKit::SubstructMatch(ref, query, match)) {
         // the whole sketch is a substructure of the reference
         for (const auto &p : match) q_to_r[p.first] = p.second;
      } else {
         // partial overlap: use the maximum common substructure
         std::vector<RDKit::ROMOL_SPTR> mols;
         mols.push_back(RDKit::ROMOL_SPTR(new RDKit::ROMol(query)));
         mols.push_back(RDKit::ROMOL_SPTR(new RDKit::ROMol(ref)));
         RDKit::MCSParameters params;
         params.BondCompareParameters.RingMatchesRingOnly = true;
         params.BondCompareParameters.CompleteRingsOnly   = true;
         RDKit::MCSResult mcs = RDKit::findMCS(mols, &params);
         if (mcs.QueryMol) {
            RDKit::MatchVectType qm, rm; // pattern -> query, pattern -> ref
            if (RDKit::SubstructMatch(query, *mcs.QueryMol, qm) &&
                RDKit::SubstructMatch(ref,   *mcs.QueryMol, rm)) {
               std::map<int, int> pat_to_r;
               for (const auto &p : rm) pat_to_r[p.first] = p.second;
               for (const auto &p : qm) {
                  std::map<int, int>::const_iterator it = pat_to_r.find(p.first);
                  if (it != pat_to_r.end()) q_to_r[p.second] = it->second;
               }
            }
         }
      }

      if (q_to_r.empty()) {
         info_dialog(parent, "Could not match " + comp_id + " to your molecule.");
         return;
      }

      // Build, indexed by query atom index: the visible heavy-atom name, and the
      // space-separated names of that heavy atom's hydrogens in the reference.
      // The sketch has no hydrogen atoms, so their names are stashed on the heavy
      // atom they belong to (the "hydrogen_names" property) for acedrg to use.
      std::vector<std::string> names(query.getNumAtoms());
      std::vector<std::string> hydrogen_names(query.getNumAtoms());
      int n_named = 0;
      int n_h_named = 0;
      for (const auto &qr : q_to_r) {
         const RDKit::Atom *ra = ref.getAtomWithIdx(qr.second);
         std::string nm;
         if (ra->getPropIfPresent("name", nm)) {
            names[qr.first] = strip_spaces(nm);
            if (! names[qr.first].empty()) n_named++;
         }
         // gather the reference atom's bonded hydrogen names
         std::string h_join;
         for (const auto &bnd : boost::make_iterator_range(ref.getAtomBonds(ra))) {
            const RDKit::Atom *other = ref[bnd]->getOtherAtom(ra);
            if (other->getSymbol() == "H") {
               std::string hn;
               if (other->getPropIfPresent("name", hn)) {
                  std::string s = strip_spaces(hn);
                  if (! s.empty()) {
                     if (! h_join.empty()) h_join += " ";
                     h_join += s;
                  }
               }
            }
         }
         if (! h_join.empty()) {
            hydrogen_names[qr.first] = h_join;
            n_h_named++;
         }
      }

      coot_ligand_editor_canvas_set_atom_names(canvas, mol_idx, names, hydrogen_names);
      coot_ligand_editor_canvas_set_display_mode(canvas,
                                                 coot::ligand_editor_canvas::DisplayMode::AtomNames);

      info_dialog(parent, "Assigned " + std::to_string(n_named) + " atom name(s) (and "
                  + std::to_string(n_h_named) + " with hydrogen names) from " + comp_id + ".");
   }

   struct assign_ctx_t {
      results_window_t *rw;
      std::string comp_id;
   };

   void on_assign_clicked(GtkButton *btn, gpointer) {
      assign_ctx_t *ctx = static_cast<assign_ctx_t *>(g_object_get_data(G_OBJECT(btn), "assign_ctx"));
      if (! ctx) return;
      results_window_t *rw = ctx->rw;
      assign_atom_names_from_monomer(GTK_WINDOW(rw->window), rw->canvas, rw->mol_idx,
                                     rw->geom_p, ctx->comp_id);
   }
}

std::string
coot::layla::cod_types_db_path() {
   const char *e = getenv("COOT_COD_TYPES_DB");
   if (e)
      return std::string(e);
   return coot::package_data_dir() + "/cod-types-db-normalized.bin";
}

void
coot::layla::search_for_similar_ligands(GtkWindow *parent,
                                        CootLigandEditorCanvas *canvas,
                                        unsigned int mol_idx,
                                        coot::protein_geometry *geom_p,
                                        const RDKit::ROMol &query_mol,
                                        double min_score) {

   // 1. Type the query molecule.
   std::set<std::string> query_types;
   try {
      query_types = cod::normalized_types_for_rdkit_mol(query_mol);
   }
   catch (const std::exception &e) {
      info_dialog(parent, std::string("Could not compute atom types for this molecule:\n") + e.what());
      return;
   }
   if (query_types.empty()) {
      info_dialog(parent, "Could not compute atom types for this molecule.");
      return;
   }

   // 2. Load the database (once).
   static cod::cod_types_db_t db;
   static bool db_loaded = false;
   if (! db_loaded) {
      std::string path = cod_types_db_path();
      try {
         db = cod::load_cod_types_db_binary(path);
         db_loaded = true;
      }
      catch (const std::exception &e) {
         info_dialog(parent, std::string("Could not load the similarity database:\n") + e.what());
         return;
      }
   }

   // 3. Search.
   std::vector<cod::cod_match_t> matches = cod::search(query_types, db, min_score);
   if (matches.empty()) {
      info_dialog(parent, "No similar monomers found (score > "
                  + std::to_string(min_score) + ").");
      return;
   }

   // 4. Build the results window.
   GtkWidget *win = gtk_window_new();
   gtk_window_set_title(GTK_WINDOW(win), "Coot Layla: Similar Monomers");
   if (parent)
      gtk_window_set_transient_for(GTK_WINDOW(win), parent);
   gtk_window_set_default_size(GTK_WINDOW(win), 720, 700);

   GtkWidget *scrolled = gtk_scrolled_window_new();
   gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrolled),
                                  GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);
   gtk_window_set_child(GTK_WINDOW(win), scrolled);

   GtkWidget *flowbox = gtk_flow_box_new();
   gtk_flow_box_set_selection_mode(GTK_FLOW_BOX(flowbox), GTK_SELECTION_NONE);
   // gtk_flow_box_set_max_children_per_line(GTK_FLOW_BOX(flowbox), 3);
   gtk_flow_box_set_homogeneous(GTK_FLOW_BOX(flowbox), TRUE);
   gtk_widget_set_margin_start(flowbox, 4);
   gtk_widget_set_margin_end(flowbox, 4);
   gtk_widget_set_margin_top(flowbox, 4);
   gtk_widget_set_margin_bottom(flowbox, 4);
   gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scrolled), flowbox);

   results_window_t *rw = new results_window_t;
   rw->window = win;
   rw->scrolled = scrolled;
   rw->canvas = canvas;
   rw->mol_idx = mol_idx;
   rw->geom_p = geom_p;

   for (const auto &m : matches) {
      GtkWidget *tile = gtk_box_new(GTK_ORIENTATION_VERTICAL, 4);

      GtkWidget *picture = gtk_picture_new();
      gtk_widget_set_size_request(picture, 200, 200);
      gtk_picture_set_can_shrink(GTK_PICTURE(picture), TRUE);
      gtk_box_append(GTK_BOX(tile), picture);

      char label_text[64];
      snprintf(label_text, sizeof(label_text), "%s   %.2f", m.comp_id.c_str(), m.score);
      GtkWidget *label = gtk_label_new(label_text);
      gtk_box_append(GTK_BOX(tile), label);

      // The explicit action button - only this assigns names (not the picture).
      GtkWidget *assign_btn = gtk_button_new_with_label("Assign atom names from this ligand");
      // Size the button to its content instead of stretching to the column width.
      gtk_widget_set_halign(assign_btn, GTK_ALIGN_CENTER);
      gtk_widget_set_hexpand(assign_btn, FALSE);
      assign_ctx_t *ctx = new assign_ctx_t{ rw, m.comp_id };
      g_object_set_data_full(G_OBJECT(assign_btn), "assign_ctx", ctx,
                             +[](gpointer p){ delete static_cast<assign_ctx_t *>(p); });
      g_signal_connect(assign_btn, "clicked", G_CALLBACK(on_assign_clicked), nullptr);
      gtk_box_append(GTK_BOX(tile), assign_btn);

      gtk_flow_box_append(GTK_FLOW_BOX(flowbox), tile);
      GtkWidget *child = gtk_widget_get_parent(tile); // the GtkFlowBoxChild

      tile_task_t *t = new tile_task_t;
      t->comp_id = m.comp_id;
      t->geom_p = geom_p;
      t->child = child;
      t->picture = GTK_PICTURE(picture);
      t->started = false;
      rw->tiles.push_back(t);
   }

   // Own the rw struct on the window; free tiles when the window is destroyed.
   g_object_set_data_full(G_OBJECT(win), "results_window_data", rw, results_window_free);

   // Lazy loading: when the user scrolls, and once after the window is shown.
   GtkAdjustment *vadj = gtk_scrolled_window_get_vadjustment(GTK_SCROLLED_WINDOW(scrolled));
   g_signal_connect(vadj, "value-changed", G_CALLBACK(on_vadjustment_value_changed), rw);

   gtk_window_present(GTK_WINDOW(win));
   g_idle_add(load_visible_idle, rw); // first batch, after allocation
}
