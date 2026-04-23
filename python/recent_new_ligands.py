#    recent_new_ligands.py
#    Copyright 2026 by Medical Research Council
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Browse new CCD (Chemical Component Dictionary) entities released in the last week.
# Displays compound metadata and 2D diagrams in a GTK4 scrolled window.

import gi
gi.require_version('Gtk', '4.0')
gi.require_version('GdkPixbuf', '2.0')
from gi.repository import Gtk, Gdk, GLib, GdkPixbuf, Gio, Pango

import json
import os
import threading
import urllib.request
import urllib.parse
import urllib.error
import socket
from datetime import datetime, timedelta

socket.setdefaulttimeout(15)

# Cache directory for ligand 2D diagrams
LIGAND_IMAGE_CACHE_DIR = os.path.join(os.path.expanduser("~"), ".cache", "Coot", "monomer-images")

def _ensure_cache_dir(d):
   if not os.path.isdir(d):
      os.makedirs(d, exist_ok=True)

def _week_range(weeks_ago=0):
   """Return (start_date, end_date) as datetime for the week `weeks_ago` weeks in the past."""
   end = datetime.today() - timedelta(weeks=weeks_ago)
   start = end - timedelta(days=7)
   return start, end

def _format_date_solr(dt):
   """Format a datetime as a Solr date string: YYYY-MM-DDT00:00:00Z"""
   return dt.strftime("%Y-%m-%dT00:00:00Z")

def _week_label(weeks_ago):
   """Human-readable label for the week range."""
   start, end = _week_range(weeks_ago)
   return "{} to {}".format(start.strftime("%d %b %Y"), end.strftime("%d %b %Y"))

def _pdbe_search_url(weeks_ago=0):
   """Build the PDBe Solr search URL for entries released in a given week,
   returning their ligands of interest."""
   base = "https://www.ebi.ac.uk/pdbe/search/pdb/select"
   fields = "pdb_id,ligand_of_interest,ligand_of_interest_name,release_date"
   start, end = _week_range(weeks_ago)
   date_fq = "release_date:[{} TO {}]".format(_format_date_solr(start), _format_date_solr(end))
   params = {
      "q": "*:*",
      "wt": "json",
      "rows": "2000",
      "group": "true",
      "group.field": "pdb_id",
      "group.ngroups": "true",
      "json.nl": "map",
      "fl": fields,
      "fq": [
         "has_bound_molecule:Y",
         date_fq,
      ],
      "sort": "release_date desc",
   }
   parts = [base + "?"]
   for key, val in params.items():
      if key == "fq":
         for fq_val in val:
            parts.append("fq=" + urllib.parse.quote(fq_val, safe="") + "&")
      else:
         parts.append(key + "=" + urllib.parse.quote(str(val), safe=",:-") + "&")
   return "".join(parts).rstrip("&")

def _fetch_json(url):
   """Fetch and parse JSON from a URL. Returns dict or None."""
   try:
      req = urllib.request.Request(url)
      req.add_header("Accept", "application/json")
      with urllib.request.urlopen(req, timeout=15) as resp:
         return json.loads(resp.read().decode("utf-8"))
   except Exception as e:
      print("recent_new_ligands: failed to fetch JSON:", e)
      return None

def _collect_new_ligands(data):
   """Extract unique 5-letter ligand codes from grouped Solr response.
   Also collects the ligand_of_interest_name map for display."""
   grouped = data.get("grouped", {}).get("pdb_id", {})
   groups = grouped.get("groups", [])
   all_tlcs = set()
   name_map = {}
   for g in groups:
      doclist = g.get("doclist", {}).get("docs", [])
      if doclist:
         doc = doclist[0]
         for tlc in doc.get("ligand_of_interest", []):
            if len(tlc) == 5:
               all_tlcs.add(tlc)
         for entry in doc.get("ligand_of_interest_name", []):
            if " : " in entry:
               tlc_part, name_part = entry.split(" : ", 1)
               tlc_part = tlc_part.strip()
               if len(tlc_part) == 5:
                  name_map[tlc_part] = name_part.strip()
   return sorted(all_tlcs), name_map

def _find_new_compounds(weeks_ago=0):
   """Find 5-letter CCD entities from entries released in a given week.
   Returns a list of (comp_id, details_dict) tuples, or None on failure."""
   url = _pdbe_search_url(weeks_ago)
   data = _fetch_json(url)
   if data is None:
      return None
   tlcs, name_map = _collect_new_ligands(data)
   if not tlcs:
      return []

   compounds = []
   for tlc in tlcs:
      details = _fetch_compound_details(tlc)
      if details is None:
         # Use the name from Solr if PDBe compound API fails
         details = {"name": name_map.get(tlc, "")}
      compounds.append((tlc, details))
   return compounds

def _fetch_compound_details(comp_id):
   """Fetch compound details from PDBe compound summary API.
   Returns a dict with name, formula, weight, creation_date, smiles, etc. or None."""
   # https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/{comp_id}
   url = "https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/{}".format(comp_id)
   try:
      req = urllib.request.Request(url)
      req.add_header("Accept", "application/json")
      with urllib.request.urlopen(req, timeout=15) as resp:
         data = json.loads(resp.read().decode("utf-8"))
      if comp_id in data and data[comp_id]:
         return data[comp_id][0]
   except Exception as e:
      print("recent_new_ligands: failed to fetch details for {}: {}".format(comp_id, e))
   return None

def _ligand_image_url(tlc):
   """URL for the RCSB CCD labelled 2D diagram (SVG)."""
   return "https://cdn.rcsb.org/images/ccd/labeled/{}/{}.svg".format(tlc[0], tlc)

def _download_ligand_image(tlc):
   """Download a ligand 2D diagram to the cache. Returns the local path or None."""
   _ensure_cache_dir(LIGAND_IMAGE_CACHE_DIR)
   local_path = os.path.join(LIGAND_IMAGE_CACHE_DIR, tlc + ".svg")
   if os.path.isfile(local_path) and os.path.getsize(local_path) > 0:
      return local_path
   url = _ligand_image_url(tlc)
   try:
      urllib.request.urlretrieve(url, local_path)
      if os.path.isfile(local_path) and os.path.getsize(local_path) > 100:
         return local_path
   except Exception as e:
      print("recent_new_ligands: image download failed for", tlc, ":", e)
   return None

def _truncate(text, max_len=80):
   if not text:
      return ""
   if len(text) <= max_len:
      return text
   return text[:max_len - 3] + "..."

def _set_image(image_widget, path, pixel_size):
   """Set the image widget from a file path (called on main thread)."""
   try:
      gfile = Gio.File.new_for_path(path)
      texture = Gdk.Texture.new_from_file(gfile)
      image_widget.set_from_paintable(texture)
      image_widget.set_pixel_size(pixel_size)
   except Exception as e:
      print("recent_new_ligands: failed to load image:", e)
   return False

def _make_ligand_row(comp_id, details):
   """Create a GTK4 Box widget for one CCD entry."""

   name = details.get("name", "") if details else ""
   formula = details.get("formula", "") if details else ""
   weight = details.get("weight", None) if details else None
   creation_date = details.get("creation_date", "") if details else ""
   smiles_list = details.get("smiles", []) if details else []
   first_observed = details.get("first_observed_in", []) if details else []

   smiles_str = ""
   if smiles_list:
      for s in smiles_list:
         if isinstance(s, dict):
            smiles_str = s.get("name", "")
            break
         elif isinstance(s, str):
            smiles_str = s
            break

   row = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=12)
   row.set_margin_start(8)
   row.set_margin_end(8)
   row.set_margin_top(6)
   row.set_margin_bottom(6)

   # 2D diagram placeholder
   image_frame = Gtk.Frame()
   image_frame.set_size_request(160, 160)
   image_widget = Gtk.Image()
   image_widget.set_pixel_size(150)
   image_widget.set_from_icon_name("image-missing")
   image_frame.set_child(image_widget)
   row.append(image_frame)

   # Text column
   text_box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=4)
   text_box.set_hexpand(True)

   # TLC + weight + date
   id_label = Gtk.Label()
   weight_str = "  {:.1f} Da".format(weight) if weight else ""
   date_str = creation_date[:10] if creation_date else ""
   id_label.set_markup(
      "<b><big>{}</big></b>{}   <small>{}</small>".format(
         GLib.markup_escape_text(comp_id),
         GLib.markup_escape_text(weight_str),
         GLib.markup_escape_text(date_str))
   )
   id_label.set_halign(Gtk.Align.START)
   id_label.set_selectable(True)
   text_box.append(id_label)

   # Name
   if name:
      name_label = Gtk.Label(label=name)
      name_label.set_halign(Gtk.Align.START)
      name_label.set_wrap(True)
      name_label.set_wrap_mode(Pango.WrapMode.WORD_CHAR)
      name_label.set_max_width_chars(80)
      text_box.append(name_label)

   # Formula
   if formula:
      formula_label = Gtk.Label()
      formula_label.set_markup("<small>Formula: {}</small>".format(
         GLib.markup_escape_text(formula)))
      formula_label.set_halign(Gtk.Align.START)
      text_box.append(formula_label)

   # PDB entry
   if first_observed:
      pdb_ids_str = ", ".join(pid.upper() for pid in first_observed[:5])
      if len(first_observed) > 5:
         pdb_ids_str += " ..."
      pdb_label = Gtk.Label()
      pdb_label.set_markup("<small>PDB: {}</small>".format(
         GLib.markup_escape_text(pdb_ids_str)))
      pdb_label.set_halign(Gtk.Align.START)
      pdb_label.set_selectable(True)
      text_box.append(pdb_label)

   # SMILES
   if smiles_str:
      smiles_label = Gtk.Label()
      smiles_label.set_markup("<small>SMILES: {}</small>".format(
         GLib.markup_escape_text(_truncate(smiles_str, 100))))
      smiles_label.set_halign(Gtk.Align.START)
      smiles_label.set_selectable(True)
      smiles_label.set_wrap(True)
      smiles_label.set_max_width_chars(80)
      text_box.append(smiles_label)

   row.append(text_box)
   return row, image_widget

def _load_image_async(comp_id, image_widget):
   """Download the 2D diagram in a thread, then update the widget on the main thread."""
   def do_download():
      path = _download_ligand_image(comp_id)
      if path:
         GLib.idle_add(_set_image, image_widget, path, 150)
   thread = threading.Thread(target=do_download, daemon=True)
   thread.start()

def recent_new_ligands_browser(weeks_ago=0):
   """Open a GTK4 window showing new CCD entities released in a given week."""

   state = {"weeks_ago": weeks_ago}

   window = Gtk.Window()
   window.set_title("Recent New Ligands (CCD)")
   window.set_default_size(750, 600)

   outer_box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=4)
   outer_box.set_margin_start(4)
   outer_box.set_margin_end(4)
   outer_box.set_margin_top(4)
   outer_box.set_margin_bottom(4)

   # Navigation bar
   nav_box = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
   nav_box.set_margin_bottom(4)

   prev_button = Gtk.Button(label="< Previous week")
   next_button = Gtk.Button(label="Next week >")
   date_label = Gtk.Label()
   date_label.set_hexpand(True)
   nav_box.append(prev_button)
   nav_box.append(date_label)
   nav_box.append(next_button)
   outer_box.append(nav_box)

   # Status label
   status_label = Gtk.Label(label="Fetching...")
   outer_box.append(status_label)

   scrolled = Gtk.ScrolledWindow()
   scrolled.set_vexpand(True)
   scrolled.set_hexpand(True)
   scrolled.set_policy(Gtk.PolicyType.AUTOMATIC, Gtk.PolicyType.AUTOMATIC)

   results_box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=0)
   scrolled.set_child(results_box)
   outer_box.append(scrolled)

   # Close button
   bottom_box = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
   bottom_box.set_halign(Gtk.Align.END)
   bottom_box.set_margin_top(4)
   close_button = Gtk.Button(label="  Close  ")
   close_button.connect("clicked", lambda _btn: window.destroy())
   bottom_box.append(close_button)
   outer_box.append(bottom_box)

   window.set_child(outer_box)

   def _load_week(weeks_ago):
      state["weeks_ago"] = weeks_ago
      date_label.set_markup("<b>{}</b>".format(_week_label(weeks_ago)))
      status_label.set_text("Fetching new CCD entries...")
      prev_button.set_sensitive(False)
      next_button.set_sensitive(False)

      # Clear previous results
      while True:
         child = results_box.get_first_child()
         if child is None:
            break
         results_box.remove(child)

      def fetch_and_populate():
         compounds = _find_new_compounds(weeks_ago)
         GLib.idle_add(_on_data_received, compounds)

      def _on_data_received(compounds):
         _populate_results(compounds, results_box, status_label, weeks_ago)
         prev_button.set_sensitive(True)
         next_button.set_sensitive(state["weeks_ago"] > 0)
         return False

      thread = threading.Thread(target=fetch_and_populate, daemon=True)
      thread.start()

   def on_prev_clicked(_btn):
      _load_week(state["weeks_ago"] + 1)

   def on_next_clicked(_btn):
      if state["weeks_ago"] > 0:
         _load_week(state["weeks_ago"] - 1)

   prev_button.connect("clicked", on_prev_clicked)
   next_button.connect("clicked", on_next_clicked)

   window.present()
   _load_week(weeks_ago)


def _populate_results(compounds, results_box, status_label, weeks_ago):
   """Populate the results box with ligand rows (called on main thread)."""

   if compounds is None:
      status_label.set_text("Failed to retrieve data from RCSB.")
      return False

   if not compounds:
      status_label.set_text("No new CCD entities found for this week.")
      return False

   status_label.set_text("{} new CCD entries for {}".format(
      len(compounds), _week_label(weeks_ago)))

   for comp_id, details in compounds:
      row, image_widget = _make_ligand_row(comp_id, details)
      sep = Gtk.Separator(orientation=Gtk.Orientation.HORIZONTAL)
      results_box.append(sep)
      results_box.append(row)
      _load_image_async(comp_id, image_widget)

   return False


if __name__ == "__main__":
   import sys

   weeks_ago = 0
   if len(sys.argv) > 1:
      try:
         weeks_ago = int(sys.argv[1])
      except ValueError:
         pass

   app = Gtk.Application(application_id="org.coot.recent_new_ligands")

   def on_activate(app):
      recent_new_ligands_browser(weeks_ago)
      for w in Gtk.Window.get_toplevels():
         if isinstance(w, Gtk.Window):
            w.set_application(app)
            break

   app.connect("activate", on_activate)
   app.run(None)
