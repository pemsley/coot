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
#
# Uses the PDBe Solr `latest_chemistry` document type to find genuinely new
# CCD codes (rather than all ligands appearing in newly-released entries).

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
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta

socket.setdefaulttimeout(15)

# Cache directory for ligand 2D diagrams
LIGAND_IMAGE_CACHE_DIR = os.path.join(os.path.expanduser("~"), ".cache", "Coot", "monomer-images")

# PDB releases on Wednesdays. weeks_ago=0 means the most recent Wednesday release.
PDB_RELEASE_WEEKDAY = 2  # Monday=0, Tuesday=1, Wednesday=2

def _ensure_cache_dir(d):
   if not os.path.isdir(d):
      os.makedirs(d, exist_ok=True)

def _most_recent_release_date():
   """Return the date of the most recent Wednesday PDB release (today if today is Wed)."""
   today = datetime.today()
   days_back = (today.weekday() - PDB_RELEASE_WEEKDAY) % 7
   return today - timedelta(days=days_back)

def _release_date_for_week(weeks_ago):
   """Return the release-date datetime for `weeks_ago` Wednesdays in the past."""
   return _most_recent_release_date() - timedelta(weeks=weeks_ago)

def _release_window(weeks_ago):
   """Return (start, end) datetimes that bracket a single Wednesday release.
   Solr release_date values are stamped at 01:00 UTC on the release day, so a
   24-hour window centred on the release day captures exactly one release."""
   release = _release_date_for_week(weeks_ago)
   start = release.replace(hour=0, minute=0, second=0, microsecond=0)
   end = start + timedelta(days=1)
   return start, end

def _format_date_solr(dt):
   return dt.strftime("%Y-%m-%dT%H:%M:%SZ")

def _release_label(weeks_ago):
   release = _release_date_for_week(weeks_ago)
   return release.strftime("%d %b %Y")

def _pdbe_search_url(weeks_ago=0):
   """Build the PDBe Solr search URL for new CCD entities released in a given week.
   Uses the latest_chemistry document type (added to PDBe Solr in 2026)."""
   base = "https://www.ebi.ac.uk/pdbe/search/pdb/select"
   start, end = _release_window(weeks_ago)
   params = {
      "q": "document_type:latest_chemistry AND latest_chemistry_entry_type:new",
      "wt": "json",
      "rows": "500",
      "json.nl": "map",
      "fl": "pdb_id,new_ligand,release_date",
      "fq": "release_date:[{} TO {}]".format(_format_date_solr(start), _format_date_solr(end)),
      "group": "true",
      "group.field": "pdb_id",
      "group.ngroups": "true",
   }
   parts = [base + "?"]
   for key, val in params.items():
      parts.append(key + "=" + urllib.parse.quote(str(val), safe=",:-[] ") + "&")
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

def _parse_new_ligand_field(value):
   """Parse a `new_ligand` field value into [(tlc, name), ...].
   The field is multi-valued (list) at the Solr level; each value is
   itself a comma-separated list of `CODE : name` pairs when an entry
   introduces multiple new CCDs simultaneously."""
   if value is None:
      return []
   if isinstance(value, str):
      raw_entries = [value]
   else:
      raw_entries = list(value)

   results = []
   for entry in raw_entries:
      # An entry may itself contain multiple "CODE : name" pairs separated by ", A1"
      # but splitting that reliably is messy because names contain commas.
      # In practice each Solr value is one CODE : name pair, so handle that primarily.
      if " : " in entry:
         tlc_part, name_part = entry.split(" : ", 1)
         tlc_part = tlc_part.strip()
         if 1 <= len(tlc_part) <= 5:
            results.append((tlc_part, name_part.strip()))
   return results

def _collect_new_ligands(data):
   """Extract unique new CCD codes from grouped Solr response.
   Returns a dict: {tlc: {"name": str, "pdb_ids": set([...])}}"""
   grouped = data.get("grouped", {}).get("pdb_id", {})
   groups = grouped.get("groups", [])
   ligands = {}
   for g in groups:
      pdb_id = g.get("groupValue", "")
      doclist = g.get("doclist", {}).get("docs", [])
      if not doclist:
         continue
      doc = doclist[0]
      for tlc, name in _parse_new_ligand_field(doc.get("new_ligand")):
         if tlc not in ligands:
            ligands[tlc] = {"name": name, "pdb_ids": set()}
         if pdb_id:
            ligands[tlc]["pdb_ids"].add(pdb_id.upper())
   return ligands

def _fetch_compound_details(comp_id):
   """Fetch compound details from PDBe compound summary API.
   Used for enrichment (formula, weight, SMILES) the name already comes
   from the latest_chemistry document."""
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

def _find_new_compounds(weeks_ago=0):
   """Find genuinely-new 5-letter CCD entities released in a given week.
   Returns a list of (comp_id, details_dict) tuples, or None on failure.
   Compound enrichment (formula, weight, SMILES) is fetched concurrently."""
   url = _pdbe_search_url(weeks_ago)
   data = _fetch_json(url)
   if data is None:
      return None
   ligands = _collect_new_ligands(data)
   if not ligands:
      return []

   # Concurrent enrichment, say 30+ ligands sequentially is sluggish, in parallel it's snappy.
   def enrich(item):
      tlc, base_info = item
      details = _fetch_compound_details(tlc) or {}
      # Authoritative name comes from latest_chemistry; only fall back to the
      # compound API name if Solr didn't give us one.
      if base_info["name"]:
         details["name"] = base_info["name"]
      details["_pdb_ids_release"] = sorted(base_info["pdb_ids"])
      return tlc, details

   compounds = []
   with ThreadPoolExecutor(max_workers=8) as pool:
      for result in pool.map(enrich, sorted(ligands.items())):
         compounds.append(result)
   return compounds

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
   smiles_list = details.get("smiles", []) if details else []
   # Prefer the entries from THIS week's release - that's what the user is browsing.
   # Fall back to first_observed_in if the release-list is somehow empty.
   pdb_ids = details.get("_pdb_ids_release", []) if details else []
   if not pdb_ids:
      pdb_ids = details.get("first_observed_in", []) if details else []

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

   # TLC + weight
   id_label = Gtk.Label()
   weight_str = "  {:.1f} Da".format(weight) if weight else ""
   id_label.set_markup(
      "<b><big>{}</big></b>{}".format(
         GLib.markup_escape_text(comp_id),
         GLib.markup_escape_text(weight_str))
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

   # PDB entries that introduced this CCD this week
   if pdb_ids:
      pdb_ids_str = ", ".join(pid.upper() for pid in pdb_ids[:5])
      if len(pdb_ids) > 5:
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

   prev_button = Gtk.Button(label="< Previous release")
   next_button = Gtk.Button(label="Next release >")
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
      date_label.set_markup("<b>Release of {}</b>".format(_release_label(weeks_ago)))
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
      status_label.set_text("Failed to retrieve data from PDBe.")
      return False

   if not compounds:
      status_label.set_text("No new CCD entities found for the {} release.".format(
         _release_label(weeks_ago)))
      return False

   status_label.set_text("{} new CCD entries in the {} release".format(
      len(compounds), _release_label(weeks_ago)))

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
