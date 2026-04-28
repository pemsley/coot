#    recent_cryo_em.py
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

# Browse cryo-EM structures released in the last week from PDBe.
# Displays metadata and entry images in a GTK4 scrolled window.

import gi
gi.require_version('Gtk', '4.0')
gi.require_version('GdkPixbuf', '2.0')
from gi.repository import Gtk, Gdk, GLib, GdkPixbuf, Gio, Pango

import gzip
import json
import os
import shutil
import threading
import urllib.request
import urllib.parse
import urllib.error
import socket
from datetime import datetime, timedelta

socket.setdefaulttimeout(15)

# Cache directories
ENTRY_IMAGE_CACHE_DIR = os.path.join(os.path.expanduser("~"), ".cache", "Coot", "structure-images")
LIGAND_IMAGE_CACHE_DIR = os.path.join(os.path.expanduser("~"), ".cache", "Coot", "monomer-images")
DOWNLOAD_DIR = os.path.join(os.path.expanduser("~"), ".cache", "Coot", "coot-download")

def _ensure_cache_dir(d):
   if not os.path.isdir(d):
      os.makedirs(d, exist_ok=True)

def _week_range(weeks_ago=0):
   """Return (start_date, end_date) as datetime for the week `weeks_ago` weeks in the past.
   weeks_ago=0 means the current week (last 7 days ending today)."""
   end = datetime.today() - timedelta(weeks=weeks_ago)
   start = end - timedelta(days=7)
   return start, end


def _format_date(dt):
   """Format a datetime as a Solr date string: YYYY-MM-DDT00:00:00Z"""
   return dt.strftime("%Y-%m-%dT00:00:00Z")


def _week_label(weeks_ago):
   """Human-readable label for the week range."""
   start, end = _week_range(weeks_ago)
   return "{} to {}".format(start.strftime("%d %b %Y"), end.strftime("%d %b %Y"))


def _pdbe_search_url(weeks_ago=0):
   """Build the PDBe Solr search URL for cryo-EM structures released in a given week."""
   base = "https://www.ebi.ac.uk/pdbe/search/pdb/select"
   fields = ("pdb_id,emdb_id,title,resolution,release_date,experimental_method,"
             "pubmed_author_list,entry_type,compound_id,"
             "ligand_of_interest,ligand_of_interest_name,"
             "organism_scientific_name,number_of_polymer_entities,status")
   start, end = _week_range(weeks_ago)
   date_fq = "release_date:[{} TO {}]".format(_format_date(start), _format_date(end))
   params = {
      "q": "*:*",
      "wt": "json",
      "rows": "500",
      "group": "true",
      "group.field": "pdb_id",
      "group.ngroups": "true",
      "json.nl": "map",
      "fl": fields,
      "fq": [
         "experimental_method:\"Electron Microscopy\"",
         date_fq,
      ],
      "sort": "release_date desc",
   }
   # Build URL manually because fq appears multiple times
   parts = [base + "?"]
   for key, val in params.items():
      if key == "fq":
         for fq_val in val:
            parts.append("fq=" + urllib.parse.quote(fq_val, safe="") + "&")
      else:
         parts.append(key + "=" + urllib.parse.quote(str(val), safe=",:-") + "&")
   return "".join(parts).rstrip("&")


def _entry_image_url(pdb_id):
   """URL for the PDBe entry front-image (200x200 PNG)."""
   return ("https://www.ebi.ac.uk/pdbe/static/entry/"
           + pdb_id.lower()
           + "_deposited_chain_front_image-200x200.png")


def _fetch_json(url):
   """Fetch and parse JSON from a URL. Returns dict or None."""
   try:
      req = urllib.request.Request(url)
      req.add_header("Accept", "application/json")
      with urllib.request.urlopen(req, timeout=15) as resp:
         return json.loads(resp.read().decode("utf-8"))
   except Exception as e:
      print("recent_cryo_em: failed to fetch JSON:", e)
      return None


def _download_image(pdb_id):
   """Download the entry image to the cache. Returns the local file path or None."""
   _ensure_cache_dir(ENTRY_IMAGE_CACHE_DIR)
   local_path = os.path.join(ENTRY_IMAGE_CACHE_DIR, pdb_id.lower() + "_front_200.png")
   if os.path.isfile(local_path) and os.path.getsize(local_path) > 0:
      return local_path
   url = _entry_image_url(pdb_id)
   try:
      urllib.request.urlretrieve(url, local_path)
      if os.path.isfile(local_path) and os.path.getsize(local_path) > 100:
         return local_path
   except Exception as e:
      print("recent_cryo_em: image download failed for", pdb_id, ":", e)
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
      print("recent_cryo_em: ligand image download failed for", tlc, ":", e)
   return None


def _fetch_emdb_summary(emdb_id):
   """Fetch reconstruction info from the EMDB API. Returns a short summary
   string, or None if unavailable."""
   url = "https://www.ebi.ac.uk/emdb/api/entry/" + emdb_id.upper()
   data = _fetch_json(url)
   if not data:
      return None
   title = data.get("admin", {}).get("title", "")
   if title == "SUPPRESSED":
      return "suppressed"
   try:
      sd = data["structure_determination_list"]["structure_determination"][0]
      method = sd.get("method", "")
      ip = sd.get("image_processing", [{}])[0]
      fr = ip.get("final_reconstruction", {}) or {}
      res = fr.get("resolution", {}).get("valueOf_")
      n = fr.get("number_images_used")
      sw = [s.get("name") for s in fr.get("software_list", {}).get("software", []) if s.get("name")]
      parts = []
      if method:
         parts.append(method)
      if n:
         parts.append("{} particles".format(n))
      if res:
         parts.append("{} Å".format(res))
      if sw:
         parts.append("(" + ", ".join(sw) + ")")
      return ", ".join(parts) if parts else None
   except (KeyError, IndexError, TypeError):
      return None


def _load_emdb_info_async(emdb_label_pairs):
   """For each (emdb_id, label_widget), fetch the EMDB summary in a thread
   and update the label on the main thread."""
   if not emdb_label_pairs:
      return

   def do_fetch():
      for emdb_id, lbl in emdb_label_pairs:
         summary = _fetch_emdb_summary(emdb_id)
         text = "{}: {}".format(emdb_id.upper(), summary) if summary else "{}: (no info)".format(emdb_id.upper())
         GLib.idle_add(_set_emdb_label, lbl, text)

   thread = threading.Thread(target=do_fetch, daemon=True)
   thread.start()


def _set_emdb_label(label_widget, text):
   label_widget.set_markup("<small>{}</small>".format(GLib.markup_escape_text(text)))
   return False


def _mmcif_url(pdb_id):
   """PDBe updated mmCIF download URL."""
   return "https://www.ebi.ac.uk/pdbe/entry-files/download/{}.cif".format(pdb_id.lower())


def _emdb_map_url(emdb_id):
   """EMDB map (gzipped CCP4) download URL. emdb_id is e.g. 'EMD-12345'."""
   num = emdb_id.upper().replace("EMD-", "")
   return ("https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-{n}"
           "/map/emd_{n}.map.gz").format(n=num)


def _download_with_size_check(url, final_path):
   """Download `url` to final_path via a .tmp file. After transfer, verify the
   on-disk size matches the server's Content-Length (if provided), then rename.
   Returns (ok, message)."""
   tmp_path = final_path + ".tmp"
   try:
      req = urllib.request.Request(url)
      with urllib.request.urlopen(req, timeout=30) as resp:
         expected = resp.headers.get("Content-Length")
         expected = int(expected) if expected is not None else None
         with open(tmp_path, "wb") as f:
            written = 0
            while True:
               chunk = resp.read(65536)
               if not chunk:
                  break
               f.write(chunk)
               written += len(chunk)
      if expected is not None and written != expected:
         os.remove(tmp_path)
         return False, "size mismatch: got {} expected {}".format(written, expected)
      if written < 100:
         os.remove(tmp_path)
         return False, "file too small ({} bytes)".format(written)
      os.rename(tmp_path, final_path)
      return True, "ok ({} bytes)".format(written)
   except Exception as e:
      try:
         if os.path.isfile(tmp_path):
            os.remove(tmp_path)
      except Exception:
         pass
      return False, str(e)


def _download_entry_files_async(pdb_id, emdb_ids, status_setter):
   """Download mmCIF for pdb_id and map for the first emdb_id (if any) into
   DOWNLOAD_DIR. status_setter(text) is called on the main thread to update UI."""
   def report(msg):
      GLib.idle_add(lambda: (status_setter(msg), False)[1])

   def do_download():
      _ensure_cache_dir(DOWNLOAD_DIR)

      cif_path = os.path.join(DOWNLOAD_DIR, pdb_id.lower() + ".cif")
      report("Downloading {}.cif...".format(pdb_id.upper()))
      ok, msg = _download_with_size_check(_mmcif_url(pdb_id), cif_path)
      if not ok:
         report("mmCIF failed: " + msg)
         return

      if not emdb_ids:
         report("Downloaded {}.cif (no EMDB id)".format(pdb_id.upper()))
         return

      last_msg = ""
      for emdb_id in emdb_ids:
         map_path = os.path.join(DOWNLOAD_DIR,
                                 "emd_" + emdb_id.upper().replace("EMD-", "") + ".map.gz")
         report("Downloading {} map...".format(emdb_id.upper()))
         ok, msg = _download_with_size_check(_emdb_map_url(emdb_id), map_path)
         if ok:
            report("Decompressing {}...".format(emdb_id.upper()))
            unzipped = map_path[:-3]  # strip .gz
            tmp_unzipped = unzipped + ".tmp"
            try:
               with gzip.open(map_path, "rb") as gz_in, open(tmp_unzipped, "wb") as out:
                  shutil.copyfileobj(gz_in, out, 1024 * 1024)
               os.rename(tmp_unzipped, unzipped)
               os.remove(map_path)
               report("Downloaded {} + {} (unzipped)".format(
                  pdb_id.upper(), emdb_id.upper()))
            except Exception as e:
               if os.path.isfile(tmp_unzipped):
                  try: os.remove(tmp_unzipped)
                  except Exception: pass
               report("gunzip failed: {}".format(e))
            return
         last_msg = "{}: {}".format(emdb_id.upper(), msg)
      report("Map failed: " + last_msg)

   thread = threading.Thread(target=do_download, daemon=True)
   thread.start()


def _truncate(text, max_len=80):
   if not text:
      return ""
   if len(text) <= max_len:
      return text
   return text[:max_len - 3] + "..."


def _make_entry_row(doc):
   """Create a GTK4 Box widget for one PDB entry."""
   pdb_id = doc.get("pdb_id", "????")
   emdb_ids = doc.get("emdb_id", [])
   title = doc.get("title", "No title")
   resolution = doc.get("resolution")
   release_date = doc.get("release_date", "")
   authors = doc.get("pubmed_author_list", [])
   organisms = doc.get("organism_scientific_name", [])
   compounds = doc.get("compound_id", [])

   # Outer row: horizontal box with image placeholder + text
   row = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=12)
   row.set_margin_start(8)
   row.set_margin_end(8)
   row.set_margin_top(6)
   row.set_margin_bottom(6)

   # Image placeholder (will be filled async)
   image_frame = Gtk.Frame()
   image_frame.set_size_request(130, 130)
   image_widget = Gtk.Image()
   image_widget.set_pixel_size(120)
   image_widget.set_from_icon_name("image-missing")
   image_frame.set_child(image_widget)
   row.append(image_frame)

   # Text column
   text_box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=4)
   text_box.set_hexpand(True)

   # PDB ID + EMDB + resolution line
   id_label = Gtk.Label()
   emdb_str = "  ({})".format(", ".join(e.upper() for e in emdb_ids)) if emdb_ids else ""
   res_str = "  {:.2f} A".format(resolution) if resolution else ""
   date_str = release_date[:10] if release_date else ""
   status = doc.get("status", "")
   if isinstance(status, list):
      status = status[0] if status else ""
   status_labels = {
      "HOLD": "On Hold",
      "HPUB": "Hold for Publication",
      "PROC": "Processing",
      "WAIT": "Waiting",
      "REFI": "Refinement",
      "WDRN": "Withdrawn",
      "OBS":  "Obsolete",
   }
   status_badge = ""
   if status and status != "REL":
      label = status_labels.get(status, status)
      status_badge = ("   <span background='#ffcc66' foreground='black'>"
                      " <b>{}</b> </span>").format(GLib.markup_escape_text(label))
   id_label.set_markup(
      "<b><big>{}</big></b>{}{}   <small>{}</small>{}".format(
         pdb_id.upper(), emdb_str, res_str, date_str, status_badge)
   )
   id_label.set_halign(Gtk.Align.START)
   id_label.set_selectable(True)
   text_box.append(id_label)

   # Title
   title_label = Gtk.Label(label=_truncate(title, 100))
   title_label.set_halign(Gtk.Align.START)
   title_label.set_wrap(True)
   title_label.set_wrap_mode(Pango.WrapMode.WORD_CHAR)
   title_label.set_max_width_chars(80)
   text_box.append(title_label)

   # EMDB reconstruction info (lazy)
   emdb_label_pairs = []
   if emdb_ids:
      for eid in emdb_ids:
         emdb_lbl = Gtk.Label()
         emdb_lbl.set_markup("<small>{}: fetching info...</small>".format(eid.upper()))
         emdb_lbl.set_halign(Gtk.Align.START)
         emdb_lbl.set_wrap(True)
         emdb_lbl.set_max_width_chars(80)
         text_box.append(emdb_lbl)
         emdb_label_pairs.append((eid, emdb_lbl))

   # Authors (short)
   if authors:
      auth_str = ", ".join(authors[:4])
      if len(authors) > 4:
         auth_str += " et al."
      auth_label = Gtk.Label()
      auth_label.set_markup("<small>{}</small>".format(GLib.markup_escape_text(auth_str)))
      auth_label.set_halign(Gtk.Align.START)
      auth_label.set_wrap(True)
      auth_label.set_max_width_chars(80)
      text_box.append(auth_label)

   # Organism
   if organisms:
      org_label = Gtk.Label()
      org_label.set_markup("<small><i>{}</i></small>".format(
         GLib.markup_escape_text(", ".join(organisms[:2]))))
      org_label.set_halign(Gtk.Align.START)
      text_box.append(org_label)

   # Ligands of interest (with images)
   ligands_of_interest = doc.get("ligand_of_interest", [])
   ligand_names = doc.get("ligand_of_interest_name", [])
   # Build a TLC -> name map from the "TLC : name" strings
   name_map = {}
   for entry in ligand_names:
      if " : " in entry:
         tlc_part, name_part = entry.split(" : ", 1)
         name_map[tlc_part.strip()] = name_part.strip()

   # Deduplicate
   seen_lig = set()
   unique_ligands = []
   for tlc in ligands_of_interest:
      if tlc not in seen_lig:
         seen_lig.add(tlc)
         unique_ligands.append(tlc)

   ligand_image_widgets = []
   if unique_ligands:
      lig_header = Gtk.Label()
      lig_header.set_markup("<small><b>Ligands of interest:</b></small>")
      lig_header.set_halign(Gtk.Align.START)
      text_box.append(lig_header)

      lig_row = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
      for tlc in unique_ligands[:6]:
         lig_box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=2)
         lig_img = Gtk.Image()
         lig_img.set_pixel_size(80)
         lig_img.set_from_icon_name("image-missing")
         lig_box.append(lig_img)
         ligand_image_widgets.append((tlc, lig_img))

         name = name_map.get(tlc, "")
         caption = tlc if not name else "{} ({})".format(tlc, _truncate(name, 30))
         cap_label = Gtk.Label()
         cap_label.set_markup("<small>{}</small>".format(
            GLib.markup_escape_text(caption)))
         cap_label.set_max_width_chars(25)
         cap_label.set_wrap(True)
         lig_box.append(cap_label)
         lig_row.append(lig_box)
      text_box.append(lig_row)

   # Download button + status label
   dl_box = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
   dl_button = Gtk.Button(label="Download mmCIF + map")
   dl_status = Gtk.Label()
   dl_status.set_halign(Gtk.Align.START)
   dl_status.set_hexpand(True)
   dl_box.append(dl_button)
   dl_box.append(dl_status)
   text_box.append(dl_box)

   def on_download_clicked(_btn):
      dl_button.set_sensitive(False)
      def set_status(msg):
         dl_status.set_markup("<small>{}</small>".format(GLib.markup_escape_text(msg)))
         # Re-enable the button when we reach a terminal state
         if msg.startswith("Downloaded") or "failed" in msg:
            dl_button.set_sensitive(True)
      _download_entry_files_async(pdb_id, emdb_ids, set_status)

   dl_button.connect("clicked", on_download_clicked)
   if status and status != "REL":
      dl_button.set_sensitive(False)
      dl_status.set_markup("<small><i>files not yet available</i></small>")

   row.append(text_box)
   return row, image_widget, pdb_id, ligand_image_widgets, emdb_label_pairs


def _load_image_async(pdb_id, image_widget):
   """Download the entry image in a thread, then update the widget on the main thread."""
   def do_download():
      path = _download_image(pdb_id)
      if path:
         GLib.idle_add(_set_image, image_widget, path, 120)

   thread = threading.Thread(target=do_download, daemon=True)
   thread.start()


def _load_ligand_images_async(ligand_image_widgets):
   """Download ligand 2D diagrams in a thread, then update widgets on the main thread."""
   def do_download():
      for tlc, image_widget in ligand_image_widgets:
         path = _download_ligand_image(tlc)
         if path:
            GLib.idle_add(_set_image, image_widget, path, 80)

   if ligand_image_widgets:
      thread = threading.Thread(target=do_download, daemon=True)
      thread.start()


def _set_image(image_widget, path, pixel_size):
   """Set the image widget from a file path (called on main thread)."""
   try:
      gfile = Gio.File.new_for_path(path)
      texture = Gdk.Texture.new_from_file(gfile)
      image_widget.set_from_paintable(texture)
      image_widget.set_pixel_size(pixel_size)
   except Exception as e:
      print("recent_cryo_em: failed to load image:", e)
   return False  # Remove from idle


def recent_cryo_em_browser(weeks_ago=0):
   """Open a GTK4 window showing cryo-EM structures released in a given week.
   weeks_ago=0 means the most recent week."""

   state = {"weeks_ago": weeks_ago}

   window = Gtk.Window()
   window.set_title("Recent Cryo-EM Structures")
   window.set_default_size(750, 600)

   outer_box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=4)
   outer_box.set_margin_start(4)
   outer_box.set_margin_end(4)
   outer_box.set_margin_top(4)
   outer_box.set_margin_bottom(4)

   # Navigation bar: Previous Week | date label | Next Week
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
      status_label.set_text("Fetching cryo-EM entries from PDBe...")
      # Disable nav buttons while loading
      prev_button.set_sensitive(False)
      next_button.set_sensitive(False)

      # Clear previous results
      while True:
         child = results_box.get_first_child()
         if child is None:
            break
         results_box.remove(child)

      def fetch_and_populate():
         url = _pdbe_search_url(weeks_ago)
         data = _fetch_json(url)
         GLib.idle_add(_on_data_received, data)

      def _on_data_received(data):
         _populate_results(data, results_box, status_label, weeks_ago)
         # Re-enable nav buttons; disable "Next week" if we're on the current week
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


def _extract_docs_from_grouped(data):
   """Extract one doc per PDB entry from the grouped Solr response."""
   grouped = data.get("grouped", {}).get("pdb_id", {})
   ngroups = grouped.get("ngroups", 0)
   groups = grouped.get("groups", [])
   docs = []
   for g in groups:
      doclist = g.get("doclist", {}).get("docs", [])
      if doclist:
         doc = doclist[0]
         # Deduplicate compound_id list
         compounds = doc.get("compound_id", [])
         if compounds:
            seen = set()
            unique = []
            for c in compounds:
               if c not in seen:
                  seen.add(c)
                  unique.append(c)
            doc["compound_id"] = unique
         docs.append(doc)
   return ngroups, docs


def _populate_results(data, results_box, status_label, weeks_ago):
   """Populate the results box with entry rows (called on main thread)."""

   if data is None:
      status_label.set_text("Failed to retrieve data from PDBe.")
      return False

   ngroups, docs = _extract_docs_from_grouped(data)

   if not docs:
      status_label.set_text("No cryo-EM structures found for this week.")
      return False

   status_label.set_text("{} cryo-EM entries for {}".format(
      ngroups, _week_label(weeks_ago)))

   for doc in docs:
      row, image_widget, pdb_id, ligand_image_widgets, emdb_label_pairs = _make_entry_row(doc)

      # Add a separator between entries
      sep = Gtk.Separator(orientation=Gtk.Orientation.HORIZONTAL)
      results_box.append(sep)
      results_box.append(row)

      # Start async image downloads
      _load_image_async(pdb_id, image_widget)
      _load_ligand_images_async(ligand_image_widgets)
      _load_emdb_info_async(emdb_label_pairs)

   return False  # Remove from idle


# Entry point: can be called from Coot's Python console
# or run as a standalone script with GTK4 main loop.
if __name__ == "__main__":
   import sys

   weeks_ago = 0
   if len(sys.argv) > 1:
      try:
         weeks_ago = int(sys.argv[1])
      except ValueError:
         pass

   app = Gtk.Application(application_id="org.coot.recent_cryo_em")

   def on_activate(app):
      recent_cryo_em_browser(weeks_ago)
      for w in Gtk.Window.get_toplevels():
         if isinstance(w, Gtk.Window):
            w.set_application(app)
            break

   app.connect("activate", on_activate)
   app.run(None)
