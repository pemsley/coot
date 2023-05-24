# Makefile.am
# 
# Copyright 2012 The University of Oxford
# Copyright 2013 by Medical Research Council
# Author: Paul Emsley
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

# desktop icons go in icons.  Application icons go here.
# 
SUBDIRS = icons

cootpixmapdir = $(pkgdatadir)/pixmaps

# no longer use xpms.

dist_cootpixmap_DATA = \
coot-0.9.8.8.png     \
add-alt-conf.png      dna.svg			   probe-clash.svg                \
add-alt-conf.svg      edit-backbone.svg		   radacina_handcuffs-small.svg   \
add-peptide-1.png     edit-chi.png		   ramachandran.png               \
add-peptide-1.svg     edit-chi.svg		   refine-1.svg                   \
add-water.svg	      emsley.jpg		   refine.png                     \
anchor.svg	      env-distances.png		   regularize-1.svg               \
aniso-atom.png	      env-distances.svg		   regularize.png                 \
aniso-atom.svg	      eugene.jpg		   renumber-residues.svg          \
appointment-new.png   flip-peptide.png		   reset-view-32.png              \
atom-at-pointer.png   flip-peptide.svg		   reset-view-32.svg              \
atom-at-pointer.svg   generic-display-objects.svg  reset-view.png                 \
auto-fit-rotamer.png  geom.svg			   reset-view.svg                 \
auto-fit-rotamer.svg  ghosts.png		   residue-info.png               \
azerbaijan2.svg       ghosts.svg		   residue-info.svg               \
azerbaijan.svg	      gln-asn-b-factors.svg	   rigid-body.png                 \
calculate-ghosts.svg  gln-asn-b-factorssvg.svg	   rigid-body.svg                 \
cell+symm.png	      glyco-structures.svg	   rotamers.png                   \
cell+symm.svg	      gnome-run.png		   rotamers.svg                   \
change-chain-id.svg   go-to-atom.svg		   rtz.png  rrz.png               \
chiral.png	      go-to-ligand.svg		   rtz.svg  rrz.svg               \
chiral.svg	      gtk-connect-16.png	   scroll-wheel.svg               \
connect-to-ccp4.svg   gtk-font.png		   secondary-structure.svg        \
connect-to-oca.svg    gtk-open-16.png		   sequence-view.png              \
connect-to-usf.svg    gtk-open-24.png		   sequence-view.svg              \
connect-to-world.svg  gtk-open.png		   side-chain-180.png             \
	    kevin.jpg			   side-chain-180.svg             \
	    kevin-orig.jpg		   smiles.png                     \
	    kleywegt.svg		   smiles.svg                     \
	    map.svg			   spin-view.png                  \
coot-favicon.png      merge-molecules.svg	   spin-view.svg                  \
coot-icon.png	      mutate-auto-fit.svg	   stereo-view.svg                \
coot-save.png	      mutate.svg		   stop.svg                       \
coot-save-symm.png    panel-screenshot.png	   temperature.png                \
crosshairs.png	      peptide-omega.png		   temperature.svg                \
crosshairs.svg	      peptide-omega.svg		   torsion-general.svg            \
delete.png	      phenixed.svg		   undo-1.png                     \
delete.svg	      pink-pointer.png		   undo-1.svg                     \
density-fit.svg       pink-pointer.svg		   unknown-blob.png               \
diff-map.png	      plane-distances.png	   unknown-blob.svg               \
diff-map.svg	      plane-distances.svg	   water-diff-map.svg             \
display-manager.png   pointer-distances.png	   water-drop.svg                 \
display-manager.svg   pointer-distances.svg \
koppi-green-tick-64.png rygle-red-cross-64.png forest51690-dwcheckno-incomplete-64.png \
  go-to-ligand.png go-to-atom.png azerbaijan.png add-water.png mutate.png \
  mutate-auto-fit.png edit-backbone.png torsion-general.png anchor.png regularize-1.png \
  refine-1.png stock-zoom-follow-window-12.svg stock-zoom-follow-window-12.png \
  rama_all.png rama_gly.png rama_pro.png rama_non_gly_pro.png \
  curlew-long.png curlew-128-filled.svg

if CLIPPER_TOP8000
dist_cootpixmap_DATA += \
  rama2_all.png 	rama2_gly.png 	rama2_pro.png   \
  rama2_non_gly_pro.png rama2_ileval.png rama2_pre_pro.png \
  rama2_non_gly_pro_pre_pro_ileval.png
endif

# removed the svg version of curlew for now because of non anti-aliasing


# This is done by the Makefile in icons now.
# dist_icons_DATA = coot.png
# iconsdir = ${datarootdir}/icons/hicolor/48x48/apps
