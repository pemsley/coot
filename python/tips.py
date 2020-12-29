# Copyright 2006 by The University of York
# Copyright 2006 by Bernhard Lohkamp 

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA

# a list of tips for Coot
def tip_list():
 return [
 """To centre on a particular atom, click it with middle mouse (if that 
  doesn't seem to work it may be because the molecule is not active).""",
 "+ and - on the keyboard change the contour level",
 """To move just one atom (in Regularize, RS Refine, or Rotate/Translate 
  mode) use Ctrl Left-mouse to pick an atom (you have to be accurate).""",
 """There is a mailing list for Coot development and discussion at
  http://www.jiscmail.ac.uk/lists/coot.html""",
 """Use Ctrl Left-mouse to drag (for example) a blob of density to the
  the pointer""",
 """Slow recentering?  Try Draw -> Smooth Recentering and reduce the 
  number of steps"""
 """Want distance of the atoms to the pointer?  
  Use Measures -> Pointer Distances...""",
 "To label an atom: Double click it. Or Shift Left-mouse on it",
 """To make atoms be insensitive to clicking, deactivate it by unclicking
  that molecule's \"Active\" button in the Display Control window""",
 """Can't label or centre on some symmetry atoms?  It's a known bug.  
  (For now, drag that atom closer to the centre of the screen).""",
 "Use function key 'F8' to make a rendered snapshot.",
 "Use the Ctrl key to rotate the view when changing Chi angles.",
 """Too many cis peptides when using dragged refinement?  Use the 
  'Planar Peptide Restraints' suggested in the Coot FAQ.""",
 "Use the Ctrl to rotate the view when using Delete.",
 """When in skeleton mode, new skeleton bones can be displayed around the
  current point using the 'S' key.""",
 """When in baton mode, the baton can be rotated independently from the 
  Guide Points by using the 'B' key (it's a toggle).""",
 """When Editing Chi Angles, switch between the angles quickly using the
 '1', '2', '3', '4' keys.""",
 "Use \"refmac-extra-params\" to pass refmac your personal parameters.",
 """Use \"poly_ala(imol)\" to turn molecule number imol into poly-ALA.  
  Use \"poly_ala(imol,\"SER\")\" to turn it into poly-SER.""",
 """Use \"fit_protein(imol)\" to rotamer search and real-space refine all
  residue of molecule number imol.""",
 "Shift Ctrl right-mouse rotates round screen Z.",
 """Ctrl  + right-mouse + horizontal (left to right) mouse movement moves 
  the view in screen Z.""",
 "Ctrl right-mouse up-down drag changes the slab.",
 """Use \"set_idle_function_rotate_angle(0.05)\" to change the 
  spin speed.""",
 """\"ligand_expert()\" enables the GUI editting of some ligand-fitting
  parameters.""",
 "Use keyboard + and - to zoom in Ramachandran and Kleywegt Plots.",
 "The 'U' key undoes last nagivation (e.g. re-centering on new pdb file).",
 "The 'D' and 'F' keys change the clipping/slabbing.",
 """\"view_matrix()\" prints the current view matrix (useful for 
  molscript, perhaps).""",
 """Baton-building low resolution maps is better down with maps that 
  have increased sampling rate (2.0 or 2.5).""",
 "Coot can read SHELXL .ins/.res files. (It can write them too.)",
 """\"Clear All Atom Labels\" can be found under Measures -> 
  Distances & Angles. Obviously.""",
 """Esc and Return are keyboard accelerators for Reject/Accept for 
  Refinement and Regularization.""",
 "To disable coot tips: add \"no_coot_tips()\" to your $COOT_HOME/.coot.py file.",
 "To validate chiral centres use Validate -> Incorrect Chiral Volumes.." 
 ]

# Function to turn off coot tips at start
def no_coot_tips():
    global do_coot_tips_flag
    do_coot_tips_flag = False
    coot.set_tip_of_the_day_flag(0)
