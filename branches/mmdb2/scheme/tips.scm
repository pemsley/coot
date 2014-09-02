
;; a list of tips for Coot
(define tip-list (list
 "To centre on a particular atom, click it with middle mouse (if that 
  doesn't seem to work it may be because the molecule is not active)."
 "+ and - on the keyboard change the contour level"
 "To move just one atom (in Regularize, RS Refine, or Rotate/Translate 
  mode) use Ctrl Left-mouse to pick an atom (you have to be accurate)."
 "There is a mailing list for Coot development and discussion at
  http://www.jiscmail.ac.uk/lists/coot.html"
 "Use Ctrl Left-mouse to drag (for example) a blob of density to the
  the pointer"
 "Slow recentering?  Try Draw -> Smooth Recentering and reduce the 
  number of steps"
 "Want distance of the atoms to the pointer?  
  Use Measures -> Pointer Distances..."
 "To label an atom: Double click it. Or Shift Left-mouse on it"
 "To make atoms be insensitive to clicking, deactivate it by unclicking
  that molecule's \"Active\" button in the Display Control window"
 "Can't label or centre on some symmetry atoms?  It's a known bug.  
  (For now, drag that atom closer to the centre of the screen)."
 "Use function key 'F8' to make a rendered snapshot."
 "Use the Ctrl key to rotate the view when changing Chi angles."
 "Too many cis peptides when using dragged refinement?  Use the 
  'Planar Peptide Restraints' suggested in the Coot FAQ."
 "Use the Ctrl to rotate the view when using Delete."
 "When in skeleton mode, new skeleton bones can be displayed around the 
  current point using the 'S' key."
 "When in baton mode, the baton can be rotated independently from the 
  Guide Points by using the 'B' key (it's a toggle)."
 "When Editing Chi Angles, switch between the angles quickly using the
 '1', '2', '3', '4' keys."
 "Use \"refmac-extra-params\" to pass refmac your personal parameters."
 "Use \"(poly-ala imol)\" to turn molecule number imol into poly-ALA.  
  Use \"(poly-ala imol 'SER)\" to turn it into poly-SER."
 "Use \"(fit-protein imol)\" to rotamer search and real-space refine all
  residue of molecule number imol."
 "Shift Ctrl right-mouse rotates round screen Z."
 "Ctrl + right-mouse + horizontal (left to right) mouse movement moves 
  the view in screen Z."
 "Ctrl right-mouse up-down drag changes the slab."
 "Use \"(set-idle-function-rotate-angle 0.05)\" to change the 
  spin speed."
 "\"(ligand-expert)\" enables the GUI editting of some ligand-fitting
  parameters."
 "Use keyboard + and - to zoom in Ramachandran and Kleywegt Plots"
 "The 'U' key undoes last nagivation (e.g. re-centering on new pdb file)."
 "The 'D' and 'F' keys change the clipping/slabbing."
 "\"(view-matrix)\" prints the current view matrix (useful for 
  molscript, perhaps)."
 "Baton-building low resolution maps is better done with maps that 
  have increased sampling rate (2.0 or 2.5)."
 "Coot can read SHELXL .res files. (It can write them too.)"
 "\"Clear All Atom Labels\" can be found under Measures -> 
  Distances & Angles. Obviously."
 "Esc and Return are keyboard accelerators for Reject/Accept for 
  Refinement and Regularization."
 "Restraints for alpha helical and beta-strand structure can in the 
  Refinement/Regularization Control Panel"
 "To disable coot tips: add \"(no-coot-tips)\" to your ~/.coot file."
 "To validate chiral centres use Validate -> Incorrect Chiral Volumes.."
))

;; Function to turn off coot tips at start
(define no-coot-tips
  (lambda ()
    (set! *do-coot-tips-flag* #f)
    (set-tip-of-the-day-flag 0)))
