---
layout: post
title:  "Coot Cryo-EM Tutorial 2: Fitting and Mutating"
date: Tue 20 Oct 16:53:52 BST 2020
---

<div style="text-align: right"> By Ana Casa&ntilde;al &amp; Paul Emsley</div>
This tutorial is designed for 0.9.1 or later.

Aim: we will fit a domain of the Cleavage and Polyadenylation Factor

1: Setting Up
-------------

### 1.1 Fetch the Files

Using a web browser, download the bundle for EMD-3908

  - [http://www.ebi.ac.uk/pdbe/entry/emdb/EMD-3908](http://www.ebi.ac.uk/pdbe/entry/emdb/EMD-3908)

Move that tar file here and extract it.

Let's download the sequence of a fragment of the structure:

  - [https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/files/CPF-X-domain.seq](https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/files/CPF-X-domain.seq)

Move that sequence file to the current directory (for easy access).

Start _Coot_:

`$ coot --map EMD-3908/map/emd_3908.map`

### 1.2 Map Manipulation

Let's see more of the map

  - **Edit** &rarr; **Map Parameters** &rarr; **Map Radius EM** &rarr; `70`
  - **OK**

Let's use a smoother map

   - **Calculate** &rarr; **Modules** &rarr; **Cryo-EM**
   - **Cryo-EM** &rarr; **Sharpen/Blur...**
   - activate the "Resample" checkbutton
   - **Make Map**
   - **Close**

You should now have an extra map ("emd_3908.map Blur 20.00").
Compare this maps with the original. You should find that the
new (smooth) map is easier to understand.

  - use the **Display Manager** to delete the first (mrc) map.

Change the contour Step for the new Maps:

**Display Manager** &rarr; **Properties** &rarr; **Change by rmsd** `0.33` &rarr; **OK**

As a rule of thumb, a good default/starting contour level is 5.5 rmsd.

### 1.3 Get the Homolog:

  - **File** &rarr; **Open Coordinates...** `coot_tutorial_2/6f9n.pdb`
  - if _Coot_ gives you a "Nomenclature Errors" dialog, just click "Cancel"
  - Use the **Display Manager** to change the representation
  - Instead of "Bonds (Colour by Atom)" use "C-alphas/Backbone"

Move back to the middle of the molecule map

  - **Cryo-EM** &rarr; **Go To Map Molecule Middle**

  - Now move the homolog to here (the centre of the map):
  - **Calculate** &rarr; **Move Molecule Here**
  - Choose the "6f9n" molecule
  - **Move It**

You can now change the colo(u)rs if you wish:
  - **Display Manager** &rarr; **Properties** &rarr; **Colour**
  - **Edit** &rarr; **Bond Colours** [slide the slider]

I like to work with a blue map and the model as orange or green (this is good colour scheme
for working but not good for screenshots used in presentations).

2: Jiggle
---------

  - **Morph** &rarr; **Jiggle Fit This Molecule with Fourier Filter** 

It should roughly fit now. If it doesn't, try jiggling again once or twice more.


3: Extract Our Fragment
-----------------------

Extract the worst-fitting (WD40) domain:

  - You can see that there are 2 chains.
  - We want to extract a wD40 domain from the larger A chain
  - Using **Jones' Rainbow**, find the domain start and end residues numbers (you are trying
    to find a doughnut-shaped molecule that fits this doughnut-shaped density)
  - Let's imagine that you think that the residues at the ends of the domain are 517 and 1011:
  - **Edit** &rarr; **Copy Fragment** &rarr; [Use 6f9n molecule:] `//A/517-1011` &rarr; **OK**

Let's work on this fragment:

  - **Display Manager** &rarr; **Last Only**
  - **Reset View**

Let's delete the sidechains from the atom selection molecule (number 4 usually):

  - **Calculate** &rarr; **Modelling** &rarr; **Delete Side-chains for Active Chain**

For the most recent model (bottom of the list), in the **Display Manager** use

  - **C-alphas/Backbone**
  - **OK** [Dismiss the Display Manager]


4: Setup Refinement
-------------------

We need to adjust the weighting of the map to the model:

  - **R/RC** &rarr; **Refinement Weight** &rarr; `60.00` &rarr; **OK**
    (don't click "Estimate" we need tighter weights at the moment)

    Let's add some local distance restraints:

  - **Calculate** &rarr; **Modules** &rarr; **Restraints**

    Usually 5.0A works well for models with no sidechains:
  - **Restaints** &rarr; **Generate All Molecule Self Restraints 5.0**
  - _{Coot displays new self-restraints as thin grey lines}_
    Briefly look at these restraints, then undisplay them:
  - **Restaints** &rarr; **Undisplay Extra Distance Restraints**

    Let's tighten up the Geman-McClure restraints a bit:

  - **Refine** &rarr; **Set Geman-McClure alpha 0.01**


5: Refinement
-------------

This is where the power of _Coot_ 0.9 becomes apparent.

_{I'd like to
note here that on a modern computer with multiple cores and threads,
this is a pleasant experience. If that is not like your computer then
the following section can be confusing and frustrating}_

  - **Refine** &rarr; **Chain Refine**
  - {wait and watch, you can turn the view if you wish (be careful not
     to click an atom)... this takes about 20 seconds on my computer
     (a PC from 2019)}
  - When the refinement dialog says "Success":
  - **Accept**
  Now let's refine again with less tight Geman-McClure restraints:
  - **Refine** &rarr; **Set Geman-McClure alpha 0.1**
  - **Refine** &rarr; **Chain Refine**
  - _{wait and watch ~10 seconds}_


_{A note on "yanking": by which I mean "smooth displacement" - not jerky shifts}_

When the refinement dialog says "Success," this time we don't yet
press Accept...  Examine the model, being careful not to inadvertently
pull on an atom. Maybe you will see that there is a region that
doesn't fit, if so, yank on the worst fitting CA and pull it to where
you think it should go.

   - Double-clicking on an atom will release the pull restraint

When you are happy, dismiss the Refinement dialog:

  - **Accept**

### 5.1 Redo

It can be non-trivial to decide what needs to move and how to move
it. It is worth undoing your modifications and refining again for
practice.

This time perhaps with the distance restraints shown:

  - Undo (left-facing arrow icon)
  - **Calculate** &rarr; **Scripting** &rarr; **Python**
  - `set_draw_moving_atoms_restraints(1)`
  - **Refine** &rarr; **Chain Refine**
  - _yank as needed_
  - **OK**

... or with different cut-off for the Geman-McClure restraints, or a
different alpha for the Geman-McClure restraints, or a different
weight for the map. Or a different blur for the map. You can delete
the current extra restraints with **Restraints** &rarr; **Delete All
Extra Restraints**.

  - Try proportial editing: with the Real Space Refinement active, use Ctrl Middle-mouse scroll
    to change the radius of the atoms affected by the atom pull displacement.

  - Test, play, refine, yank until satisfied.

Reset Geman-McClure alpha to 1:

  - **Refine** &rarr; **Set Geman-McClure alpha 1**



6: Review and Trim
------------------

Upon review, you will notice that there are parts of the model that
don't fit the map. Try yanking them around with Tandem Refine. Other
parts of the model don't have density - so delete the residue range -
this may help the alignment we are about to do.

Maybe the density fit validation dialog will be useful? You might need
to reset the weight: `1.5` seems like a good number

  - **Validate** &rarr; **Density Fit Analysis** [Choose the "atom selection from pdb6f9n"
    molecule]


7: Mutate
---------

This is _Coot's_ version of "Homology Modelling" - except that the model geometry
optimization occurs in the context of the experiemental data:

  - We have finished with the "Self Restraints" - let's delete them:
  - **Restraints** &rarr; **Delete Extra Restraints**
  - Coot displays a molecule chooser dialog:
  - Select the "atom selection from pdb6f9n.ent" molecule
  - **OK**

  - **Calculate &rarr; Use ClustalW for Alignment, Then Mutate**
    
    The chain for mutation is the "atom selection from pdb6f9n.ent" and the
    target sequence is in the file `CPF-X-domain.seq`
  - _{wait}_
  - **Refine** &rarr; **Chain Refine**
  - _{wait}_
  - **OK**


8: Check & Edit
---------------

Go through the structure residue by residue looking for things to fix (this can take
some time, so you might just do a few, for a flavour, before moving onto Section 8.2)

### 8.1 Do some Fix-ups

  - Out of register errors
  - Missing side-chains
  - Missing loops
  - Loops with too many residues in the model
  - Bad rotamers
  - Bad peptides
  - Clashes
  - Use: **Validate** &rarr; **Alignment vs. PIR...** to help you adjust the residue numbers

There will be places where you need to close (or open) a loop by renumbering residues.

  - Use **Validate &rarr; Validation Outliers**

    or
  - **Validate &rarr; Overlaps, Peptides, CBeta, Rama and Rota Outliers**

Side chains atoms can be added the residue at the centre of the screen using "K"
and deleted using "Shift K".

### 8.2 Unclash

Part of the work on the refinement has been to change the way non-bonded contacts
are minimized. To reduce/remove atom overlaps (or "Clashes" as Molprobity would
call them)

  - First determine the model number of the domain we are fixing - for me it's `2`
  - and now add hydrogen atoms to that model:
  - **Calculate** &rarr; **Scripting** &rarr; **Python**
  - `coot_reduce(2)`    
  - **Refine** &rarr; **Chain Refine**
  - _{wait - this takes a lot longer than before - there are lots of interactions now}_
  - **Accept**

_{For your notes, don't act on this now, to delete the hydrogen atoms, use the Python
  scripting function `delete_hydrogens(4)` (or whatever your model number is)}_

### 8.3 Ramachandran Outliers

It can be convenient to use Ramachandran Restraints to improve the model validation. This
is how I do it:

   - Click the **R/RC** button
   - Tick **Use Torsion Restraints**
   - Tick **Ramachandran Restraints**
   - **OK**
   - **Validate** &rarr; **Ramachandran Plot** [and chose the "atom selection" molecule]
   - **Outliers Only**

Now click on a red spot - they may be many of them (maybe 50 or so):

   - Look at the model. Sometimes the outlier is problematic because there is a model-building
     error - you should not attempt to fix such problems using Ramachandran restraints.
   - But sometimes, things are good, it's just that the main-chain torsion angles need
     tweaking:
   - **Tandem Refine**
   - Now flipping and back-flipping on This Peptide or the Next Peptide often fixes the
     problem.
   - If needed you can change the Ramachandran restraints weight using the menu items
     in the Refine menu. (But note however, that _Coot_ refinement often gets upset when the
     Ramachandran weights are high.)

9: Done
-------

You can check how well you did by comparing against the reference model - the accession code
for that is 6oej `coot_tutorial_2/6eoj.pdb`
