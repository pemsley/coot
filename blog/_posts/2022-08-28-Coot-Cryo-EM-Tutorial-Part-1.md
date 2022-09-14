---
layout: post
title:  "Coot Cryo-EM Tutorial 2022: Part 1"
date: Sat 27 Aug 2022 10:43:23 BST
---

_Coot_ Cryo-EM Tutorial: Part 1: Using An AlphaFold2 Model
=================================================

The idea is to complete a partial model of the parathyroid+PTH1R recpeptor complex using models provided by AlphaFold2.

This tutorial is designed for 0.9.8.4 or any later version in the 0.9.x series. _Coot_ 1 works differently.


1: Start
--------

  - Open a Terminal
  - `mkdir coot_tutorial_2`
  - `cd coot_tutorial_2`
  - if you have an NVIDIA graphics card:
  - `export __GL_FSAA_MODE=5`
  - Start _Coot_:
  - `coot`
  - Keep the terminal window where _Coot_ is running always visible. You will
    see that as many of the algorithms run they will be writing information.

Activate the Cryo-EM Module:

  - **Calculate** &rarr; **Modules** &rarr; **Cryo-EM**

2: Load the Map and Model for this Tutorial
------------------------------------------

  - **File** &rarr; **Open Coordinates** &rarr; ``32143-partial-model.pdb`` &rarr; **Open**
  - **File** &rarr; **Open Map** &rarr; ``EMD-32143.map`` &rarr; **Open**

  Expand the map radius to cover the molecule:
  - **Edit** &rarr; **Map Parameters** &rarr; **Map Radius** &rarr; `77` &rarr; **Apply**

  It is frequently the case that to make the map more interpretable we need to unsharpen (_i.e._ blur) and resample it:
  - **Cryo-EM** &rarr; **Sharpen/Blur...** &rarr; [check **Resample Factor**] &rarr; **OK**

Now use the Display Manager to Delete the original `EMD-32143.map` map.

If you have a slow computer, dragging the view while displaying EM maps can be annoying. If that is the case for you,
turn off dynamic recontouring:

   - **Edit** &rarr; **Preferences** &rarr; **Maps** &rarr; **Dragged Map** &rarr; **OK**

Make a copy of the partial model to which we will add various components:

   - **Edit** &rarr; **Copy Molecule** &rarr; **OK**

Let's change the name of the copy:

   - **Calculate** &rarr; **Scripting** &rarr; **Python**
   - `set_molecule_name(2, 'working.pdb')`
   - Press "Enter" to execute.

3: Rigid-Body Fitting
----------

To bring the parts of the map that do not yet have an atomic interpretation to our attention, we need
to mask the map by the current model:

   - **Calculate** &rarr; **Map Tools** &rarr; **Mask Map by Atom Selection** &rarr; **Atom Selection** &rarr; `//` &rarr; **OK**
   - Using the Display Manager, undisplay the starting map (for me, that's map number 1).

In the residual map you should now you should see a WD40 domain, an immunoglobulin domain next to it, and a helix
somewwhat separate.
The domain with the best density is the WD40 domain. Let's start by modelling that with an AlphaFold2 model:

   - Using middle-mouse click and drag, bring the centre of the WD40 domain to the centre of the screen

[Note: it is only this part of the tutorial (_i.e._ "where is the middle of the domain?") that needs any skill,
everything else is just following the script, and thus is scriptable.]

   - **File** &rarr; **Fetch AlphaFold2 Model using UniProt ID...** &rarr; `P54311` &rarr; **OK**
   - Press "U" on the keyboard to bring the view back to the map for the WD40 domain

If you do not have AlphaFold2 utilities in your copy of _Coot_ then use the UniProt web site
to download the AlphaFold2 model for ``P54311``

https://www.uniprot.org/

   - **Calculate** &rarr; **Move Molecule Here** &rarr; `3: AF-P54311-F1-model_v3.pdb` &rarr; **OK**
   - **Cryo-EM** &rarr; **Jiggle Fit with Fourier Filtering**

[wait a couple of seconds - you can think of this step as local molecular placement. We could have used MOLREP to achieve the same result with no human intervention, but considerably more palaver]

   - _{The AlphaFold2 model should jump into the density}_

You can check that this is the correct solution because it finds some density for the helical extension.

4: Real Space Refinement Morph Fitting
--------------------------------------

RSR morphing in _Coot_ requires that we add some additional restraints:

   - **Calculate** &rarr; **Modules** &rarr; **Refine**
   - **Calculate** &rarr; **Modules** &rarr; **Restraints**
   - **Restraints** &rarr; **Generate All-Molecule Self Restraints 5.0**
   - **Undisplay Extra Restraints**
   - **R/RC** &rarr; **Refinement Weight** &rarr; **Estimate**
   - **More Control** Geman-McClure weight to `0.01` &rarr; **OK**
   - **Refine** &rarr; **Chain Refine**
   - _{watch and wait}_
   - **Accept**

Now let's reduce the weight on the GM restraints:

   - **More Control** Geman-McClure weight to `0.1` &rarr; **OK**
   - **Refine** &rarr; **Chain Refine**

This should allow the model to make local/smaller scale adjustments to better fit the map

   - _{watch and wait}_
   - **Accept**

5: Trim
-------

See that not all of the model is well-fitting into the reconstruction. We should trim off the bits we don't need.

There are two types of trimming that we can do:

either:

   - Look at the model and use the Delete function to remove misfitting residue ranges

or (using a recent version of _Coot_)
   - use the Python scripting function `alphafold_pLDDT_colours(3)` to colour by pLDDT and remove the parts of the model that you don't want (the pLDDT replaces the B-factors in AlphaFold2 models). Unlike B-factors, the lower the pLDDT, the lower our confidence in the correctness of that part of the model. A typical value to chop is the 70% level.

6: More Refinement
------------------

 ### 6.1 Unclash

Part of the work on the refinement has been to change the way non-bonded contacts
are minimized. To reduce/remove atom overlaps (or "Clashes" as Molprobity would
call them)

  - First determine the model number of the domain we are refining - for me it's `3` -
  - and now add hydrogen atoms to that model:
  - **Calculate** &rarr; **Scripting** &rarr; **Python** &rarr; `coot_reduce(3)`

The Extra Distance Restraints ar now out of date, we should delete them:

  - **Restraints** &rarr; **Delete Extra Distance Restraints**

Now refine:

  - **Refine** &rarr; **Chain Refine**
  - _{wait - this takes a lot longer than before - there are lots of interactions now}_
  - **Accept**

_{should we wish to delete the hydrogen atoms (we don't right now), we would use the Python
  scripting function `delete_hydrogens(4)` (4 being the model number)}_

### 6.2 Ramachandran Plot Outliers

It can be convenient to use Ramachandran Restraints to improve the Ramachandran Plot of the model:

   - Click the **R/RC** button
   - Tick **Use Torsion Restraints**
   - Tick **Ramachandran Restraints**
   - **OK**
   - **Validate** &rarr; **Ramachandran Plot** [and chose the `atom selection` molecule]
   - **Outliers Only**

   Sometimes, the main-chain  &phi; and  &psi; torsion angles need tweaking:

   - Click on a red spot in the Ramachandran Plot
   - _{Coot recentres the view to that residue}_
   - **Tandem Refine**
   - (if you don't have a **Tandem Refine** button, enable it by clicking right-mouse on the right-hand side
     of the vertical divider in the horizontal toolbar and selecting **Manage Buttons** in the pop-up.)
   - Flipping and back-flipping by clicking on **This Peptide** or the **Next Peptide** often fixes the problem.
     (The greener the Rama-balls, the higher the probability of the  &phi;, &psi; pair for that residue)
   - If needed, you can change the Ramachandran restraints weight using the menu items
     in the **Refine** menu or the **R/RC** dialog. (But note however, that the _Coot_ minimizer often gets
     upset when the Ramachandran weights are high.)

7: Merge
--------

Merge this into the `working.pdb` molecule:

   - **Edit** &rarr; **Merge Molecules**
   - Check the check-button for `AF-P54311-F1-model_v3.pdb`
   - Change the molecule selector combobox to `working.pdb`
   - **Merge**


8: How did we do? [Optional]
-----------------

If you wish to continue to Part 2 (as would typically be the case), then now is not the time for this.

If however, you think "That's enough" then you can compare your model with the published model, ``7vvl``.

   - Download the model using **File** &rarr; **Fetch Model using Accession Code**&rarr; `7vvl` &rarr; **OK**
   - Compare and contrast.

