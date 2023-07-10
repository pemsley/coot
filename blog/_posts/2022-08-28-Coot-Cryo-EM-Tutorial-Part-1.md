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
    * `export __GL_FSAA_MODE=5`
    * (this will make your graphics more smooth on some PCs, nice to have, not important)
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

  Here are the links if you don't have the files:

    - https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/files/32143-partial-model.pdb
    - https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-32143/map/emd_32143.map.gz

  Now that we have read in the map, let's expand the map radius to cover the molecule:
  - **Edit** &rarr; **Map Parameters** &rarr; **Map Radius** &rarr; `77` &rarr; **Apply**

  It is frequently the case that to make the map more interpretable we need to unsharpen (_i.e._ blur) and resample it:
  - **Cryo-EM** &rarr; **Sharpen/Blur...** &rarr; [check **Resample Factor**] &rarr; **OK**

If you like, compare and contrast these two maps, you should see that the second map is more easy to interpret.

Now use the Display Manager to Delete the original `EMD-32143.map` map.

If you have a slow computer, dragging the view while displaying EM maps can be annoying. If that is the case for you,
turn off dynamic recontouring:

   - **Edit** &rarr; **Preferences** &rarr; **Maps** &rarr; **Dragged Map** &rarr; No &rarr; **OK**

Make a copy of the partial model to which we will add various components:

   - **Edit** &rarr; **Copy Molecule** &rarr; **OK** (no need to change the selection)

Let's change the name of the copy:

   - **Calculate** &rarr; **Scripting** &rarr; **Python**

 in the "Command:" entry:
   - `set_molecule_name(3, 'working.pdb')`
   - Press "Enter" to execute (you want to be setting the name of the newly created copy of the model (not a map!))

_Coot_ will echo your command in the green box below.

To see that that worked, look at the molecule name for the model in the Display Manager.
The model number is 2 if you have followed the instructions exactly.

3: Rigid-Body Fitting
----------

To bring the parts of the map that do not yet have an atomic interpretation to our attention, we need
to mask the map by the current model:

   - **Calculate** &rarr; **Map Tools** &rarr; **Mask Map by Atom Selection** &rarr; **Atom Selection** &rarr; `//` &rarr; **OK**
   - Using the Display Manager, undisplay the starting map (for me, that's map number 1).

Sometimes, Coot will set the contour level of the new (residual) map to zero! So let's scroll the contour level of the new map
to about 8.0 r.m.s.d.

In the residual map you should now you should see a WD40 domain (beta propeller), an immunoglobulin domain next to it, and a helix
somewhat separate.
The domain with the best density is the WD40 domain. Let's start by modelling that with an AlphaFold2 model:

   - Using middle-mouse click and drag, bring the centre of the WD40 domain to the centre of the screen

[Note: it is only this part of the tutorial (_i.e._ "where is the middle of the domain?") that needs any skill,
everything else is just following the script (and thus is scriptable).]

   - **File** &rarr; **Fetch AlphaFold2 Model using UniProt ID...** &rarr; `P54311` &rarr; **Get it**

If you do not have AlphaFold2 utilities in your copy of _Coot_ then use the UniProt web site
to download the AlphaFold2 model for ``P54311``

  - https://www.uniprot.org/

Check that you have the "Morph" menu item in the main menu bar. If you do not, then you will need to get the "Morph" tools from Curlew

  - **File** &rarr; **Curlew**
  - Look for the **Morph** tool and click its **Install** button

  - Use the combobox in the Display Manager to change the representation style of the AlphaFold2 ("AF-P54311-F1-model_v3.pdb") model
    to "CA + Ligands".
    Doing so makes the overall fit of the model to the map more clear.

 - **Morph** &rarr; **Jiggle Fit with Fourier Filtering**
   - _{Coot adds the "Morph" menu to the main menu bar}_
   - Close the Curlew dialog.

[wait a couple of seconds - you can think of this step as local molecular placement. We could have used MOLREP or PHASER
to achieve the same result with no human intervention, but with considerably more palaver]

   - _{The AlphaFold2 model should jump into the density}_

You can check that this is the correct solution because it finds some density for the N-terminal helical extension (the
appropriate density level for this helix is at a rather lowver level than the rest of the molecule, so some scrolling will
be in order).

4: Real Space Refinement Morph Fitting
--------------------------------------

RSR morphing in _Coot_ requires that we add some additional restraints:

   - **Calculate** &rarr; **Modules** &rarr; **Refine**
   - _{Coot adds the "Refine" menu to the main menu bar}_
   - **Calculate** &rarr; **Modules** &rarr; **Restraints**
   - _{Coot adds the "Restraints" menu to the main menu bar}_
   - **Restraints** &rarr; **Generate All-Molecule Self Restraints 5.0**
   - _{Coot displayes a grey mesh of additional restraints}_
   - **Restraints** &rarr; **Undisplay Extra Restraints**
   - **R/RC** &rarr; **Refinement Weight** &rarr; **Estimate**
         o "R/RC.." is labelled "Refine/Regularization Control..." if you have labels on your vertical toolbar
   - **More Control** Geman-McClure weight to `0.01` &rarr; **OK**
   - **Refine** &rarr; **Chain Refine**
   - _{watch and wait}_
   - **Accept**

Now let's reduce the weight on the GM restraints:

   - **More Control** Geman-McClure weight to `0.1` &rarr; **OK**
   - **Refine** &rarr; **Chain Refine**

This should allow the model to make local/smaller scale adjustments to better fit the map

   - _{watch and wait}_
   - **Accept** (in the refinement dialog)

5: Trim
-------

See that not all of the model is well-fitting into the reconstruction. We should trim off the bits we don't need.

There are two types of trimming that we can do:

either:

   - Look at the model and use the Delete function to remove residues for which there is little to no density

or (using a recent version of _Coot_)
   - use the Python scripting function `alphafold_pLDDT_colours(3)` to colour by pLDDT and remove the parts of the model that you don't want (the pLDDT replaces the B-factors in AlphaFold2 models). Unlike B-factors, the lower the pLDDT, the lower our confidence in the correctness of that part of the model. A typical value to chop is the 70% level.

6: More Refinement
------------------

 ### 6.1 Unclash

Part of the work on the Coot's refinement has been to change the way non-bonded contacts
are minimized. To reduce/remove atom overlaps (or "Clashes" as Molprobity would
call them)

  - First determine the model number of the domain we are refining - for me it's `4` -
  - and now add hydrogen atoms to that model:
  - **Calculate** &rarr; **Modelling** &rarr; **Add Hydrogen Atoms**

The Extra Distance Restraints are now out of date, we should delete them:

  - **Restraints** &rarr; **Delete Extra Restraints**
  - Choose "AF-54311-F1-model_v3.pdb" in the combobox
  - **OK**

Now refine:

  - **Refine** &rarr; **Chain Refine**
  - _{wait - this takes a lot longer than before - there are lots of interactions now, the atoms should not move
      much if we have done well in the previous fitting and refinement steps}_
  - **Accept**

_{should we wish to delete the hydrogen atoms (we don't right now), we would use the Python
  scripting function `delete_hydrogens(4)` (4 being the model number)}_

### 6.2 Ramachandran Plot Outliers

It can be convenient to use Ramachandran Restraints to improve the Ramachandran Plot of the model (not everyone approves of this):

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

Merge this molecule into the `working.pdb` molecule:

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

