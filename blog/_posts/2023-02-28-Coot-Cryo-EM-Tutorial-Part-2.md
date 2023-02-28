---
layout: post
title:  "Coot Cryo-EM Tutorial 2022: Part 2"
date: Tue 28 Feb 2023 10:43:23 BST
---

_Coot_ Cryo-EM Tutorial: Part 2: Molecular Extraction
=====================================================

Part 2 of this tutorial is similar to Part 1, but the setup is a little different.

1: Find the Nanobody
--------------------

Now we want to model the nanobody domain.

We have done a BLAST search for the seqence of our nanobody and know that the structure was published already in a different complex - it's the "E" chain in the entry ``6GDG``. So let's fetch that structure, extract the "E" chain and then fit it to our reconstruction.

   - First centre the view on the centre of the nanobody/immunoglobulin domain. You should
   be able to do this by eyeballing it and using middle-mouse click and drag - but if not, then:
   - **Calculate** &rarr; **Scripting** &rarr; **Python**
   - `set_rotation_centre(121, 71, 102)`
   - Press "Enter" to execute.

2: Fetch the Reference Structure
----------------------

Now let's download the reference structure:

   - **File** &rarr; **Fetch PDB using Accession Code** &rarr; `6gdg` &rarr; **Get it**
   - _{Coot downloads 6gdg and moves the centre of the screen it its centre}_
   - **Draw** &rarr; **Undo Last Navigation**
   - _{Coot changes the centre of the screen to were we had been before loading 6dgd}_
   - **Edit** &rarr; **Copy Fragment**
   - Change the molecule selector combobox to the most recently loaded molecule: _i.e._ `coot-download/6gdg.cif`
   - **Atom Selection for fragment** &rarr; `//E`
   - Check the **Move new molecule here** checkbutton
   - **OK**
   - _{Coot places the "E" chain at the centre of the screen}_

We no longer need to see `6dgd.cif` so use the Display Manager to undisplay it.

3: Jiggle-Fit
-------------

Let's make the "E" chain a bit more clear by changing the molecular representation:

   - Ctrl Shift mouse-forward-scroll
   - _{Coot changes to CA-bonds mode}_
   - **Cryo-EM** &rarr; **Jiggle-Fit with Fourier Filter**
   - _{Coot changes the "E" chain to the correct orientation to fit the map}_

4: RSR Morph Refinement
-----------------------

Now we proceed in a similar way to Part 1:

   - **Restraints** &rarr; **Generate All-Molecule Self Restraints 5.0**
   - **Undisplay Extra Restraints**
   - **R/RC** &rarr; **Refinement Weight** &rarr; **Estimate**
   - **More Control** Geman-McClure weight to `0.01` &rarr; **OK**
   - **Refine** &rarr; **Chain Refine**
   - _{watch and wait}_
   - **Accept**

Now let's reduce the weight on the GM restraints:

   - **More Control** Geman-McClure weight to `0.1` &rarr; **OK**

This should allow the model to make local/smaller scale adjustments to better fit the map

   - **Refine** &rarr; **Chain Refine**
   - _{watch and wait}_
   - **Accept**

5: Merge
--------

Merge this into the `working.pdb` molecule:

   - **Edit** &rarr; **Merge Molecules**
   - Check the button for `Atom selection from 6dgd.cif`
   - Change the molecule selector combobox to `working.pdb`
   - **Merge**

