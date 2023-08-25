---
layout: post
title:  "Coot Cryo-EM Tutorial 2022: Part 3"
date: Tue 28 Feb 2023 11:00:00 GMT
---


_Coot_ Cryo-EM Tutorial: Part 3: De Novo Model-Building
=======================================================

Now it's time to build the parathyroid hormone (or at least, the part of it that we can see in the reconstruction)

1: Build a Poly-Ala Helix
----------------

   - Centre on the middle of the helix that is the part of the map for the parathyroid hormone. You should be able to eyeball this, but if not:
   - **Calculate** &rarr; **Scripting** &rarr; **Python**
   - `set_rotation_centre(121, 71, 102)`
   - Press "Enter" to execute.
   - **Calculate** &rarr; **Other Modelling Tools** &rarr; **Place Helix Here**
   - _{Coot creates a helix}_
   - **Refine** &rarr; **Chain Refine**
   - **Accept**


2: Add Sidechains
----------------

Now let's put on some sidechains:

   - **Calculate** &rarr; **Assign Sequence...** &rarr; **1: Associate Sequence File...**
   - Change the molecule selector combobox to the most recently created molecule: _i.e._ `Helix-5`
    (your molecule number may be different)
   - Use the **File...** button to find and open the file `parathyroid-hormone.pir`
   - **OK**
   - **Cryo-EM** &rarr; **Assign Sequence Based on Associated Sequence**
   - _{Coot adds sidechains}_
   - **Refine** &rarr; **Chain Refine**
   - For me, the side chain for `LYS 58` is in the wrong direction, so pull it into the correct tube (the density is not great)
   - Bring `TRP 54` to the centre of the screen and
   - **Backrub Rotamer Fit**
   - _{Coot flips the sidechain into its density}_
   - **Accept**


3: Merge
--------

Merge this into the `working.pdb` molecule:

   - **Edit** &rarr; **Merge Molecules**
   - Check the button for `Helix-5` (or whatever your molecule number was)
   - Change the molecule selector combobox to `working.pdb`
   - **Merge**


4: How did we do?
-----------------

Now we can compare your model with the published model, ``7vvl``.

   - Download the model using **File** &rarr; **Fetch Model using Accession Code**&rarr; `7vvl` &rarr; **OK**
   - Compare and contrast.

