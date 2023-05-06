---
layout: post
title:  "Coot Tutorial 2023 Version"
date: Fri 5 May 13:10:32 BST 2023
---



# Coot Tutorial

Let's do the Coot Tutorial using modern methods.

 - **Edit** &rarr; **Settings** &rarr; **Install Template Keybindings**

Load the tutorial:

 - **Calculate** &rarr; **Load Tutorial Model and Data**

The map is quite nice looking - the resolution of the data is 1.8 A. Let's turn on "Updating Maps"

  - **Calculate** &rarr; **Updating Maps**
  - Check the **Auto-Update** button
  - **OK**

Set the default termperature factor:
 - **Calculate** &rarr; **Scripting...** &rarr; **Python...**
 - `median_temperature_factor(0)`
 - You should get an answer close to 16.
 - `set_default_temperature_factor_for_new_atoms(16)`
 - **Close** the Python window

Let's view the structure anomolies:

 - **Validate** &rarr; **Overlaps, Peptides, CBeta Ramas & Rotas**

_[A dialog is presented with Interesting/Outliers listed by chain-id and residue number]_

- Click on the first button (for A 2)
- Take a look at what is interacting with this VAL
    - It's actually the PHE at A 89 that's wrong.
    - Recentre on that residue (middle-mouse click on the CB)

Fix the rotamer:
 - `J` (on the keyboard)

If things have gone to plan, the sidechain will jump into the density and the maps will update and the difference maps will disappear in this region.

 - `Shift R` to tidy up the geometry
 - **Accept**

Next problem is at A 32. The density here is not very good and there may be a number of differrent conformers/rotamers. Let's, for now just choose one of them that has higher probability than the current model

 - `J`
 - `Shift R`

Notice that the number of Outliers/Problems has decreased.

Now look at the rest of the labels in that dialog. Many of them are around the same residues, 40, 41 and 42. Let's look at the Ramachandran Plot for this chain

 - **Validate** &rarr; **Ramachandran Plot** then pull across to select the protein molecule

_[A Ramachandran Plot appears]_

In the new Ramachandran daialog activate the **Selection** check-button and make the selection text read `//A`. Now click **Apply**

_[The number of Rama spots decreases]_

If you mouse over the red spots they will tell you the residue number.

 - Let's click on the spot for A 41

_[Coot recentres onto the specified residue]_

 - **OK** on the Ramachandran dialog to undisplay it for now.
 - See that there is red density, green density and a red non-pro _cis_ peptide and a twisted trans (yellow ribbon)
 - Use **Sphere Refine** to tidy this up

_[That didn't do the job, there are still red Rama balls]_

 - **Flip This Peptide**

_[The Rama balls go green]_

You may notice that the side-chain orientation is incorrect

 - Use **Backrub Rotamer Fit** to fix it
 - **Accept**

Back to the validation dialog:

  - Click on the button for "A 72"
  - Here we have a missing main-chain atoms and missing side-chain atoms
  - Delete the residue using `Ctrl D`
  - Add in a replacement ALA using `Y`
  - `Shift R` to refine
  - **Accept**
  - **Mutate & Auto-fit**,
  - click on an atom in the new ALA - change it to a "CYS"
  - **Add Alt Conf..** and split the whole residue
  - Rotamer number 1 is the correct one
  - **OK**

One problem left!

- Click on the rotamer outlier button for A91
- `J`

Now the **Interesting/Outliers/Problems** dialog is empty

Let's find parts of the map that are as yet not modelled, _i.e._ find unmodelled blobs:

 - **Validate** &rarr; **Unmodelled blobs**
 - Defaults are fine
 - **Find Blobs**

Blobs are sorted by size.

 - Click on the button **Blob 1**

What are we looking at?

Blob 4 (if present in the dialog) and Blob 3 are the same thing. Let's do Blob 3 first because the density is better. Notice the shape of the blobs and the environement. The precipitant for this structure was ammonium sulfate. So let's add a SO4 here:

 - **Calculate** &rarr; **Modelling** &rarr; **Add Other Solvent Molecules**
 - **SO4: SULPHATE ION**

_[Coot adds a sulfate ion using "Jiggle Fit"]_

 - Now you can do the same thing for **Blob 4** if you wish
 - **Close** the **Solvent Ligands** dialog

Blob 2 is the ligand. Let's build a ligand into that:

 - File &rarr; Get Monomer
 - `3GP` (this is 3'-guanosine monophosphate)
 - **OK**
 - Check the chemistry (yes the phosphate is on the 3' end)
 - Undisplay the new 3GP using the Display Manager
 - **Ligand** &rarr; **Find Ligands..**
 - In the following dialog, check the button for the `3GP_from_dict` ligand (no need to make it flexible)
 - In **Where to Search?** choose "**Right Here**"
 - **Find Ligands**

The ligand fitting function doesn't fit the ligand as well as it used to
now that we are using the new monomer library. We will need to do some
interactive changes to the ligand, using eigen-flipping (`E` key) and real space refinement (`Shift R`).

When the ligand is corrected placed, let's merge it into the main molecule

 - **Edit** &rarr;  **Merge Molecules...**
 - Choose the fitting ligand and add it to molecule 0 (the protein)

**Blob 1** is a piece of missing protein.

 - **Draw** &rarr; **Sequence View** and pull acrooss to the protein model

You can see on the right hand side of the sequce view that the A chain is mssing the residues "QTC"

 - Right-mouse over the sequence view to **Close** it.
 - Middle-mouse click the "C" atom of residue A 93.
 - **Calculate** &rarr; **Fit Loop** &rarr; **Fit Loop by Rama Search...**
 - Choose Residue numbers 94 and 96 in the A chain
 - Add the sequence "QTC"
 - **Fit Loop**

_[Coot builds the loop]_

 - **Calculate** &rarr; **Other Modelling Tools** &rarr; **Add OXT to Residue...**
 - **Add it**

Is there still work to be done?

Maybe. You can use Ctrl Shift Arrow keys to reorient the CYS A 96 if needed (Use Ctrl Arrow keys to
rotate the residue (about an axis that's into the screen)).

---

This tutorial is based on the structure 2sar.


