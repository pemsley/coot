---
layout: post
title:  "Coot Cryo-EM Tutorial 1: Basics"
date: Tue 20 Oct 16:40:39 BST 2020
categories: x
---
<div style="text-align: right"> By Ana Casa&ntilde;al &amp; Paul Emsley</div>
This tutorial is designed for 0.9.1

1: Start
--------

Download this file:
[https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/files/coot-cryo-em-tutorial-basics-files.tar.gz](https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/files/coot-cryo-em-tutorial-basics-files.tar.gz)

Untar it and start Coot from that directory.

  - Open Coot
  - Keep the terminal window where Coot is running always visible. You will see that as many of the algorithms run they will be writing useful information. 

### 1.1 Load the Tutorial Data

  - **Calculate** &rarr; **Load Tutorial Model and Data**
  - Open the **Display Manager**
  - Delete the difference map

### 1.2 Representation

  - **Edit** &rarr; **Map Parameters...**
  - A Map Radius of 70&Aring;  is a lot of map, useful for a cryo-EM reconstruction, but for now we are focussing on details:
  - **Map Radius** &rarr; `16` _{sensible for crystallographic map}_
  - **Edit** &rarr; **Map Colour...**
  - Change the map colour _{see that dark pastel colours are the best}_
  - **Display Manager** &rarr; **Properties**  &rarr; **Solid/Transparent**
  - Tweak the **Opacity %age**
  - Useful? Maybe. (About 30% seems best)
  - Back to **Standard Lines**

Note: if you have a speedy computer, then you can make a copy of the map (Calculate &rarr; Map Tools &rarr; Copy Map Molecule...) and use Standard Lines for one and Solid/Transparent for the copied map - that looks quite good and can be helpful. It can be laggy though with a big radius, and that is the opposite of useful.

Now let's try with a white background:
  - **Edit** &rarr; **Background Colour** &rarr; **White**
  - _{this is good for making figures, but you will need a very pastel map colour (not so dark)}_
  - To "highlight" a residue:
  - **Draw** &rarr; **Additional Representation...**
  - **Ball & Stick** &rarr; **Add Representation**
  - Use the Display manager's check buttons to undisplay the ligand when done.


2: Navigation
---------------

Moving around the molecule:
  - Middle-mouse to centre on an atom
  - Double-click to toggle the atom label
  - "Space" is forwards
  - "Shift Space" is backwards
  - Middle-mouse click and drag to "pan" the view _{just like PyMOL}_

To change the "residue to residue" speed:

 - **Edit** &rarr; **Preferences...**
 - **Smooth Recentering** &rarr; ``10`` steps _{smaller numbers are quicker}_

**Add Here Dragged Map Off**

You might like the residue reorienting mode:
   - **Calculate** &rarr; **Scripting**  &rarr;  **Python**
   - `set_reorienting_next_residue_mode(1)`
   - Now the main-chain orientation of the next residue will rotate to match the current one - that might take a bit of getting used to, so you can turn it off:
   - `set_reorienting_next_residue_mode(0)`

3: Short-cuts and Interface
---------------------------

Go To Residue:

  - "Ctrl G" _{see the little dialog?}_
  - Type the residue number (and the chain id, if needed)
  - Press Enter _{Coot recentres on the specified residue}_
  - The recentering also works with 3 single-letter-code letters, e.g.``HEY``

You can toggle spin and rock the view with:
	- "I", "Ctrl R"

### 3.1 Install a Few Buttons

  - Right-mouse click on the right-hand side of the horizontal tool-button bar
  - Add "**Sphere Refine +**", "**Tandem Refine**" and "**Backrub Rotamers**"
  - Toggle on the "**Backrub Rotamers**" button _{the button will turn pink}_

### 3.2 Install Key-bindings

  - **Edit** &rarr; **Settings** &rarr; **Install Template Key-bindings**
  - **Edit** &rarr; **Settings** &rarr; **Key-Bindings** _{see what the keybindings are}_

### 3.3 Install a Few Extensions

  - **File** &rarr; **Curlew**
  - From there, select and install "**Black Box Morph and Fit**", "**Chain Refine**" and "**Morph**"

### 3.4 Try it out

Use the cursor to navigate using the map
 
  - "G" _{brings the blob under the cursor to the centre of the screen}_

To undo last navigation/screen recentre:

  - "U"

Also:
- "+" and "-" Change the contour level
- "M" and "N" for zoom in and out
- "D" and "F" for depth of field/clipping planes, increase and decrease

### 3.5 Quick Save-As

The new version of coot is a lot more advanced and a bit more crashy than the standard version. It's a good idea to save every now and then - making it easier to recover.

  - "Ctrl S" _{for quick-save-as}_
  
  Quick Save as will save all models that have not been saved and will save the session too.


4: Model-Building
-------------------

### 4.1 Find and Fix a Bad Rotamer

- **Validate** &rarr; **Density Fit Analysis...**
- Click on the bar for residue A89

Fix the Rotamer, First the Traditional usage
   - click on **Auto-Fit Rotamer** in the vertical tool-bar

But let's use the Key-bindings to speeds up, so undo the previous fit:
  - "Ctrl Z"
  - "J" _{auto-fit rotamer}_

And to tidy up after that: 
  - "X"  _{does real-space refinement}_

### 4.2 Let's Find Something Else to Fix

  **Validate** &rarr; **Ramachandran Plot**
   - click on "**Outliers Only**"
    - click on Residue A41
 
Note: Peptide Flags:
    - Green: _cis_-peptide in front of a PRO
    - Yellow: twisted _trans_-peptide
    - Red: _cis_-peptide in front of a non-PRO

When it comes to fixing up model-building/validation problems,
often times "Neighbour Refine" is the first tool to try. If we try that, we see that it doesn't work in this case. We have to look and think. We see that the one of the peptides is facing the wrong way and we need to flip it:

So, let's fix this with Real Space Refine and flipping:
  - "H" _{does a neighbour Refine}_
  - "Q"  _{does a pep-flip}_
  - "H"  _{tidy up with a RSR}_

5: Now To the Demo Box!
----------------------------

The point of this exercise is for you to see/learn how Coot behaves with lower resolution data.

For these refinements, it might be useful to use the Coot atom interaction dots:

**Calculate** &rarr; **Scripting** &rarr; **Python**
   - `set_do_coot_probe_dots_during_refine(1)`

OK, let's find the script and run it:
   - **Calculate** &rarr; **Run Script** &rarr; `demo-box-of-buttons-madrid.scm`
   - **Crankshaft Peptide**  _{brings the problem to the centre of the screen}_
   - Take a look at the problem...
	   - _{Yikes!}_
   - **Morph** &rarr; **Crankshaft**
   - Does it fit better now?
   - **Sphere Refine**
     - _{things improve, balls are green}_
     - **OK**
   - **Wonky N-terminus**
   - **Rotamer Markup On**
   - Centre on Residue A4
   - **Tandem Refine**
   - _{Play with it}_ - make the balls go green and remove the clashes

Now try that with a blurred map:

 - Undo
 - **Calculate** &rarr; **Map Sharpen and Blur...**
 - Blur to about 60 or 70&Aring;&sup2;  _{to resemble cryo-EM map}_
 - **Tandem Refine** again


<!--stackedit_data:
eyJoaXN0b3J5IjpbLTE5OTM1ODM3MzAsNjUyMTIzNzQ5LDEyNj
A4NTYyNzEsMTUzNjg4Mzg2MCwxMDUyOTg2MDM2LDEyNjg3NjQ4
NjgsLTIwMTU2MDM5MDZdfQ==
-->
