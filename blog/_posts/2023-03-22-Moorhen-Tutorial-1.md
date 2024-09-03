---
layout: post
title: "Moorhen Tutorial 1: Fix up the Cyclin-Dependent Kinase"
date: Tue 3 Sep 14:41:45 GMT 2024
---

Welcome to Moorhen ("Coot on the Web").

## Getting Started

[layout](https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/web/tutorial/images/moorhen-gui-items.png)


  - Open a web browser window and point it at [moorhen.org](https://moorhen.org).

Let's load some data:

  - **File** &rarr; **Load Tutorial Data ...** &rarr; **Tutorial 1** &rarr; **OK**

_Moorhen displays a protein model, a blue 2Fo-Fc-style map and an Fo-Fc-style map in green (positive) and red (negative)._

 - Use Left-Mouse click and drag to rotate the view.
   (just click and drag, when using at trackpad)

 - Use scroll-wheel scroll to zoom in and out
   (use 2-finger drag on a trackpad)

 - Use middle-mouse click to centre on an atom
    (use Option-click on a trackpad on a Mac)

 - To pan the view, use middle-mouse click and drag
   (use Shift Option click and drag on a trackpad on a Mac)

You can change the speed that moving the mouse spins the view:

  - **Preferences** &rarr; "Mouse Sensitivity" &rarr; 0.4 (# for example)
  - Click off the Preferences dialog to make the dialog disappear.

  - Likewise one can change the sensitivity of the zoom using the "Mouse wheel zoom sensitivity" slider.
    I find a value of around 7 to 8 to be more useful than the default.

  - Similarly, you can change the thickness of the map lines if you wish using "Map contour settings..."

  - Use "[" and "]" on the keyboard to adjust the radius of the density.
  - Ctrl middle-mouse scroll to change the contour level (one step at time seems good to me).

## Let's Go!

Our job is to fix and amend the protein model in a way that is consistent with the data. Let's first look
at the Ramachandran plot:

## Models

   - Click on the **Models** button
You will see buttons for various representation styles for the model.
   - Click on the **Sequences** accordian button
You will see a sequence viewer - let's use that to move around the structure. Click on a few letters (that represent the residues in the protein). Notice that the map density mesh is redrawn around the new centre.

See the grey rectangle over the sequence numbers?

  - Move your cursor inside the box and use left-mouse click and drag to move around the sequence.

You can resize the rectangle to display more residue letters (by clicking and dragging on
the _edge_ of the box), but if you make it too wide it will not display any.

  - Click on the **Ligands** tab

Notice that there are "No ligands." We will add one later.

For higher-end computers we can use a smoother representation of the bonds and atoms.

 - Click on the Settings icon in the vertical toolbar
 - Change the "Bond Smoothness" to **Nice**.

 - Use the **X** icon to close the **Models** dialog (as you normally would, to close a window).

## Maps
   - Click on the **Maps** button

The 2Fo-Fc-style map has a blue icon and the difference map has an icon with red and green. You can use the sliders there to adjust the contour level. Click or slide the slider for the 2Fo-Fc map so that the level is about 0.42 (1.40 rmsd) or thereabouts. Set the difference map to a level of around 0.53 (3.6 rmsd) (if it isn't already).

  - Close the dialog and then you can also use _Ctrl scroll-wheel scroll_ to change the contour level.

## Validation Tools

  - Open the **Validation** window

  - Choose **Ramachandran Plot**

_Moorhen shows the Ramachandran Plot for the "A" chain of this protein_

   - You can resize the dialog for a better view of the plot (pull slowly on the bottom right triangle).

You will see that there are several interesting red spots.

  - Let's click on the spot at the middle top

_Moorhen will put residue A180 at the centre of the screen_

Take a look at the region...

Hmmm... the carbonyl oxygen atoms are a bit close to each other. Are there any other Ramachandran outliers in the area?

   - To measure the distance between atoms, hold M and click on an atom, then (with "M" still pressed) click on another atom.

2.26 A is too close.

  - To undisplay the distance, press the "C" key.

   - Click on **Models**,
   - then click on the **Rama Balls** button for this protein model.

_Moorhen displays Ramachandran balls that uses colour to represent the probability of those phi, psi angles for that residue_

Residue A180 has a red ball and just upstream at A178 the LYS has a orange ball. Maybe both of these can be improved.

Notice that there is a green blob close to the N of A180. Maybe it would be better if that blob were fitted by the carbonyl of the peptide. So let's try that:

## Flipping a Peptide

 - Click middle-mouse the N atom of residue A 180
 - Click right-mouse on an atom in the peptide (O atom of 179)

_Moorhen displays a grid with icons for modelling (with which you may already be familiar if you have used Coot)_

As you move the mouse over the icons in the toolbar, you will see a tooltip for that icon.

  - To flip the peptide you want to use the "Flip Peptide" button. Click it.


_Moorhen flips the Peptide and the Ramachandran ball for that residues turns green_

Yay. Progress. Let's see if we can do the same for peptide bond bettwen A177 and A178.

- Bring A177 CA to the centre of the screen
- Click right-mouse over a peptide and then
- Click on the "Flip Peptide" button

_Moorhen flips the peptide and the Ramachandran ball for 177 turrns green_

More progress. Good stuff.

  - Look at the Ramachandran Plot again. Notice that the red spots for the problematic residues have disappeared.

**Navigation tip 1**: Use middle mouse button click on an atom to put the centre of the residue at the centre of the screen. Use Alt left mouse to put the clicked atom at the centre of the screen.

**Navigate tip 2**: To go to a specific residue, use **Edit** &rarr; **Go to...** and in the "Atom selection" entry
type the residue selection, e.g. "//A/156" to mean residue 156 in the "A" chain.

## Real Space Refinement

  - Right click on an atom in residue A179
  - Click on the "Refine Residues" button (the icon is a target).

_Moorhen refines residues along the chain_

Now the local backbone fits quite nicely into the blue map (but the Rama balls may not be completely green).

## Connect the Maps: Updating Maps

Wouldn't it be nice though, if the difference map updated when you add atoms to green blobs or removed them from red blobs. A good correction of the model would then make those blobs disappear.

Let's try that:

  - **File** &rarr; **Connect molecule and map for updating...** &rarr; **OK**.  (No need to change the values in the options menus because they are already setup to be correct.)

_Moorhen will display a "toast" top left informing you of the current R-factor and the number of Moorhen points that you have collected (so far none, because we have just started)_

**Note**: Using updating/connected maps will slow down the model-building process somewhat but we now have the advantage of collection Moorhen points and watching the R-factor go down as we make changes. Moorhen points indicate progress in flattening the difference map.

## Difference Map Peaks

In the **Validation** dialog, the active tool in the Validation option menu is the **Ramachandran Plot**

 - Let's change that to **Difference Map Peaks**.

_Moorhen displays the difference map peaks in a waterfall plot_

On the left of the waterfall plot are the most positive peaks (and if there were any the most negative peaks would be displayed on the far right).

  - Here I find it useful to adjust the contour level to 0.64 (and 0.47 for the difference map).
  - Let's open the **Validation Tools** card again and click on the biggest/leftmost peak.

What are we looking at? An orange Ramachandran ball?

## Flipping... flipping

It's a flipped peptide.

  - So let's flip it to the correct orientation.

_As we do so, Moorhen makes several updates. It moves the model, it updates the maps in the light of the new model, it updates the Ramachandran balls so that they are both green now. And, in the toast, it updates the R-factor (a tiny amount) and gives us some Moorhen points._

Flipping a peptide to the correct orientation generally gives you 15-20 Moorhen points.

Have a look around... At the end of the helix, residues A198 has an orange Ramachandran ball and a green difference map peak along the peptide. Let's flip that too while we are here.

  - Click the "Flip Peptide" button and pick the "C" atom of A197.
    You can use the key-bding "Shift Q" to flip the peptide when you hover of an atom of the residue and see golden balls.
    

_As before, the Rama balls update, the maps update (the green blob disappears), the R-factor updates and we get more Moorhen points._

You will notice that that the Difference Map Peaks graph has been updated too - the leftmost peak has gone. Have a look at the next 12 or so peaks. What do you notice?

## Adding Waters

You will notice that they are mostly peaks of waters. We could add waters one by one, but a more
automated method is to do many at the same time.

  - **Calculate** &rarr; **Add waters...** &rarr; **OK**

This will add around 100 waters. And as above, the maps and the R-factors will update and we will get many Moorhen points. The map should improve a bit so that the ligand is more easy to make out.

## Contact Dots and Clashes

  - Use the sequence viewer to navigate to residue A194.

What this? It's a flipped peptide - let's flip it back to where it should be. But having done that, what do you notice? Let's use Moorhen's clash analysis:

  - In the Model window for the tutorial structure you will see a box labelled **Cont. dots** - click it.

_Moorhen display contact dots and clash interactions_

Wooh! Pink sticks. Bad news! So let's also flip the pepide on the next rung of the helix - residue A197.

_The contact dots between A194 and A197 disappear_

  - OK, you can turn off the contact dots for now by unclicking the **Cont. dots** box in the Model window.

## Change the Residue Type

Now navigate to residues A193. What do the maps tell you is going on here? What is the residue type? What does the model say? What does the map say?

OK, so first let's fill the side-chain with the atoms of the type from the main-chain atoms: "TYR" - use right-mouse
to get the context-sensitive menu and click on the "Auto-fit Rotamer" button (top left).

What do we see? What does that suggest?

It suggests that the sequence of the model doesn't match the sequence of the protein from which the data were collected. OK, so let's mutate it.

  - In the context-sensitive menu, click the "Mutate Residue" button (it's the radtion symbol), then from the residues type chooser currently "ALA (A)" choose "PHE (F)"

_Moorhen updates the map so that the red blob goes away_

   - Navigate to residue A168.
   
What do we see? What should it be instead?

## Mutate

OK, so let's mutate it:

  - Use the (right-mouse menu) "Mutate Residue" button change the type to "TYR (Y")

(More Moorhen Points - yay).

  - Let's go back to the **Difference Map Peaks** in the **Validation tools**

Now you can see negative peaks. Let's have a look at those.

Can you find a negative difference map peak that is close to resiude A187? Have a look at the model? What needs to be changed?

## Rotamers

In the Model window, click the box for "Rota. dodec"

_Moorhen displays rotamer dodecamers coloured by rotamer probability_

Yellow is not so high probability. Can we find something with higher probability than the current rotamer?

 - Let's change the rotamer using the "Auto-fit Rotamer" button.

_On improvement of the rotamer probability, Moorhen will change the colour of the dodecahedron to be more green_

Can you find a negative difference map peak close to residue A141? Again examine. What do we need to do? Let's fix the rotamer then using "Auto-fit Rotamer" as before.

## Fit the Ligand

OK, now it's time to fit the ligand!

  - Use the **Difference Map Peaks** to navigate to the ligand.

Several of the top 5 peaks should now correspond to the ligand.

[ligand](https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/web/tutorial/images/LZA-coot-render-v2.png)

  - **Ligand** &rarr; **Get Monomer** &rarr;
  - Add "LZA" in the "Monomer identifier" entry &rarr; **OK**

_Moorhen imports the LZA ligand_

OK, fine. Now let's undisplay it:

  - In the Model window, the bottom card should now be the card for the newly imported ligand ("#3 Mol LZA"). Click on the eye icon to undisplay the ligand

_Moorhen changes the icon to an uncrossed eye and the ligand disappears_

  - **Ligand** &rarr; **Find ligand...**

  - Check the options then click **Find**.

_Moorhen presents a new dialog saying that 1 possible ligand location_

To visualize the ligand, click on the "boxed circle" icon and the ligand will appear in the density.

It should be a reasonably close fit but not exact because the algorithm didn't use conformational variation.

To refine ligand click the "target" icon in that dialog.

To merge the ligand, use the "arrow" icon in that dialog.

_Moorhen updates the maps so that the difference map blobs change_

## Add a Water

  - Navigate to the nearby water peak using middle-mouse click and drag (or Shift-Alt Left-Mouse on a PC) to drag
    the view to the water blob.
  - You can also use the arrow-keys to pan the view

Let's add a water here

  - In **Edit** menu click on the **Add simple** menu item
  - Change the option menu to read "HOH" (if needed)
  - **OK**

  _Moorhen adds a water at the centre of the screen_

There are several water peaks in the map similar to this.

At some stage, when you add a water, you will see a the contours of negative density over part or all the water peak. What does that mean?

## More Validation Tools

  - Click the **Validation** button
  - Click on **Validation Plot**

_Moorhen displays interactive validation graphs._

## Over to You!

 - Click on the validation bars to navigate to interesting parts of the structure and make some fixes.

You should be able to collect about 1800 Moorhen points. Maybe more!

## Make a Pretty Picture

  - Using the cards in the drawer, undisplay the maps using the eye icon
  - Click on "Bonds" to undisplay the "Bonds" representation of the model
  - Likewise undisplay the Rama ball and Rota dodecs if you have the displayed
  - Click on **Ligands** to display the ligand in the model
  - Click on **Ribbons** to display the model in Ribbon mode

To navigate to the ligand:

  - **Ligand** &rarr; **Centre on ligand...**
  - Click "mol-1"
  - Click "Chain A"
  - Click the "/1/A/301(LZA)" label

  - **View** &rarr; **Scene Setting** &rarr; **Background Colour** - change it if you wish
  - Adjust the sliders for clipping and fogging to make the ligand more clearly visible
  - **OK**

Use keyboard "S" to activate the "in application" screen capture.


## Export Your Molecule

  - In the for the Model window, click the "download" icon

Then it's time to think about Reciprocal Space Refinement.

## Notes:

 The tutorial is based on 2vtq "Identification of N-(4-piperidinyl)-4-(2,6-dichlorobenzoylamino)-1H- pyrazole-3-carboxamide (AT7519), a Novel Cyclin Dependent Kinase Inhibitor Using Fragment-Based X-Ray Crystallography and Structure Based Drug Design" Wyatt, P.G. _et al._ (2008), _J. Med. Chem_. **51**, 4986.
