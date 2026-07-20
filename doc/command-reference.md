# Coot command reference

The commands below are typed into the **Command** tab of the Python/AI terminal. Text is matched case-insensitively and extra whitespace is ignored. Where a command takes a model or map number, omitting it acts on the *active* molecule.

> This file is generated from the command definitions (`python/coot_commands/`). Do not edit by hand - run `python3 -m coot_commands.docs` to regenerate.

## Building

### `add solvent SO4 here`

Add a monomer (by 3-letter code) at the screen centre.

Examples:

- `add solvent SO4 here`
- `add ligand GOL`

Fetches the named monomer from the dictionary and moves it to the screen centre. Use "add water" for a single water atom.

### `add water here`

Add a water at the screen centre.

Examples:

- `add water here`
- `add water`

Places a water atom at the screen centre (the pointer position), like the 'Place Atom At Pointer' tool.


## Delete

### `delete all`

Ask to confirm closing every loaded map and model.

Examples:

- `delete all`
- `close everything`

Asks for confirmation before closing every loaded map and model. Type 'delete all confirm' to go ahead.

### `delete all confirm`

Close every loaded map and model (confirmed).

Confirms and carries out 'delete all', closing every loaded map and model.

### `delete map 1`

Delete (close) a map.

Examples:

- `delete map 1`
- `close map 1`

Closes the given map, freeing its molecule number.

### `delete model 0`

Delete (close) a model.

Examples:

- `delete model 0`
- `close model 0`

Closes the given model, freeing its molecule number.

### `delete molecule 2`

Delete (close) a molecule of either kind.

Examples:

- `delete molecule 2`
- `close mol 2`

Closes the given molecule (map or model) by number.


## Display

### `show only active`

Show only the active model, hiding the others.

Hides every model except the active one.

### `hide residue environment`

Hide the residue environment distances.

Examples:

- `hide residue environment`
- `hide environment`

Hides the residue environment distances.

### `hide map 1`

Hide (undisplay) a map.

Examples:

- `hide map 1`
- `undisplay map`

With no map number, acts on the active map.

### `hide model 0`

Hide (undisplay) a model.

Examples:

- `hide model 0`
- `undisplay model`

With no model number, acts on the active model.

### `show residue environment`

Show the environment distances around a residue.

Examples:

- `show residue environment`
- `show environment A/89`

Shows the environment distances (contacts and H-bonds) around a residue, centring on it. With no residue named, uses the active residue (the one at the centre of the screen).

### `show map 1`

Show (display) a map.

Examples:

- `show map 1`
- `display map`

With no map number, acts on the active map.

### `show model 0`

Show (display) a model.

Examples:

- `show model 0`
- `display model`

With no model number, acts on the active model.


## Fetch

### `fetch 3GP`

Fetch a monomer from the CCD dictionary.

Examples:

- `fetch 3GP`
- `fetch monomer ATP`

Fetches a monomer (by 3-letter code) from the CCD dictionary. Four-or-more character codes are treated as PDB accessions.

### `fetch 1abc`

Fetch coordinates from the PDBe.

Examples:

- `fetch 1abc`
- `fetch pdb 4hhb`

Fetches coordinates from the PDBe. Add "and map" to also fetch maps, or "from pdb-redo" for the re-refined version.

### `fetch 1abc and map`

Fetch coordinates and maps (via EDS).

Examples:

- `fetch 1abc and map`
- `fetch map for 1abc`

Fetches coordinates and the 2Fo-Fc / Fo-Fc maps from the Electron Density Server (EDS).

### `fetch 1abc from pdb-redo`

Fetch a re-refined model and maps from PDB-REDO.

Examples:

- `fetch 1abc from pdb-redo`
- `fetch pdb-redo 1abc`

Fetches the re-refined model and maps from PDB-REDO.


## Help

### `help`

List the available commands, grouped by category.


## Labels

### `clear labels`

Remove all atom labels.

Examples:

- `clear labels`
- `remove all labels`


## Ligand

### `fit ligand LIG`

Fit a ligand into density by its three-letter code.

Examples:

- `fit ligand LIG`
- `fit ligand ATP here`
- `find ligand NAG into model 0 map 1`

Fits a ligand (by three-letter code) into density, the same as Ligand > Find Ligands. Uses the active model and the refinement map by default. Add 'here' to search only at the current view centre - handy after 'go to blob N'. Returns the molecule numbers of the fitted solutions.


## Maps

### `colour map 1 blue`

Set a map's colour.

Examples:

- `colour map 1 blue`
- `colour map cyan`

With no map number, acts on the active map.

### `contour map 1 to 0.35`

Set a map's absolute contour level.

Examples:

- `contour map 1 to 0.35`
- `contour to 0.3`

Sets the absolute contour level (map units). With no map number, acts on the active map. Add "sigma" to contour in RMSD units instead.

### `contour map 1 to 1.5 sigma`

Set a map's contour level in sigma.

Examples:

- `contour map 1 to 1.5 sigma`
- `contour to 1.2 rmsd`

Sets the contour level in sigma (map RMSD units). With no map number, acts on the active map.

### `map 1 is a difference map`

Mark a map as a difference map.

Examples:

- `map 1 is a difference map`
- `make map 2 a difference map`

Marks the map as a difference map (green/red, contoured either side of zero). With no map number, acts on the active map.


## Model editing

### `add altconf A/72`

Add an alternate conformation to a residue.

Examples:

- `add altconf A/72`
- `add alt conf A 72`
- `add altconf`

Adds an alternate conformer to the named residue. With no residue named, acts on the active residue (the one at the centre of the screen).

### `add OXT to A/89`

Add a terminal OXT atom to a residue.

Examples:

- `add OXT to A/89`
- `add OXT to A 89`
- `add OXT`

Adds a terminal OXT oxygen to a residue (usually a chain's C-terminus). With no residue named, acts on the active residue (the one at the centre of the screen); with no model number, the active model.

### `add terminal residue to A/89`

Add a terminal residue onto the end of a chain.

Examples:

- `add terminal residue to A/89`
- `add terminal residue`
- `add terminal residue to A 89 as ALA`

Builds a new residue onto the end of a chain, attached to the named terminal residue and fitted against the refinement map. The type defaults to 'auto' (guessed from any sequence); add 'as ALA' to force one. With no residue named, acts on the active residue; with no model number, the active model.

### `autofit A/72`

Auto-fit the best rotamer for a residue.

Examples:

- `autofit A/72`
- `auto-fit rotamer A 72`
- `autofit`

Auto-fits the best rotamer for the residue against the refinement map (open a map first). With no residue named, acts on the active residue (the one at the centre of the screen).

### `backrub rotamer A/89`

Apply a backrub rotamer fit to a residue.

Examples:

- `backrub rotamer A/89`
- `backrub A 89`
- `backrub rotamer`

Applies a backrub rotamer fit to the named residue. With no residue named, acts on the active residue (the one at the centre of the screen).

### `delete residue A/72`

Delete a single residue.

Examples:

- `delete residue A/72`
- `delete residue A 72`
- `delete residue`

Deletes the named residue. With no residue named, acts on the active residue (the one at the centre of the screen).

### `pepflip A/89`

Flip the peptide following a residue.

Examples:

- `pepflip A/89`
- `pep flip A 89`
- `pepflip`

Flips the peptide bond following the named residue by 180. With no residue named, acts on the active residue (the one at the centre of the screen).

### `replace residue A/72 with ALA`

Mutate a residue to another type.

Examples:

- `replace residue A/72 with ALA`
- `mutate residue A 72 to GLY`
- `mutate residue to ALA`

Mutates the residue to the given (1- or 3-letter) type. With no residue named, acts on the active residue (the one at the centre of the screen).


## Models

### `merge model 1 and model 2`

Merge one model into another.

Merges the second model into the first, so the combined model keeps the first model's number.

### `superpose model 0 onto model 1`

Superpose one model onto another (SSM).

Examples:

- `superpose model 0 onto model 1`
- `superpose model 0 onto model 1 in place`

Superposes the source model onto the target (reference) model by secondary-structure matching (SSM). By default a superposed copy is made, leaving the source model untouched; add 'in place' to move the source model itself instead.


## Navigation

### `centre at 12.0 4.5 -3.2`

Centre the view on an x, y, z position.

Centres the view on the given orthogonal coordinates.

### `go to A 45`

Centre the view on a chain/residue.

Examples:

- `go to A 45`
- `centre on A/45`
- `go to model 0 B 12`

Centres on the given chain/residue, picking a sensible atom (CA for protein, P/C1' for nucleotides, and so on) so it works for waters, ligands and nucleic acids too. With no model number, acts on the active model. Chain and residue may be separated by a space or a slash.

### `next residue`

Centre on the next residue.

Examples:

- `next residue`
- `next`
- `forward residue`

Centres on the next residue after the active one (like the space bar). The active residue is the one at the centre of the screen.

### `previous residue`

Centre on the previous residue.

Examples:

- `previous residue`
- `prev`
- `back residue`

Centres on the residue before the active one (like shift-space). The active residue is the one at the centre of the screen.


## Refinement

### `refine`

Real-space refine the active residue.

Examples:

- `refine`
- `refine active residue`
- `refine here`

Real-space refines the active residue (the one at the centre of the screen) against the refinement map.

### `refine b factors`

Refine the atomic B-factors (ADPs) of a model.

Examples:

- `refine b factors`
- `refine b-factors of model 0`
- `refine adps`

Refines the atomic B-factors (ADPs) of the whole model using the shiftfield method, against the reflection data. Needs a map with reflection data set as the refinement map (unlike real-space refinement, this is a whole-molecule reciprocal-space step). With no model number, acts on the active model.

### `refine chain A`

Real-space refine an entire chain.

Examples:

- `refine chain A`
- `refine chain B`

Real-space refines a whole chain (its full residue range) against the refinement map. With no model number, acts on the active model.

### `refine A 45 to 50`

Real-space refine a range of residues.

Examples:

- `refine A 45 to 50`
- `refine A/45-50`
- `real space refine A 100 to 105`

Real-space refines the given residue range against the refinement map. With no model number, acts on the active model.

### `refine A 45`

Real-space refine a single residue.

Examples:

- `refine A 45`
- `refine A/45`
- `refine residue B 12`

Real-space refines a single residue against the refinement map. With no model number, acts on the active model.

### `refine sphere A/89`

Real-space refine the sphere around a residue.

Examples:

- `refine sphere A/89`
- `refine sphere A 89`
- `refine sphere`

Real-space refines the sphere of atoms around the named residue against the refinement map. With no residue named, acts on the active residue (the one at the centre of the screen).


## Representation

### `colour carbons coloured`

Use element-coloured carbons for a model.

Examples:

- `colour carbons coloured`
- `coloured carbons`
- `coloured carbons for model 0`

Uses per-element carbon colouring. With no model number, acts on the active model.

### `colour carbons grey`

Use grey carbon colours for a model.

Examples:

- `colour carbons grey`
- `grey carbons`
- `grey carbons for model 0`

Uses grey for carbon atoms. With no model number, acts on the active model.

### `hide ribbons`

Hide the ribbon representation of a model.

Examples:

- `hide ribbons`
- `hide ribbons for model 0`

Removes the ribbon representation added by 'show ribbons'. With no model number, acts on the active model.

### `hide surface`

Hide the molecular surface of a model.

Examples:

- `hide surface`
- `hide surface for model 0`

Removes the molecular surface added by 'show surface'. With no model number, acts on the active model.

### `hide symmetry`

Hide symmetry-related molecules.

### `show ribbons`

Show a ribbon representation of a model.

Examples:

- `show ribbons`
- `show ribbons for model 0`

Draws a ribbon (cartoon) representation of the whole model, coloured by chain. With no model number, acts on the active model.

### `show surface`

Show a molecular surface of a model.

Examples:

- `show surface`
- `show surface for model 0`

Draws a molecular surface for the whole model, coloured by chain. With no model number, acts on the active model.

### `show symmetry`

Show symmetry-related molecules.


## Session

### `close sequence`

Close the sequence view for a model.

Examples:

- `close sequence`
- `close sequence of model 0`

Closes the sequence view opened by 'open sequence'. With no model number, acts on the active model.

### `list maps`

List the loaded maps.

Lists the loaded maps with their molecule number and name, marking which are difference maps.

### `list models`

List the loaded models.

Examples:

- `list models`
- `list molecules`

Lists the loaded models with their molecule number and name.

### `load tutorial`

Load the tutorial model and data.

Examples:

- `load tutorial`
- `load tutorial model and data`

Loads the bundled tutorial model and its data (map coefficients), the same as File > Open Tutorial.

### `open sequence`

Open the sequence view for a model.

Examples:

- `open sequence`
- `open sequence of model 0`

Opens the sequence view for the model. Close it again with 'close sequence'. With no model number, acts on the active model.


## Settings

### `get bond thickness`

Report the current default bond thickness.

Examples:

- `get bond thickness`
- `what is the bond thickness`

### `get contour level`

Report the current contour level of a map.

Examples:

- `get contour level`
- `what is the contour level of map 1`

With no map number, acts on the active map.

### `get contour step`

Report the contour-level scroll step for normal maps.

Examples:

- `get contour step`
- `what is the contour step`

### `get difference map contour step`

Report the contour-level scroll step for difference maps.

Examples:

- `get difference map contour step`
- `what is the difference map contour step`

### `get font size`

Report the current on-screen label font size.

Examples:

- `get font size`
- `what is the font size`

### `get map radius`

Report the current map display radius.

Examples:

- `get map radius`
- `what is the map radius`

### `get map sampling rate`

Report the current map sampling rate.

Examples:

- `get map sampling rate`
- `what is the map sampling rate`

### `set bond thickness to 4`

Set the default bond (stick) thickness.

Examples:

- `set bond thickness to 4`
- `set bond thickness 3`

Line/stick thickness for model bonds, in pixels. A whole number; the default is 5.

### `set contour level to 1.5 sigma`

Set the contour level of a map (absolute, or in sigma).

Examples:

- `set contour level to 1.5 sigma`
- `set contour level of map 1 to 0.3`
- `set contour level of map 1 to 2 sigma`

Sets the contour level of a map, in absolute units by default or in sigma when the value ends with 'sigma'. With no map number, acts on the active map.

### `set contour step to 0.1`

Set the contour-level scroll step for normal maps.

Examples:

- `set contour step to 0.1`
- `set contour step 0.05`

How much a scroll changes the contour level of a normal (non-difference) map. See 'set difference map contour step' for difference maps.

### `set difference map contour step to 0.1`

Set the contour-level scroll step for difference maps.

Examples:

- `set difference map contour step to 0.1`
- `set diff map contour step 0.05`

How much a scroll changes the contour level of a difference map.

### `set font size to 2`

Set the on-screen label font size.

Examples:

- `set font size to 2`
- `set font size 3`

Size of on-screen labels: 1 (small), 2 (medium) or 3 (large).

### `set map radius to 20`

Set the map display radius (Angstroms).

Examples:

- `set map radius to 20`
- `set map radius 15`

Radius, in Angstroms, of the sphere of density drawn around the screen centre. Larger radii are slower to contour.

### `set map sampling rate to 1.8`

Set the map sampling rate for maps read from now on.

Examples:

- `set map sampling rate to 1.8`
- `set map sampling 2`

Finer sampling (higher rate) makes smoother maps at the cost of memory. Applies to maps read after it is set; typical values are 1.5-2.5.

### `set updating maps on`

Turn Coot's auto-updating (sfcalc) maps on or off.

Examples:

- `set updating maps on`
- `set updating maps off`

Turns Coot's auto-updating maps on or off. Coot's auto-updating maps recompute the 2Fo-Fc and difference maps from the reflection data as you edit the model, so density follows the atoms. Needs the refinement map (a map with reflection data attached) and a difference map to be loaded - 'load tutorial' or opening an MTZ provides both. Only one updating-maps session can run at a time.

### `stop updating maps`

Stop Coot's auto-updating (sfcalc) maps.

Turns off auto-updating maps (same as 'set updating maps off'). Coot's auto-updating maps recompute the 2Fo-Fc and difference maps from the reflection data as you edit the model, so density follows the atoms. Needs the refinement map (a map with reflection data attached) and a difference map to be loaded - 'load tutorial' or opening an MTZ provides both. Only one updating-maps session can run at a time.


## State

### `undo`

Undo the last action.

Undoes the last modification.


## Validation

### `check cis peptides`

Count cis peptide bonds in a model.

Examples:

- `check cis peptides`
- `cis peptides`

Counts cis peptide bonds. Cis peptides before proline are common and usually fine; other cis peptides are worth checking. With no model number, acts on the active model.

### `check clashes`

List steric clashes (atom overlaps) for a model.

Examples:

- `check clashes`
- `atom overlaps`

Lists atom overlaps (steric clashes) with a clash volume above 2 A^3. With no model number, acts on the active model.

### `check gln and asn`

Open the Gln/Asn side-chain flip (B-factor outlier) analysis.

Examples:

- `check gln and asn`
- `gln asn outliers`

Opens the Gln/Asn B-factor outlier analysis, flagging glutamine and asparagine side chains that may need a 180 degree flip. With no model number, acts on the active model.

### `check missing atoms`

List residues with missing atoms in a model.

Examples:

- `check missing atoms`
- `missing atoms`

Lists residues that are missing modelled atoms. With no model number, acts on the active model.

### `check non-standard residues`

List non-standard residue types in a model.

Examples:

- `check non-standard residues`
- `non-standard residues`

Lists residue types that are not standard amino acids or water - ligands, modified residues and the like. With no model number, acts on the active model.

### `check ramachandran`

List Ramachandran outliers for a model.

Examples:

- `check ramachandran`
- `ramachandran outliers`

Lists residues in improbable regions of the Ramachandran plot (probability below 0.02). With no model number, acts on the active model.

### `check rotamers`

List the least probable rotamers for a model.

Examples:

- `check rotamers`
- `rotamer outliers`

Reports the least probable side-chain rotamers (low probability = unusual). With no model number, acts on the active model.

### `check waters`

List highly-coordinated waters (possible ions) in a model.

Examples:

- `check waters`
- `highly coordinated waters`

Lists waters with an unusually high number of close contacts (coordination number 5 or more within 3.2 A), which may be misassigned ions. With no model number, acts on the active model.

### `go to blob 1`

Centre the view on a blob from the last 'find blobs' result.

Examples:

- `go to blob 1`
- `centre on blob 2`

Centres the view on one of the blobs from the most recent 'find blobs' command, numbered from 1 (largest first). Run 'find blobs' first to populate the list.

### `open validation`

Open the interactive validation overlay for a model.

Examples:

- `open validation`
- `validation overlay`

Opens Coot's interactive validation overlay (Ramachandran, rotamers, density fit and more) for a model against the active map. With no model number, acts on the active model. With no map number, acts on the active map.

### `difference map peaks`

Mark peaks in the difference map (missing/wrong density).

Examples:

- `difference map peaks`
- `difference map peaks above 4 sigma`

Marks peaks in the Fo-Fc difference map - candidate sites for missing atoms, waters or ligands, and for parts of the model in wrong density. Uses a 4 sigma cut-off by default. Needs a difference map to be loaded. With no model number, acts on the active model.

### `validate anomalies`

Summarise model outliers (Ramachandran, clashes, C-beta, chirals).

Examples:

- `validate anomalies`
- `find outliers`

Summarises geometry outliers - Ramachandran improbables, atom overlaps (clashes), C-beta deviations and chiral volume errors - as text. With no model number, acts on the active model.

### `find blobs`

Summarise unmodelled density blobs (candidate build sites).

Examples:

- `find blobs`
- `check unmodelled blobs`
- `find blobs above 2 sigma`

Lists peaks of density not accounted for by the model - candidate sites for waters, ligands or unbuilt residues. The search masks out density within 1.9 A of the model and finds the peaks left over, so it belongs on a difference (mFo-DFc) map: by default it uses the loaded difference map at a 3.0 sigma cut-off. On a 2mFo-DFc map almost everything is above a low cut-off, so blobs would appear all over - name a map with 'using map N' only if you mean to. Add e.g. 'above 2 sigma' to change the threshold.


## View

### `rock`

Toggle rocking the view.

Toggles idle rocking; issue again to stop.

### `background black`

Set the background colour.

Examples:

- `background black`
- `set background colour to white`

Colour names: black, white, grey, and the other named colours accepted by colour commands.

### `orthographic`

Switch between orthographic and perspective projection.

Examples:

- `orthographic`
- `perspective view`

### `zoom to 30`

Set the view zoom factor.

Examples:

- `zoom to 30`
- `set zoom 50`

Larger numbers zoom out. Typical range ~10-100.

### `spin`

Toggle spinning the view.

Toggles idle spinning; issue again to stop.

### `fullscreen`

Toggle fullscreen mode.

Toggles fullscreen; issue again to leave fullscreen.
