---
layout: post
title:  "Gtk4 Coot Update"
date: Sun 23 Apr 07:36:06 BST 2023
---

After something of a hiatus, I have refreshed my efforts on GTK4-based Coot!

To that end Jakub Smulski has been recruited to the team, funded by
Global Phasing. Jakub has been working on the the two major stumbling
blocks in the conversion of Coot to use GTK4.  Those are (or were!)
the validation graphics and the ligand-building gui.

Those two sub-projects, as well as the Ramachandran plot viewer and
the sequence viewer are/were/had been based on GooCanvas. Goocanvas
though is not available in GTK4, so these sub-projects all need a
rewrite from scratch.

We have now done a sizable chunk of that work. Hence this update.
Here for example is the validation graphs in action:

![Coot with Validation Graph]({{"../../../images/gtk4-coot-with-validation-graph.png"}})

We are working on the ligand-building GUI and the Ramachandran plot
now (the plan for the Ramachandran plot is to use OpenGL rather than
GTK).

You can follow the progress of this project by tracking the `gtk4`
branch of coot on GitHub.

