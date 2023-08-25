---
layout: post
title:  "GTK4 Coot Progress"
date: Sat 27 Aug 04:13:20 BST 2022
---

So the way forward for _Coot_ on MacOS is to use GTK4 (instead of Gtk3) for the GUI. And
that means a whole nother rewrite.  Oh boy. Oh well, in we go then (I had known that GTK4
was in Coot's future because it's the route to Vulkan and all those new graphical goodies
you can do on the graphics card now, but had not planned to start on it for a few years).

Turns out that GTK4 is new technology and there is not much documentation or help
especially in the areas of the Python API, OpenGL and deployment on MacOS. At times it has
been frustrating and slow-going.

Until very recently there was no way to edit the XML for the gui description for GTK4
applications other than "by hand" - i.e. there was no RAD Tool (because Glade doesn't work
for GTK4).  Had you asked me three months ago if it was possible to edit the GUI by
editing its XML file I would have said "No way! Are you crazy?".

But after a couple of months, it turns out that it is nowhere near as hideous as I had
thought though.

Maybe Cambalache can help. I am looking forward to trying out Cambalache - But at the
moment, I can't work how to install the new version without reinstalling the OS. So that
has stalled until I sort that out.

So... what has happen in the last couple of months?

 - I have got Python booting, and got PyGObject working (it was extremely tricky to
   compile) and they are big steps forward.  Now the Python-based gui tools are in play.
   
 - The Basic framebuffer works and we can draw and see molecules.
 
 - (However the effects framebuffers (_i.e._ "Fancy" style) are proving to be
   frustratingly difficult to get working.)
 
 - Martin Noble's MoleculesToTriangles has been incorporated.
 
 - The gestures and event controllers have been rewritten from scratch (mostly complete
   now).  That means we can recentre, rotate, zoom and click on atoms.

 - The gui buttons are about 80% done. Refinement works and Delete works. Rotation,
   Translation and Chi angles do not (yet).
   
 - I am trying to reduce the number of dialogs as I go. New interfaces will be revealers or
   overlays as much as possible. If there is going to be a dialog, it will be a transient
   (no more dialog windows hidden under the main window).
 
 - The interface has changed to "select first" - _i.e._ bring the residue of interest to the
   centre of the screen and act on it with the tool of choice. No residue picking needed
   (except for residue range operations). You heard it here first.
   
 - While the menu system has been rewritten from scratch (there are 170 menus), only a few
   of the actions have been connected up. There are around 1500 dialog widget event controller
   to either edit or rewrite.
 
![Select First Functions]({{"../../../images/SelectFirstRefinement.png"}})


 You can follow the progress of this project - it's the **gtk4** branch of _Coot_ on github.
 
 You can build binaries on MacOS using Homebrew quite easily using
 cootgtk4.rb. See **Issue 33** on _Coot_ github page if you want to try
 that.
