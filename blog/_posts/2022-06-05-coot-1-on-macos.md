---
layout: post
title:  "Coot 1 on MacOS"
date: Sun 5 Jun 2022 10:33:23 BST
---

The release of _Coot 1_ was met with some enthusiasm from MacOS users.
The Homebrew formula for _Coot_, written by Yoshitaka Moriwaki was very helpful to that end.

![Coot 1 on MacOS]({{"../../../images/coot-1-macos-screeny.png"}})


As you may know, MacOS is not my main platform, and it was only (relatively)
recently that I got a Mac laptop that was able to compile _Coot 1_. Once I had done so,
it became apparent was that the performance (i.e. frames/second) of _Coot 1_
on MacOS was terrible - close to unusable - if not actually unusable.

Googling around, it seemed that other GTK-using applications that were
targeting MacOS (such as the GIMP or Inkscape) were also very
slow. This [Phoronix
article](https://www.phoronix.com/scan.php?page=news_item&px=GTK4-macOS-Improved)
was very interesting. It seems that Christian Hergert is the main
developer of [the MacOS
backend](https://blogs.gnome.org/chergert/2020/12/15/gtk-4-got-a-new-macos-backend-now-with-opengl/).

Should I care that Coot 1 on MacOS is slow?

I decided that I should care.

I thought it would be interesting to see if moving to GTK4 would solve the performance problems for _Coot_.

So I extracted the new OpenGL code in _Coot_, i.e. the "Coot Render Engine" and built that into a new
application with a rudimentary gui that would compile with Gtk+3 and, with some changes, also compile with GTK4.

I made a [video about the results](https://www.youtube.com/watch?v=_c7NO3_8KNc).

It was enough to show me that it was the GTK backend that was the problem for MacOS and a solution
was available.

So... how to move _Coot_ to GTK4? A quick change to configure.ac to enable the `--with-gtk4` flag
and let's see what happens.

Many, many compilation errors, that's what.

OK, so a 30,000 line patch and 5 days later, _Coot_ be compiled with GTK4.

So what's different?
 - packing widgets has changed,
 - handling the mouse and keyboard input events has changed,
 - color usage has changed,
 - window placement has changed,
 - there is no canvas,
 - menus have changed,
 - there is no longer a toolbar or menubar,
 - radio buttons have gone,
 - drag and drop has changed,
 - Glade doesn't work with GTK4 - so the GUI has to be developed by editing an XML
   file by hand - _ugh!_ (the Glade XML file for Gtk+3 Coot is 30,000 lines or more)
 - Python deadlocks on startup (and so was commented out).

I should emphasise that this is not a useful build of _Coot_. I just
hacked on it until it compiled and displayed something. To make that
happen I used `#ifdef` for many hundreds of blocks of code. To get _Coot_
back into a useful state, those code blocks will need to be converted to the new
system and entered back into the compilation path.  Python will need to be restored

So the work to
make the update is substantial (several months?). It's worth doing. I had in
mind that I would need to change to GTK4 at some stage, but it seems
that that is going to happen sooner than I expected.

Lucrezia Catapano has made [a new tutorial video using Coot
1](https://www.youtube.com/watch?v=Xhonm4K1y0c) in its current
state. The graphics are (to my eye) pretty, but slow.


