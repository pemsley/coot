---
layout: post
title:  "Moorhen"
date: Sun 16 Apr 20:20:01 GMT 2023
---

Moorhen is "Coot on the Web"!

Q: OK, so why "Moorhen"?

A: Althought "WebCoot" may have been a more obvious name for the project, we decided that "Moorhen" would
more fairly reflect that this project is the work of more than just me. The "meta" of the name is mean to be
interpretted as "quite like (a) Coot internally, but a bit different on the outside."

The Moorhen project is being developed, at the moment, by Filomeno
Rodriguez Sanchez, Stuart McNichollas, Martin Noble and me.

Q: So... what is Moorhen?

A: Moorhen is, at least at the moment, Coot in the web browser. There is no server!
(Other that the trivial server that sends the app to your web browser.) All the computation is done in
the browser.

Q: OK, go on.

A: Moorhen is made possible (or tractable) by [WebAssembly](https://webassembly.org/). WebAssembly is a
binary format designed to work on a virtual machine (inside your web browser).

Q: OK, how does that help?

A: From the beginning, Coot has had architecture that separated the comutational libraries from the
display. In 2022, this was concretized into a full-Coot like interface (in that the molecules are refered to
by index) but without using the GUI (which is complex and heavyweight, requiring Gtk+ and OpenGL).
This interface is called [libcootapi](https://pemsley.github.io/coot/blog/2023/02/26/libcootapi.html).

In retrospect this division was not as clean and precise as it might have been - in particular the
real-space refinement has an intimate and intricate interface between the computational components and the
part of the program that actually draws the molecules. The use of multi-threading made this interface more
complex still.

Nevertheless, I was able to create an interface that encapsulated a good deal of _Coot_ functionality.

Because it has been divorced from the standard GUI, **libcootapi** is available for additional/alternative
interfaces. Naturally, one of those is Python. With additional interface functions, another is Blender
(Coot+Blender will have its own blog post later).

Q: OK, I want to try it.

A: You can try Moorhen at [moorhen.org](https://moorhen.org).

