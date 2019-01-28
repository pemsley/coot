
*Coot*
----

*Coot* is a toolkit for Macromolecular Crystallography and
model-building.  *Coot* uses widgets (with the gui builder glade),
mmdb, clipper, and OpenGL, together with a new approach to map
contouring and importing/creation and other modelling and building
operations.

Blog
----

[Coot Development Blog](https://pemsley.github.io/coot/ "Coot Development Blog")

Redistribution
--------------

This software package (coot and its relatives) is Free Software.

This software is substantially dependent on Free and Open Source
Software is distributed under the GNU General Public Licence version 3.
You are able to use and redistribute this software with *your* 
software (provided you also distribute this source code).

You may find *Coot* distributed with CCP4 Software.  Be aware that Coot
is not covered under the CCP4 Software Licence, either part (i) or
part (ii).

This is Free Software - anyone referring to it as "Open Source" shall
be eaten by a gnu.

Prerequisites:
-------------

 To build this software from scratch, you will need:

 * a C++ compiler that knows about the STL,
 * gtk+ version 2, 
 * gtkglext version 1.2,
 * mmdb-2,
 * fftw-2.1.5, 
 * clipper and its dependencies (its dependencies are ccp4 libs, mmdb, fftw)
 * OpenGL: (GL, glu, glut)

 And optionally (coot is better if you have these):
 
 * gsl-1.3 or later (http://www.mirror.ac.uk/sites/ftp.gnu.org/gnu/gsl)
   [This is quite important, you can't regularize or refine without this] 
 * goocanvas
 * guile-1.8.x
   (http://www.mirror.ac.uk/sites/ftp.gnu.org/gnu/guile/) [scripting]
 * guile-gui-0.2.tar.gz (http://www.ossau.uklinux.net/guile/) [interactive scripting]
 * guile-www-1.1.5 or later
   [OCA Server interface]
 * guile-gtk (ftp://ftp.gnu.org/gnu/guile-gtk) [interactive scripting]
 * goosh-1.3 or later (http://arglist.com/guile/) 
   [external program (e.g. refmac, lsqman) interface]
 * SSMlib [superposition]

 (you can also get many of these from http://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/dependencies)

Binaries
--------

Binaries for Coot are availble from here:

http://www2.mrc-lmb.cam.ac.uk/Personal/pemsley/coot/binaries/

Compiling
---------

*Coot* is tricky to compile by hand and has several dependencies. The
build-it script (in this directory) is how I build binaries.  You are
encouraged to use it if you want to compile *Coot*.


Idea Thanks
-----------

   I have a memory of a conversation with Mike Hartshorn - which
   planted the seed of the idea of writing a free model-building
   program in about 1991.  It took many years for this idea to bear
   fruit however. [Mike went on to write AstexView more rapidly].

   Many of the ideas in Coot are of course influenced by Alwyn Jones'
   O (and to a lesser extent Frodo).  I am grateful to him also
   because it in the late 80s that was sitting in front of a PS300
   running Frodo in York that I though "wow, this is what I want to
   do", without him and Frodo, I wouldn't be a crystallographer. 
   Thanks Alwyn.

   Many evening sessions with Adrian Lapthorn undoubtedly put some
   flesh on the bones of the ideas - certainly in those days we formed
   many opinions about The Right Way to do things.  Some of them are
   now present in Coot.

   Some credit for the way the interactive refinement works should go
   to Warren DeLano.  At a CCP4 Study Weekend (2003?) he presented
   "sculpt" mode in PyMol.  The audience (and I!) was really wowed.
   Until that point, I'd rejected the idea of such interactivity
   because it was not reproducable in the history (not easily,
   anyway).  But then I saw Warren's demo and I thought, "we've go to
   have that something like that too!"  It then took me a couple of
   days to bolt it into the existing real space refinement code.
   Later on I was demonstrated and impressed by PyMol's
   "view"s.  18 months later I added a view system to Coot.


Code Thanks...
---------

  particularly go to:

  * Kevin Cowtan for many things, including clipper, and explaining to me
    how it works :)
  * Eugene Krissinel for mmdb, 
  * Eugene Krissinel and Kim Henrick for libssm,
  * Stuart McNicholas for MGTree code,
  * Alexi Vagin for Refmac cif dictionaries,
  * And Mr Behind-the-scenes-man (who allowed the whole thing to happen):  
    Keith Wilson,

  Finally, thanks to those who have emailed advice[1], the Coot
  testers[2] and particularly those who tried to compile the source
  code[3] and sent patches[4]: they suffered so that you don't have
  to.

  Thank you all - I have very much enjoyed doing Coot.

[1] Ethan Merritt, Gerard Kleywegt, George Sheldrick

[2] Eleanor Dodson, Miguel Ortiz Lombardia, Charlie Bond,
    Jan Dohnalek, Garib Murshudov, Jean Wittingham, Florence Vincent,
    Tracy Gloucester, Constantina Fotinou, Roberto Steiner, Adrian
    Lapthorn, Claudia Schnick, Rosa Grenha, Ezra Peisach, Ben Luisi, 
    Frank von Delft, Karen McLuskey, Marcin Cymborowski, Stephen Graham.

[3] William G. Scott, Bernhard Lohkamp, Alex Schuettlekopf, Luca Jovine, 
    Bob Nolte.

[4] Ezra Peisach, Charlie Bond, Mike Hartshorn

---

Paul Emsley
Now in Cambridge: pemsley at mrc dash lmb dot cam dot ac dot uk
