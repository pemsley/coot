RCrane, version 0.11

OVERVIEW

RCrane allows for semi-automated building of RNA structure within Coot.  This
is done using the methodology described in:
    Keating KS and Pyle AM.  Semiautomated model building for RNA
      crystallography using a directed rotameric approach.  Proc Natl Acad
      Sci USA, 107: 8177-8182 (2010).
Please note that all publications resulting from the use of RCrane must
acknowledge this manuscript.

RCrane is copyright 2010-2012, Kevin Keating, and is licensed under the
Educational Community License, Version 2.0.


REQUIREMENTS

RCrane requires Coot 0.6.2 or newer, and Coot must have been built with GTK2
and with Python scripting enabled.

If you do not already have an appropriate version of Coot installed, it may be
downloaded from http://lmb.bioch.ox.ac.uk/coot/software/binaries/releases/.  If
installing on Linux, you will need a build of Coot with "-python-gtk2" in the
file name (i.e. coot-0.6.2-binary-Linux-i686-ubuntu-10.04.2-python-
gtk2.tar.gz).  If installing on Windows, use the WinCoot build (i.e.
WinCoot-0.6.2.exe).  Versions of Coot for OS X may be found at
http://sage.ucsc.edu/~wgscott/xtal/wiki/index.php/Stand-Alone_Coot.


INSTALLATION

RCrane is included with Coot 0.7 or newer.  To launch RCrane, click on
Extensions -> RCrane launch.

If you are using Coot 0.6.2, you must install RCrane separately.  To do this,
first unzip rcrane.zip.  If you are installing RCrane for all users of the
computers, it is recommended that you put the RCrane directory within the Coot
installation directory (e.g. /opt/coot/rcrane or C:\WinCoot\rcrane).  If you
are installing RCrane for only yourself, you may put the rcrane directory
anywhere within your home directory (i.e. /home/name/Documents/rcrane or
C:\Users\name\Documents\rcrane).


TO RUN RCRANE FOR THE CURRENT COOT SESSION

RCrane is included with Coot 0.7 or newer.  To launch RCrane, click on
Extensions -> RCrane launch.

If you are using Coot 0.6.2, click on Calculate -> Run Script...  Then choose
launch.py within the RCrane directory.  This will create an RCrane menu.


TO AUTOMATICALLY LAUNCH RCRANE EVERY TIME COOT IS STARTED

If you are using Coot 0.7, edit the rcrane.py file within your .coot-
preferences directory.  In Linux and OS X, this file is located at ~/.coot-
preferences/rcrane.py.  In Windows, this file is located within the WinCoot
directory at WinCoot/.coot-preferences/rcrane.py.  If the file does not exist,
create it.  Add the line
rcane.createRCraneMenuWithCheck()
to the file.
If you are using Coot 0.6.2, click on Calculate -> Run Script...  Then choose
installRCrane.py within the RCrane directory.


USING RCRANE - TRACING A NEW STRUCTURE

To start a new RCrane trace, first move the center of the screen near where you
would like to build the first phosphate, then select New Trace 5'->3'... from
the RCrane menu.  An RCrane window will appear and allow you to select the
starting phosphate location.  Once the desired phosphate is selected, click on
the Accept button.  RCrane will then attempt to trace the first nucleotide of
the chain.  If this is not the nucleotide you wished to build, then you may
change the phosphate and base coordinates using the appropriate buttons.  Once
you are happy with the current nucleotide, click Accept Nt and RCrane will
attempt to trace another nucleotide.

After you have finished tracing a chain or chain fragment, accept the final
nucleotide using the Accept Nt button and then click on Build Backbone.  RCrane
will then calculate the backbone atomic coordinates for the traced nucleotides.
After coordinate calculation is complete, you will be able to review each suite
and select alternate conformers where appropriate.


USING RCRANE – ROTAMERIZING AN EXISTING STRUCTURE

To rebuild a portion of an existing structure, select Rotamerize Existing
Structure... from the RCrane menu, then click on atoms in the first and last
suite you would like to rebuild.  RCrane will calculate new backbone
coordinates for the specified suites, and then allow you to review each suite
and select alternate conformers where appropriate.

Note that when you specify the region to rebuild, you are selecting suites and
not nucleotides.  Therefore, if you click on the phosphate of nucleotide 18,
you have selected suite 18, which spans from the sugar of nucleotide 17 to the
sugar of nucleotide 18.

Also note that rotamerize will not work on portions of a molecule with
insertion codes or alternate locations.  Additionally, rotamerize currently
clears all fixed atoms and custom restraints on the specified molecule when
run.

If you would like to rebuild an existing structure without taking electron
density into account, select Rotamerize Without Density... instead.  This is
primarily intended for use when building theoretical models (i.e. when no
electron density is available).  Note that Rotamerize Without Density... is
only available when running a version of Coot newer than r3728 (0.7-pre).


KEYBOARD SHORTCUTS

Keyboard shortcuts can be assigned for starting a new trace or rotamerizing an
existing structure.  To set the keyboard shortcuts, see CUSTOMIZING RCRANE
below.  By default, no keyboard shortcuts are set.


CUSTOMIZING RCRANE

There are a small number of options that control how RCrane operates.  Note
that you do not need to set any of these options to use RCrane with the default
settings.  To change the options, first install RCrane (see TO AUTOMATICALLY
LAUNCH RCRANE EVERY TIME COOT IS STARTED above).  Next, edit the rcrane.py file
within your .coot-preferences directory.  In Linux and OS X, this file is
located at ~/.coot-preferences/rcrane.py.  In Windows, this file is located
within the WinCoot directory at WinCoot/.coot-preferences/rcrane.py.  If
rcrane.py contains a run_python_script line, all setting commands must be added
*after* that line.

If you would like to add a keyboard shortcut for a new 5' to 3' trace, then add
the line
add_key_binding("New RCrane trace 5'->3'", "N", rcrane.newTrace5to3)
and replace the "N " with the desired keyboard shortcut.  Note that using an
uppercase letter will create a shortcut for Shift+letter.

If you would like to define a custom keyboard shortcut for a new 3' to 5'
trace, then add the line
add_key_binding("New RCrane trace 3'->5'", "M", rcrane.newTrace3to5)
and replace the "M " with the desired keyboard shortcut.  Note that using an
uppercase letter will create a shortcut for Shift+letter.

If you would like to create a keyboard shortcut for the rotamerization of
existing structures, then add the line
add_key_binding("RCrane rotamerize", "J", rcrane.newRotamerize)
and replace the "J " with the desired keyboard shortcut.  Note that using an
uppercase letter will create a shortcut for Shift+letter.

If you would like to create a keyboard shortcut for the rotamerization of
existing structures without using the electron density, then add the line
add_key_binding("RCrane rotamerize without density", "K",
rcrane.newRotamerizeWithoutDensity)
and replace the "K " with the desired keyboard shortcut.  Note that using an
uppercase letter will create a shortcut for Shift+letter.

If you are using Coot 0.7 or newer and would like RCrane to automatically
launch with Coot, add the line
rcane.createRCraneMenuWithCheck()

By  default, RCrane limits the number of nucleotides that can be rotamerized at
once (default: 20).  This is done to prevent excessively long run times due to
an erroneous mouse click.  To increase the limit, add the line
rcrane.setRotamerizeMaxNucleotides(40)
and replace "40" with the desired limit.  To remove the limit, add the line
rcrane.setRotamerizeMaxNucleotides(-1)

By default, RCrane uses REFMAC5 parameters for ideal bond lengths and angles
during minimizations.  These are the same parameters used by Coot during
refinement or regularizement.  To instead use the Phenix/MolProbity pucker-
specific parameters in RCrane, add the line
rcrane.enablePhenixRestraints()
Note that this requires Coot 0.7-pre r3926 or later.  To enable the
Phenix/MolProbity restraints for the current session only, open Coot’s Python
scripting window (Calculate -> Scripting... -> Python...) and type the command
rcrane.enablePhenixRestraints()
To revert to the REFMAC5 parameters, type the command
rcrane.disablePhenixRestraints()
To determine the current parameter set, type the command
rcrane.usingPhenixRestraints()
A return value of True indicates that the Phenix parameters are being used.  A
return value of False indicates that the REFMAC5 parameters are being used.


INSTALLING RCRANE FOR ALL USERS

RCrane is included with Coot 0.7 and newer.  If you are using Coot 0.6.2 on
Linux and OS X, installing RCrane will only install it for the current user.
Additional users must also follow the steps in TO AUTOMATICALLY LAUNCH RCRANE
EVERY TIME COOT IS STARTED above.  Alternatively, RCrane can be installed for
all users of a computer by:
1. Install RCrane for the current user by following the instructions in TO
   AUTOMATICALLY LAUNCH RCRANE EVERY TIME COOT IS STARTED above.
2. Create a pyextras directory within the Coot installation directory (e.g.
   /opt/coot/pyextras).
3. Move ~/.coot-preferences/rcrane.py to the newly created pyextras directory.
4. Create the file /etc/profile.d/rcrane.sh containing the line:
       export COOT_PYTHON_EXTRAS_DIR=/opt/coot/pyextras
   (using the appropriate path for the pyextras directory).
5. Any users who are currently logged in must log out and log in again.

RCrane will now be automatically launched with Coot for all users of the
computer.


UNINSTALLING RCRANE

RCrane is included with Coot 0.7 and newer and cannot be uninstalled.  If you
are using Coot 0.6.2, you may prevent RCrane from running when Coot is started
by deleting the rcrane.py file within your .coot-preferences directory.  In
Linux and OS X, this file is located at ~/.coot-preferences/rcrane.py.  In
Windows, this file is located within the WinCoot directory at WinCoot/.coot-
preferences/rcrane.py.  To completely uninstall RCrane from Coot 0.6.2, delete
the RCrane directory as well.
