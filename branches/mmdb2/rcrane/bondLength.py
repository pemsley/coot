#!/usr/bin/env python
"""List the maximum length of a bond that Coot will draw as connected."""

# Copyright 2010, 2011, 2012 Kevin Keating
# 
# Licensed under the Educational Community License, Version 2.0 (the
# "License"); you may not use this file except in compliance with the
# License. You may obtain a copy of the License at
# 
# http://www.osedu.org/licenses/ECL-2.0
# 
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an "AS IS"
# BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
# or implied. See the License for the specific language governing
# permissions and limitations under the License.

from coot import coot_version

LONG_BOND_DIST_CUTOFF = 2.16

#Coot 0.6.2 and 0.7 have different distance cut-offs for deciding whether or not to draw a bond
if coot_version()[0:3] == "0.6":
    BOND_DIST_CUTOFF  = 1.64 #value for the Coot's 0.6.2 branch
else:
    BOND_DIST_CUTOFF  = 1.71 #value for the Coot trunk (i.e. 0.7 and newer)

BOND_LIST = [(" P  ", " OP1", LONG_BOND_DIST_CUTOFF),   #the P-OP1 and P-OP2 bonds must be the first two elements of this list
             (" P  ", " OP2", LONG_BOND_DIST_CUTOFF),   #since we use BOND_LIST[0:2] as the bond list for the last nucleotide
             (" P  ", " O5'", LONG_BOND_DIST_CUTOFF),
             (" C1'", " C2'", BOND_DIST_CUTOFF),
             (" C2'", " C3'", BOND_DIST_CUTOFF),
             (" C3'", " C4'", BOND_DIST_CUTOFF),
             (" C3'", " O3'", BOND_DIST_CUTOFF),
             (" C4'", " C5'", BOND_DIST_CUTOFF),
             (" C4'", " O4'", BOND_DIST_CUTOFF),
             (" O4'", " C1'", BOND_DIST_CUTOFF),
             (" C5'", " O5'", BOND_DIST_CUTOFF),
             (" C2'", " O2'", BOND_DIST_CUTOFF)]

# Note that when using BOND_DIST_CUTOFF, bonds greater than or equal to the distance are not drawn.
# When using LONG_BOND_DIST_CUTOFF, bond greater than the distance are not drawn,
# but bonds equal to the distance *are* drawn.
# For now, we ignore that distinction, since getting a bond length of *exactly* 1.64 (or 1.71) should be fairly uncommon

# Also note that the P-O3' bond uses the BOND_DIST_CUTOFF, not LONG_BOND_DIST_CUTOFF
# This is presumably a bug/oversight in Coot due to the fact that it's an inter-residue bond.


#a list of all bonds (excluding hydrogens) expected for backbone and nucleosides
BOND_LIST_FULL = {"A": [], "G": [], "C": [], "U": [], "backbone": []}

for (atom1, atom2, bondLength) in BOND_LIST:
    atom1 = atom1.strip()
    atom2 = atom2.strip()
    BOND_LIST_FULL["backbone"].append((atom1, atom2))


BOND_LIST_FULL["A"] = [("C1'", "N9"),
                       ("N9",  "C8"),
                       ("C8",  "N7"),
                       ("N7",  "C5"),
                       ("C5",  "C4"),
                       ("C4",  "N9"),
                       ("C5",  "C6"),
                       ("C6",  "N1"),
                       ("N1",  "C2"),
                       ("C2",  "N3"),
                       ("N3",  "C4"),
                       ("C6",  "N6")]


BOND_LIST_FULL["G"] = [("C1'", "N9"),
                       ("N9",  "C8"),
                       ("C8",  "N7"),
                       ("N7",  "C5"),
                       ("C5",  "C4"),
                       ("C4",  "N9"),
                       ("C5",  "C6"),
                       ("C6",  "N1"),
                       ("N1",  "C2"),
                       ("C2",  "N3"),
                       ("N3",  "C4"),
                       ("C6",  "O6"),
                       ("C2",  "N2")]

BOND_LIST_FULL["C"] = [("C1'", "N1"),
                       ("N1",  "C2"),
                       ("C2",  "N3"),
                       ("N3",  "C4"),
                       ("C4",  "C5"),
                       ("C5",  "C6"),
                       ("C6",  "N1"),
                       ("C4",  "N4"),
                       ("C2",  "O2")]

BOND_LIST_FULL["U"] = [("C1'", "N1"),
                       ("N1",  "C2"),
                       ("C2",  "N3"),
                       ("N3",  "C4"),
                       ("C4",  "C5"),
                       ("C5",  "C6"),
                       ("C6",  "N1"),
                       ("C4",  "O4"),
                       ("C2",  "O2")]