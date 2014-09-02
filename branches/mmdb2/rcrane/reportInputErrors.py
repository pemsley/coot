#!/usr/bin/env python
"""Common functions for reporting errors in rotamerize and extendChain input."""

# Copyright 2012 Kevin Keating
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

STANDARD_BASES = frozenset(["A", "G", "C", "U"]) #the names of the four standard RNA bases
    #used to make sure that the user isn't trying to rotamerize a modified nucleotide

def reportPDB2Error(resNum):
    """Report that a residue contains PDB2-formatted atom names
    
    ARGUMENTS:
        resNum - the number of the residue containing PDB2 atom names
    RETURNS:
        None
    """
    
    errorMsg = ["This molecule appears to use PDB2 naming (nucleotide " + str(resNum) + " contains a C1* atom).",
                    "RCrane requires PDB3 naming (e.g. C1' instead of C1*).",
                    "This PDB file can be updated using Remediator",
                    "(http://kinemage.biochem.duke.edu/software/remediator.php)."]
    for curline in errorMsg:
        print curline
    #add_status_bar_text("This molecule appears to use PDB2 naming.  RCrane requires PDB3 naming.  Cannot rotamerize.")
    
    #Since this isn't an intuitive problem (i.e. it's far less obvious for the user than "you clicked on
    #two different chains"), we create a dialog box instead of a status bar message
    errorMsg[0] += " " #we want two spaces between sentences
    errorMsg[1] += " "
    createRCraneErrorDialog(" ".join(errorMsg) + "\n\nRotamerization canceled.")


def reportInsCodeError(resNum):
    """Report that a residue contains an insertion code
    
    ARGUMENTS:
        resNum - the number of the residue containing the insertion code
    RETURNS:
        None
    """
    
    print "Residue " + str(resNum) + " contains an insertion code.  RCrane currently cannot rotamerize nucleotides with insertion codes."
    add_status_bar_text("Residue " + str(resNum) + " contains an insertion code.  RCrane currently cannot rotamerize nucleotides with insertion codes.")
    
    
def reportModifiedNucError(resNum, resType):
    """Report that a residue contains an insertion code
    
    ARGUMENTS:
        resNum - the number of the residue containing the modified nucleotide
    RETURNS:
        None
    """
    
    print "Nucleotide " + str(resNum) + " is a modified nucleotide (" + str(resType) + ").  RCrane does not currently support rotamerizing modified nucleotides."
    add_status_bar_text("Nucleotide " + str(resNum) + " is a modified nucleotide (" + str(resType) + ").  RCrane does not currently support rotamerizing modified nucleotides.")
