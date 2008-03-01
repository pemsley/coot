# Copyright 2008 by The University of York
# Author: Bernhard Lohkamp
# Copyright 2007, 2008 by The University of Oxford
# Author: Paul Emsley

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA

import unittest
import os

class NcsTestFunctions(unittest.TestCase):

    def test01_0(self):
        """NCS chains info"""

        # should return False
        ncs_chain_info = ncs_chain_ids(-1)
        self.failIf(ncs_chain_info, "   Fail: ncs-chains returns %s, should be False" %ncs_chain_info)

        # a normal case
        make_ncs_ghosts_maybe(imol_rnase)
        ncs_chain_info = ncs_chain_ids(imol_rnase)
        self.failUnless(ncs_chain_info, "   Fail: ncs-chain-ids returns False")
        self.failUnless(len(ncs_chain_info) > 0, "   Fail: ncs-chains returns %s" %ncs_chain_info)
        first_ghost = ncs_chain_info[0]
        self.failUnless(len(first_ghost) > 1, "    Fail: first-ghost %s" %first_ghost)
        self.failUnless(type(first_ghost[0]) is StringType, "    Fail: not strings first-ghost %s" %first_ghost[0])
        self.failUnless(type(first_ghost[1]) is StringType, "    Fail: not strings first-ghost %s" %first_ghost[1])
        print "   NCS info: ", ncs_chain_info

    def test02_0(self):
        """NCS deviation info"""

        # should return False
        ncs_chain_info = ncs_chain_differences(-1, "XX")
        self.failIf(ncs_chain_info, "   Fail: ncs-chains returns %s, should be False" %ncs_chain_info)

        # should return False for insulin
        ncs_chain_info = ncs_chain_info(imol_insulin, "A")
        self.failIf(ncs_chain_info, "   Fail: ncs-chains for insulin returns %s, should be False" %ncs_chain_info)

        # a normal case
        make_ncs_ghosts_maybe(imol_rnase)
        ncs_chain_info = ncs_chain_differences(imol_rnase, "A")
        self.failUnless(ncs_chain_info, "   Fail: ncs-chain-differences returns False")
        self.failUnlessEqual(len(ncs_chain_info), 3,
                             """   Fail on length: length ncs-chain-differences should be 3 is %s\n
                                     ncs-chain-differences returns %s""" %(len(ncs_chain_info), ncs_chain_info))
        
