# Copyright 2007, 2008 by The University of York
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

rnase_pir     = os.path.join(unittest_data_dir, "rnase.pir")
poly_ala_frag = os.path.join(unittest_data_dir, "crashes_on_cootaneering.pdb")

class CootaneerTestFunctions(unittest.TestCase):

    def test01_0(self):
        """Cootaneer Beta Strand"""
        imol_model = read_pdb(poly_ala_frag)
        imol_map = make_and_draw_map(rnase_mtz, "FWT", "PHWT", "", 0, 0)

        self.failUnless(os.path.isfile(rnase_pir), "missing rnase pir file")

        assign_sequence_from_file(imol_model, rnase_pir)

        set_rotation_centre(64.271, 7.036, 14.42)

        n_atom = closest_atom(imol_model)
        self.failUnless(n_atom, "missing closest atom")
        imol     = n_atom[0]
        chain_id = n_atom[1]
        resno    = n_atom[2]
        inscode  = n_atom[3]
        at_name  = n_atom[4]
        alt_conf = n_atom[5]
        print "   Cootaneering: imol %s chain-id %s resno %s inscode %s at-name %s alt-conf %s" \
              %(imol, chain_id, resno, inscode, at_name, alt_conf)
        cootaneer(imol_map, imol, [chain_id, resno, inscode, at_name, alt_conf])
