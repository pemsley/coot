# Copyright 2006 Joel Bard
# Copyright 2006 by Bernhard Lohkamp
# Copyright 2006 by Paul Emsley, The University of York

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
# Foundation, Inc.,  51 Franklin Street, Fifth Floor, 
# Boston, MA 02110-1301, USA

# Read in cns coeff-data (given filenames) and a pdb molecule filename to make
# maps.
#
def cns2coot(twofofc_coeffs, fofc_coeffs, model_pdb):
    import os, string

    # remove any trailing spaces of string
    # BL says: maybe somthing for coot-utils?
    # ..and has already a buikld-in function in python
    def trim_trailing_spaces(s):
        f = string.rstrip(s)
        return f

    # return [] on bad, else return a 6 memberes list.
    def get_cell_from_pdb(pdb_file):
        s = []
        fin = False
        try:
           fin = open(pdb_file,'r')
        except IOError:
           print("BL WARNING:: Cannot read ", pdb_file)
        if (fin):
           lines = fin.readlines()
           for line in lines:
               if "CRYST1" in line:
                  s = string.split(line)
                  s = s[1:7]
                  fin.close()
                  break
        else:
           print("BL WARNING:: IOError in ", pdb_file) 
        return s

    #  return "" on bad, else return a string that is the space group.
    #   
    #  have to go through all sorts of hoops because sometimes we might
    #  get P 21 and sometimes we might get P 1 21 1 for example
    #  sftools needs P21, can't handle P1211
    #  that being said - there must be a better way
    # set SYMM1=`awk 'BEGIN { FIELDWIDTHS="55 11 4" } /CRYST1/ {gsub(" ","",$2); coot_utils.printf "%s\n", $2}' $3`
    #
    def get_symm_from_pdb_file(pdb_file):

        s = ""
        fin = False
        try:
           fin = open(pdb_file,'r')
        except IOError:
           print("BL WARNING:: Cannot read ", pdb_file)
        if (fin):
           lines = fin.readlines()
           for line in lines: 
               if "CRYST1" in line:
                  s = string.split(line)
                  s = s[7:]
                  s = string.join(s,"")
                  fin.close()
                  break
        else:
           print("BL WARNING:: IOError in ", pdb_file)
        return s

    def get_symm_from_pdb_file_2(pdb_file):

        s = ""
        fin = False
        try:
           fin = open(pdb_file,'r')
        except IOError:
           print("BL WARNING:: Cannot read ", pdb_file)
        if (fin):
           lines = fin.readlines()
           for line in lines:
               if "CRYST1" in line:
                  s = string.split(line)
                  s = s[7:]
                  s = string.join(s)
                  fin.close()
                  break
        else:
           print("BL WARNING:: IOError in ", pdb_file)
        return s

    # Return a symmetry string
    # on problem, return False. Use that to exit (elsewhere).
    #
    def get_symm_from_symop_lib(symm2):
        
        try: 
           e = os.environ['CCP4']
        except:
           print("CCP4 not set up.")
           return False
        if (e):
           s = False
           e = os.path.join(e,"lib","data","symop.lib")
           fin = False
           try:
              fin=open(e,'r')
           except:
              print("Problems opening ", e)
           if (fin):
              lines = fin.readlines()
              for line in lines: 
                 if symm2 in line:
                    t = string.split(line)
                    s = t[3]
                    break
              fin.close()
           return s
        else:
           print("Problems with CCP4 set up")
        
    # main body
    map1_prefix = coot_utils.strip_extension(twofofc_coeffs)
    map2_prefix = coot_utils.strip_extension(fofc_coeffs)
    map1_tmp = map1_prefix + "_tmp.pdb"
    map2_tmp = map2_prefix + "_tmp.pdb"
    cell = get_cell_from_pdb(model_pdb)

    print("map1_prefix: ", map1_prefix)
    print("map2_prefix: ", map2_prefix)
    print("cell: ", cell)

    symm1 = get_symm_from_pdb_file(model_pdb)
    pdbset_log = "cns2coot-pdbset-tmp.log"

    print("symm1: ", symm1)
    coot_utils.popen_command("pdbset",["XYZIN",model_pdb,"XYZOUT",map1_tmp],["CELL " + string.join(cell), "SPACEGROUP " + symm1],pdbset_log,0)
    symm2 = get_symm_from_pdb_file_2(map1_tmp)
    print("symm2: ", symm2)
    symm = get_symm_from_symop_lib(symm2)

    if not(symm):
       print("Failed to find symm in symop.lib!")
    else:
       print("INFO:: SYMM is ", symm)
       map1_mtz = map1_prefix + ".mtz"
       map2_mtz = map2_prefix + ".mtz"
       map_coot_mtz = map1_prefix + "-coot.mtz"
       cad_log = "cns2coot-cad-tmp.log"
       sftools_1_log = "cns2coot-sftools-1-tmp.log"
       sftools_2_log = "cns2coot-sftools-2-tmp.log"

       if os.path.isfile(map1_mtz):
          os.remove(map1_mtz)
       if os.path.isfile(map2_mtz):
          os.remove(map2_mtz)
   
       coot_utils.popen_command("sftools", [] ,["read " + twofofc_coeffs, "cns", string.join(cell), symm, "END", "W", "P", "R", "SET LABELS", "FOM", "PHIC", "SCALE", "FWT", "PHWT", "WRITE " + map1_mtz] , sftools_1_log, 0)

       coot_utils.popen_command("sftools", [] ,["read " + fofc_coeffs, "cns", string.join(cell), symm, "END", "W", "P", "R", "SET LABELS", "FOM", "PHIC", "SCALE", "DELFWT", "PHDELWT", "WRITE " + map2_mtz] , sftools_2_log, 0)

       coot_utils.popen_command("cad", ["HKLIN1", map1_mtz, "HKLIN2", map2_mtz, "HKLOUT", map_coot_mtz],["LABIN FILE_NUMBER 1 E1=FOM E2=PHIC E3=FWT E4=PHWT", "LABIN FILE_NUMBER 2 E1=DELFWT E2=PHDELWT", "END"], cad_log, 0)

       if os.path.isfile(map1_mtz):
          os.remove(map1_mtz)
       if os.path.isfile(map2_mtz):
          os.remove(map2_mtz)

       # now load them in Coot
       #
       coot.read_pdb(model_pdb)
       coot.make_and_draw_map_with_reso_with_refmac_params(map_coot_mtz, "FWT", "PHWT", "", 0, 0, 0, "Fobs:None-specified", "SigF:None-specified", "RFree:None-specified", 0, 0, 0, -1.00, -1.00)
       coot.make_and_draw_map_with_reso_with_refmac_params(map_coot_mtz, "DELFWT", "DELPHWT", "", 0, 1, 0, "Fobs:None-specified", "SigF:None-specified", "RFree:None-specified", 0, 0, 0, -1.00, -1.00)

#cns_coot("2fo-fc_Fum480.coeff","fo-fc_Fum480.coeff","bindividual_Fum480.pdb")
