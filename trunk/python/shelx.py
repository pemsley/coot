# Copyright 2004, 2005 by The University of York
# Copyright 2005 by Bernhard Lohkamp
 
#;;; This program is free software; you can redistribute it and/or modify
#;;; it under the terms of the GNU General Public License as published by
#;;; the Free Software Foundation; either version 2 of the License, or (at
#;;; your option) any later version.
 
#;;; This program is distributed in the hope that it will be useful, but
#;;; WITHOUT ANY WARRANTY; without even the implied warranty of
#;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#;;; General Public License for more details.
 
#;;; You should have received a copy of the GNU General Public License
#;;; along with this program; if not, write to the Free Software
#;;; Foundation, Inc.,  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


#;; direction is either 'forwards or 'backwards
#;; 
#;; start-resno is higher than stop-resno if we are building backwards
#;; 
#;; (fit-gap 0 "A"  23 26)   ; we'll build forwards
#;; (fit-gap 0 "A"  26 23)   ; we'll build backwards
#;;

# BL python adaption of shelx.scm by Paul Emsley

#;; Provide functions to handle shelx fcf file (which I believe to be old-style
#;; cif rather than mmCIF).
#;; 
#;; So, in Coot, if it sees the extention of .fcf then it calls the
#;; fcf handler and the handler does the conversion, writes a new
#;; mmCIF file and then reads it.
#;; 

# BL says: Where do we get SG from? Otherwise clipper can't cif file!
# what shall we do???
# otherwise fine, if we include SG manually, e.g
# _symmetry.space_group_name_H-M  'C 2 2 21'


def handle_shelx_fcf_file(filename):
    import os

    filename = os.path.normpath(filename)

    #; OK, so we have been passed a shelx fcf file.
    #; We need to run this through the awk filter.
    # BL says: I guess in python it's easier to do it within python

    output_mmCIF_file_name = filename + ".cif"

#    print "BL DEBUG:: filename is", filename
#    print "BL DEBUG:: fnd cif filename is", output_mmCIF_file_name
    convert_shelx_fcf_to_cif(filename,output_mmCIF_file_name)

    read_cif_data_with_phases(output_mmCIF_file_name)



def convert_shelx_fcf_to_cif(fcf_filename,cif_filename):

   import string, math

   intable = 0

   try:
    fin = open(fcf_filename,'r')
   except IOError:
    print "BL INFO:: Cannot read ", fcf_filename
   try:
    fout = file(cif_filename,"w")
   except IOError:
    print "BL INFO:: Cannot write ", cif_filename
   if (fin and fout):
    lines = fin.readlines()
    for line in lines:
       if "_shelx_title" in line:
          # dont write
          pass
       elif "_shelx_refln_list_code" in line:
          # dont write
          pass
       elif "_shelx_F_calc_maximum" in line:
          # dont write
          pass
       elif "_exptl_crystal_F_000" in line:
          # dont write
          pass
       elif "_reflns_d_resolution_high" in line:
          # dont write
          pass
       elif "_shelx_refln_list_code" in line:
          # dont write
          pass
       elif "_shelx_F_calc_maximum" in line:
          # dont write
          pass
       elif "_exptl_crystal_F_000" in line:
          # dont write
          pass
       elif "_reflns_d_resolution_high" in line:
          # dont write
          pass
       else:
        if "_refln_F_squared_meas" in line:
         line = line.replace("_refln_F_squared_meas"," _refln.F_meas")
        if "_refln_F_squared_sigma" in line:
         line = line.replace("_refln_F_squared_sigma"," _refln.F_meas_sigma")
        # symm equivs:
        if " _symmetry_equiv_pos_as_xyz" in line:
          line = line.replace(" _symmetry_equiv_pos_as_xyz","_symmetry_equiv.pos_as_xyz")
        if "_cell_length_a" in line:
          line = line.replace("_cell_length_a","_cell.length_a     ")
        if "_refln_F_squared_meas" in line:
         line = line.replace("_refln_F_squared_meas"," _refln.F_meas")
        if "_refln_F_squared_sigma" in line:
         line = line.replace("_refln_F_squared_sigma"," _refln.F_meas_sigma")
        # symm equivs:
        if "_symmetry_equiv_pos_as_xyz" in line:
          line = line.replace("_symmetry_equiv_pos_as_xyz","_symmetry_equiv.pos_as_xyz")
        if "_cell_length_a" in line:
          line = line.replace("_cell_length_a","_cell.length_a     ")
        if "_cell_length_b" in line:
          line = line.replace("_cell_length_b","_cell.length_b     ")
        if "_cell_length_c" in line:
          line = line.replace("_cell_length_c","_cell.length_c     ")
        if "_cell_angle_alpha" in line:
          line = line.replace("_cell_angle_alpha","_cell.angle_alpha  ")
        if "_cell_angle_beta" in line:
          line = line.replace("_cell_angle_beta","_cell.angle_beta   ")
        if "_cell_angle_gamma" in line:
          line = line.replace("_cell_angle_gamma","_cell.angle_gamma  ")
        if "_refln_index_h" in line:
          line = line.replace("_refln_index_h"," _refln.index_h")
        if "_cell_length_b" in line:
          line = line.replace("_cell_length_b","_cell.length_b     ")
        if "_cell_length_c" in line:
          line = line.replace("_cell_length_c","_cell.length_c     ")
        if "_cell_angle_alpha" in line:
          line = line.replace("_cell_angle_alpha","_cell.angle_alpha  ")
        if "_cell_angle_beta" in line:
          line = line.replace("_cell_angle_beta","_cell.angle_beta   ")
        if "_cell_angle_gamma" in line:
          line = line.replace("_cell_angle_gamma","_cell.angle_gamma  ")
        if "_refln_index_h" in line:
          line = line.replace("_refln_index_h"," _refln.index_h")
        if "_refln_index_k" in line:
          line = line.replace("_refln_index_k"," _refln.index_k")
        if "_refln_index_l" in line:
          line = line.replace("_refln_index_l"," _refln.index_l")
        if "_refln_F_calc" in line:
          line = line.replace("_refln_F_calc"," _refln.F_calc")
        if "_refln_phase_calc" in line:
          line = line.replace("_refln_phase_calc"," _refln.phase_calc")
          intable = 1
        if (intable == 1 and len(string.split(line))>5):
          col = string.split(line)
          col[3] = str(math.sqrt(float(col[3])))
#BL says: dunno if that is conversion Paul intends, but otherwise division by 0
          if (math.sqrt(float(col[3]))) == 0.0 :
             col[4] = str(0.5 * float(col[4]))
          else:
             col[4] = str(0.5 * (float(col[4]) + 0.0) / math.sqrt(float(col[3])))
          line = string.join(col) + "\n"
        fout.write(line)

    fin.close()
    fout.close()
   else:
    print "BL INFO:: IO error in conversion!"


#handle_shelx_fcf_file("test.fcf")
