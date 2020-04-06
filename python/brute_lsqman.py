#   brute_lsqman.py
#   Copyright (C) 2008 by Bernhard Lohkamp, The University of York
#   Copyright (C) 2003 by Charlie Bond, The University of Dundee
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


# brute-lsqman - run lsqman on chains A of two pdb-files and read in
# the result to coot. Charlie Bond 2003.
# Can keep a count of the number of successful runs if necessary

global lsqman_count
lsqman_count = 0

def brute_lsqman(pdb1_imol, pdb2_imol):

    import subprocess
    global lsqman_count

    pdbout_filename = "coot-lsqman-" + str(lsqman_count) + ".pdb"
    
    lsqman_command = "lsqman"
    command_line_args = [" -b "]
    data_lines = ["re m1 coot-tmp1.pdb",
                  "re m2 coot-tmp2.pdb",
                  "brute m1 a m2 a 50 25 100",
                  "imp m1 * m2 *", 
                  "app m1 m2",
                  "wr m2 ", pdbout_filename,
                  "quit"]
    lsqman_log = "coot-lsqman" + str(lsqman_count) + ".log"

    write_pdb_file(pdb1_imol, "coot-tmp1.pdb")
    write_pdb_file(pdb2_imol, "coot-tmp2.pdb")
    status = popen_command(lsqman_command, command_line_args, data_lines, 
              lsqman_log, True)
    if (status == 0): # lsqman ran OK
        handle_read_draw_molecule(pdbout_filename)
        lsqman_count += 1
    else:
        print("lsqman failed - sorry. I don't know what to say.")
