import glob
import os
import subprocess

def get_top_lines(file_name):
    f = open(file_name)
    lines = f.readlines()
    return lines[:18]

def get_needs_whole_block(lines):
    for line in lines:
        if '//' in line:
            if not '://' in line:
                return True
        if '#define' in line:
            return True
        if '#include' in line:
            return True
        if '#shader' in line:
            return True
        if '/* ' in line:
            return False
    return False

def get_copyright_year(fn):

    year = "1900"
    result = subprocess.run(['git', 'log', fn], capture_output=True, text=True)
    lines = result.stdout.split("\n")
    for line in lines:
        if 'Date: ' in line:
            # print(line)
            parts = line.split()
            year = parts[5]
    print(fn, year)
    return year

def get_copyright_license_lines(fn):
    cry = get_copyright_year(fn)
    author = "Paul Emsley"
    author_line = " * Author: " + author + "\n"
    file_line = " * " + fn + "\n"
    crh = "University of York"
    if cry > "2012":
        crh = "Medical Research Council"
    copyright_line = " * Copyright " + cry + " by " + crh + "\n"
    new_lines = [ "/*\n",
                  file_line,
                  " *\n",
                  copyright_line,
                  author_line,
                  " *\n",
                  " * This file is part of Coot\n",
                  " *\n",
                  " * This program is free software; you can redistribute it and/or modify\n",
                  " * it under the terms of the GNU Lesser General Public License as published\n",
                  " * by the Free Software Foundation; either version 3 of the License, or (at\n",
                  " * your option) any later version.\n",
                  " *\n",
                  " * This program is distributed in the hope that it will be useful, but\n",
                  " * WITHOUT ANY WARRANTY; without even the implied warranty of\n",
                  " * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU\n",
                  " * Lesser General Public License for more details.\n",
                  " *\n",
                  " * You should have received a copies of the GNU General Public License and\n",
                  " * the GNU Lesser General Public License along with this program; if not,\n",
                  " * write to the Free Software Foundation, Inc., 51 Franklin Street,\n",
                  " * Fifth Floor, Boston, MA, 02110-1301, USA.\n",
                  " * See http://www.gnu.org/licenses/\n",
                  " *\n",
                  " */\n" ]
    return new_lines


# this f is a file descriptor
def write_block_to_file(f, fn):
    lines = get_copyright_license_lines(fn)
    for line in lines:
        f.write(line)

def add_block_to_file(fn):
    tmp_file = fn + ".tmp"
    f = open(tmp_file, 'w')
    write_block_to_file(f, fn)
    fr = open(fn)
    lines = fr.readlines()
    for line in lines:
        f.write(line)
    fr.close()
    if os.path.exists(tmp_file):
        os.rename(tmp_file, fn)
        # pass

source_dir = 'lidia-core'
source_dir = 'ligand'
source_dir = 'mini-mol'
source_dir = 'pli'
source_dir = 'skeleton'
source_dir = 'validation-graphs'
source_dir = 'src'

for fn in glob.glob(source_dir + "/*"):
    if fn == 'coot-utils/json.hpp': continue
    if '.tmp' in fn: continue
    if 'ligands-2016.db' in fn: continue
    if 'side-chain-data.tar.gz' in fn: continue
    if os.path.isdir(fn): continue
    print("processing", fn)
    top_lines = get_top_lines(fn)
    needs_whole_block = get_needs_whole_block(top_lines)
    if needs_whole_block:
        add_block_to_file(fn)

for fn in glob.glob(source_dir + "/*"):
    if '.tmp' in fn: continue
    if 'ligands-2016.db' in fn: continue
    if 'side-chain-data.tar.gz' in fn: continue
    if os.path.isdir(fn): continue
    top_lines = get_top_lines(fn)
    for line in top_lines:
        if "ll rights reserved" in line:
            print("WARNING::", line.strip(), fn)
