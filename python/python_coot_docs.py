#    python_coot_docs.py
#    Copyright (C) 2008  Bernhard Lohkamp, The University of York
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

# BL script to make texi files from coot python modules

def is_empty_file(file_name):

    print "to be done?"

def make_texi_file(python_file, doc_file):

    import os

    def extracted_function(line_no):
        fcn = lines[line_no][4:-1]
        lb = fcn.find("(")
        # to get a space between the funcn and '('
        # will result in better readability of html output
        if (fcn[lb-1] != " "):
            fcn = fcn.replace("(", " (")
        # check for multiple line functions
        count = 0
        while ("(" in fcn and not ")" in fcn):
            count += 1
            # remove white space too (maybe - FIXME - check)
            lb = lines[line_no+count].lstrip()
            # remove EOL '\n' if exists from fcn
            fcn = fcn[0:-1]
            fcn += lb
        return fcn
    
    def extract_info(line_no):
        extracted_function = lines[line_no][4:-1]
        for i in range(line_no - 1, 0, -1):
            if '#' in lines[i][0]:
                pass
            else:
                write_info(i, line_no)
                break

    def write_info(start_line, stop_line):
        docout.write("\n")
        docout.write("@c BL python scripted from " + python_file + ":" + str(start_line) + "\n")
        docout.write("@deffn procedure " + extracted_function(stop_line) + "\n")
        for i in range(start_line, stop_line):
            line = lines[i][2:]
	    if line == "":
                docout.write("\n")
            else:
                docout.write(line)
        docout.write("\n")
        docout.write("@end deffn\n")

    pyin = open(python_file, 'r')
    docout = open(doc_file, 'w')
    lines = pyin.readlines()
    for i in range(len(lines)):
        line = str(lines[i])
        if "def" in line[0:3]:
            extract_info(i)
        # This is for xtra documentations!
	elif "# deftexi" in line[0:9]:
            lines[i] = "def" + line[9:]
            extract_info(i)
        else:
            pass

    pyin.close()
    docout.close()

def make_section_file(section_file, base_name):

    import os

    so = open(section_file, 'w')
    so.write("@node    " + base_name + "\n")
    so.write("@section " + base_name + "\n")
    so.write("@cindex  " + base_name + "\n")
    so.close()

def append_file(master_file_name, append_file_name):

    master = open(master_file_name, 'a')
    append = open(append_file_name, 'r')
    lines = append.readlines()
    master.writelines(lines)
    master.close()
    append.close()

def remove_file(file_name):

    import os

    os.remove(file_name)

import glob, os

source_files = glob.glob("*.py")
# take out coot.py if there
try:
	source_files.remove("coot.py")
	source_files.remove("python_coot_docs.py")
except: 
	pass

all_doc_file_name = "coot-python-functions.texi"

# first make the menu
def make_menu(file_name):

    doco = open(file_name, 'w')
    doco.write("@menu\n")  
    
    for source_file in source_files:
        base_name, ext = os.path.splitext(source_file)
        line = "* " + base_name + "::\n"
        doco.write(line)

    doco.write("@end menu\n")

    doco.close()

make_menu(all_doc_file_name)

# make the texi files and sections
for source_file in source_files:

    base_name, ext = os.path.splitext(source_file)
    doc_file =  base_name + ".texi"
    section_file = base_name + "-section.texi"
    make_texi_file(source_file, doc_file)
    make_section_file(section_file, base_name)

    append_file(all_doc_file_name, section_file)
    append_file(all_doc_file_name, doc_file)

    remove_file(section_file)
    remove_file(doc_file)
    
