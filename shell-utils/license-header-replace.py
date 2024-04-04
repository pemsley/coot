
import sys
import os

def write_lines(lines, file_name):

    f = open(file_name, "w")
    try:
        for l in lines:
            f.write(l)
    except TypeError as e:
        print(e, "with l", l)
    f.close()

def make_copyright_line(person, copyright_year):

    pp = person
    if person == "Martin": pp = "Martin Noble, University of Oxford"
    if person == "Kevin":  pp = "Kevin Cowtan & University York"

    line = " * Copyright " + copyright_year + " by " + pp + "\n"
    return line

def insert_new_header(person, file_name, copyright_year, lines):

    copyright_line = make_copyright_line(person, copyright_year)
    file_line = " * " + file_name + "\n"
    author_line = " * Author: " + person + "\n"
    if person == "Martin": author_line = " * Author: Martin Noble\n"

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
                  " * General Public License for more details.\n",
                  " *\n",
                  " * You should have received a copy of the GNU General Public License\n",
                  " * along with this program; if not, write to the Free Software\n",
                  " * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA\n",
                  " * 02110-1301, USA\n",
                  " */\n" ]

    el = enumerate(new_lines)
    for idx,new_line in el:
        lines.insert(idx, new_line)

def get_copyright_year(line):

    y = ""
    parts = line.split()
    # now find the first part after "Copyright" that contains "20"
    idx_copyright = -1
    for idx,part in enumerate(parts):
        if part == "Copyright":
            idx_copyright = idx
        if part == "(C)":
            idx_copyright = idx
    if idx_copyright != -1:
        for idx,part in enumerate(parts):
            if idx > idx_copyright:
                if "20" in part:
                    y = part
    if "Created on" in line: y = parts[-1]

    return y

def get_author(line):
    parts = line.split()
    author = parts[-2] + " " + parts[-1]
    return author

def replace_lines(file_name, lines):

    done = False

    # ---------------- Martin's code style

    if True:
        raw_top_do_it = False
        if '#include ' in lines[0]: raw_top_do_it = True
        if '#define '  in lines[1]: raw_top_do_it = True
        if '#ifndef '  in lines[1]: raw_top_do_it = True
        if '#ifndef '  in lines[0]: raw_top_do_it = True

        if raw_top_do_it:
            insert_new_header("Martin", file_name, "2009", lines)
            done = True

        for idx,line in enumerate(lines):
            if done: continue
            if idx < 10:
                if line == " */\n":
                    for i in range(idx+1):
                        lines.pop(0)
                    insert_new_header("Martin", file_name, "2009", lines)
                    done = True

        fparts = 0
        for idx,line in enumerate(lines):
            if done: continue
            if idx > 10: continue
            if "//  MoleculesToTriangles"    in line: fparts += 1
            if "//  Created by Martin Noble" in line: fparts += 1
            if "//  Copyright "              in line: fparts += 1
            if "Copyright" in line: copyright_year = get_copyright_year(line)
            if line == "//\n":
                if fparts == 3:
                    for i in range (idx+1):
                        lines.pop(0)
                    insert_new_header("Martin", file_name, copyright_year, lines)
                    done = True

    # ---------------- Kevin's code style

    if False:
        copyright_year = ""
        for idx,line in enumerate(lines):
            if done: continue
            if " (C) "        in line: copyright_year = get_copyright_year(line)
            if " Copyright "  in line: copyright_year = get_copyright_year(line)
            if " Created on:" in line: copyright_year = get_copyright_year(line)
            if " Author: "    in line: author = get_author(line)
            if "York all rights reserved" in line:
                lines.pop(idx)
                lines.pop(idx-1)
                insert_new_header("Kevin Cowtan", file_name, copyright_year, lines)
                done = True

            if idx < 10:
                if " */\n" in line:
                    if len(line) == 4:
                        for i in range(idx+1):
                            lines.pop(0)
                        insert_new_header(author, file_name, copyright_year, lines)
                        done = True

if len(sys.argv) > 1:
    fn = sys.argv[1]
    new_file_name = fn + ".new"
    print("Converting header:", fn)
    f = open(fn)
    lines = f.readlines()
    if len(lines) > 6:
        replace_lines(fn, lines)
        write_lines(lines, new_file_name)
