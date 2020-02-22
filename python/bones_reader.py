#    bones_reader.py
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

mapman_exe = os.path.join(os.getenv('HOME'), "bin/lx_mapman")

def generic_object_from_bones(bones_file):

    import operator
    from types import IntType
    
    global atom_xyz_list, bone_list, conn_list
    atom_xyz_list = []
    bone_list = []
    conn_list = []
    #save_conn = False

    def add_ordinate(s):
        global atom_xyz_list
        obj = map(float, s.split())
        atom_xyz_list.append(obj)

    def add_bone(s):
        global bone_list
        obj = map(int, s.split())
        bone_list += obj

    #def add_connection_pairs(s):
    #    conn_list += s.split()

    def add_connection(s):
        global conn_list
        obj = map(int, s.split())
        conn_list += obj     

    def my_indexer(n):
        # un-fortran indexing
        if (n % 2 == 0):
            # even
            ret = n // 2
        else:
            ret = (n - 1) // 2
        return ret - 1

    # main body
    if (os.path.isfile(bones_file)):
        fin = open(bones_file, 'r')
        lines = fin.readlines()
        fin.close()
        reading_atoms_flag = False
        reading_atom_bone_flag = False
        reading_connectivity_flag = False
        reading_atom_colour_flag = False
        reading_residue_name_flag = False
        reading_residue_type_flag = False
        reading_residue_pointers_flag = False
        
        for line in lines:
            
            # atom colour
            if ("SKL_ATOM_COLOUR" in line):
                print "found atom_colour."
                reading_connectivity_flag = False
                #saved_conn = False

            # connectivity
            if (reading_connectivity_flag):
                add_connection(line)
            if ("SKL_CONNECTIVITY" in line):
                print "found atom_connectivity."
                reading_connectivity_flag = True
                reading_atom_bone_flag = False

            # bones
            if (reading_atom_bone_flag):
                add_bone(line)
            if ("SKL_ATOM_BONE" in line):
                print "found atom_bone."
                reading_atoms_flag = False
                reading_atom_bone_flag = True

            # atoms
            if (reading_atoms_flag):
                add_ordinate(line)
            if ("SKL_ATOM_XYZ" in line):
                print "found atom_xyz"
                reading_atoms_flag = True

        print "length atom_xyz_list: %s length bone_list: %s length conn_list: %s" \
              %(len(atom_xyz_list), len(bone_list), len(conn_list))

        lines_obj = new_generic_object_number("Lines")

        current_position = False
        for connections in conn_list:
            if (not current_position):
                #print "(from %s) init move to : %s" %(connections,
                #                                      my_indexer(connections))
                current_position = my_indexer(connections)
            else:
                xyz_index = my_indexer(connections)
                if (connections % 2 == 0):
                    # even connections
                    # a move to:
                    #print "(from %s) move to: %s" %(connections, xyz_index)
                    current_position = xyz_index
                else:
                    # a draw to:
                    #print "(from %s) line: %s %s" %(connections,
                    #                                current_position,
                    #                                xyz_index)
                    to_generic_object_add_line(lines_obj, "green", 3,
                                               atom_xyz_list[current_position][0],
                                               atom_xyz_list[current_position][1],
                                               atom_xyz_list[current_position][2],
                                               atom_xyz_list[xyz_index][0],
                                               atom_xyz_list[xyz_index][1],
                                               atom_xyz_list[xyz_index][2])
                    current_position = my_indexer(connections)

        set_display_generic_object(lines_obj, 1)
            
    else:
        print "BL INFO:: no bones file %s found" %bones_file


def bones_it(map_file_name):
    bones_file = "my.bones"
    data_lines = ["read m1 " + map_file_name + " ccp4",
                  "bo sk m1 0.5 0.15 1",
                  "bones connect",
                  bones_file,
                  "skl",
                  "5",
                  "quit"]

    popen_command(mapman_exe, [], data_lines, "coot-mapman.log", True)
    # now read in bones_file
    generic_object_from_bones(bones_file)

# add a menu item:
menu = coot_menubar_menu("Mapman")

def bonesing_func(imol):
    print "bonesing", imol
    export_map(imol, "tmp.map")
    bones_it("tmp.map")
    # remove tmp map?
    if (os.path.isfile("tmp.map")):
        print "BL INFO:: removing temporary map file tmp.map"
        os.remove("tmp.map")

add_simple_coot_menu_menuitem(menu, "Mapman Bones...",
                              lambda func: map_molecule_chooser_gui("Map to Bonesify:",
                                                       lambda imol: bonesing_func(imol)))

#generic_object_from_bones("my.bones")
