
import coot_utils
import coot_gui_api
from gi.repository import Gtk
import coot
import gi
gi.require_version('Gtk', '3.0')
import coot_gui


def res_name_from_atom_spec(atom_spec):
    imol = atom_spec[1]
    res_name = coot.residue_name(imol,
                                 atom_spec[2],
                                 atom_spec[3],
                                 atom_spec[4])
    return res_name


def user_defined_add_single_bond_restraint():

    m = "Click on 2 atoms to define the additional bond restraint"
    coot.add_status_bar_text(m)

    def make_restr(*args):
        atom_spec_1 = args[0]
        atom_spec_2 = args[1]
        imol = atom_spec_1[1]
        print("BL DEBUG:: imol: %s spec 1: %s and 2: %s" %
              (imol, atom_spec_1, atom_spec_2))
        coot.add_extra_bond_restraint(imol,
                                      atom_spec_1[2],
                                      atom_spec_1[3],
                                      atom_spec_1[4],
                                      atom_spec_1[5],
                                      atom_spec_1[6],
                                      atom_spec_2[2],
                                      atom_spec_2[3],
                                      atom_spec_2[4],
                                      atom_spec_2[5],
                                      atom_spec_2[6],
                                      1.54, 0.02)  # generic distance?!

    coot.user_defined_click_py(2, make_restr)


def user_defined_add_arbitrary_length_bond_restraint(bond_length=2.0):

    # maybe make a generic one....
    def make_restr(text_list, continue_qm):
        s = "Now click on 2 atoms to define the additional bond restraint"
        coot.add_status_bar_text(s)
        dist = text_list[0]
        try:
            bl = float(dist)
        except:
            bl = False
            coot.add_status_bar_text(
                "Must define a number for the bond length")
        if bl:
            # save distance for future use?!
            def make_restr_dist(*args):
                atom_spec_1 = args[0]
                atom_spec_2 = args[1]
                imol = atom_spec_1[1]
                print("BL DEBUG:: imol: %s spec 1: %s and 2: %s" %
                      (imol, atom_spec_1, atom_spec_2))
                coot.add_extra_bond_restraint(imol,
                                              atom_spec_1[2],
                                              atom_spec_1[3],
                                              atom_spec_1[4],
                                              atom_spec_1[5],
                                              atom_spec_1[6],
                                              atom_spec_2[2],
                                              atom_spec_2[3],
                                              atom_spec_2[4],
                                              atom_spec_2[5],
                                              atom_spec_2[6],
                                              bl, 0.035)
            coot.user_defined_click(2, make_restr_dist)
            if continue_qm:
                user_defined_add_arbitrary_length_bond_restraint(bl)

    def stay_open(*args):
        pass
    # generic_single_entry("Add a User-defined extra distance restraint",
    #                     "2.0",
    #                     "OK...",
    #                     lambda text: make_restr(text))
    coot_gui.generic_multiple_entries_with_check_button(
        [["Add a User-defined extra distance restraint",
          str(bond_length)]],
        ["Stay open?", lambda active_state: stay_open(active_state)],
        "OK...",
        lambda text, stay_open_qm: make_restr(text, stay_open_qm))

# spec_1 and spec_2 are 7-element coot_utils.atom_specs
#


def add_base_restraint(imol, spec_1, spec_2, atom_name_1, atom_name_2, dist):
    """
    spec_1 and spec_2 are 7-element coot_utils.atom_specs
    """

    print("add_base_restraint", imol, spec_1,
          spec_2, atom_name_1, atom_name_2, dist)

    coot.add_extra_bond_restraint(imol,
                                  spec_1[2],
                                  spec_1[3],
                                  spec_1[4],
                                  atom_name_1,
                                  spec_1[6],
                                  spec_2[2],
                                  spec_2[3],
                                  spec_2[4],
                                  atom_name_2,
                                  spec_2[6],
                                  dist, 0.035)


def a_u_restraints(spec_1, spec_2):

    imol = spec_1[1]
    print("BL DEBUG:: add_base_restraint a u", imol,
          spec_1, spec_2, " N6 ", " O4 ", 3.12)
    add_base_restraint(imol, spec_1, spec_2, " N6 ", " O4 ", 3.12)
    add_base_restraint(imol, spec_1, spec_2, " N1 ", " N3 ", 3.05)
    add_base_restraint(imol, spec_1, spec_2, " C2 ", " O2 ", 3.90)
    add_base_restraint(imol, spec_1, spec_2, " N3 ", " O2 ", 5.12)
    add_base_restraint(imol, spec_1, spec_2, " C6 ", " O4 ", 3.92)
    add_base_restraint(imol, spec_1, spec_2, " C4 ", " C6 ", 8.38)


def g_c_restraints(spec_1, spec_2):

    imol = spec_1[1]
    print("BL DEBUG:: add_base_restraint gc", imol,
          spec_1, spec_2, " O6 ", " N4 ", 3.08)
    add_base_restraint(imol, spec_1, spec_2, " N6 ", " O4 ", 3.12)
    add_base_restraint(imol, spec_1, spec_2, " N1 ", " N3 ", 3.04)
    add_base_restraint(imol, spec_1, spec_2, " N2 ", " O2 ", 3.14)
    add_base_restraint(imol, spec_1, spec_2, " C4 ", " N1 ", 7.73)
    add_base_restraint(imol, spec_1, spec_2, " C5 ", " C5 ", 7.21)


def dna_a_t_restraints(spec_1, spec_2):

    imol = spec_1[1]
    add_base_restraint(imol, spec_1, spec_2, " C2 ", " O2 ", 3.49)
    add_base_restraint(imol, spec_1, spec_2, " N1 ", " N3 ", 2.85)
    add_base_restraint(imol, spec_1, spec_2, " N6 ", " O4 ", 3.23)
    add_base_restraint(imol, spec_1, spec_2, " C6 ", " C8 ", 9.94)
    add_base_restraint(imol, spec_1, spec_2, " N1 ", " O2 ", 3.52)


def dna_g_c_restraints(spec_1, spec_2):

    imol = spec_1[1]
    add_base_restraint(imol, spec_1, spec_2, " O6 ", " N4 ", 2.72)
    add_base_restraint(imol, spec_1, spec_2, " N1 ", " N3 ", 2.81)
    add_base_restraint(imol, spec_1, spec_2, " N2 ", " O2 ", 2.83)
    add_base_restraint(imol, spec_1, spec_2, " N9 ", " N1 ", 8.83)


def user_defined_RNA_A_form():
    def make_restr(*args):
        spec_1 = args[0]
        spec_2 = args[1]
        print("BL DEBUG:: have specs", spec_1, spec_2)
        res_name_1 = res_name_from_atom_spec(spec_1)
        res_name_2 = res_name_from_atom_spec(spec_2)
        print("BL DEBUG:: have resnames", res_name_1, res_name_2)
        # just check the first letter, should be save
        if (res_name_1[0] == "G" and
                res_name_2[0] == "C"):
            g_c_restraints(spec_1, spec_2)

        if (res_name_1[0] == "C" and
                res_name_2[0] == "G"):
            g_c_restraints(spec_2, spec_1)

        if (res_name_1[0] == "A" and
                res_name_2[0] == "U"):
            a_u_restraints(spec_1, spec_2)

        if (res_name_1[0] == "U" and
                res_name_2[0] == "A"):
            a_u_restraints(spec_2, spec_1)

    coot.user_defined_click_py(2, make_restr)


def user_defined_DNA_B_form():
    def make_restr(*args):
        spec_1 = args[0]
        spec_2 = args[1]
        res_name_1 = res_name_from_atom_spec(spec_1)
        res_name_2 = res_name_from_atom_spec(spec_2)
        if (res_name_1 == "DG" and
                res_name_2 == "DC"):
            dna_g_c_restraints(spec_1, spec_2)

        if (res_name_1 == "DC" and
                res_name_2 == "DG"):
            dna_g_c_restraints(spec_2, spec_1)

        if (res_name_1 == "DA" and
                res_name_2 == "DT"):
            dna_a_t_restraints(spec_1, spec_2)

        if (res_name_1 == "DT" and
                res_name_2 == "DA"):
            dna_a_t_restraints(spec_2, spec_1)

    coot.user_defined_click_py(2, make_restr)


def user_defined_add_helix_restraints():
    def make_restr(*args):
        spec_1 = args[0]
        spec_2 = args[1]
        chain_id_1 = spec_1[2]
        chain_id_2 = spec_2[2]
        res_no_1 = spec_1[3]
        res_no_2 = spec_2[3]
        imol = spec_1[1]

        if (chain_id_1 == chain_id_2):
            # if backwards, swap 'em
            if res_no_2 < res_no_1:
                tmp = res_no_1
                res_no_1 = res_no_2
                res_no_2 = tmp
            for rn in range(res_no_1, res_no_2 - 2):
                coot.add_extra_bond_restraint(imol,
                                              chain_id_1, rn, "", " O  ", "",
                                              chain_id_1, rn + 3, "", " N  ", "",
                                              3.18, 0.035)
                if (rn + 4 <= res_no_2):
                    coot.add_extra_bond_restraint(imol,
                                                  chain_id_1, rn, "", " O  ", "",
                                                  chain_id_1, rn + 4, "", " N  ", "",
                                                  2.91, 0.035)

    coot.user_defined_click_py(2, make_restr)


def user_defined_delete_restraint():
    def del_restr(*args):
        spec_1 = args[0]
        spec_2 = args[1]
        imol = spec_1[1]
        coot.delete_extra_restraint_py(imol, ["bond", spec_1[2:], spec_2[2:]])

    coot.user_defined_click_py(2, del_restr)


# exte dist first chain A resi 19 ins . atom  N   second chain A resi 19 ins . atom  OG  value 2.70618 sigma 0.4
#
def extra_restraints2refmac_restraints_file(imol, file_name):
    restraints = list_extra_restraints(imol)

    if restraints:
        fin = open(file_name, 'w')
        for restraint in restraints:
            if (restraint[0] == 'bond'):
                chain_id_1 = restraint[1][1]
                resno_1 = restraint[1][2]
                inscode_1 = restraint[1][3]
                atom_1 = restraint[1][4]
                chain_id_2 = restraint[2][1]
                resno_2 = restraint[2][2]
                inscode_2 = restraint[2][3]
                atom_2 = restraint[2][4]
                value = restraint[3]
                esd = restraint[4]
                fin.write("EXTE DIST FIRST CHAIN %s RESI %i INS %s ATOM %s "
                          % (chain_id_1 if (chain_id_1 != "" and chain_id_1 != " ") else ".",
                             resno_1,
                             inscode_1 if (
                                 inscode_1 != "" and inscode_1 != " ") else ".",
                             atom_1))
                fin.write(" SECOND CHAIN %s RESI %i INS %s ATOM %s "
                          % (chain_id_2 if (chain_id_2 != "" and chain_id_2 != " ") else ".",
                             resno_2,
                             inscode_2 if (
                                 inscode_2 != "" and inscode_2 != " ") else ".",
                             atom_2))
                fin.write("VALUE %f SIGMA %f\n" % (value, esd))
        fin.close()


def res_name2plane_atom_name_list(res_name):

    if not (isinstance(res_name, str)):
        return False
    else:
        if (res_name == "DG"):
            return ["N1", "C6", "O6", "C2", "N2", "N3", "C5", "C4", "N9", "N7", "C8", "C1'"]
        elif (res_name == "DA"):
            return ["N1", "C6", "C2", "N6", "N3", "C5", "C4", "N9", "N7", "C8", "C1'"]
        elif (res_name == "DT"):
            return ["N3", "C2", "O2", "C4", "O4", "C5", "C7", "C6", "N1", "C1'"]
        elif (res_name == "DC"):
            return ["N3", "C2", "O2", "C4", "N4", "C5", "C6", "N1", "C1'"]

        elif (res_name == "G"):
            return ["N1", "C6", "O6", "C2", "N2", "N3", "C5", "C4", "N9", "N7", "C8", "C1'"]
        elif (res_name == "A"):
            return ["N1", "C6", "C2", "N6", "N3", "C5", "C4", "N9", "N7", "C8", "C1'"]
        elif (res_name == "T"):
            return ["N3", "C2", "O2", "C4", "O4", "C5", "C7", "C6", "N1", "C1'"]
        elif (res_name == "U"):
            return ["N3", "C2", "O2", "C4", "O4", "C5", "C6", "N1", "C1'"]
        elif (res_name == "C"):
            return ["N3", "C2", "O2", "C4", "N4", "C5", "C6", "N1", "C1'"]
        else:
            return []


# example?
# exte stac plan 1  firs resi 99 chai A atoms { CB CG CD1 CD2 CE1 CE2 CZ OH }  plan 2  firs resi 61 chai B atoms { CB CG CD1 CD2 CE1 CE2 CZ }  dist 3.4 sddi 0.2  sdan 6.0 type 1
#
def write_refmac_parallel_plane_restraint(file_name,
                                          res_spec_0, res_spec_1,
                                          atom_list_0, atom_list_1):

    fout = open(file_name, 'w')
    fout.write("EXTE STACK PLAN 1 FIRST RESIDUE ")
    fout.write(str(coot_utils.residue_spec_to_res_no(res_spec_0)))
    fout.write(" INS . ")  # hack
    fout.write(" CHAIN ")
    fout.write(coot_utils.residue_spec_to_chain_id(res_spec_0))
    fout.write(" ATOMS { ")
    for atom_name in atom_list_0:
        fout.write(" " + atom_name + " ")
    fout.write(" } PLAN 2 FIRST RESIDUE ")
    fout.write(str(coot_utils.residue_spec_to_res_no(res_spec_1)))
    fout.write(" INS . ")  # hack
    fout.write(" CHAIN ")
    fout.write(coot_utils.residue_spec_to_chain_id(res_spec_1))
    fout.write(" ATOMS { ")
    for atom_name in atom_list_1:
        fout.write(" " + atom_name + " ")
    fout.write(" } DIST 3.4 SDDI 0.2 SDAN 6.0 TYPE 1 \n")
    fout.close()
    # maybe some return value that it worked to write at some point


def add_parallel_planes_restraint(imol, rs_0, rs_1):

    print("in add_parallel_planes_restraint: rs_0: %s rs_1 %s" % (rs_0, rs_1))

    rn_0 = coot_utils.residue_spec_to_residue_name(imol, rs_0)
    rn_1 = coot_utils.residue_spec_to_residue_name(imol, rs_1)
    atom_ls_0 = res_name2plane_atom_name_list(rn_0)
    atom_ls_1 = res_name2plane_atom_name_list(rn_1)

    write_refmac_parallel_plane_restraint(
        "tmp.rst", rs_0, rs_1, atom_ls_0, atom_ls_1)

    coot.add_refmac_extra_restraints(imol, "tmp.rst")


def user_defined_add_planes_restraint():

    coot.add_status_bar_text(
        "Click on 2 atoms to define the additional parallel planes restraint")

    def make_restr(*args):
        atom_0 = args[0]
        atom_1 = args[1]
        rs_0 = coot_utils.atom_spec_to_residue_spec(atom_0)
        rs_1 = coot_utils.atom_spec_to_residue_spec(atom_1)
        imol = coot_utils.atom_spec_to_imol(atom_0)

        rn_0 = coot.residue_name(imol,
                                 res_spec_utils.residue_spec_to_chain_id(
                                     coot_utils.atom_spec_to_residue_spec(atom_0)),
                                 res_spec_utils.residue_spec_to_res_no(
                                     coot_utils.atom_spec_to_residue_spec(atom_0)),
                                 coot_utils.residue_spec_to_ins_code(coot_utils.atom_spec_to_residue_spec(atom_0)))
        rn_1 = coot.residue_name(imol,
                                 res_spec_utils.residue_spec_to_chain_id(
                                     coot_utils.atom_spec_to_residue_spec(atom_1)),
                                 res_spec_utils.residue_spec_to_res_no(
                                     coot_utils.atom_spec_to_residue_spec(atom_1)),
                                 coot_utils.residue_spec_to_ins_code(coot_utils.atom_spec_to_residue_spec(atom_1)))

        print("BL DEBUG:: got resname 0", rn_0)
        print("BL DEBUG:: got resname 1", rn_1)

        atom_ls_0 = res_name2plane_atom_name_list(rn_0)
        atom_ls_1 = res_name2plane_atom_name_list(rn_1)

        write_refmac_parallel_plane_restraint("tmp.rst",
                                              rs_0, rs_1,
                                              atom_ls_0, atom_ls_1)
        coot.add_refmac_extra_restraints(imol, "tmp.rst")

    coot.user_defined_click_py(2, make_restr)


if True:
    if coot_gui_api.main_menubar():

        menu = coot_gui.coot_menubar_menu("Restraints")

        coot_gui.add_simple_coot_menu_menuitem(
            menu,
            "Add Simple C-C Single Bond Restraint...",
            lambda func: user_defined_add_single_bond_restraint())

        coot_gui.add_simple_coot_menu_menuitem(
            menu,
            "Add Distance Restraint...",
            lambda func: user_defined_add_arbitrary_length_bond_restraint())

        coot_gui.add_simple_coot_menu_menuitem(
            menu,
            "Add Helix Restraints...",
            lambda func: user_defined_add_helix_restraints())

        coot_gui.add_simple_coot_menu_menuitem(
            menu,
            "RNA A form bond restraints...",
            lambda func: user_defined_RNA_A_form())

        coot_gui.add_simple_coot_menu_menuitem(
            menu,
            "DNA B form bond restraints...",
            lambda func: user_defined_DNA_B_form())

        def combobox_to_molecule_number(combobox):
            imol = -1
            tree_iter = combobox.get_active_iter()
            if tree_iter is not None:
                model = combobox.get_model()
                it = model[tree_iter]
                imol = it[0]
            return imol

        coot_gui.add_simple_coot_menu_menuitem(
            menu,
            "Read Refmac Extra Restraints...",
            lambda func:
            coot_gui.generic_chooser_and_file_selector("Apply restraints to molecule",
                                                       coot_utils.valid_model_molecule_qm,
                                                       "File:", "",
                                                       lambda imol, file_name:
                                                       coot.add_refmac_extra_restraints(imol, file_name)))

        coot_gui.add_simple_coot_menu_menuitem(
            menu,
            "Delete Extra Restraints...",
            lambda func:
            coot_gui.molecule_chooser_gui("Delete Extra Restraints for Molecule:",
                                          lambda imol:
                                          coot.delete_all_extra_restraints(imol)))

        def set_prosmart_sigma_limit_func(low, high):
            with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                                  aa_ins_code, aa_atom_name, aa_alt_conf]:
                coot.set_extra_restraints_prosmart_sigma_limits(
                    aa_imol, low, high)

        # coot_gui.add_simple_coot_menu_menuitem(
        #     menu,
        #     "ProSMART restraints interesting limit to 0.5...",
        #     lambda func: set_prosmart_sigma_limit_func(-0.5, 0.5)
        # )

        # coot_gui.add_simple_coot_menu_menuitem(
        #     menu,
        #     "ProSMART restraints interesting limit to 2.5...",
        #     lambda func: set_prosmart_sigma_limit_func(-2.5, 2.5)
        # )

        # def set_extra_restraints_display_func(state):
        #     with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
        #                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
        #         coot.set_show_extra_restraints(aa_imol, state)

        # coot_gui.add_simple_coot_menu_menuitem(
        #     menu,
        #     "Undisplay Extra Restraints",
        #     lambda func: set_extra_restraints_display_func(0))

        # coot_gui.add_simple_coot_menu_menuitem(
        #     menu,
        #     "Display Extra Restraints",
        #     lambda func: set_extra_restraints_display_func(1))

        # def set_prosmart_display_CA_func(state):
        #     with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
        #                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
        #         coot.set_extra_restraints_representation_for_bonds_go_to_CA(
        #             aa_imol, state)

        # ProSMART menu item - but who needs it?
        # coot_gui.add_simple_coot_menu_menuitem(
        #     menu,
        #     "Extra Restraints to CA",
        #     lambda func: set_prosmart_display_CA_func(1))

        coot_gui.add_simple_coot_menu_menuitem(
            menu,
            "Extra Restraints Standard Representation",
            lambda func: set_prosmart_display_CA_func(0))

        coot_gui.add_simple_coot_menu_menuitem(
            menu,
            "Delete an Extra Restraint...",
            lambda func: user_defined_delete_restraint())

        def delete_restraints_func():
            with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code, aa_atom_name, aa_alt_conf]:
                coot.delete_extra_restraints_for_residue(
                    aa_imol, aa_chain_id, aa_res_no, aa_ins_code)

        coot_gui.add_simple_coot_menu_menuitem(
            menu,
            "Delete Restraints for this residue",
            lambda func: delete_restraints_func()
        )

        def del_deviant_restr_func(text):
            try:
                n = float(text)
                with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code, aa_atom_name, aa_alt_conf]:
                    coot.delete_extra_restraints_worse_than(aa_imol, n)
            except:
                print("BL WARNING:: no float given")

        coot_gui.add_simple_coot_menu_menuitem(
            menu,
            "Delete Deviant Extra Restraints...",
            lambda func: coot_gui.generic_single_entry("Delete Restraints worse than ",
                                                       "4.0", " Delete Outlying Restraints ",
                                                       lambda text: del_deviant_restr_func(text))
        )

        coot_gui.add_simple_coot_menu_menuitem(
            menu,
            "Save as REFMAC restraints...",
            lambda func:
            coot_gui.generic_chooser_and_file_selector("Save REFMAC restraints for molecule",
                                                       coot_utils.valid_model_molecule_qm,
                                                       " Restraints file name:  ",
                                                       "refmac-restraints.txt",
                                                       lambda imol, file_name:
                                                       extra_restraints2refmac_restraints_file(imol, file_name)))

        coot_gui.add_simple_coot_menu_menuitem(
            menu,
            "Add Parallel Planes Restraint...",
            lambda func:
            user_defined_add_planes_restraint()
        )
