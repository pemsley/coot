
import coot
import coot_gui
import coot_utils
import os

global pisa_command
pisa_command = "pisa"
global pisa_min_version
pisa_min_version = "v1.13"

def pisa_molecule_chooser_gui(mode):

    if mode == "interfaces":
        coot_gui.molecule_chooser_gui("Choose molecule for PISA assembly analysis",
                                      lambda imol: pisa_interfaces(imol))

    if mode == "assemblies":
        coot_gui.molecule_chooser_gui("Choose molecule for PISA assembly analysis",
                                      lambda imol: pisa_assemblies(imol))


def pisa_assemblies(imol):

    global pisa_min_version
    pisa_exe = pisa_new_enough_qm()
    if not pisa_exe:
        msg = "Your pisa version it too old.  Need at least " + pisa_min_version + "."
        coot.info_dialog(msg)
    else:
        #
        # main line
        pdb_file_name, pisa_config, pisa_xml_file_name = prep_for_pisa('assemblies', imol)

        status_1 = coot_utils.popen_command(pisa_exe,
                                 ["pisa", "-analyse", pdb_file_name],
                                 [], "pisa.log", False)

        if (status_1 != 0):
            coot.info_dialog("Ooops PISA failed to deliver the goods!\n\n(Go for a curry instead?)")
        else:
            # good
            # used to be print "BL DEBUG:: 2nd pisa args", [pisa_project_name, "-xml", pisa_config]
            status_2 = coot_utils.popen_command(pisa_exe,
                                     ["pisa", "-xml", "assemblies"],
                                     [], pisa_xml_file_name, False)
            if (status_2 == 0):
                pisa_assemblies_xml(imol, pisa_xml_file_name)


# called by pisa_assemblies, which is the entry point from the main
# program gui.
#
#  it calls parse_pisa_assemblies which does the work

def pisa_assemblies_xml(imol, file_name):

    from xml.etree.ElementTree import ElementTree
    import os

    if (os.path.isfile(file_name)):
        print("opened", file_name)
        sm = ElementTree()
        xml_file = clean_xml_file(file_name)
        if xml_file:
            sm.parse(xml_file)
            if (xml_file != file_name):
                os.remove(xml_file)
            parse_pisa_assemblies(imol, sm)
            # rm xml_file ?!

# return a filename which is clean (pisa) xml or False if there was problems
# not sure if needed any more if v1.13
def clean_xml_file(filename):

    if (os.path.isfile(filename)):
        fin = open(filename, 'r')
        lines = fin.readlines()
        fin.close()
        if ("<pisa_" == lines[0][0:6] or
            "<pdb_" == lines[0][0:5]):
            # have a clean file already, return filename
            return filename
        else:
            # not a clean xml file, so clean it up
            new_lines = []
            start_write = 0
            for line in lines:
                if ("<pisa_" == line[0:6] or
                    "<pdb_" == line[0:5]):
                    start_write += 1
                if (start_write == 1):
                    new_lines.append(line)
                if (start_write == 2):
                    new_lines.append(line)
                    break
            if new_lines:
                tmp_file = "xml-tmp.xml"
                fout = open(tmp_file, 'w')
                fout.writelines(new_lines)
                fout.close()
                return tmp_file
            else:
                return False

    else:
        return False

    
# Exported to the main level.  A molecule is common to assemblies and interfaces.
# 
# If pisa-result-type is 'assemblies', return a molecule number or False
# 
# If pisa-result-type is 'interfaces', return a interface-molecule record or False
#
# If pisa-result-type is neither of the above, return False
#
# a interface-molecule record contains information about pvalue and residues.
#
def pisa_handle_xml_molecule(imol, molecule, pisa_results_type):

    # was it a chain? (or residue selection disguised as
    # a chain?)
    #
    def pisa_make_atom_selection_string(chain_id_raw):

        #  first try to split the chain-id-raw on a "]".  If there was no
        # "]" then we had a simple chain-id.  If there was a "]" then we
        # have something like "[CL]D:32", or "[ZN]-:2" from which we need
        # to extract a residue number and a chain id.  Split on ":" and
        # construct the left and right hand sides of s.  Then use those
        # together to make an mmdb selection string. If chaind-id is "-"
        # then reset it to blank. 
        #
        if ("]" in chain_id_raw):
            # print "found a residue selection chain_id", chain_id_raw
            s = chain_id_raw[chain_id_raw.find("]") + 1:] # e.g. "D:32"
            if (":" in s):
                residue_string = s[s.find(":") + 1:]
                chain_string   = s[0:s.find(":")]
                element_string = chain_id_raw[1:chain_id_raw.find("]")]
                if (chain_string == "-"):
                    atom_selection_string = "// /" + residue_string + "/" + element_string
                else:
                    atom_selection_string = "//" + chain_string + "/" + residue_string
                # print "BL coot_utils.debug:: atom_selection_string: %s from %s" \
                #      %(atom_selection_string, chain_id_raw)
                return atom_selection_string
        else:
            #print "found a simple chain_id", chain_id_raw
            return "//" + chain_id_raw

    # return a list of residue dictionaries
    # filter out the ones with bsa < 0.1
    # residues_xml are parsed as a list of xml elements
    #
    def handle_residues(residues_xml):
        residues = []
        for residue in residues_xml:
            res_dic = dict((ele.tag, ele.text) for ele in residue)
            if not res_dic["ins_code"]:
                res_dic["ins_code"] = ""
            if (float(res_dic["bsa"]) > 0.1):
                residues.append(res_dic)
        return residues


    ( STRING,
      STDOUT) = list(range(2))
    # record in an element in xml element tree (molecule)
    # port is STRING or STDOUT
    #
    def print_molecule(record_xml, port):
        # first make a dictionary
        rec_dic = dict((ele.tag, ele.text) for ele in record_xml)
        # rec_dic = record
        # now either print to stdout or write to string (which will
        # be returned - otherwise return None, False on error)
        ret = None
        if (port == STRING):
            import io
            out = io.StringIO()
        elif (port == STDOUT):
            import sys
            out = sys.stdout
        else:
            return False
        
        out.write("[molecule: id %s\n"       %rec_dic['id'])
        out.write(" molecule: chain-id %s\n" %rec_dic['chain-id'])
        out.write(" molecule: class %s\n"    %rec_dic['class'])
        out.write(" molecule: symop %s\n"    %rec_dic['symop'])
        out.write(" molecule: symop-no %s\n" %rec_dic['symop-no'])
        out.write(" molecule: natoms %s\n"   %rec_dic['natoms'])
        out.write(" molecule: nres %s\n"     %rec_dic['nres'])
        out.write(" molecule: area %s\n"     %rec_dic['area'])
        out.write(" molecule: solv-en %s\n"  %rec_dic['solv-en'])
        out.write(" molecule: pvalue %s\n"   %rec_dic['pvalue'])
        out.write(" molecule: with %s residues]\n" %len(record_xml.getiterator("residues")))
        if (port == STRING):
            ret = out.getvalue()
            out.close()
        return ret

    # record in an element in xml element tree (residue)
    def print_residue(record_xml, port):
        # first make a dictionary
        rec_dic = dict((ele.tag, ele.text) for ele in record_xml)
        # rec_dic = record
        # now either print to stdout or write to string (which will
        # be returned - otherwise return None, False on error)
        ret = None
        if (port == STRING):
            import io
            out = io.StringIO()
        elif (port == STDOUT):
            import sys
            out = sys.stdout
        else:
            return False
        
        out.write("[residue: ser-no %s,\n" %rec_dic['ser-no'])
        out.write(" name %s,\n"            %rec_dic['nam'])
        out.write(" seq-num %s,\n"         %rec_dic['seq-num'])
        out.write(" ins-code %s,\n"        %rec_dic['ins-code'])
        out.write(" asa %s,\n"             %rec_dic['asa'])
        out.write(" bsa %s,\n"             %rec_dic['bsa'])
        out.write(" solv-en %s]\n"         %rec_dic['solv-en'])
        if (port == STRING):
            ret = out.getvalue()
            out.close()
        return ret

    #
    def process_molecule_internal(molecule_xml):
        
        # molecule is an element in the xml tree
        # lets make a dictionary out of it for easier handling
        mol_dic = dict((ele.tag, ele.text) for ele in molecule)
        #
        # rtop-symbols are common to interfaces and assemblies
        rtop_symbols = ['rxx', 'rxy', 'rxz', 'ryx', 'ryy', 'ryz', 'rzx', 'rzy', 'rzz',
                        'tx', 'ty', 'tz']
        # matrix elements
        # these symbols only interfaces have
        extra_symbols = ['pvalue', 'residues']

        ass_rtop_symbols = []        # association
        atom_selection_string = "//" # default, everything
        symm_name_part = ""          # the atom selection info without
                                     # "/" because that would chop the
                                     # name in the display manager.


        chain_id_raw = mol_dic["chain_id"]
        atom_selection_string = pisa_make_atom_selection_string(chain_id_raw)
        if (chain_id_raw == atom_selection_string):
            symm_name_part = "chain " + chain_id_raw
        else:
            symm_name_part = chain_id_raw

        # was it on of the rotation or coot_utils.translation symbols?
        for symbol in rtop_symbols:
            ass_rtop_symbols.append([symbol, float(mol_dic[symbol])])

        # do something with residues (if interface)
        if "residues" in mol_dic:
            residues = molecule.getiterator("residue")
            mol_dic["residues"] = handle_residues(residues)
            
        if not (len(ass_rtop_symbols) == 12):
            return False # ass_rtop_symbols were not all set
        else:
            mat = [sym[1] for sym in ass_rtop_symbols]
            #print "== atom-selection string %s   mat:::: %s" %(atom_selection_string, mat)
            #print "currently %s molecules" %(coot.graphics_n_molecules())
            new_molecule_name = "Symmetry copy of " + \
                                str(imol) +\
                                " using " + symm_name_part
            # new-molecule-by-symop-with-atom-selection,
            # perhaps? (20100222 doesn't seem needed because
            # the transformation is contained in mat.)
            new_mol_no = coot.new_molecule_by_symmetry_with_atom_selection(
                imol,
                new_molecule_name,
                atom_selection_string,
                *(mat + [0,0,0]))

            # the return value depends on pisa-results-type
            if pisa_results_type == 'assemblies':
                # assemblies:
                return new_mol_no
            elif pisa_results_type == 'interfaces':
                # interfaces:
                # in python list of 3 (2 would be enough?!)
                return [new_mol_no, mol_dic["symop"],
                        mol_dic]
            else:
                return False
                
    # main line
    #
    if (pisa_results_type == 'assemblies' or
        pisa_results_type == 'interfaces'):
        return process_molecule_internal(molecule)
    else:
        return False


# ----------------------------------------------------
#                      pisa assemblies:
# ----------------------------------------------------
#
# pisa results
#    name
#    status
#    total_asm
#    asm_set
#       ser_no
#       assembly
#          id
#          size
#          mmsize
#          diss_energy
#          asa
#          bas
#          entropy
#          diss_area
#          int_energy
#          n_uc
#          n_diss
#          symNumber
#          molecule
#             chain_id
#             rxx
#             rxy
#             rxz
#             tx
#             ryx
#             ryy
#             ryz
#             ty
#             rzx
#             rzy
#             rzz
#             tz
#             rxx-f
#             rxy-f
#             rxz-f
#             tx-f
#             ryx-f
#             ryy-f
#             ryz-f
#             ty-f
#             rzx-f
#             rzy-f
#             rzz-f
#             tz-f

def parse_pisa_assemblies(imol, entity):

    # Return the model number of the new assembly molecule
    #
    def create_assembly_set_molecule(assembly_set_molecule_numbers,
                                     assembly_set_number = ""):
        if not assembly_set_molecule_numbers:
            return False
        else:
            first_copy = coot.copy_molecule(assembly_set_molecule_numbers[0])
            if not coot_utils.valid_model_molecule_qm(first_copy):
                return False
            else:
                rest = assembly_set_molecule_numbers[1:]
                merge_molecules(rest, first_copy)
                coot.set_molecule_name(first_copy, "Assembly Set " + assembly_set_number)
                return first_copy

    ( STRING,
      STDOUT) = list(range(2))
    # record in an element in xml element tree
    # port is STRING or STDOUT
    #
    def print_assembly(record, port):
        # first make a dictionary
        #rec_dic = dict((ele.tag, ele.text) for ele in record)
        rec_dic = record
        # now either print to stdout or write to string (which will
        # be returned - otherwise return None, False on error)
        ret = None
        if (port == STRING):
            import io
            out = io.StringIO()
        elif (port == STDOUT):
            import sys
            out = sys.stdout
        else:
            return False
        
        out.write("assembly id: %s\n" %rec_dic['id'])
        out.write("assembly size: %s\n" %rec_dic['size'])
        out.write("assembly symm-number: %s\n" %rec_dic['symNumber'])
        out.write("assembly asa: %s\n" %rec_dic['asa'])
        out.write("assembly bsa: %s\n" %rec_dic['bsa'])
        out.write("assembly diss_energy: %s\n" %rec_dic['diss_energy'])
        out.write("assembly entropy: %s\n" %rec_dic['entropy'])
        out.write("assembly diss_area: %s\n" %rec_dic['diss_area'])
        out.write("assembly int_energy: %s\n" %rec_dic['int_energy'])
        # maybe split next line when more than 2 or so molecules?!
        out.write("assembly components: %s\n" %rec_dic['molecule'])
        if (port == STRING):
            ret = out.getvalue()
            out.close()
        return ret

    # Return a list of model numbers
    #
    def create_assembly_molecule(assembly_molecule_numbers):
        return assembly_molecule_numbers

    # return an assembly record:
    # this is now a dictionary with all assembly properties
    # molecule is replaced with a list of symmetry created molecules
    # arg assembly is assembly element from xml tree
    #
    def handle_assembly(assembly):
        assembly_molecule_numbers = []
        assembly_dic = dict((ele.tag, ele.text) for ele in assembly)
        molecules = assembly.getiterator("molecule")
        for molecule in molecules:
            mol_no = pisa_handle_xml_molecule(imol, molecule, 'assemblies')
            assembly_molecule_numbers.append(mol_no)
        assembly_dic["molecule"] = assembly_molecule_numbers
        for mol_no in assembly_molecule_numbers:
            if (coot_utils.valid_model_molecule_qm(mol_no)):
                coot.set_mol_displayed(mol_no, 0)
        return assembly_dic


    # handle assembly-set (the xml tree element for asm_set.)
    #
    # return values [False or the molecule number of the
    # assembly-set] and the assembly-record-set (in dictionary
    # format - which includes the assembly set molecules list!!
    #
    def handle_assembly_set(assembly_set):
        assembly_records = assembly_set.find("assembly")
        assembly_record = handle_assembly(assembly_records)
        assembly_set_molecules = assembly_record["molecule"]
        new_mol = create_assembly_set_molecule(assembly_set_molecules,
                                               assembly_record["id"])

        return new_mol, assembly_record

    # main line
    #

    #
    first_assembly_set_is_displayed_already_qm = False
    top_assembly_set = False

    assemblies = entity.getiterator("asm_set")
    for ass in assemblies:
        molecule_number, assembly_record_set = handle_assembly_set(ass)
        if not top_assembly_set:
            # make a string out of the dictionary
            top_assembly_set = print_assembly(assembly_record_set, STRING)
        if first_assembly_set_is_displayed_already_qm:
            coot.set_mol_displayed(molecule_number, 0)
        else:
            first_assembly_set_is_displayed_already_qm = True


    if top_assembly_set:
        print("top assembly-set:\n", top_assembly_set)
        s = "top assembly-set: \n\n" + top_assembly_set
    else:
        s = "no assembly-sets found"
    coot.info_dialog(s)
            

def make_pisa_config(pisa_coot_dir, config_file_name):
    s = os.getenv("CCP4")
    if (os.path.isdir(s)):
        ls = [["DATA_ROOT", os.path.join(s, "share", "pisa")],
              ["SRS_DIR", os.path.join(s, "share", "ccp4srs")],
              ["PISTORE_DIR", os.path.join(s, "share", "pisa")],
              ["PISA_WORK_ROOT", pisa_coot_dir],
              # according to Paule (and I agree):
              # that we need these next 3 is ridiculous
              # BL:: hope not really needed (as for Win exe is missing)
              # maybe already fixed!?
              ["MOLREF_DIR", os.path.join(s, "share", "pisa")],
              ["RASMOL_COM", os.path.join(s, "bin", "rasmol")],
              ["CCP4MG_COM", os.path.join(s, "bin", "ccp4mg")],
              # the parent dir of this needs to be writable, usually isnt
              # since it's all in a ccp4 install dir (which is DATA_ROOT?!?!)
              ["SESSION_PREFIX", "pisa_"],
              ]
        fin = open(config_file_name, 'w')
        for item_pair in ls:
            fin.write(item_pair[0] + "\n")
            fin.write(item_pair[1] + "\n")
        fin.close()


# 20100213 prep-for-pisa needs to make directory, config file, write
# the pdb file and the return value should be #f if there was a
# problem or some value where we can check that pisa -analyse ran
# (probably a directory).  It is not clear right now where the output
# is going.  config files has PISA_WORK_ROOT coot-pisa but things
# seems to be written to DATA_ROOT
# /home/emsley/ccp4/ccp4-6.1.3/share/pisa which seems like a bug (or
# something like it) in pisa.  Needs rechecking
#
# maybe santisation of xml fiel can go here too?!?!
def prep_for_pisa(mode, imol):

    #
    def make_stubbed_name(imol):
        return coot_utils.strip_extension(os.path.basename(coot.molecule_name(imol)))

    if coot_utils.valid_model_molecule_qm(imol):
        pisa_coot_dir = "coot-pisa"
        stubbed_name = make_stubbed_name(imol)
        pdb_file_name      = os.path.join(pisa_coot_dir, stubbed_name + ".pdb")
        pisa_config        = os.path.join(pisa_coot_dir, stubbed_name + "-pisa.cfg")
        pisa_xml_file_name = os.path.join(pisa_coot_dir, stubbed_name + "-" + str(mode) + ".xml")
        #pisa_project_name  = os.path.join(stubbed_name)

        coot.make_directory_maybe(pisa_coot_dir)
        make_pisa_config(pisa_coot_dir, pisa_config)
        coot.write_pdb_file(imol, pdb_file_name)
        if (os.path.isfile(pdb_file_name)):
            return pdb_file_name, pisa_config, pisa_xml_file_name
        else:
            return False, False, False

# needs fleshing out (see notes for prep-for-pisa).
#
def cached_pisa_analysis(dir):
    return False
    
#
# return pisa_command_exe or False
#
def pisa_new_enough_qm():
    global pisa_command
    global pisa_min_version
    if not pisa_command:
        pisa_command="pisa"
    pisa_exe = coot_utils.find_exe(pisa_command, "CBIN", "CCP4_BIN", "PATH")
    if pisa_exe:
        tmp_file = "pisa-version.log"
        process = coot_utils.popen_command(pisa_exe, [], [], tmp_file)
        fin = open(tmp_file, 'r')
        lines = fin.readlines()
        fin.close()
        os.remove(tmp_file)
        for line in lines:
            if " v" in line:
                if line.split()[0] >= pisa_min_version:
                    return pisa_exe
                else:
                    return False
    else:
        return False
    

# -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------
#                                   interfaces
# -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------

def pisa_interfaces(imol):
    global pisa_min_version
    pisa_exe = pisa_new_enough_qm()
    if not pisa_exe:
        msg = "Your pisa version it too old.  Need at least" \
              + pisa_min_version + "."
        coot.info_dialog(msg)
    else:
        pdb_file_name, pisa_config, pisa_xml_file_name = prep_for_pisa("interfaces", imol)

        if pisa_config:
            if not cached_pisa_analysis(pisa_config):
                # pisa analysis
                status_1 = coot_utils.popen_command(pisa_exe,
                                         ["pisa", "-analyse", pdb_file_name],
                                         [],
                                         "pisa-analysis.log", False)
                if (status_1 == 0):
                    status_2 = coot_utils.popen_command(pisa_exe,
                                             ["pisa", "-xml", "interfaces"],
                                             [],
                                             pisa_xml_file_name, False)
                    if (status_2 == 0):  # good
                        pisa_interfaces_xml(imol, pisa_xml_file_name)


def pisa_interfaces_xml(imol, file_name):
    if not os.path.isfile(file_name):
        print("WARNING:: in pisa_interfaces_xml: %s does not exist" %file_name)
    else:
        from xml.etree.ElementTree import ElementTree
        sm = ElementTree()
        xml_file = clean_xml_file(file_name)
        if xml_file:
            sm.parse(xml_file)
            if (xml_file != file_name):
                os.remove(xml_file)
            parse_pisa_interfaces(imol, sm)

# pdb_entry
#    pdb_code
#    status
#    n_interfaces
#    interface
#       id
#       type
#       n_occ
#       int_area
#       int_solv_en
#       pvalue
#       stab_en
#       css
#       overlap
#       x-rel
#       fixed
#       h-bonds
#          n_bonds
#          bond
#             chain-1
#             res-1
#             seqnum-1
#             inscode-1
#             atname-1
#             chain-2
#             res-2
#             seqnum-2
#             inscode-2
#             atname-2
#             dist
#       salt-bridges
#          n-bonds
#          bond
#             chain-1
#             res-1
#             seqnum-1
#             inscode-1
#             atname-1
#             chain-2
#             res-2
#             seqnum-2
#             inscode-2
#             atname-2
#             dist
#       ss-bonds
#          n-bonds
#          bond
#             chain-1
#             res-1
#             seqnum-1
#             inscode-1
#             atname-1
#             chain-2
#             res-2
#             seqnum-2
#             inscode-2
#             atname-2
#             dist
#       cov-bonds
#          n-bonds
#          bond
#             chain-1
#             res-1
#             seqnum-1
#             inscode-1
#             atname-1
#             chain-2
#             res-2
#             seqnum-2
#             inscode-2
#             atname-2
#             dist
#       molecule
#          id
#          chain_id
#          class
#          symop
#          symop_no
#          cell_i
#          cell_j
#          cell_k
#          rxx
#          rxy
#          rxz
#          tx
#          ryx
#          ryy
#          ryz
#          ty
#          rzx
#          rzy
#          rzz
#          tz
#          int_natoms
#          int_nres
#          int_area
#          int_solv_en
#          pvalue
#          residues
#             residue
#                ser_no
#                name
#                seq_num
#                ins_code
#                bonds
#                asa
#                bsa
#                solv_en
#                
def parse_pisa_interfaces(imol, xml_entity):
    
    # return a list of bonds:
    #
    # [bond-type, atom-spec-1, atom-spec-2]
    # 
    # where bond-type is 'hbond', 'salt-bridge', 'ss-bond', or 'cov-bond'.
    #
    def molecule_bonds(entity):

        # return a list of 2 atom specs given something like: 
        #
        # xml entity bond
        # return False if the atom specs are not fully set
        #
        def parse_bond(bond_description):

            bond_dic = dict((ele.tag, ele.text) for ele in bond_description)
            chain_id_1  = bond_dic["chain-1"]
            chain_id_2  = bond_dic["chain-2"]
            res_no_1    = bond_dic["seqnum-1"]
            res_no_2    = bond_dic["seqnum-2"]
            ins_code_1  = bond_dic["inscode-1"]
            ins_code_2  = bond_dic["inscode-2"]
            atom_name_1 = bond_dic["atname-1"]
            atom_name_2 = bond_dic["atname-2"]

            if not (chain_id_1 and chain_id_2 and
                    res_no_1 and res_no_2 and
                    atom_name_1 and atom_name_2):
                return False
            else:
                if not ins_code_1:
                    ins_code_1 = ""
                if not ins_code_2:
                    ins_code_2 = ""
                return [[chain_id_1, int(res_no_1), ins_code_1, atom_name_1, ""],
                        [chain_id_2, int(res_no_2), ins_code_2, atom_name_2, ""]]
            
                    
        ret_bonds = []
        bond_type = entity.tag
        bonds = entity.getiterator("bond")
        for bond in bonds:
            coot_utils.atom_specs = parse_bond(bond)
            if coot_utils.atom_specs:
                coot_utils.atom_specs.insert(0, bond_type)
                ret_bonds.append(atom_specs)
        return ret_bonds
            

    def pisa_handle_xml_interface(interface_entity):
        molecules = []
        bonds = []

        # molecule info is either False or a pair of a new
        # molecule number and a pisa molecule record (which contains
        # things like pvalue, residues, id, class).
        for molecule in interface_entity.getiterator("molecule"):
            molecule_info = pisa_handle_xml_molecule(imol, molecule, 'interfaces')
            if molecule_info:
                molecules.append(molecule_info)
        for bond_type in ['h-bonds', 'salt-bridges', 'ss-bonds', 'cov-bonds']:
            bonds += (molecule_bonds(interface_entity.find(bond_type)))
        pvalue  = interface_entity.find("pvalue").text
        area    = interface_entity.find("int_area").text
        solv_en = interface_entity.find("int_solv_en").text
        stab_en = interface_entity.find("stab_en").text

        return [molecules, bonds, area, solv_en, pvalue, stab_en]

    def pisa_handle_pdb_entry(pdb_entry_entity):
        interfaces = []
        for interface_xml in pdb_entry_entity.getiterator("interface"):
            interface = pisa_handle_xml_interface(interface_xml)
            if interface:
                interfaces.append(interface)
        return interfaces
    
    # main line
    #
    pdb_entries = xml_entity.getiterator("pdb_entry")
    for pdb_entry_xml in pdb_entries:
        pdb_entry = pisa_handle_pdb_entry(pdb_entry_xml)
        handle_pisa_interfaces(pdb_entry)
        
        
