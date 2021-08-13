

# simply read the file and make a molecule
#
def generate_molecule_from_mmcif(comp_id, mmcif_file_name):

   read_cif_dictionary(mmcif_file_name)
   return get_monomer(comp_id)

# return a molecule number. Return -1 on fail.
#
def generate_molecule_from_mmcif_by_dict_gen(comp_id, mmcif_file_name):
    """return a molecule number. Return -1 on fail."""

    global use_mogul
    import os

    # return a new molecule number, generator is "pyrogen" or "acedrg"
    #
    def dict_gen(generator, args, working_dir):

        stub = comp_id + "-" + generator
        log_file_name = os.path.join(working_dir, stub + ".log")

        status = popen_command(generator, args, [], log_file_name, True)
        if (status != 0):
            return -1 # bad molecule
        else:
            pdb_name = os.path.join(working_dir, stub + ".pdb")
            cif_name = os.path.join(working_dir, stub + ".cif")

            imol = read_pdb(pdb_name)
            read_cif_dictionary(cif_name)
            return imol

    if not os.path.isfile(mmcif_file_name):
        return -1 # fail
    
    if enhanced_ligand_coot_p():
        # Use pyrogen if we have mogul
        #
        if use_mogul:
            # pyrogen
            #
            working_dir = get_directory("coot-pyrogen")
            args = ["-r", comp_id, "-d", working_dir, "-c", mmcif_file_name]
            dict_gen("pyrogen", args, working_dir)
        else:
            # acedrg
            #
            working_dir = get_directory("coot-acedrg")
            stub = os.path.join(working_dir, comp_id + "-acedrg")
            args = ["-r", comp_id, "-c", mmcif_file_name, "-o", stub]
            dict_gen("acedrg", args, working_dir)
    else:
        # acedrg
        #
        working_dir = get_directory("coot-acedrg")
        args = ["-M", "-r", comp_id, "-c", mmcif_file_name]
        # BL says:: different args tp above?! Wonder why!? Need to check.
        dict_gen("acedrg", args, working_dir)
                
            

