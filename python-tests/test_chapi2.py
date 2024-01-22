import  chapi 
mc = chapi.molecules_container_t(True)
imol = mc.read_pdb('data/tutorial-modern.pdb')
residue_list = mc.get_residues_near_residue(imol, '//A/20', 6)
for r in residue_list:
    print(r.chain_id, r.res_no)
