a_molecule = [
				# a list of models
		[		# a list of chains
			["A",[
				[1081,"","ASN",[
					[[" O  ",""],
					[1.0,22.87," O"],
					[0.287,20.822,40.742]],
					[[" C  ",""],
					[1.0,22.17," C"],
					[0.303,20.8649,39.5379]],
					[[" ND2",""],
					[1.0,29.17," N"],
					[1.363,16.5869,37.908]],
					[[" OD1",""],
					[1.0,32.810001373291," O"],
					[2.877,18.049,38.597]],
					[[" CG ",""],
					[1.0,27.4," C"],
					[1.711,17.676,38.529]],
					[[" CB ",""],
					[1.0,24.58," C"],
					[0.547,18.528,39.102]],
					[[" CA ",""],
					[1.0,23.18," C"],
					[0.657,19.99,38.7]],
					[[" N  ",""],
					[1.0,22.87," N"],
					[0.531,20.0909,37.2929]]]
				],
				[1082,"","GLN",[
					[[" O  ",""],
					[1.0,19.223," O"],
					[0.419,24.068,39.868]],
					[[" C  ",""],
					[1.0,20.54," C"],
					[1.431,23.549,40.3401]],
					[[" NE2",""],
					[1.0,14.77," N"],
					[5.637,22.1359,39.458]],
					[[" OE1",""],
					[1.0,25.2299," O"],
					[5.2859,23.3729,41.34]],
					[[" CD ",""],
					[1.0,25.7," C"],
					[5.08,23.136,40.01]],
					[[" CG ",""],
					[1.0,19.98," C"],
					[4.176,24.007,39.253]],
					[[" CB ",""],
					[1.0,21.63," C"],
					[3.038,23.178,38.604]],
				  	[[" CA ",""],
					[1.0,20.56," C"],
					[2.118,22.47,39.593]],
				  	[[" N  ",""],
					[1.0,21.99," N"],
					[1.156,21.65,38.894]]]
				]
			    ]	# end of chain A content
			]	# end of chain A
		]		# end of list of chains
	]			# end of list of models

# Return a python_mol
#
#
#
def jiggled_mol(reference_mol, current_mol, traj_frac):

	def jiggle_random():
		import random
		return (random.random() - 0.5)

	def jiggled_pos(ref_pos, current_pos):
		if (traj_frac < 0):	# magic value
			# make a starting set of coords
			return map(lambda x1: 
				jiggle_random() * 1.0 + x1,
				ref_pos)
		else:
			q = 1 - traj_frac 
			return map(lambda x, x_ref:
				traj_frac * x_ref + q * x + 
				0.4 * q * jiggle_random(),
				current_pos, ref_pos)

	def jiggled_atom(ref_atom, current_atom):
		ref_pos = ref_atom[2]
		cur_pos = current_atom[2]
		ret = [ref_atom[0], ref_atom[1], jiggled_pos(ref_pos, cur_pos)]
		return ret

	def jiggled_residue(ref_res, cur_res):
		ret = [ref_res[0], ref_res[1], ref_res[2], map(jiggled_atom, ref_res[3], cur_res[3])]
		return ret

	def jiggled_chain(ref_chain, cur_chain):
		ret = [ref_chain[0],
			       map(jiggled_residue, 
				   ref_chain[1],
				   cur_chain[1])]
		return ret

	def jiggled_model(ref_model, cur_model):
		ret = map(jiggled_chain, ref_model, cur_model)
		return ret

        j_mod = map (jiggled_model, reference_mol, current_mol)
	return j_mod


#
def disrupt(reference_mol, biggness):

	ret = jiggled_mol(reference_mol, reference_mol, -1)
	return ret

#
max_count = 5000.0

mol_no = add_molecule(a_molecule, "test molecule")
if (not mol_no == -1):
	set_rotation_centre(*centre_of_mass(mol_no))
	for count in range(int(max_count) + 1):
		current_mol = disrupt(a_molecule,0.8)
		new_mol = jiggled_mol(a_molecule, current_mol, 
			count / max_count)
		#print "cycle ", count, count/max_count, max_count
		clear_and_update_molecule(mol_no, new_mol)
		# gtk stuff
		if gtk.events_pending():
			gtk.main_iteration(False)
