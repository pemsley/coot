# lidia-core/test-fei-vs-mine.py
# 
# Copyright 2013 by Medical Research Council
# Author: Paul Emsley
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA
#

import os
import subprocess
import glob

def get_fei_types(comp_id):
    file_name = "../src/cod/AtomTypeTests2/"
    file_name += comp_id
    file_name += "_CodAtomType.txt"

    try:
        f = file(file_name)
        lines = f.readlines()
        atom_lines = lines[1:]
        types = []
        for atom_line in atom_lines:
            part = atom_line.split('\t')
            type = part[2][:-1]
            types.append(type)
        return types
    except:
        print 'problem with file', file_name


def get_pe_types(comp_id):
    something = subprocess.check_output(['./test-cod-atom-types', comp_id])
    lines = something.split('\n')
    hit_types = False
    types = []
    for line in lines:
        if hit_types == True:
            parts = line.split()
            if len(parts) == 2:
                types.append(parts[1])
        if line[:10] == 'PE-TYPES::':
            hit_types = True
    return types

    
# return success stats:
def diff_for_comp_id(comp_id):
    fei_types = get_fei_types(comp_id)
    # print 'fei_types:\n', fei_types
    pe_types = get_pe_types(comp_id)
    # print 'pe_types:\n', pe_types

    if len(fei_types) != len(pe_types):
        print 'Problem with mismatched type lengths for', comp_id, ' pe had', len(pe_types), 'types'
        for i in range(len(pe_types)):
            print '  ',i,' ',pe_types[i]
        return [10, 0] # fake all failures

    else:
        # happy path
        n_match = 0
        fail_match = []
        for i in range(len(fei_types)):
            if pe_types[i] == fei_types[i]:
                n_match += 1
            else:
                fail_match.append([i, pe_types[i], fei_types[i]])
        print '::: comp-id:',comp_id,'matched', n_match, 'of', len(fei_types)
        for failure in fail_match:
            print '   ', failure[0], " ", failure[1],  '\t', failure[2]
        # return success stats:
        return [len(fei_types), len(fei_types)-len(fail_match)]

def test_0s():

    test_set = ["0AD", "024", "090", "099", "008", "074", "00A", "023", "057", "039", "0CP",
                "0PJ", "032", "002", "059", "0E4", "047", "055", "064", "033", "01G", "0G6",
                "094", "0PQ", "066", "098", "003", "061", "0Z6", "096", "088", "0UN",
                "062", "0PN", "071", "0BD", "0JZ", "084", "0LI", "0CE", "000", "073", "017",
                "012", "065", "0CO", "0EZ", "0E3", "031", "0BI", "0E2", "0MA", "072", "009",
                "001", "006", "0AS", "01A", "0CL", "01K", "075", "0A5", "0PY", "0D5", "00C",
                "0CS", "0NC", "0A0", "0AH", "0AJ", "069", "0A6", "0AK", "0AZ", "004", "0CR",
                "0AB", "0A7", "0A8", "0PA", "00B", "0FA", "005", "0A9", "0AG", "0A1", "0BN",
                "0YG", "028", "0AF", "01W", "0AO", "0AI", "0AY", "018", "007", "020", "0A2",
                "0DC", "041", "00G", "0IN", "0DT", "01I", "0A4", "077", "0AP",
                "0DG", "093", "068", "0MO", "097", "0A3", "0HG", "06C",
                "042"]

#    test_set = ['06C']
#    test_set = ['024']
#    test_set = ['039']
#    test_set = ['0CP']  # bridged cyclohexane
#    test_set = ['03R']  # coords of ideal are 0, hence null molecule
#    test_set = ['017']  # bridge atom problem (I think ref types are wrong :)
#    test_set = ['0CE']  # same as 017
#    test_set = ['0JZ']  # macrocycle including Se, problem with Se order
#    test_set = ['093']  # 

    sum_trials = 0
    sum_success = 0
    for test_id in test_set:
        stats = diff_for_comp_id(test_id)
        sum_trials += stats[0]
        sum_success += stats[1]
    print '      overall: {0}/{1} = {2:.2f}%'.format(sum_success, sum_trials,
                                                     float(sum_success*100)/float(sum_trials))

def gen_all():

    output_dir = 'pe-types'
    dict_dir = os.path.expandvars('$HOME/ccp4/ccp4-6.3.0/lib/data/monomers')
    dirs = glob.glob(dict_dir + '/*')
    dirs = [ os.path.join(dict_dir, "8"),
	     os.path.join(dict_dir, "9") ]
    for dir in dirs:
	if (os.path.isdir(dir)):
	    dir_basename = os.path.basename(dir)
	    # print dir
	    results_dir = os.path.join(output_dir, dir_basename)
	    if (not (os.path.isdir(results_dir))):
		os.makedirs(results_dir)
	    print dir_basename, results_dir
	    files = glob.glob(dir + "/*.cif")
	    for file in files:
		# print file
		# get the comp-id from the file name
		p = os.path.splitext(os.path.basename(file))
		comp_id = p[0]
		# A8D is a genuine library error.
		# Note to self: what are the SMILES strings for these things!?
		# B1M now passes
		exclude = [        'BVA', 'B51', 'B13', 'BCB', 'BF4', 'BCL', 'B12', 'BEF', 'BFD',
			    '7HE', 'ITM', 'ICA', 'JM1', 'AG1', 'A8D', 'AC9', 'APW', 'ALB', 'CN1',
			    'CHL', 'DVT', 'DW1', 'DWC', 'DAQ', 'DAE', 'DW2', 'CO3', 'CLZ', 'CL7',
			    'E52', 'C2C', 'CL1', 'FCI', 'FEM', 'FDC', 'FLL', 'FS2', 'FNE', 'CLA',
			    'HB1', 'HC0', 'HEG', 'HF5', 'HC1', 'HCN', 'CFC', 'CMO', 'CNF', 'KEG',
			    'KYT', 'KYS', 'LCO', 'MGF', 'NFV', 'NMQ', 'NFC', 'MAP', 'ONP', 'OXX',
			    'OEC', 'MO7', 'PFC', 'PTE', 'PEJ', 'PNQ', 'PMR', 'PCD', 'PHF', 'ME3',
			    'MNQ', 'MF4', 'RU7', 'REO', 'REP', 'RUC', 'RTC', 'TBR', 'TL2', 'V7O',
			    'VEA', 'WO3', 'WO2', 'YBT', 'ZRC', '202', '34B', '39B', '39E' ]
		do_it = True
		for c in exclude:
		    if (comp_id == c):
			do_it = False
		if (do_it):
		    out_file = os.path.join(results_dir, comp_id + '.types')
		    print "    ", out_file
		    atom_types = get_pe_types(comp_id)
		    f = open(out_file, 'w')
		    for type in atom_types:
			f.write(type)
			f.write('\n')
		    f.close()

def types_hist():

    atom_types = {}
    output_dir = 'pe-types'
    dirs = glob.glob(output_dir + '/*')
    for dir in dirs:
        files = glob.glob(dir + "/*.types")
	for file in files:
	    f = open(file)
	    for line in f:
		type = line.rstrip()
		try:
		    atom_types[type] += 1
		except:
		    atom_types[type] = 1
	    f.close()
		    
    sum = 0

    if True:
	for key in atom_types.keys():
	    sum += atom_types[key]
	    print atom_types[key], "   ", key

    if False:
	l = []
	for value in atom_types:
	    t = [atom_types[value], value]
	    l.append(t)
	for i in l:
	    print i
    

if __name__ == "__main__":

    # test_0s()
    # gen_all()

    types_hist()

    
