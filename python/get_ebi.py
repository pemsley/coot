# get-ebi.py
# Copyright 2005, 2006 by Bernhard Lohkamp
# Copyright 2005, 2006 by Paul Emsley, The University of York
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


# this is get-ebi stuff
global oca_server
#oca_server = "http://bip.weizmann.ac.il"
#oca_server = "http://oca.ebi.ac.uk"
#oca_server = "http://structure.embl-hamburg.de"
oca_server = "http://www.ebi.ac.uk/msd-srv/oca"

global oca_pdb_request_stub
oca_pdb_request_stub = "oca-bin/save-pdb?id="
#oca_pdb_request_stub = "oca-bin/send-pdb?id="
#oca_pdb_request_stub = "oca-bin/send-x-pdb?id="

global pdb_file_name_tail
pdb_file_name_tail = ""

global oca_sfs_request_stub 
oca_sfs_request_stub = "oca-bin/send-sf?id="
global oca_sfs_request_tail 
oca_sfs_request_tail = ""
#oca_sfs_request_tail = "sf.ent.Z"

global coot_tmp_dir
coot_tmp_dir = "coot-download"

# e.g. (ebi-get-pdb "1crn")
# 
# no useful return value
#
# Note that that for sf data, we need to construct something like the
# string: http://oca.ebi.ac.uk/oca-bin/send-sf?r2acesf.ent.Z and we
# don't need to strip any html (thank goodness). Also not that the
# accession code now is lower case.
#
# data-type(here string) can be 'pdb' or 'sfs' (structure factors). 
# We might like to use
# 'coordinates rather than 'pdb in the future.
# 
# The optional argument imol-coords-arg-list is necessary for
# ouptutting sfs, because we need coordinates from which we can
# calculate phases.
#

# we dont need something like net-get-url in python 
# since we have build in functions like urlretrieve (in module urllib)

# check the directory and get url url_string.
#
def check_dir_and_get_url(dir,file_name,url_string):
    import os,urllib

    if (os.path.isfile(dir) or os.path.isdir(dir)):
       if (os.path.isfile(dir)):
          print dir," is atually a file and not a dir, so we can't write to it"
       else:
          if (os.path.isdir(dir)):
             try:
                urllib.urlretrieve(url_string, file_name)
             except IOError:
                print "ERROR:: We can't get ",url_string
          else:
             print "ERROR:: Oops - Can't write to ",dir," directory!"
    else:
       os.mkdir(dir)
       if (os.path.isdir(dir)):
          try:
             urllib.urlretrieve(url_string, file_name)
          except IOError:
             print "ERROR:: We can't get ",url_string
       else:
         print "ERROR:: Oops - create-directory ",dir," failed!"

# get url_string for data type (string actually) 'pdb' or 'sfs'
#
def get_url_str(id, url_string, data_type, imol_coords_arg_list):
    import operator

    #print "DEBUG:: in get_url_string:", id, url_string, data_type

    if (data_type == "pdb"):
       pdb_file_name = coot_tmp_dir + "/" + id + ".pdb" + pdb_file_name_tail
       check_dir_and_get_url(coot_tmp_dir,pdb_file_name,url_string)
       imol_coords = handle_read_draw_molecule(pdb_file_name)
       return imol_coords

    if (data_type == "sfs"):
       sfs_file_name = coot_tmp_dir + "/" + id + ".cif"
#       print "BL DEBUG:: cif output file is: ",sfs_file_name
       imol_coords = imol_coords_arg_list
       if (operator.isNumberType(imol_coords) and imol_coords>=-1):
         check_dir_and_get_url(coot_tmp_dir,sfs_file_name,url_string)
         read_cif_data(sfs_file_name,imol_coords_arg_list) 

# Get the pdb and sfs. @var{id} is the accession code
def get_ebi_pdb_and_sfs(id):
    import operator,string

    imol_coords = get_ebi_pdb(id)
#    print "BL DEBUG:: imol_coords is (3)", imol_coords
    if (not operator.isNumberType(imol_coords)):
       print "Failed at reading coordinates. imol-coords was ",imol_coords

    if (imol_coords < 0):	# -1 is coot code for failed read.
       print "failed to read coordinates."
    else:
       down_id = string.lower(id)
       url_str = oca_server + "/" + oca_sfs_request_stub + down_id + oca_sfs_request_tail
#       print "BL DEBUG:: get_url_str with",id,url_str,"sfs",imol_coords
       get_url_str(id,url_str,"sfs",imol_coords)

# Return a molecule number on success
# or not a number (False) or -1 on error.
#
def get_ebi_pdb(id):
    import urllib, string

    up_id = string.upper(id)
    url_str = oca_server + "/" + oca_pdb_request_stub + up_id
    imol_coords = get_url_str(id,url_str,"pdb",None)
    return imol_coords


# Get data and pdb for accession code id from the Electron Density
# Server.
#
# @var{id} is the accession code.
#
# returns imol of read pdb or False on error.
#
# 20050725 EDS code
#
def get_eds_pdb_and_mtz(id):
    import string
    import urllib

    # Gerard DVD Kleywegt says we can find the coords/mtz thusly:
    #
    # - model = http://eds.bmc.uu.se/eds/sfd/1cbs/pdb1cbs.ent
    # - mtz   = http://eds.bmc.uu.se/eds/sfd/1cbs/1cbs_sigmaa.mtz
    #
    # 20091222
    # - newprefix: http://eds.bmc.uu.se/eds/dfs/cb/1cbs/
    # 
    # URL:: "http://eds.bmc.uu.se/eds/sfd/sa/2sar/pdb2sar.ent"
    # URL:: "http://eds.bmc.uu.se/eds/sfd/sa/2sar/2sar_sigmaa.mtz"

    eds_site = "http://eds.bmc.uu.se/eds"

    # "1cbds" -> "cb/"
    def mid_chars(id_code):
        if not id_code:  # check for string?
            return "//fail//"
        if not (len(id_code) == 4):
            return "/FAIL/"
        else:
            return id_code[1:3] + "/"

    r = coot_mkdir(coot_tmp_dir)
  
    if (r):
      down_id = string.lower(id)
      eds_url = eds_site + "/dfs/"
      target_pdb_file = "pdb" + down_id + ".ent"
      dir_target_pdb_file = coot_tmp_dir + "/" + target_pdb_file
      mc = mid_chars(down_id)
      model_url = eds_url + mc + down_id + "/" + target_pdb_file
      target_mtz_file = down_id + "_sigmaa.mtz"
      dir_target_mtz_file = coot_tmp_dir + "/" + target_mtz_file
      mtz_url = eds_url + mc +down_id + "/" + target_mtz_file

      try:
        s1 = urllib.urlretrieve(model_url, dir_target_pdb_file)
        print "INFO:: read model status: ",s1
      except IOError:
        print "BL ERROR:: We can't open ", model_url
      try:
        s2 = urllib.urlretrieve(mtz_url, dir_target_mtz_file)
        print "INFO:: read mtz   status: ",s2
      except IOError:
        print "BL ERROR:: We can't open ", mtz_url 


      r_imol = handle_read_draw_molecule(dir_target_pdb_file)
      sc_map = make_and_draw_map(dir_target_mtz_file,"2FOFCWT","PH2FOFCWT","",0,0)
      make_and_draw_map(dir_target_mtz_file,"FOFCWT","PHFOFCWT","",0,1)
      set_scrollable_map(sc_map)
      if (valid_model_molecule_qm(r_imol)):
          return r_imol
      else:
          return False

    else:
      print "Can't make directory ",coot_tmp_dir



# BL says: to test, some examples

#id = "2BSX"
#get_ebi_pdb_and_sfs(id)
#get_eds_pdb_and_mtz(id)
