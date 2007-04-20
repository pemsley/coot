# Copyright 2004, 2005 by The University of York
# Copyright 2005 by Bernhard Lohkamp
 
#;;; This program is free software; you can redistribute it and/or modify
#;;; it under the terms of the GNU General Public License as published by
#;;; the Free Software Foundation; either version 2 of the License, or (at
#;;; your option) any later version.
 
#;;; This program is distributed in the hope that it will be useful, but
#;;; WITHOUT ANY WARRANTY; without even the implied warranty of
#;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#;;; General Public License for more details.
 
#;;; You should have received a copy of the GNU General Public License
#;;; along with this program; if not, write to the Free Software
#;;; Foundation, Inc.,  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


#;; direction is either 'forwards or 'backwards
#;; 
#;; start-resno is higher than stop-resno if we are building backwards
#;; 
#;; (fit-gap 0 "A"  23 26)   ; we'll build forwards
#;; (fit-gap 0 "A"  26 23)   ; we'll build backwards
#;;

global coot_tmp_dir
coot_tmp_dir = "coot-download"

def coot_mkdir(dir_name):
  import os
  if (os.path.isfile(dir_name)):
     return False
  else:
     if (os.path.isdir(dir_name)):
       return True
     else:
       os.mkdir(dir_name)
       return True


# this is get-ebi stuff
global oca_server
oca_server = "http://bip.weizmann.ac.il"
#oca_server = "http://oca.ebi.ac.uk"
#oca_server = "http://structure.embl-hamburg.de"

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

#now the proper stuff

def check_dir_and_get_url(dir,file_name,url_string):
    import os,urllib

#    print "BL DEBUG:: dir, filename, url_string: ",dir,file_name,url_string
  
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


def get_url_str(id,url_string,data_type,imol_coords_arg_list):
    import operator

    if (data_type == "pdb"):
       pdb_file_name = coot_tmp_dir + "/" + id + ".pdb" + pdb_file_name_tail
       check_dir_and_get_url(coot_tmp_dir,pdb_file_name,url_string)
#       handle_read_draw_molecule(pdb_file_name)
#       pdb_file_name = "bollocks.pdb"
       imol_coords = handle_read_draw_molecule(pdb_file_name)
#       print "BL DEBUG:: imol_coords is (1)", imol_coords
       return imol_coords

    if (data_type == "sfs"):
       sfs_file_name = coot_tmp_dir + "/" + id + ".cif"
       print "BL DEBUG:: cif output file is: ",sfs_file_name
       imol_coords = imol_coords_arg_list
#       print "BL DEBUG:: imol_coords is (2)", imol_coords
       if (operator.isNumberType(imol_coords) and imol_coords>=-1):
         check_dir_and_get_url(coot_tmp_dir,sfs_file_name,url_string)
         read_cif_data(sfs_file_name,imol_coords_arg_list) 


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


def get_ebi_pdb(id):
    import urllib, string

    up_id = string.upper(id)
    url_str = oca_server + "/" + oca_pdb_request_stub + up_id
    imol_coords = get_url_str(id,url_str,"pdb",None)
#    print "BL DEBUG:: imol_coords is (4)", imol_coords
    return imol_coords


def get_eds_pdb_and_mtz(id):
    import string
    import urllib

    #; Gerard DVD Kleywegt says we can find the coords/mtz thusly:
    #;
    #; - model = http://eds.bmc.uu.se/eds/sfd/1cbs/pdb1cbs.ent
    #; - mtz   = http://eds.bmc.uu.se/eds/sfd/1cbs/1cbs_sigmaa.mtz

    eds_site = "http://eds.bmc.uu.se/eds"

    r = coot_mkdir(coot_tmp_dir)
  
    if (r):
      down_id = string.lower(id)
      eds_url = eds_site + "/sfd/"
      target_pdb_file = "pdb" + down_id + ".ent"
      dir_target_pdb_file = coot_tmp_dir + "/" + target_pdb_file
      model_url = eds_url + down_id + "/" + target_pdb_file
      target_mtz_file = down_id + "_sigmaa.mtz"
      dir_target_mtz_file = coot_tmp_dir + "/" + target_mtz_file
      mtz_url = eds_url + down_id + "/" + target_mtz_file

      try:
        s1 = urllib.urlretrieve(model_url, dir_target_pdb_file)
        print "read model status: ",s1
      except IOError:
        print "BL INFO:: We can't open ", model_url
      try:
        s2 = urllib.urlretrieve(mtz_url, dir_target_mtz_file)
        print "read mtz   status: ",s2
      except IOError:
        print "BL INFO:: We can't open ", mtz_url 


      handle_read_draw_molecule(dir_target_pdb_file)
      sc_map = make_and_draw_map(dir_target_mtz_file,"2FOFCWT","PH2FOFCWT","",0,0)
      make_and_draw_map(dir_target_mtz_file,"FOFCWT","PHFOFCWT","",0,1)
      set_scrollable_map(sc_map,0)

    else:
      print "Can't make directory ",coot_tmp_dir



# BL says: to test, some examples

#id = "2BSX"
#get_ebi_pdb_and_sfs(id)
#get_eds_pdb_and_mtz(id)
