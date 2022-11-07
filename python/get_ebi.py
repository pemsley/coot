
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

import os
import numbers
import coot
import coot_utils
import pdbe_validation_data


pdbe_server = "https://www.ebi.ac.uk"
pdbe_pdb_file_dir = "pdbe/entry-files/download"

pdbe_file_name_tail = "ent"

# sf example http://www.ebi.ac.uk/pdbe-srv/view/files/r4hrhsf.ent

# 20151126-PE No, we can't have coot-download created on coot-startup, it must be
#             made only when we need it.
# global coot_tmp_dir
# coot_tmp_dir = coot_utils.get_directory("coot-download")

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
# outputting sfs, because we need coordinates from which we can
# calculate phases.
#

# helper function to avoid downloading empty files
# returns download filename upon success or False when fail
# should we return the response header too?
#
def coot_urlretrieve(url, file_name, reporthook=None):

    """Helper function to avoid downloading empty files
    returns download filename upon success or False when fail."""

    import urllib.request, urllib.parse, urllib.error
    local_filename = False
    class CootURLopener(urllib.request.FancyURLopener):
        def http_error_default(self, url, fp, errcode, errmsg, headers):
            # handle errors the way you'd like to
            # we just pass
            pass

    # opener = CootURLopener(context=ssl_context)
    opener = CootURLopener()
    try:
        local_filename, header = opener.retrieve(url, file_name, reporthook)
    except:
        # we could catch more here, but dont bother for now
        print("BL WARNING:: retrieve of url %s failed" %url)

    return local_filename

# we dont need something like net-get-url in python 
# since we have build in functions like urlretrieve (in module urllib)

# check the directory and get url url_string.
#
def check_dir_and_get_url(dir, file_name, url_string):
    import os,urllib.request,urllib.parse,urllib.error

    # FIXME logic, can be done better
    if (os.path.isfile(dir) or os.path.isdir(dir)):
       if (os.path.isfile(dir)):
          print(dir, " is atually a file and not a dir, so we can't write to it")
       else:
          if (os.path.isdir(dir)):
              coot_urlretrieve(url_string, file_name)
          else:
              print("ERROR:: Oops - Can't write to ", dir, " directory!")
    else:
       os.makedirs(dir)
       if (os.path.isdir(dir)):
           coot_urlretrieve(url_string, file_name)
       else:
         print("ERROR:: Oops - create-directory ",dir," failed!")

# get url_string for data type (string actually) 'pdb' or 'sfs'
#
def get_url_str(id, url_string, data_type, imol_coords_arg_list):
    import operator

    #print "DEBUG:: in get_url_string:", id, url_string, data_type

    coot_download_dir = get_directory("coot-download")
    if data_type == "pdb":
       pdb_file_name = coot_download_dir + "/" + id + ".pdb." + \
           pdbe_file_name_tail
       check_dir_and_get_url(coot_download_dir, pdb_file_name, url_string)
       imol_coords = handle_read_draw_molecule(pdb_file_name)
       return imol_coords

    if data_type == "cif":
       pdb_file_name = coot_download_dir + "/" + id + ".cif"
       check_dir_and_get_url(coot_download_dir, pdb_file_name, url_string)
       imol_coords = handle_read_draw_molecule(pdb_file_name)
       return imol_coords

    if data_type == "sfs":
       sfs_file_name = coot_download_dir + "/" + id + ".cif"
       #       print "BL DEBUG:: cif output file is: ",sfs_file_name
       imol_coords = imol_coords_arg_list
       if (isinstance(imol_coords, numbers.Number) and imol_coords>=-1):
         check_dir_and_get_url(coot_tmp_dir, sfs_file_name, url_string)
         coot.read_cif_data(sfs_file_name, imol_coords_arg_list)
         # do we need to return something here too?!

# Get the pdb and sfs. @var{id} is the accession code
#
def get_ebi_pdb_and_sfs(id):

    import operator,string

    imol_coords = get_ebi_pdb(id)
    if (not isinstance(imol_coords, numbers.Number)):
       print("Failed at reading coordinates. imol-coords was ",imol_coords)

    if (imol_coords < 0):	# -1 is coot code for failed read.
       print("failed to read coordinates.")
    else:
       down_id = id.lower()
       url_str = pdbe_server + "/" + pdbe_pdb_file_dir + "/" + \
                 "r" + down_id + "sf." + \
                 pdbe_file_name_tail
       get_url_str(id, url_str, "sfs", imol_coords)

# Return a molecule number on success
# or not a number (False) or -1 on error.
#
def get_ebi_pdb(id):
    import urllib.request, urllib.parse, urllib.error, string

    # print "======= id:", id
    down_id = string.lower(id)
    pdb_url_str = pdbe_server + "/" + pdbe_pdb_file_dir + "/" + down_id + ".ent"
    cif_url_str = pdbe_server + "/" + pdbe_pdb_file_dir + "/" + down_id + ".cif"
    url_status = get_url_str(id, pdb_url_str, "pdb", None)
    # e.g. http://ftp.ebi.ac.uk/pub/databases/pdb +
    #      /validation_reports/cb/1cbs/1cbs_validation.xml.gz
    # print "BL DEBUG:: get-ebi-pdb ======= url-status", url_status
    if valid_model_molecule_qm(url_status):
        pdb_validate(down_id, url_status)
        return url_status
    else:
        cif_url_status = get_url_str(id, cif_url_str, "cif", None)
        if valid_model_molecule_qm(cif_url_status):
            # print "BL DEBUG:: get-ebi-pdb ======= cif_url_status", cif_url_status
            pdb_validate(down_id, cif_url_status)
            return cif_url_status

    return False


# Return a list of molecules (i.e. the model molecule and the 2 maps).
# or, if it didn't work then return False
#
# @var{id} is the accession code.
#
# returns imol of read pdb or False on error.
#
# 20050725 EDS code
#
# return a list of 3 molecule numbers [imol, map, diff_map] or False
#
#
def get_eds_pdb_and_mtz(id):
    import string
    import urllib.request, urllib.parse, urllib.error

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
    #
    # 20161010 new prefix
    # http://www.ebi.ac.uk/pdbe/coordinates/
    # http://www.ebi.ac.uk/pdbe/coordinates/files/1cbs_map.mtz

    def get_cached_eds_files(accession_code):
        down_code = accession_code.lower()
        dir_name = coot_utils.get_directory("coot-download")
        pdb_file_name = os.path.join(dir_name, "pdb" + down_code + ".ent")
        mtz_file_name = os.path.join(dir_name, down_code + "_map.mtz")

        print( "::::::::: pdb_file_name:", pdb_file_name)
        print( "::::::::: mtz_file_name:", mtz_file_name)
        if not os.path.isfile(pdb_file_name):
            return False
        else:
            if not os.path.isfile(mtz_file_name):
                return False
            else:
                imol = coot.read_pdb(pdb_file_name)
                imol_map = coot.make_and_draw_map(mtz_file_name, "FWT", "PHWT", "", 0, 0)
                imol_map_d = coot.make_and_draw_map(mtz_file_name, "DELFWT", "PHDELWT", "", 0, 1)
                if not (coot_utils.valid_model_molecule_qm(imol) and
                        coot_utils.valid_map_molecule_qm(imol_map) and
                        coot_utils.valid_map_molecule_qm(imol_map_d)):
                    coot.close_molecule(imol)
                    coot.close_molecule(imol_map)
                    coot.close_molecule(imol_map_d)
                    return False
                else:
                    return [imol, imol_map, imol_map_d]

    eds_site = "https://www.ebi.ac.uk/pdbe/coordinates"
    # https://www.ebi.ac.uk/pdbe/entry/pdb/6tje
    # eds_core = "some://thing" ;; for web pages
    eds_core = "https://www.ebi.ac.uk/pdbe/entry/pdb"
    eds_coords_site = "https://www.ebi.ac.uk/pdbe/entry-files/download"
    # now the map mtz files are like this:
    # https://www.ebi.ac.uk/pdbe/coordinates/files/zz/4zzn/4zzn_map.mtz

    # "1cbds" -> "cb/"
    #
    def mid_chars(id_code):
        if not id_code:  # check for string?
            return "//fail//"
        if not (len(id_code) == 4):
            return "/FAIL/"
        else:
            return id_code[1:3] + "/"

    # main line
    #
    cached_status = get_cached_eds_files(id)
    if isinstance(cached_status, list):
        return cached_status
    else:
        coot_tmp_dir = coot_utils.get_directory("coot-download")	
        r = coot_utils.coot_mkdir(coot_tmp_dir)

        if (r):
            down_id = id.lower()
            target_pdb_file = "pdb" + down_id + ".ent"
            target_cif_file = down_id + ".cif"
            dir_target_pdb_file = coot_tmp_dir + "/" + target_pdb_file
            dir_target_cif_file = coot_tmp_dir + "/" + target_cif_file
            model_url = eds_coords_site + "/" + target_pdb_file
            model_cif_url = eds_coords_site + "/" + target_cif_file
            target_mtz_file = down_id + "_map.mtz"
            dir_target_mtz_file = coot_tmp_dir + "/" + target_mtz_file
            # mtz_url = eds_site  + "/files/" + target_mtz_file
            #mtz_url = eds_site + "/files/" + mid_chars(down_id) + "/" + \
            #          down_id + "/" + down_id + "_map.mtz"
            mtz_url = eds_coords_site + "/" + down_id + "_map.mtz"
            eds_info_page = eds_core + "/" + down_id
            bad_map_status = False

            print("model_url:", model_url)
            print("  mtz_url:", mtz_url)

            # I am not sure this exists now
            # print("eds_info_page:", eds_info_page)

            try:
                # pre_download_info = coot.coot_get_url_as_string(eds_info_page)
                # print "INFO:: --------------- pre-download-info:", pre_download_info
                #  bad_map_status = "No reliable map available" in pre_download_info
                # if "There is no structure factor entry" in pre_download_info:
                #     print("WARNING:: no sfs available for entry %s, so wont download." %id)
                #     # no pdb and no mtz
                #    return False
                pass
            except:
                print("ERROR:: could not get pre_download_info from", eds_core)
                # we probably wont get anything else, so bail out.
                return False

            s1 = coot_urlretrieve(model_url, dir_target_pdb_file)
            s2 = coot_urlretrieve(mtz_url, dir_target_mtz_file)

            if bad_map_status:
                s = "This map (" + down_id + ") is marked by the EDS as \"not a reliable map\""
                coot.info_dialog(s)

            # maybe should then not load the map!?

            print("INFO:: read pdb model status: ",s1)
            print("INFO:: read mtz data  status: ",s2)

            if s1 and os.path.isfile(s1):
                r_imol = coot.handle_read_draw_molecule(dir_target_pdb_file)
                if not coot_utils.valid_model_molecule_qm(r_imol):
                    s1_cif = coot_urlretrieve(model_cif_url, dir_target_cif_file)
                    print("INFO:: read cif model status: ",s1_cif)
                    if (s1_cif == 0):
                        r_imol = coot.handle_read_draw_molecule(dir_target_pdb_file)
                        if not coot_utils.valid_model_molecule_qm(r_imol):
                            return False
                        else:
                            return r_imol
                    else:
                        return False
            if s2 and os.path.isfile(s2):
                map_1 = coot.make_and_draw_map(dir_target_mtz_file, "FWT", "PHWT","",0,0)
                map_2 = coot.make_and_draw_map(dir_target_mtz_file, "DELFWT", "PHDELWT",
                              "", 0, 1)
                coot.set_scrollable_map(map_1)
                return [map_1, map_2]  # r_imol not in this scope (at the moment)
            else:
                return False

        else:
            print("Can't make directory ",coot_tmp_dir)

# not sure if coot function is better or python script function coot_urlretrieve
# return 0 on success
def net_get_url(my_url, file_name):
    coot.coot_get_url(my_url, file_name)

def get_pdb_redo(text):

    if not isinstance(text, str):
        print("WARNING:: No string. No accession code.")
    else:
        if not (len(text) == 4):
            print("WARNING:: Accession code not 4 chars.")
        else:
            text = text.lower()
            stub = "https://pdb-redo.eu/db/" + \
                   text + "/" + text + "_final"
            pdb_file_name = text + "_final.pdb"
            mtz_file_name = text + "_final.mtz"
            py_file_name = text + ".py"
            url_pdb = stub + ".pdb"
            url_mtz = stub + ".mtz"
            url_py = stub + ".py"

            print("DEBUG:: getting", url_pdb)
            status = net_get_url(url_pdb, pdb_file_name)
            if not status == 0:
                print("Failed to get %s %s status %s" %(url_pdb, pdb_file_name, status))
            print("DEBUG:: getting", url_mtz)
            status = net_get_url(url_mtz, mtz_file_name)
            if not status == 0:
                print("Failed to get %s %s status %s" %(url_mtz, mtz_file_name, status))
            print("DEBUG:: getting", url_py)
            status = net_get_url(url_py, py_file_name)
            if not status == 0:
                print("Failed to get %s %s status %s" %(url_py, py_file_name, status))

            status_imol = read_pdb(pdb_file_name)
            if status_imol < 0:
                print("INFO:: problem opening pdb file. Most likely something went wrong in the download")
            else:
                print("DEBUG:: make-and-draw-map with", mtz_file_name)
                coot.make_and_draw_map(mtz_file_name, "FWT", "PHWT", "", 0, 0)
                coot.make_and_draw_map(mtz_file_name, "DELFWT", "PHDELWT", "", 0, 1)
                anom_map = coot.make_and_draw_map(mtz_file_name, "FAN", "PHAN", "", 0, 1)
                if anom_map > -1:
                    coot.set_map_colour(anom_map, 0.5, 0.5, 0)
                exec(compile(open(py_file_name, "rb").read(), py_file_name, 'exec'))
            

# BL says: to test, some examples

#id = "2BSX"
#get_ebi_pdb_and_sfs(id)
#get_eds_pdb_and_mtz(id)
