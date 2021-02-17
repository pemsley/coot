
import os
import sys
import re
import urllib2
from xml.etree import ElementTree

# "[[Rosuvastatin]]"

# return a dictionary, which might contain key 'REDIRECT'
#
# data is a tuple of db-id e.g. "DrugBank", "ChemSpider" and
#                    an entry_id (string) "e.g. DB01234" or "236" or False
#
def parse_wiki_drug_xml(tree, was_redirected=False):

    def unbracket(redirect_str):
        # print "unbracket this:", redirect_str
        bracket_open_match  = redirect_str.find('[[')
        bracket_close_match = redirect_str.find(']]')
        # print "bracket_open_match",  bracket_open_match
        # print "bracket_close_match", bracket_close_match
        s = redirect_str[bracket_open_match+2:bracket_close_match]
        # print "unbracketted ->", s
        return s

    drug_bank_id = False
    drugbank_str = "DrugBank *="
    chemspid_str = "ChemSpiderID *="
    pubchem_str  = "PubChem *="
    drugbank_re = re.compile(drugbank_str)
    chemspid_re = re.compile(chemspid_str)
    pubchem_re  = re.compile(chemspid_str)
    redirect_re = re.compile("#REDIRECT")
    query_ele = tree.find("query")
    db_dict = {}
    try:
        rev_iter = query_ele.getiterator("rev")
        chemspider_id = False
        if len(rev_iter) > 0:
            rev_text = rev_iter[0].text
            decode_text = rev_text.encode('ascii', 'ignore')
            lines = decode_text.split("\n")
            for line in lines:
                # print line
                if drugbank_re.search(line):
                    drugbank_id = line.split('=')[-1].strip()
                    db_dict["DrugBank"] = drugbank_id
                if chemspid_re.search(line):
                    chemspider_id = line.split('=')[-1].strip()
                    db_dict["ChemSpiderID"] = chemspider_id
                if pubchem_re.search(line):
                    pubchem_id = line.split('=')[-1].strip()
                    db_dict["PubChem"] = pubchem_id
                if redirect_re.search(line):
                    if was_redirected:
                        print("oops - found #REDIRECT but was redirected already!")
                        return db_dict
                    else:
                        inn = unbracket(line)
                        db_dict['REDIRECT'] = inn
                        return db_dict
            print("returning db_dict: {}".format(db_dict))
            return db_dict # happy case
    except:
        print("wikipedia didn't understand the query, bad format?")
        return {}

    return {}


# as above:
# return a dictionary, which might contain key 'REDIRECT'
#
# data is a tuple of db-id e.g. "DrugBank", "ChemSpider" and
#                    an entry_id (string) "e.g. DB01234" or "236" or False
#
def name_to_web_chem_db_id(drug_name, was_redirected=False):

    def unspace(name):
        match = name.find(" ")
        if match != -1:
            new_name = name[:match] + "%20" + name[match+1:]
            return unspace(new_name)
        else:
            return name

    debug = False
    dn = drug_name;
    # don't lower-case redirections: e.g. "mdma" -> "MDMA"
    if not was_redirected:
        dn = drug_name.lower()
    dn = unspace(dn)
    url = "http://en.wikipedia.org/w/api.php?format=xml&action=query&titles=" + \
           dn + "&prop=revisions&rvprop=content"
    if debug:
        print("fetching wikipedia url:", url)

    response = urllib2.urlopen(url)
    xml_string = response.read()
    xml_tree = ElementTree.fromstring(xml_string)
    db_dict = parse_wiki_drug_xml(xml_tree, was_redirected)
    if debug:
        print("parse_wiki_drug_xml db_dict:", db_dict)
        fn = drug_name + '.xml'
        # DeepCode hates this:
        # f = open(fn, 'w+')
        # f.write(xml_string)
        # f.close()
    if 'REDIRECT' in db_dict:
        return name_to_web_chem_db_id(db_dict['REDIRECT'], True)
    return db_dict
    

def fetch_molecule(drug_name):

    db_dict = name_to_web_chem_db_id(drug_name)

    # print "in fetch_molecule(), db_dict:", db_dict
    try_chemspider = False # usually not needed
    if "DrugBank" in db_dict:
       try:
          db_mol_uri = "https://www.drugbank.ca/structures/small_molecule_drugs/" + \
                       db_dict["DrugBank"] + ".mol"
          file_name = db_dict["DrugBank"] + ".mol"
          fetch_it = True
          if os.path.exists(file_name):
              if os.stat(file_name).st_size > 0:
                  fetch_it = False

          if fetch_it:
             response = urllib2.urlopen(db_mol_uri)
             xml_string = response.read()
             f = open(file_name, "w")
             f.write(xml_string)
             f.close()
          # print "returning", file_name
          return file_name
       except urllib2.HTTPError as e:
           print(e)
           try_chemspider = True
    else:
        try_chemspider = True

    if try_chemspider:
       if "ChemSpiderID" in db_dict:
          db_mol_uri = "http://www.chemspider.com/FilesHandler.ashx?type=str&striph=yes&id="
          db_mol_uri += db_dict['ChemSpiderID']
          file_name = "ChemSpider-" + db_dict['ChemSpiderID'] + ".mol"
          fetch_it = True
          if os.path.exists(file_name):
              if os.stat(file_name).st_size > 0:
                  fetch_it = False

          if fetch_it:
              response = urllib2.urlopen(db_mol_uri)
              mol_string = response.read()
              f = open(file_name, "w")
              f.write(mol_string)
              f.close()
              
          # print "returning", file_name
          return file_name
       else:
          print("no ChemSpider in dictionary")
          return False

    #if we got here, bad news
    return False


if __name__ == "__main__":

   mol_name = "crestor"
   if len(sys.argv) > 1:
       mol_name = sys.argv[1]
   fetch_molecule(mol_name)
   # fetch_molecule("Rosuvastatin")
   # fetch_molecule("benzene") # ChemSpider
