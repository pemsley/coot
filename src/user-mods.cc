
#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#include "python-3-interface.hh"
#endif

#include <fstream>
#include <stdexcept>
#include <sstream>

#include "utils/coot-utils.hh"
#include "coot-utils/coot-coord-utils.hh"

#include "user-mods.hh"
#include "guile-fixups.h"



// cc-interface.hh has lots of dependenciees

#include <gtk/gtk.h> // needed for cc-interface.hh
#include "cc-interface.hh"


#undef  USER_MODS_STANDALONE

// throw an exception when the conversion fails
// int str_to_int(const std::string &s);
// throw an exception when the conversion fails
// float str_to_float(const std::string &s);

coot::flips_container::flips_container(const std::string &file_name) {

   handle_user_mods(file_name);
}
   

void
coot::flips_container::parse_set_or_single(const std::string &line) {

   if (line.length() > 56) {
      // It seems at some point the format got changes, so we need to offset by -1
      // for older versions of reduce. We check for position of "sc:" to find the
      // offset
      std::size_t sc_pos = line.find("sc=");
      int offset = 0;
      if (sc_pos == 45) {
         offset = -1;
      }
      std::string set_string = line.substr(10+offset,7);
      std::string chain_id = line.substr(19+offset,1);
      std::string resno_string = line.substr(20+offset,4);
      std::string res_type =     line.substr(25+offset,3);
      std::string atom_name =    line.substr(28+offset,4);
      std::string score_string = line.substr(49+offset,8);
      std::string info_string  = line.substr(57+offset,line.length()-57);
      if (atom_name == "    ")
	 atom_name = " CA "; // e.g. for amide.
      try { 
	 int resno_int = coot::util::string_to_int(resno_string);
	 if (false) { // debug
	    std::cout <<  "chain_id :" << chain_id << ": ";
	    std::cout << "resno_string :" << resno_string << ": ";
	    std::cout << "resno_int :" << resno_int << ": ";
	    std::cout << "res_type :" << res_type << ": ";
	    std::cout << "atom_name :" << atom_name << ": ";
	    std::cout << "score_string :" << score_string << ": ";
	    std::cout << "info_string :" << info_string << ": ";
	 }
	 try {
	    float score = coot::util::string_to_float(score_string);
	    if (0)
	       std::cout << "score :" << score << ":" << std::endl;
	    std::string ins_code = "";
	    std::string alt_conf = "";
       std::string user_mod_rest = line.substr(57+offset);
	    coot::atom_spec_t as(chain_id, resno_int, ins_code,
				 atom_name, alt_conf);
	    flips_container::flip fl(as, set_string, user_mod_rest,
				     res_type, score);
	    if (user_mod_rest.length() > 0)
	       if (user_mod_rest.substr(0,1) == "!")
		  std::cout << "WARNING:: Best is a close contact"
			    << std::endl;
	    flips.push_back(fl);
	 }
	 catch (const std::runtime_error &score_err) {
	    // std::cout << "Ooops " << score_err.what() << std::endl;
	 } 
      }
      catch (const std::runtime_error &rte) {
	 // std::cout << "Ooops " << rte.what() << std::endl;
      }
      // std::cout << "\n";
   }
}

void
coot::flips_container::parse_no_adj(const std::string &line) {

   if (line.length() > 66) {

      // It seems at some point the format got changes, so we need to offset by -1
      // for older versions of reduce. We check for position of "sc:" to find the
      // offset

      std::size_t first = line.find(":");
      std::size_t second = line.find(":", first+1);
      int offset[3];
      std::vector<coot::atom_spec_t> specs;
      for (int ii=0; ii<3; ii++) { 
	 offset[ii] = ii*16;
    if (second == 32) {
       offset[ii] += -1;
    }
    std::string chain_id_1     = line.substr(offset[ii] + 19, 1);
	 std::string resno_string_1 = line.substr(offset[ii] + 20, 4);
	 std::string res_type_1 =     line.substr(offset[ii] + 25, 3);
	 std::string atom_name_1 =    line.substr(offset[ii] + 28, 4);
	 std::cout <<  " ii: " << ii << " ";
	 std::cout <<  "NoAdj: chain_id :" << chain_id_1 << ": ";
	 std::cout << "resno_string :" << resno_string_1 << ": ";
	 try {
	    int resno_int_1 = coot::util::string_to_int(resno_string_1);
	    std::cout << " resno_int_1 :" << resno_int_1 << ": ";
	    std::cout << "res_type_1 :" << res_type_1 << ": ";
	    std::cout << "atom_name_1 :" << atom_name_1 << ": ";
	    std::cout << std::endl;
	    coot::atom_spec_t spec_1(chain_id_1, resno_int_1, "", atom_name_1, "");
	    specs.push_back(spec_1);
	 }
	 catch (const std::runtime_error &rte) {
	    std::cout << rte.what() << std::endl;
	 }
      }
      std::string info_string = line.substr(66, line.length()-66);
      std::cout << "info_string :" << info_string << ": " << std::endl;
      no_adjust na(specs, info_string);
      no_adjustments.push_back(na);
   }
}


void
coot::flips_container::store(const std::vector<std::string> &lines) {

   for(unsigned int iline=0; iline<lines.size(); iline++) {
      // std::cout << "parsing line: " << lines[iline] << std::endl;
      if (lines[iline].length() > 17) { 
	 if (lines[iline].substr(0,6) == "ATOM  ")
	    break;
	 if (lines[iline].substr(10,6) == "Single")
	    parse_set_or_single(lines[iline]);
	 if (lines[iline].substr(10,3) == "Set")
	    parse_set_or_single(lines[iline]);
	 if (lines[iline].substr(10,5) == "NoAdj")
	    parse_no_adj(lines[iline]);
      }
   }
}

std::vector<std::string>
coot::flips_container::get_user_mods(const std::string &filename) const {

   std::vector<std::string> v;

   char line[10000];

   std::ifstream s(filename.c_str());

   if (s) {
      while (! s.eof()) {
	 s.getline(line, 9999);
	 v.push_back(line);
      }
   }
   return v;
}

void 
coot::flips_container::handle_user_mods(const std::string &filename) {

   if (coot::file_exists(filename)) { 
      std::vector<std::string> user_mod_strings = 
	 get_user_mods(filename);
      store(user_mod_strings);
   } else {
      std::cout << "File does not exist: " << filename << std::endl;
   }
}


#ifdef USER_MODS_STANDALONE
int main(int argc, char **argv) {

   if (argc < 2) {
      std::cout << "Usage: " << argv[0] << " pdb-file-name" << std::endl;
   } else {
      coot::flips_container f;
      f.handle_user_mods(argv[1]);
   }
   return 0;
}
#endif 

#ifdef USER_MODS_STANDALONE
// throw an exception when the conversion fails
int str_to_int(const std::string &s) {

   std::istringstream iss(s);
   int i;

   if (not (iss >> i)) {
      std::string mess("Fail to convert '");
      mess += s;
      mess += "' to a string";
      throw std::runtime_error(mess);
   }
   return i;
}
#endif


#ifdef USER_MODS_STANDALONE

// throw an exception when the conversion fails
float str_to_float(const std::string &s) {

   std::istringstream iss(s);
   float f;

   if (not (iss >> f)) {
      std::string mess("Fail to convert '");
      mess += s;
      mess += "' to a string";
      throw std::runtime_error(mess);
   }
   return f;
} 
#endif


#ifdef USE_GUILE
SCM
coot::flips_container::user_mods() const {

   SCM r = SCM_EOL;
   SCM f_l = SCM_EOL;
   SCM na_l = SCM_EOL;
   for(unsigned int i=0; i<flips.size(); i++) {

      // make a list e.g.
      // (list atom_spec residue-type info-string set-string score)
      // 
      SCM flip_scm = SCM_EOL;
      flip_scm = scm_cons(scm_double2num(flips[i].score), flip_scm);
      flip_scm = scm_cons(scm_from_locale_string(flips[i].set_string.c_str()), flip_scm);
      flip_scm = scm_cons(scm_from_locale_string(flips[i].info_string.c_str()), flip_scm);
      flip_scm = scm_cons(scm_from_locale_string(flips[i].residue_type.c_str()), flip_scm);
      SCM atom_spec_scm = atom_spec_to_scm(flips[i].atom_spec);
      flip_scm = scm_cons(atom_spec_scm, flip_scm);
      f_l = scm_cons(flip_scm, f_l);
   }
   for(unsigned int i=0; i<no_adjustments.size(); i++) {
      SCM no_adj_scm = SCM_EOL;
      no_adj_scm = scm_cons(scm_from_locale_string(no_adjustments[i].info_string().c_str()), no_adj_scm);
      SCM specs_scm = SCM_EOL;
      for (unsigned int ispec=0; ispec<no_adjustments[i].atom_specs.size(); ispec++) {
	 specs_scm = scm_cons(atom_spec_to_scm(no_adjustments[i].atom_specs[ispec]), specs_scm);
      }
      specs_scm = scm_reverse(specs_scm);
      no_adj_scm = scm_cons(specs_scm, no_adj_scm);
      na_l = scm_cons(no_adj_scm, na_l);
   }
   r = scm_cons(scm_reverse(na_l), r);
   r = scm_cons(scm_reverse(f_l), r);
   return r;
}
#endif // USE_GUILE


#ifdef USE_PYTHON
PyObject *
coot::flips_container::user_mods_py() const {

   PyObject *r = PyList_New(2);
   PyObject *flips_list = PyList_New(0);
   PyObject *no_adj_list = PyList_New(0);
   for (unsigned int iflip=0; iflip<flips.size(); iflip++) {

      // make a list e.g.
      // (list atom_spec residue-type info-string set-string score)
      // 
      PyObject *flip_py = PyList_New(5);
      PyObject *atom_spec_py = atom_spec_to_py(flips[iflip].atom_spec);
      PyList_SetItem(flip_py, 0, atom_spec_py);
      PyList_SetItem(flip_py, 1, myPyString_FromString(flips[iflip].residue_type.c_str()));
      PyList_SetItem(flip_py, 2, myPyString_FromString(flips[iflip].info_string.c_str()));
      PyList_SetItem(flip_py, 3, myPyString_FromString(flips[iflip].set_string.c_str()));
      PyList_SetItem(flip_py, 4, PyFloat_FromDouble(flips[iflip].score));
      PyList_Append(flips_list, flip_py);
      Py_XDECREF(flip_py);
   }
   // An adjustment is 2 items: first is a list of atom specs, second is a info-string
   for (unsigned int ina=0; ina<no_adjustments.size(); ina++) {
      PyObject *no_adjust_py = PyList_New(2);
      PyObject *info_string_py = myPyString_FromString(no_adjustments[ina].info_string().c_str());
      PyObject *no_adjust_atom_spec_list_py = PyList_New(no_adjustments[ina].atom_specs.size());
      for (unsigned int ispec=0; ispec<no_adjustments[ina].atom_specs.size(); ispec++) {
	 PyObject *atom_spec_py = atom_spec_to_py(no_adjustments[ina].atom_specs[ispec]);
	 PyList_SetItem(no_adjust_atom_spec_list_py, ispec, atom_spec_py);
      }
      PyList_SetItem(no_adjust_py, 0, no_adjust_atom_spec_list_py);
      PyList_SetItem(no_adjust_py, 1, info_string_py);
      PyList_Append(no_adj_list, no_adjust_py);
      Py_XDECREF(no_adjust_py);
   }
   PyList_SetItem(r, 0, flips_list);
   PyList_SetItem(r, 1, no_adj_list);
   return r;
}
#endif // USE_PYTHON
