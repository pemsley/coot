/* coot-utils/coot-utils.cc
 *
 * Copyright 2004, 2005, 2006 by The University of York
 * Author: Paul Emsley
 * Copyright 2007 by Paul Emsley
 * Copyright 2014 by Medical Research Council
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

// Portability (to Windows, particularly) functions go here.
//

#include <iostream>
#include <algorithm>
#include <cstring>

#include <stdexcept> // for string_to_int.
#include <sstream>   // ditto.
#include <cstdio>    // 20090806 Justin Lecher says we need this on Gentoo
#include <iomanip>

#include <math.h>  // for fabs

#include "compat/coot-sysdep.h"
#if defined _MSC_VER
#include <direct.h>
#include <windows.h>
#include <initguid.h>
#include <lm.h>
#else
#if !defined(WINDOWS_MINGW)
#include <unistd.h>
#include <pwd.h>
#endif // MINGW
#include <glob.h>
#endif

#include <ctype.h>  // for toupper

#include <sys/types.h>  // for getpwnam
#include <sys/stat.h>   // for mkdir
#ifndef _MSC_VER
#include <unistd.h> // Also needed for mkdir on Fedora?
#endif

#include "coot-utils.hh"

#include "utils/logging.hh"
extern logging logger;


// below function sets this:
static std::string real_path_for_coot_executable;

//! do this on startup
void coot::set_realpath_for_coot_executable(const std::string &argv0) {

   // 20240902-PE patch from Charles - compiling on Windows
#ifdef _MSC_VER
   char *exec_path = _fullpath(NULL, argv0.c_str(), MAX_PATH);
#else
   char *exec_path = realpath(argv0.c_str(), NULL);
#endif

   if (exec_path) {
      // std::cout << "set_realpath_for_coot_executable(): got exec_path " << exec_path << std::endl;
   } else {
      std::cout << "ERROR::  set_realpath_for_coot_executable(): null exec_path " << std::endl;
      logger.log(log_t::ERROR, "set_realpath_for_coot_executable(): null exec_path",
                 std::string(strerror(errno)));
   }

   if (exec_path) {
      real_path_for_coot_executable = exec_path;
   }
}


std::string
coot::util::append_dir_dir (const std::string &s1, const std::string &dir) {

   std::string s;

   s = s1;
   s += "/";
   s += dir;

   return s;

}

std::string
coot::util::append_dir_file(const std::string &s1, const std::string &file) {

   std::string s;

   s = s1;
   s += "/";
   s += file;

   return s;
}

// Return the userid and name (e.g ("paule", "Paul Emsley")) for use
// as a label in the database.  Not important to get right.
//
std::pair<std::string, std::string> coot::get_userid_name_pair() {

   std::pair<std::string, std::string> p("unknown","unknown");
// BL says:: in windows we dont have USER and no getpwnam, so do somethign else
// and we avoid the ugly code below!!!
#if defined WINDOWS_MINGW
   const char *u = getenv("USERNAME");
   if (u) {
      p.first  = u;
      p.second = u;
   }
#else
   const char *u = getenv("USER");
#ifdef _MSC_VER
   // Man this is ugly windows code...
   LPUSER_INFO_10 pBuf = NULL;
   NET_API_STATUS nStatus;

   // Call the NetUserGetInfo function with level 10
   nStatus = NetUserGetInfo(NULL, (LPCWSTR) u, 10, (LPBYTE *)&pBuf);
   if (nStatus == NERR_Success) {
      if (pBuf) {
	 p.first  = (char *) pBuf->usri10_name;
	 p.second = (char *) pBuf->usri10_full_name;
      }
   }

#else
   if (u) {
      struct passwd *pwbits = getpwnam(u);
      std::string uid;
      std::string nam;
      p.first  = pwbits->pw_name;
      p.second = pwbits->pw_gecos;
   }
#endif // MSC
#endif // MINGW
   return p;
}


// Return like mkdir: mkdir returns zero on success, or -1 if an error
// occurred.
//
int
coot::util::create_directory(const std::string &dir_name_in) {

   int istat = -1;
   struct stat s;
   // on Windows stat works only properly if we remove the last / (if it exists)
   // everything else following seems to be fine with the /

   std::string dir_name = remove_trailing_slash(dir_name_in);

   int fstat = stat(dir_name.c_str(), &s);

   // 20060411 Totally bizarre pathology!  (FC4) If I comment out the
   // following line, we fail to create the directory, presumably
   // because the S_ISDIR returns true (so we don't make the
   // directory).

   if ( fstat == -1 ) { // file not exist
      // not exist
      // std::cout << "INFO:: Creating directory " << dir_name << std::endl;
      logger.log(log_t::INFO, "Creating directory", dir_name);

      bool change_permission = true;

#if defined(WINDOWS_MINGW) || defined(_MSC_VER)
      istat = mkdir(dir_name.c_str());
#else

#ifdef EMSCRIPTEN
      mode_t mkdir_mode = 0755;
      istat = mkdir(dir_name.c_str(), mkdir_mode);
#else
      if (change_permission) {
         mode_t mode = S_IRUSR|S_IWUSR|S_IXUSR; // over-ridden
         mode = 511; // octal 777
         mode_t mode_o = umask(0);
         mode_t mkdir_mode = mode - mode_o;
         istat = mkdir(dir_name.c_str(), mkdir_mode);
         umask(mode_o); // oh yes we do, crazy.
      }

   } else {
      if ( ! S_ISDIR(s.st_mode) ) {
	 // exists but is not a directory
	 istat = -1;
      } else {
	 // was a directory already
	 istat = 0; // return as if we made it
      }
#endif
#endif
   }

   // and finally some output: was the directory actually created?

   struct stat buf;
   int err = stat(dir_name.c_str(), &buf);
   if (err == 0) {
      istat = 0; // all OK then.
      // std::cout << "INFO:: in create_directory() " << dir_name << " confirmed as existing" << std::endl;
   } else {
      std::cout << "ERROR:: in create_directory() \"" << dir_name << "\" does not exist!" << std::endl;
   }

   return istat;
}


std::string
coot::util::intelligent_debackslash(const std::string &s) {

   std::string filename_str = s;

#ifdef WINDOWS_MINGW
   int slen = s.length();
   for (int i=0; i<slen; i++) {
       if (filename_str[i] == '\\') {
          filename_str.replace (i,1,"/");
       }
   }
#endif
  return filename_str;
}

std::string
coot::util::remove_trailing_slash(const std::string &s) {

   // BL says:: On Windows the null termination doesnt seem to work
   //    by subsituting null.
   // That was ugly anyway, let's erase the last character.
   std::string scratch = s;

   if (s.length() > 0) {
#ifdef HAVE_CXX11
      if (s.back() == '/')
         scratch.erase(scratch.end()-1);
      if (s.back() == '\\')
         scratch.erase(scratch.end()-1);
#else
      std::string g = scratch.substr(scratch.length()-1);
      if (scratch.substr(scratch.length()-1) == "/") {
	 std::string::iterator it = scratch.end();
         // scratch.erase(it-1);
	 std::string::size_type l = scratch.length();
	 scratch=scratch.substr(0,l-1);
      } else {
	 // we need an else, because above makes "/" -> ""
	 // so scratch.length would be 0
	 if (scratch.substr(scratch.length()-1) == "\\")
	    scratch.erase(scratch.end()-1);
      }
#endif
   }
   return scratch;
}

bool
coot::util::is_number(char c) {

   return ((c >= 48) && (c<=57));
}

bool
coot::util::is_letter(char c) {

   return (((c >= 65) && (c<=90)) || ((c >= 97) && (c<=122))) ;
}


bool
coot::is_member_p(const std::vector<std::string> &v, const std::string &a) {

   bool ir = 0;
   unsigned int vsize = v.size();

   for (unsigned int i=0; i<vsize; i++) {
      if (v[i] == a) {
	 ir = 1;
	 break;
      }
   }
   return ir;
}

// Some sort of polymorphism would be appropriate here, perhaps?
bool
coot::is_member_p(const std::vector<int> &v, const int &a) {

   bool ir = 0;
   unsigned int vsize = v.size();

   for (unsigned int i=0; i<vsize; i++) {
      if (v[i] == a) {
	 ir = 1;
	 break;
      }
   }
   return ir;
}

void
coot::remove_member(std::vector<int> *v_p, const int &a) {

   int vsize = v_p->size();

   for (int i=0; i<vsize; i++) {
      if ((*v_p)[i] == a) {
	 // shift down the others and reduct the vector length by one:
	 for (int j=(i); j<(vsize-1); j++) {
	    (*v_p)[j] =  (*v_p)[j+1];
	 }
	 v_p->resize(vsize-1);
      }
   }
}

bool
coot::util::even_p(int ii) {

   bool r = ((ii%2)==0);
   return r;
}

bool
coot::util::close_double_p(const double &d1, const double &d2, const double &diff_crit) {

   double d = fabs(d1-d2);
   return (d < diff_crit);
}



std::string
coot::util::Upper(const std::string &s) {

   std::string r = s;
   int nchars = s.length();
   for (int i=0; i<nchars; i++) {
      r[i] = toupper(s[i]);
   }
   return r;
}

std::string
coot::util::remove_whitespace(const std::string &s) {

   std::string r;
   int nchars = s.length();
   for (int i=0; i<nchars; i++) {
      if (s[i] != ' ') {
	 if (s[i] != '\n') {
	    if (s[i] != '\t') {
	       r += s[i];
	    }
	 }
      }
   }
   return r;
}

// "ALA X  " -> "ALA X";
std::string
coot::util::remove_trailing_whitespace(const std::string &s) {

   int sl = s.length();
   int cutpoint = 0;
   if (sl > 0) {
      for (int i=sl-1; i>=0; i--) {
	 if (s[i] != '\n' && s[i] != '\t' && s[i] != ' ') {
	    return s.substr(0,i+1);
	 }
      }
      return ""; // got to the end of a string of whitespaces
   } else {
      return ""; // input was empty string.
   }
}




std::string
coot::util::remove_leading_spaces(const std::string &s) {

   int sl = s.length();
   int cutpoint = 0;
   if (sl > 0) {
      for (int i=0; i<sl; i++) {
	 if (!(s[i] == ' ')) {
	    cutpoint = i;
	    break;
	 }
      }
   }
   return s.substr(cutpoint, sl);
}

// Remove the first bit from long
// e.g. remove_string("Cottage", "tag") -> "Cote";
//
std::string
coot::util::remove_string(const std::string &long_string, const std::string &bit) {

   std::string r = long_string;

   std::string::size_type ipos = long_string.find(bit);
   if (ipos != std::string::npos) {
      // OK, we found a match
      if ((ipos + bit.length()) < long_string.length())
	 r = r.substr(0,ipos) + r.substr((ipos+bit.length()));
      else
	 r = r.substr(0,ipos);
   }
   return r;

}





std::string
coot::util::int_to_string(int i) {
   char s[100];
   snprintf(s,99,"%d",i);
   return std::string(s);
}

std::string
coot::util::long_int_to_string(long int i) {
   char s[100];
   snprintf(s,99,"%ld",i);
   return std::string(s);
}

std::string
coot::util::float_to_string(float f) {

#if 0
   char s[100];
   snprintf(s,99,"%5.2f",f);
   return std::string(s);
#endif

   return float_to_string_using_dec_pl(f, 2);

}

std::string
coot::util::float_to_string_using_dec_pl(float f, unsigned short int n_dec_pl) {

#if 0 // previous
   char s[100];
   std::string prec="%7.";
   prec += coot::util::int_to_string(n_dec_pl);
   prec += "f";
   // snprintf(s,99,"%7.4f",f); // haha, FIXME. (use n_dec_pl, not 4)
   snprintf(s, 99, prec.c_str() ,f);
   return std::string(s);
#endif

   // a valgrind error here means that the f that you passed to this function was not
   // initialized.

   std::stringstream s;
   s << std::right << std::fixed;
   s << std::setprecision(n_dec_pl);
   s << f;
   std::string ss = s.str();
   return ss;

}

std::string
coot::util::float_to_unspaced_string_using_dec_pl(float f, unsigned short int n_dec_pl) {
   return float_to_string_using_dec_pl(f, n_dec_pl);
}

// throw an exception on unable to convert
int
coot::util::string_to_int(const std::string &s) {

   int i;
   std::istringstream myStream(s);

   if (myStream>>i) {
      return i;
   } else {
      std::string mess = "Cannot convert \"";
      mess += s;
      mess += "\" to an integer";
      throw std::runtime_error(mess);
   }
}

// throw an exception on unable to convert
float
coot::util::string_to_float(const std::string &s) {

   float f;
   std::istringstream myStream(s);

   if (myStream>>f) {
      return f;
   } else {
      std::string mess = "Cannot convert \"";
      mess += s;
      mess += "\" to a float";
      throw std::runtime_error(mess);
   }
}

// throw an exception on unable to convert
double
coot::util::string_to_double(const std::string &s) {

   double f;
   std::istringstream ss(s);

   if (ss>>f) {
      return f;
   } else {
      std::string mess = "Cannot convert \"";
      mess += s;
      mess += "\" to a double";
      throw std::runtime_error(mess);
   }
}


std::string
coot::util::plain_text_to_sequence(const std::string &s) {

   std::string r;

   for (unsigned int i=0; i<s.length(); i++) {
      // std::cout << "testing :" << s.substr(i, 1) << ":" << std::endl;
      if (coot::util::is_fasta_aa((s.substr(i, 1))))
	 r += toupper(s[i]);
   }

   return r;
}

std::string
coot::util::plain_text_to_pir(const std::string &title, const std::string &sequence, short int il) {

   std::string r = "> ";
   r += title;
   if (il == 2) {
     // for python
     r += "\\n";
     r += "\\n";
   } else {
     r += "\n";
     r += "\n";
   }
   r += sequence;
   r += "*";
   return r;
}



short int
coot::util::is_fasta_aa(const std::string &a_in) {

   short int r = 0;

   std::string a(upcase(a_in));

   if (a == "A" || a == "G" ) {
      r = 1;
   } else {
      if (a == "B"
	  || a == "C" || a == "D" || a == "E" || a == "F" || a == "H" || a == "I"
	  || a == "K" || a == "L" || a == "M" || a == "N" || a == "P" || a == "Q"
	  || a == "R" || a == "S" || a == "T" || a == "U" || a == "V" || a == "W"
	  || a == "Y" || a == "Z" || a == "X" || a == "*" || a == "-") {
	 r = 1;
      }
   }
   return r;
}

std::string
coot::util::single_quote(const std::string &s, const std::string &quote_char) {

   std::string r = quote_char;
   r += s;
   r += quote_char;
   return r;
}


coot::sequence::fasta::fasta(const std::string &seq_in) {

   std::string seq;

   std::cout << "debug:: coot::sequence::fasta::fasta seq_in: " << seq_in << std::endl;

   int nchars = seq_in.length();
   bool found_greater = false;
   bool found_newline = false;

   std::string t;

   for (int i=0; i<nchars; i++) {

      // std::cout << "checking character: " << seq_in[i] << std::endl;

      if (found_newline && found_greater) {
         t = toupper(seq_in[i]);
         if (is_fasta_aa(t)) {
            std::cout << "adding character: " << seq_in[i] << std::endl;
            seq += t;
         }
      }
      if (seq_in[i] == '>') {
         std::cout << "DEBUG:: " << seq_in[i] << " is > (greater than)\n";
         found_greater = true;
      }
      if (seq_in[i] == '\n') {
         if (found_greater) {
            std::cout << "DEBUG:: " << seq_in[i] << " is carriage return\n";
            found_newline = true;
         }
      }
   }

   if (seq.empty()) {
      std::cout << "WARNING:: coot::utils fasta::fasta() no sequence found or improper fasta sequence format\n";
   }
}

bool
coot::sequence::fasta::is_fasta_aa(const std::string &a) const {

   bool r = 0;

   if (a == "A" || a == "G" ) {
      r = 1;
   } else {
      if (a == "B"
          || a == "C" || a == "D" || a == "E" || a == "F" || a == "H" || a == "I"
          || a == "K" || a == "L" || a == "M" || a == "N" || a == "P" || a == "Q"
          || a == "R" || a == "S" || a == "T" || a == "U" || a == "V" || a == "W"
          || a == "Y" || a == "Z" || a == "X" || a == "*" || a == "-") {
         r = 1;
      }
   }
   return r;
}

// This should be a util function return the directory of this
// filename: /d/a -> /d/      "a" -> ""
std::string coot::util::file_name_directory(const std::string &file_name) {

   int end_char = -1;
   std::string rstring = "";

   if (file_name.length() == 0)
      return rstring;

   for (int i=file_name.length()-1; i>=0; i--) {
      // std::cout << file_name.substr(0, i) << std::endl;

      // BL says:: in windows we should check for \ too. Too much pain to get
      // everything converted to / for file_chooser, e.g. with debackslash!

      // Windows specific #ifdef removed 20081010, makes indenting lower done
      // this file work again.
      //
      if (file_name[i] == '/' || file_name[i] == '\\') {
	 if (i < int(file_name.length())) {
	    end_char = i;
	 } else {
	    // never get here?
	    std::cout << "cannont happen. end_char = " << end_char
		      << " file_name.length(): " << file_name.length() << std::endl;
	    end_char = i-1;
	 }
	 break;
      }
   }
   if (end_char != -1)
      rstring = file_name.substr(0, end_char+1);

   return rstring;
}

std::string
coot::util::current_working_dir() {
   std::string s = "";
   const unsigned long l = 2480;
   char b[l];
   char *x = getcwd(b,l);
   if (x)
      s = std::string(b);
   return s;
}

// If cwd is a substring of f (starting at 0), then return the
// basename of f (i.e. cwd stripped from f).  If cwd is not a
// substring of f, then return f;
//
std::string
coot::util::relativise_file_name(const std::string &f, const std::string &cwd) {

   std::string r = f;

   std::string::size_type pos = f.find(cwd);
   if (pos == 0) {
      // found it
      if (f.length() > cwd.length()) // sanity check
	 r = f.substr(cwd.length()+1);
   }
   return r;
}

// return absolute path for filename (can include dirs)
// oupon error return input filename and throw error.
//
// I don't know how to do this in Unix - maybe getcwd(), but that
// looks wrong.
//
std::string
coot::util::absolutise_file_name(const std::string &file_name) {

   std::string ret = file_name;
   // first check if already absolute file, i.e. starts with '/' or
   // has a ':' in second position (windows)

   if (file_name.empty()) return "";

   if (file_name.substr(0,1) != "/" &&
       file_name.substr(1,1) != ":" ) {
     std::string s = current_working_dir();
#ifdef WINDOWS_MINGW
     ret = intelligent_debackslash(s + "\\" + file_name);
#else
     ret = s + "/" + file_name;
#endif
   }

   return ret;
}

std::string
coot::util::name_sans_extension(const std::string &f) {

   std::string r = f;
   std::string::size_type iext = f.find_last_of(".");
   if (iext != std::string::npos)
      r = f.substr(0, iext);
   return r;
}



std::string
coot::util::file_name_non_directory(const std::string &file_name) {

   int slash_char = -1;
   std::string rstring = "";

   for (int i=file_name.length()-1; i>=0; i--) {
#ifdef WINDOWS_MINGW
      if (file_name[i] == '/' || file_name[i] == '\\') {
	 slash_char = i;
	 break;
      }
#else
      if (file_name[i] == '/') {
	 slash_char = i;
	 break;
      }
#endif // MINGW
   }

   if (slash_char != -1)
      rstring = file_name.substr(slash_char+1);
   else
      rstring = file_name;

   // std::cout << "DEBUG:: non-directory of " << file_name << " is " << rstring << std::endl;
   return rstring;
}

// return "" on no extension /usr/blogs/thign/other
std::string
coot::util::file_name_extension(const std::string &file_name) {

   int dot_char = -1;
   std::string rstring = "";

   for (int i=file_name.length()-1; i>=0; i--) {
      if (file_name[i] == '.') {
	 dot_char = i;
	 break;
      }
   }

   if (dot_char != -1)
      rstring = file_name.substr(dot_char);

   //    std::cout << "DEBUG:: extension of " << file_name << " is " << rstring << std::endl;
   return rstring;
}


bool
coot::util::extension_is_for_shelx_coords(const std::string &ext) {

   bool r = 0;
   if ((ext == ".INS") ||
       (ext == ".ins") ||
       (ext == ".RES") ||
       (ext == ".res") ||
       (ext == ".hat") ||
       (ext == ".HAT"))
      r = 1;

   return r;
}


bool
coot::util::extension_is_for_mdl_mol_or_mol2_coords(const std::string &ext) {

   bool r = 0;
   if ((ext == ".mdl")   ||
       (ext == ".MDL")   ||
       (ext == ".mol")   ||
       (ext == ".mol2")  || // tripos mol2
       (ext == ".mol3d") || // tripos mol2
       (ext == ".MOL2")  || // tripos mol2
       (ext == ".MOL")   ||
       (ext == ".sdf")   ||
       (ext == ".SDF"))
      r = 1;
   return r;
}


bool
coot::util::extension_is_for_coords(const std::string &ext) {

   bool r = false;
   if ((ext == ".pdb") ||
       (ext == ".pdb.gz") ||
       (ext == ".ent") ||
       (ext == ".ent.gz") ||
       (ext == ".PDB"))
      r = true;
   return r;
}


// mtz only I think?!
bool
coot::util::extension_is_for_auto_datasets(const std::string &ext) {

   bool r = false;
   if (ext == ".mtz")
      r = true;
   return r;
}

bool
coot::util::extension_is_for_scripts(const std::string &ext) {

   bool r = false;
   if ((ext == ".py") ||
       (ext == ".scm"))
      r = true;
   return r;
}


short int
coot::is_mmcif_filename(const std::string &filename) {

   short int i=0;

   std::string::size_type idot = filename.find_last_of(".");
   if (idot != std::string::npos) {
      std::string t = filename.substr(idot);

      std::string::size_type icif   = t.rfind(".cif");
      std::string::size_type immcif = t.rfind(".mmcif");
      std::string::size_type immCIF = t.rfind(".mmCIF");

      if ( (icif   != std::string::npos) ||
	   (immcif != std::string::npos) ||
	   (immCIF != std::string::npos) ) {
	 i = 1;
      }
   }
   return i;
}

// base, i.e. $HOME/coot-build or /usr
std::string
coot::prefix_dir() {

   std::string s;
   char *env = getenv("COOT_PREFIX");
   if (env) {
      s = env;
   } else {
      std::string dds = package_data_dir();
      if (! dds.empty())
         if (dds.back() == '/')
            dds.erase(dds.size() - 1);    // or pop_back

      // 20230607-PE Merge conflict
      //
      // maybe I should have kept Jakub's version:
      // 20230616-PE OK, let's just replace it then
      std::string jds = util::append_dir_dir(dds,"..");
      std::string js  = util::append_dir_dir(jds,"..");
      s = js;
   }
   return s;
}



std::string
coot::get_home_dir() {
   const char *s = getenv("HOME");
   if (s) {
      return std::string(s);
   } else {
      s = getenv("COOT_HOME");
      if (s)
         return std::string(s);
   }
   return ""; //empty
}

#include <filesystem>

// The user can set COOT_DATA_DIR (in fact this is the usual case
// when using binaries) and that should over-ride the built-in
// PKGDATADIR.
//
// Use this to find things in $prefix/share/coot
std::string
coot::package_data_dir() {

   // std::string xdatadir = XDATADIR; CMake
   std::string pkgdatadir = PKGDATADIR; // CMake does this too, it seems.
   // For binary installers, they use the environment variable:

   // std::cout << "debug:: in coot::package_data_dir() xdatadir: :" << xdatadir << ":" << std::endl;
   // std::cout << "debug:: in coot::package_data_dir() pkgdatadir: :" << pkgdatadir << ":" << std::endl;

   char *env = getenv("COOT_DATA_DIR");
   if (env) {
      pkgdatadir = std::string(env);
   } else {
      char *env = getenv("COOT_PREFIX");
      if (env)
         pkgdatadir = std::string(env) + std::string("/share/coot");
   }
   if (std::filesystem::exists(pkgdatadir)) {
      // good - we are (probably) not using relocated binaries
   } else {
      // let set pkgdatadir relative to the binary we are running (which was set
      // using set_realpath_for_coot_executable())
      // std::cout << "................. real_path_for_coot_executable " << real_path_for_coot_executable << std::endl;
      if (real_path_for_coot_executable.empty()) {
	 // there is no oops if this is chapi or moorhen - so removee the messaage
         // std::cout << "OOPS:: real_path_for_coot_executable is empty " << std::endl;
      } else {
         std::filesystem::path p(real_path_for_coot_executable);
         std::filesystem::path p_1 =   p.parent_path();
         std::filesystem::path p_2 = p_1.parent_path();
         // std::cout << "here with p_2 " << p_2.string() << std::endl;
         std::filesystem::path p_3 = p_2 / "share";
         std::filesystem::path p_4 = p_3 / "coot";
         pkgdatadir = p_4.string();
      }
   }
   return pkgdatadir;
}

// Use this to find things in $prefix/share/RDKit
std::string
coot::rdkit_package_data_dir() {

   std::string p = package_data_dir();
   std::string d = util::file_name_directory(p);
   std::string r = d + "/RDKit"; // maybe append_dir_file need to be moved to utils?
   return r;
}

// if you can try to get the directoy dir in this directory.
// if not, try to make it in this directory.
// if not, try to find it in $HOME
// if not try to make it in $HOME
// if not, return the empty string
std::string
coot::get_directory(const std::string &dir) {

   struct stat s;
   int fstat = stat(dir.c_str(), &s);
   if (fstat == -1 ) { // file not exist
      int status = util::create_directory(dir);
      if (status == 0) { // success
	 return dir;
      } else {
	 // try to create in $HOME
	 const char *e = getenv("HOME");
	 if (e) {
	    std::string home(e);
	    const std::string d = util::append_dir_dir(home, dir);
	    fstat = stat(d.c_str(), &s);
	    if (fstat == -1) {
	       int status = util::create_directory(d);
	       if (status == 0) { // fine
		  return d;
	       } else {
		  // couldn't create in $HOME either
		  std::string empty;
		  return empty;
	       }
	    } else {
	       return d;
	    }
	 } else {
	    // no $HOME
	    std::string empty;
	    return empty;
	 }
      }
   } else {
      return dir;
   }
}

// Let's use C++-17 filesystem rather than stat()
#include "xdg-base.hh"

//!  Use XDG Base Directory to get the download directory
std::string
coot::get_download_directory() {

   xdg_t xdg;
   std::filesystem::path p = xdg.get_cache_home();
   return p.string();
}




std::pair<std::string, std::string>
coot::util::split_string_on_last_slash(const std::string &string_in) {

   std::string::size_type islash = string_in.find_last_of("/");
   std::string first;
   std::string second("");

   if (islash != std::string::npos) {
      first  = string_in.substr(0, islash);
      second = string_in.substr(islash);
      if (second.length() > 0)
	 second = second.substr(1);
   } else {
      first = string_in;
   }
   return std::pair<std::string, std::string> (first, second);
}

std::vector<std::string>
coot::util::split_string(const std::string &string_in,
			 const std::string &splitter) {

   std::vector<std::string> v;
   std::string s=string_in;

   while (1) {
      std::string::size_type isplit=s.find_first_of(splitter);
      if (isplit != std::string::npos) {
	 std::string f = s.substr(0, isplit);
	 v.push_back(f);
	 if (s.length() > (isplit+splitter.length())) {
	    s = s.substr(isplit+splitter.length());
	 } else {
	    break;
	 }
      } else {
	 if (! s.empty()) {
	    v.push_back(s);
	    break;
	 }
      }
   }
   return v;
}

// by default splitter is " "
std::vector<std::string>
coot::util::split_string_no_blanks(const std::string &string_in,
                                   const std::string &splitter) {

   std::vector<std::string> v;
   std::string s=string_in;

   while (true) {
      std::string::size_type isplit=s.find_first_of(splitter);
      if (isplit != std::string::npos) {
         std::string f = s.substr(0, isplit);
         if (f.length() > 0) {
            v.push_back(f);
         }
         if (s.length() >= (isplit+splitter.length())) {
            s = s.substr(isplit+splitter.length());
         } else {
            break;
         }
      } else {
         if (! s.empty()) {
            v.push_back(s);
         }
         break;
      }
   }
   return v;
}

std::vector<std::string>
coot::util::split_string_on_whitespace_no_blanks(const std::string &string_in) {

   std::vector<std::string> v;
   std::string s=string_in;

   while (true) {
      std::string::size_type isplit_s=s.find_first_of(" ");
      std::string::size_type isplit_t=s.find_first_of("\t");
      if (isplit_s != std::string::npos) {
         std::string f = s.substr(0, isplit_s);
         if (f.length() > 0) {
            v.push_back(f);
         }
         if (s.length() >= (isplit_s+1)) {
            s = s.substr(isplit_s+1);
         } else {
            break;
         }
      } else {
         if (isplit_t != std::string::npos) {
            std::string f = s.substr(0, isplit_t);
            if (f.length() > 0) {
               v.push_back(f);
            }
            if (s.length() >= (isplit_t+1)) {
               s = s.substr(isplit_t+1);
            } else {
               break;
            }
         } else {
            if (! s.empty()) {
               v.push_back(s);
            }
            break;
         }
      }
   }

   if (false) {
      for(unsigned int i=0; i<v.size(); i++) {
         std::cout << v[i] << ":";
      }
      std::cout << std::endl;
   }
   return v;
}



std::string
coot::util::downcase(const std::string &s) {

   std::string r = s;
   std::string::iterator it=r.begin();

   while ( (it!=r.end()) ) {
      *it = ::tolower(*it);
      it++;
   }
   return r;
}

std::string
coot::util::upcase(const std::string &s) {

   std::string r = s;
   std::string::iterator it=r.begin();

   while (it!=r.end()) {
      *it = ::toupper(*it);
      it++;
   }
   return r;
}

//
std::string
coot::util::capitalise(const std::string &s) {

   if (s.length() == 0) {
      return "";
   } else {
      std::string rt = s.substr(0,1);
      rt += coot::util::downcase(s.substr(1));
      return rt;
   }
}



long int
coot::util::random() {

#if defined(WINDOWS_MINGW) || defined(_MSC_VER)
          long int r = rand();
          return r;
#else
          long int r = ::random();
          return r;
#endif
}


float
coot::util::random_f() {

   float r = static_cast<float>(random());
   float irm = 1.0/static_cast<float>(RAND_MAX);

   return r * irm;
}



bool
coot::file_exists(const std::string &filename) {

   struct stat s;
   int fstat = stat(filename.c_str(), &s);
   if ( fstat == -1 ) { // file not exist
      return false;
   } else {
      return true;
   }
}


bool
coot::file_is_empty(const std::string &filename) {

   struct stat s;
   int fstat = stat(filename.c_str(), &s);
   if ( fstat == -1 ) { // file not exist
      return false;
   } else {
      off_t ss = s.st_size;
      if (ss == 0)
         return true;
   }
   return false;
}

bool
coot::file_exists_and_non_empty(const std::string &file_name) {
   bool status = false;
   if (coot::file_exists(file_name)) {
      struct stat buf;
      int istat = stat(file_name.c_str(), &buf);
      if (istat == 0) { // success
         if (buf.st_size > 0) {
            status = true;
         }
      }
   }
   return status;
}

bool
coot::file_exists_and_non_tiny(const std::string &file_name, unsigned int tiny_size_max) {

   bool status = false;
   if (coot::file_exists(file_name)) {
      struct stat buf;
      int istat = stat(file_name.c_str(), &buf);
      if (istat == 0) { // success
         if (buf.st_size > tiny_size_max) {
            status = true;
         }
      }
   }
   return status;
}


bool coot::is_directory_p(const std::string &filename) {

   bool st = 0;
   struct stat s;
   int fstat = stat(filename.c_str(), &s);
   if ( fstat == -1 ) { // file not exist
      return 0;
   } else {
      if (S_ISDIR(s.st_mode)) {
	 return 1;
      } else {
	 return 0;
      }
   }
   return st;
}




// is ALA, GLY, TRP, MET, MSE (RNA DNA too)...?
bool
coot::util::is_standard_residue_name(const std::string &residue_name) {

   if (residue_name == "ALA")
      return 1;
   if (residue_name == "ARG")
      return 1;
   if (residue_name == "ASN")
      return 1;
   if (residue_name == "ASP")
      return 1;
   if (residue_name == "CYS")
      return 1;
   if (residue_name == "GLN")
      return 1;
   if (residue_name == "GLU")
      return 1;
   if (residue_name == "GLY")
      return 1;
   if (residue_name == "HIS")
      return 1;
   if (residue_name == "ILE")
      return 1;
   if (residue_name == "LEU")
      return 1;
   if (residue_name == "LYS")
      return 1;
   if (residue_name == "MET")
      return 1;
   if (residue_name == "MSE")
      return 1;
   if (residue_name == "PHE")
      return 1;
   if (residue_name == "PRO")
      return 1;
   if (residue_name == "SER")
      return 1;
   if (residue_name == "THR")
      return 1;
   if (residue_name == "TRP")
      return 1;
   if (residue_name == "TYR")
      return 1;
   if (residue_name == "VAL")
      return 1;
   if (residue_name == "G")
      return 1;
   if (residue_name == "A")
      return 1;
   if (residue_name == "T")
      return 1;
   if (residue_name == "C")
      return 1;
   if (residue_name == "U")
      return 1;
   if (residue_name == "DG")
      return 1;
   if (residue_name == "DA")
      return 1;
   if (residue_name == "DT")
      return 1;
   if (residue_name == "DC")
      return 1;
   if (residue_name == "DU")
      return 1;
   if (residue_name == "GR")
      return 1;
   if (residue_name == "AR")
      return 1;
   if (residue_name == "UR")
      return 1;
   if (residue_name == "TR")
      return 1;
   if (residue_name == "UR")
      return 1;
   if (residue_name == "Gr")
      return 1;
   if (residue_name == "Ar")
      return 1;
   if (residue_name == "Ur")
      return 1;
   if (residue_name == "Tr")
      return 1;
   if (residue_name == "Ur")
      return 1;
   if (residue_name == "Gd")
      return 1;
   if (residue_name == "Ad")
      return 1;
   if (residue_name == "Ud")
      return 1;
   if (residue_name == "Td")
      return 1;
   if (residue_name == "Ud")
      return 1;

   return 0;

}

bool
coot::util::is_standard_amino_acid_name(const std::string &residue_name) {

   if (residue_name == "ALA")
      return 1;
   if (residue_name == "ARG")
      return 1;
   if (residue_name == "ASN")
      return 1;
   if (residue_name == "ASP")
      return 1;
   if (residue_name == "CYS")
      return 1;
   if (residue_name == "GLN")
      return 1;
   if (residue_name == "GLU")
      return 1;
   if (residue_name == "GLY")
      return 1;
   if (residue_name == "HIS")
      return 1;
   if (residue_name == "ILE")
      return 1;
   if (residue_name == "LEU")
      return 1;
   if (residue_name == "LYS")
      return 1;
   if (residue_name == "MET")
      return 1;
   if (residue_name == "MSE")
      return 1;
   if (residue_name == "PHE")
      return 1;
   if (residue_name == "PRO")
      return 1;
   if (residue_name == "SER")
      return 1;
   if (residue_name == "THR")
      return 1;
   if (residue_name == "TRP")
      return 1;
   if (residue_name == "TYR")
      return 1;
   if (residue_name == "VAL")
      return 1;

   return 0;
}

bool
coot::util::is_standard_nucleotide_name(const std::string &residue_name) {

   if (residue_name == "G")
      return 1;
   if (residue_name == "A")
      return 1;
   if (residue_name == "T")
      return 1;
   if (residue_name == "C")
      return 1;
   if (residue_name == "U")
      return 1;
   if (residue_name == "DG")
      return 1;
   if (residue_name == "DA")
      return 1;
   if (residue_name == "DT")
      return 1;
   if (residue_name == "DC")
      return 1;
   if (residue_name == "DU")
      return 1;
   if (residue_name == "GR")
      return 1;
   if (residue_name == "AR")
      return 1;
   if (residue_name == "UR")
      return 1;
   if (residue_name == "TR")
      return 1;
   if (residue_name == "UR")
      return 1;
   if (residue_name == "Gr")
      return 1;
   if (residue_name == "Ar")
      return 1;
   if (residue_name == "Ur")
      return 1;
   if (residue_name == "Tr")
      return 1;
   if (residue_name == "Ur")
      return 1;
   if (residue_name == "Gd")
      return 1;
   if (residue_name == "Ad")
      return 1;
   if (residue_name == "Ud")
      return 1;
   if (residue_name == "Td")
      return 1;
   if (residue_name == "Ud")
      return 1;
   return 0;
}




bool
coot::sequence::is_sequence_triplet(const std::string &s) {

   bool r = false;

   if (s.length() == 3) {
      std::string ss = util::upcase(s);
      // A-Z check
      if (ss[0] < 91 && ss[0] >= 65)
	 if (ss[1] < 91 && ss[1] >= 65)
	    if (ss[2] < 91 && ss[2] >= 65)
	       r = true;
   }
   return r;
}


// return a set of string that match the glob, with the directory name pre-appended
std::vector<std::string>
coot::util::glob_files(const std::string &dir, const std::string &glob_pattern) {
#ifndef _MSC_VER
   std::vector<std::string> r;
   glob_t myglob;
   std::string glob_files = append_dir_file(dir, glob_pattern);
   int flags = 0;
   glob(glob_files.c_str(), flags, 0, &myglob);
   size_t count = myglob.gl_pathc;
   for (char **p = myglob.gl_pathv; count ; p++, count--) {
      char *file(*p);
      r.push_back(file);
   }
   globfree(&myglob);
   return r;
#else
   std::vector<std::string> r;
   WIN32_FIND_DATA ffd;
   HANDLE hFind;
   std::string glob_files = append_dir_file(dir, glob_pattern);
   hFind = FindFirstFile(glob_files.c_str(), &ffd);
   if (hFind != INVALID_HANDLE_VALUE) {
     do {
       r.push_back(ffd.cFileName);
     } while (FindNextFile(hFind, &ffd) != 0);
     FindClose(hFind);
   }
   return r;
#endif
}


coot::gauss_legendre_t::gauss_legendre_t() {
   fill_weight_abscicca(16);
}

void
coot::gauss_legendre_t::fill_weight_abscicca(int N) {

   weight_abscissa_.resize(16+1);

   weight_abscissa_[1]  = std::pair<double, double>(0.1894506104550685, -0.0950125098376374);
   weight_abscissa_[2]  = std::pair<double, double>(0.1894506104550685,  0.0950125098376374);
   weight_abscissa_[3]  = std::pair<double, double>(0.1826034150449236, -0.2816035507792589);
   weight_abscissa_[4]  = std::pair<double, double>(0.1826034150449236,  0.2816035507792589);
   weight_abscissa_[5]  = std::pair<double, double>(0.1691565193950025, -0.4580167776572274);
   weight_abscissa_[6]  = std::pair<double, double>(0.1691565193950025,  0.4580167776572274);
   weight_abscissa_[7]  = std::pair<double, double>(0.1495959888165767, -0.6178762444026438);
   weight_abscissa_[8]  = std::pair<double, double>(0.1495959888165767,  0.6178762444026438);
   weight_abscissa_[9]  = std::pair<double, double>(0.1246289712555339, -0.7554044083550030);
   weight_abscissa_[10] = std::pair<double, double>(0.1246289712555339,  0.7554044083550030);
   weight_abscissa_[11] = std::pair<double, double>(0.0951585116824928, -0.8656312023878318);
   weight_abscissa_[12] = std::pair<double, double>(0.0951585116824928,  0.8656312023878318);
   weight_abscissa_[13] = std::pair<double, double>(0.0622535239386479, -0.9445750230732326);
   weight_abscissa_[14] = std::pair<double, double>(0.0622535239386479,  0.9445750230732326);
   weight_abscissa_[15] = std::pair<double, double>(0.0271524594117541, -0.9894009349916499);
   weight_abscissa_[16] = std::pair<double, double>(0.0271524594117541,  0.9894009349916499);

}


// return empty string on failure
//
std::string
coot::suggest_new_comp_id(const std::string &comp_id_in) {

   std::string r;
   if (comp_id_in.length() != 3) {
      r = comp_id_in + "1";
   } else {
      char s[3];
      for (unsigned int i=0; i<3; i++)
	 s[i] = comp_id_in[i];
      std::vector<bool> is_numeral(3, false);
      std::vector<int> n(3, -1);
      for (unsigned int i=0; i<3; i++) {
	 if (s[i] >= 48 && s[i] <58) {
	    is_numeral[i] = true;
	    n[i] = s[i] - 48;
	 }
      }
      if (is_numeral[0] && is_numeral[1] && is_numeral[2]) {
	 if (n[2]<9) {
	    r =  s[0];
	    r += s[1];
	    r += s[2]+1;
	 } else {
	    if (n[1]<9) {
	       r =  s[0];
	       r += s[1]+1;
	       r += "0";
	    }
	 }
      }
      if (is_numeral[0] && is_numeral[1] && !is_numeral[2]) {
	    r =  s[0];
	    r += s[1];
	    r += "2";
      }
      if (is_numeral[0] && !is_numeral[1] && is_numeral[2]) {
	 r =  s[0];
	 r += s[1];
	 r += s[2]+1;
      }

      if (is_numeral[0] && !is_numeral[1] && !is_numeral[2]) {
	 r =  s[0];
	 r += s[1];
	 r += "2";
      }
      if (!is_numeral[0] && is_numeral[1] && is_numeral[2]) {
	 if (n[2]<9) {
	    r =  s[0];
	    r += s[1];
	    r += s[2]+1;
	 } else {
	    if (n[1]<9) {
	       r =  s[0];
	       r += s[1]+1;
	       r += "0";
	    }
	 }
      }
      if (!is_numeral[0] && is_numeral[1] && !is_numeral[2]) {
	 r =  s[0];
	 r += s[1];
	 r += "2";
      }
      if (!is_numeral[0] && !is_numeral[1] && is_numeral[2]) {
	    r =  s[0];
	    r += s[1];
	    r += s[2]+1;
      }
      if (!is_numeral[0] && !is_numeral[1] && !is_numeral[2]) {
	 r =  s[0];
	 r += s[1];
	 r += "2";
      }
   }
   return r;
}


// can throw a std::runtime_error exception.  If this returns, it guarantees a useful result.
//
// This doesn't deal with large numbers or negative numbers.
//
std::pair<std::string, long>
coot::util::extract_number_string(const std::string &s) {

   std::pair<std::string, long> r("", 0);

   for (std::string::size_type i=0; i<s.size(); i++) {
      if (is_number(s[i])) {
	 r.first += s[i];
      } else {
	 break;
      }
   }
   if (r.first.length() == 0) {
      std::string mess("No number");
      throw(std::runtime_error(mess));
   } else {
      // string_to_long?
      r.second = string_to_int(r.first);
   }
   return r;
}

int
coot::util::round_up_by_hundreds(int num) {

   // 123 should return 200
   // 100 should return 100

   float a = static_cast<float>(num+99) * 0.01;
   float f = floorf(a);
   int ii = static_cast<int>(f) * 100;

   return ii;
}

#include <fstream>

// return true for success
bool
coot::copy_file(const std::string &from_file, const std::string &to_file) {

   auto file_to_string = [] (const std::string &file_name) {
                            std::string s;
                            std::string line;
                            std::ifstream f(file_name.c_str());
                            if (!f) {
                               std::cout << "ERROR:: copy_file() failed to open " << file_name << std::endl;
                            } else {
                               while (std::getline(f, line)) { 
                                  s += line;
                                  s += "\n"; // may be different for windows.
                               }
                            }
                            return s;
                         };

   auto string_to_file = [] (const std::string &s, const std::string &to_file) {
                            bool success = false;
                            std::ofstream f(to_file);
                            if (f) {
                               f << s;
                               f.close();
                               success = true;
                            }
                            return success;
                   };

   bool success = false;
   if (file_exists(from_file)) {
      std::string s = file_to_string(from_file);
      success = string_to_file(s, to_file);
   }

   return success;
}
