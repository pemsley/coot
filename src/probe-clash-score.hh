
#ifndef PROBE_CLASH_SCORE_HH
#define PROBE_CLASH_SCORE_HH

#include "coot-utils/coot-coord-utils.hh"

namespace coot {

   // " B  72 CYS  HG A: B  53 HIS  H   "
   // 
   class probe_atom_spec_t : public atom_spec_t {
   public: 
      probe_atom_spec_t(const std::string &s) : atom_spec_t() {
	 if (s.length() > 14) { 
	    std::string chain_local = s.substr(0,2);
	    std::string res_no_str = s.substr(2, 4);
	    std::string atom_name_local = s.substr(11, 4);
	    try {
	       int resno_local = coot::util::string_to_int(res_no_str);
	       if (chain_local[0] == ' ') { 
		  if (chain_local.length() > 1)
		     chain = std::string(chain_local.substr(1));
	       } else { 
		  chain = chain_local;
	       }
	       resno = resno_local;
	       atom_name = atom_name_local;
	    }
	    catch (const std::exception &e) {
	       std::cout << "WARNING:: " << e.what() << std::endl;
	    }
	 }
      }
      probe_atom_spec_t() : atom_spec_t() {}
   };

   class one_way_probe_contact_t {
   public:
      probe_atom_spec_t from_atom;
      std::vector<probe_atom_spec_t> to_atoms;
      one_way_probe_contact_t(const probe_atom_spec_t &spec) {
	 from_atom = spec;
      }
      unsigned int size() const { return to_atoms.size(); }
      void add(const probe_atom_spec_t &spec) {
	 std::vector<probe_atom_spec_t>::const_iterator it;
	 it = std::find(to_atoms.begin(), to_atoms.end(), spec);
	 if (it == to_atoms.end()) {
	    to_atoms.push_back(spec);
	 }
      }
   };

   class one_way_probe_contact_container_t {
   public:
      std::vector<one_way_probe_contact_t> contacts;
      void add(const probe_atom_spec_t &from_atom,
	       const probe_atom_spec_t &to_atom) {
	 bool found = false;
	 for (unsigned int i=0; i<contacts.size(); i++) { 
	    if (contacts[i].from_atom == from_atom) {
	       contacts[i].add(to_atom);
	       found = true;
	    }
	 }
	 if (! found) {
	    one_way_probe_contact_t new_contact(from_atom);
	    new_contact.add(to_atom);
	    contacts.push_back(new_contact);
	 }
      }
      unsigned int size() const {
	 unsigned int s = 0;
	 for (unsigned int i=0; i<contacts.size(); i++)
	    s += contacts[i].size();
	 return s;
      }
   };

   class probe_clash_score_t {
   public:
      bool filled;
      int n_bad_overlaps;
      int n_hydrogen_bonds;
      int n_small_overlaps;
      int n_close_contacts;
      int n_wide_contacts;
      probe_clash_score_t() {
	 filled = false;
      }
      probe_clash_score_t(const std::string &dots_file_name);
   }; 

   
   // couldn't get this to work - so I did by another method.
   // deleteable cruft.
   class spec_eraser {
   public:
      std::map<std::pair<probe_atom_spec_t, probe_atom_spec_t>, bool> ref_specs;
      spec_eraser(const std::map<std::pair<probe_atom_spec_t, probe_atom_spec_t>, bool> &ref_specs_in) {
	 ref_specs = ref_specs_in;
      }
      bool operator() (const std::pair<probe_atom_spec_t, probe_atom_spec_t> &s) const {
	 return true;
      } 
   };
}


#ifdef USE_GUILE
//! \brief return scheme false on failure or a scheme list
// (n_bad_overlaps n_hydrogen_bonds n_small_overlaps n_close_contacts
// n_wide_contacts)
// 
SCM probe_clash_score_scm(const std::string &dots_file_name);
#endif // USE_GUILE

#ifdef USE_PYTHON
//! \brief return scheme false on failure or a scheme list
// (n_bad_overlaps n_hydrogen_bonds n_small_overlaps n_close_contacts
// n_wide_contacts)
// 
PyObject *probe_clash_score_py(const std::string &dots_file_name);
#endif // USE_PYTHON


#endif // PROBE_CLASH_SCORE_HH
