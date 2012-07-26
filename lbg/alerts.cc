
#ifdef MAKE_ENTERPRISE_TOOLS

#ifdef USE_PYTHON
#include <Python.h>
#endif

#include "lbg.hh"


void
lbg_info_t::setup_user_defined_alert_smarts() {

   PyObject *m = PyImport_AddModule("__main__");
   user_defined_alerts_smarts_py = PyObject_GetAttrString(m,"user_defined_alert_smarts");

}


// return a vector of alert matches
// 
std::vector<lbg_info_t::alert_info_t>
lbg_info_t::alerts(const RDKit::ROMol &mol) const {

   std::vector<lbg_info_t::alert_info_t> v;

   std::vector<std::pair<std::string, std::string> > patterns = alert_smarts();
   std::vector<std::pair<std::string, std::string> > uda = user_defined_alert_smarts();
   for (unsigned int i=0; i<uda.size(); i++)
      patterns.push_back(uda[i]);

   
   for (unsigned int ipat=0; ipat<patterns.size(); ipat++) { 
      // std::cout << "checking pattern " << patterns[ipat] << std::endl;
      RDKit::ROMol *query = RDKit::SmartsToMol(patterns[ipat].first);
      std::vector<RDKit::MatchVectType>  matches;
      bool recursionPossible = true;
      bool useChirality = true;
      bool uniquify = true;
      int matched = RDKit::SubstructMatch(mol, *query, matches, uniquify,
					  recursionPossible, useChirality);
      if (matched) { 
	 // std::cout << "...... ALERT!" << std::endl;
	 for (unsigned int im=0; im<matches.size(); im++) {
	    alert_info_t alert(patterns[ipat].first,
			       patterns[ipat].second,
			       matches[im]);
	    v.push_back(alert);
	 }
      }
      delete query;
   }
   return v;
}

std::vector<std::pair<std::string, std::string> >
lbg_info_t::user_defined_alert_smarts() const {

   std::vector<std::pair<std::string, std::string> > v;

#ifdef USE_PYTHON

   if (user_defined_alerts_smarts_py) { 
      if (PyList_Check(user_defined_alerts_smarts_py)) {
	 int len = PyList_Size(user_defined_alerts_smarts_py);
	 for (unsigned int i=0; i<len; i++) { 
	    PyObject *item_pair = PyList_GetItem(user_defined_alerts_smarts_py, i);
	    if (PyList_Check(item_pair)) {
	       int item_len = PyList_Size(item_pair);
	       if (item_len == 2) {
		  PyObject *alert = PyList_GetItem(item_pair, 0);
		  PyObject *description = PyList_GetItem(item_pair, 1);
		  if (PyString_Check(alert)) { 
		     if (PyString_Check(description)) {
			std::string a = PyString_AsString(alert);
			std::string d = PyString_AsString(description);
			std::pair<std::string, std::string> p(a,d);
			v.push_back(p);
		     }
		  }
	       }
	    }
	 }
      }
   } 
#endif
   return v;
} 


#endif // MAKE_ENTERPRISE_TOOLS
