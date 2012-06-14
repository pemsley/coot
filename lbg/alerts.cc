
#ifdef MAKE_ENTERPRISE_TOOLS

#include "lbg.hh"

// return a vector of alert matches
// 
std::vector<lbg_info_t::alert_info_t>
lbg_info_t::alerts(const RDKit::ROMol &mol) const {

   std::vector<lbg_info_t::alert_info_t> v;

   std::vector<std::pair<std::string, std::string> > patterns = alert_smarts();
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

#endif // MAKE_ENTERPRISE_TOOLS
