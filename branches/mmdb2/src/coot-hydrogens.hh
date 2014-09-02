
namespace coot {
   // generic interface - using CCP4 tools.  Maybe this should be
   // declared in src somewhere then?
   // 
   std::pair<bool, std::string>
   add_hydrogens(CResidue *residue_p,
		 const dictionary_residue_restraints_t &restraints);


}
