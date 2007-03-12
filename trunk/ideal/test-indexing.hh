
#include <string>
#include <vector>
#include <map>


namespace coot {
   
   class testclass {
   public:
      
      // std::map <std::vector<std::map <std::string, int> > > big_index;

      std::vector<std::map <std::string, int> > atom_name_resno_to_index;
      void add_residue_atom_map(int iresno, const std::map<std::string, int> &atom_map) {
	 if (iresno > atom_name_resno_to_index.size() ) {
	    atom_name_resno_to_index.resize(iresno + 1);
	 }
	 atom_name_resno_to_index[iresno] = atom_map;
      }
      void add_atom(int iresno, const std::string &at_name, int atom_index) {
	 if (iresno > atom_name_resno_to_index.size() ) {
	    atom_name_resno_to_index.resize(iresno + 1);
	 }
	 atom_name_resno_to_index[iresno][at_name] = atom_index;
      }

      void set_big_index(const std::string &chain,
			 int iresno,
			 const std::string &at_name, int atom_index) {
	 
      }
   };

}
