
namespace coot {

   class smcif {
      clipper::Cell get_cell(PCMMCIFData data) const;
      std::pair<bool,clipper::Spacegroup> get_space_group(const std::vector<std::string> &symm_strings) const;
      std::vector<CAtom *> read_coordinates(PCMMCIFData data, const clipper::Cell &cell, const clipper::Spacegroup &spg) const;
      // e.g. "O"    -> " O"
      //      "V5+"  -> " V"
      //      "Zn2+" -> "ZN"
      //  return the oxidation state also in second - if possible. 0 if not.
      std::pair<std::string, int> symbol_to_element(const std::string &symbol) const;
   public:
      smcif() {};
      CMMDBManager *read_sm_cif(const std::string &file_name) const;

   };

   class simple_sm_u {
   public:
      std::string label; // atom name
      realtype u11, u22, u33, u12, u13, u23;
      simple_sm_u() {
	 u11 = 0;
	 u22 = 0;
	 u33 = 0;
	 u12 = 0;
	 u13 = 0;
	 u23 = 0;
      }
      simple_sm_u(const std::string label_in,
		  realtype u11_in, realtype u22_in, realtype u33_in,
		  realtype u12_in, realtype u13_in, realtype u23_in) {
	 label = label_in;
	 u11 = u11_in;
	 u22 = u22_in;
	 u33 = u33_in;
	 u12 = u12_in;
	 u13 = u13_in;
	 u23 = u23_in;
      } 
   };
   
}
