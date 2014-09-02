
#include "clipper/cif/cif_data_io.h"
#include "clipper/core/xmap.h"

namespace coot {

   class smcif {
      clipper::Cell get_cell(PCMMCIFData data) const;
      clipper::Cell get_cell_for_data(PCMMCIFData data) const;
      std::pair<bool,clipper::Spacegroup> get_space_group(const std::vector<std::string> &symm_strings) const;
      std::vector<CAtom *> read_coordinates(PCMMCIFData data, const clipper::Cell &cell, const clipper::Spacegroup &spg) const;
      std::pair<bool,clipper::Spacegroup> get_space_group(PCMMCIFData data) const;
      std::pair<bool,clipper::Spacegroup> get_space_group(PCMMCIFData data, const std::string &symm_tag) const;


      // e.g. "O"    -> " O"
      //      "V5+"  -> " V"
      //      "Zn2+" -> "ZN"
      //  return the oxidation state also in second - if possible. 0 if not.
      std::pair<std::string, int> symbol_to_element(const std::string &symbol) const;


      clipper::HKL_info mydata;
      clipper::Cell data_cell;
      clipper::Spacegroup data_spacegroup;
      clipper::Resolution data_resolution;

      // fill this
      clipper::HKL_data<clipper::datatypes::F_sigF<float> > myfsigf;
      // and this (from the real and imaginary components)
      clipper::HKL_data<clipper::datatypes::F_phi<float> >  my_fphi;
      
      // c.f. get_cell() from a coords file
      clipper::Cell get_cell_for_data(const std::string &file_name) const;
      
      std::pair<bool,clipper::Spacegroup> get_space_group(const std::string &file_name) const;

      clipper::Resolution get_resolution(const clipper::Cell &cell,
					 const std::string &file_name) const;
      void setup_hkls(const std::string &file_name);


   public:
      smcif() {};
      smcif(const std::string &file_name) {
	 read_sm_cif(file_name);
      }
      CMMDBManager *read_sm_cif(const std::string &file_name) const;
      // return success status, true is good
      bool read_data_sm_cif(const std::string &file_name);
      // return an empty map if not possible
      clipper::Xmap<float> map() const;
      // return an empty map in first if not possible
      std::pair<clipper::Xmap<float>, clipper::Xmap<float> > sigmaa_maps() const;

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
