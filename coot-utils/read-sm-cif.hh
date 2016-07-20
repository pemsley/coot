
#include "clipper/cif/cif_data_io.h"
#include "clipper/core/xmap.h"

namespace coot {

   class smcif {
      clipper::Cell get_cell(mmdb::mmcif::PData data) const;
      clipper::Cell get_cell_for_data(mmdb::mmcif::PData data) const;
      std::pair<bool,clipper::Spacegroup> get_space_group(const std::vector<std::string> &symm_strings) const;
      std::vector<mmdb::Atom *> read_coordinates(mmdb::mmcif::PData data, const clipper::Cell &cell, const clipper::Spacegroup &spg) const;
      std::pair<bool,clipper::Spacegroup> get_space_group(mmdb::mmcif::Data *data) const;
      std::pair<bool,clipper::Spacegroup> get_space_group(mmdb::mmcif::Data *data, const std::string &symm_tag) const;
      std::pair<bool,clipper::Spacegroup> get_space_group_from_loop(mmdb::mmcif::Data *data, const std::string &symm_tag) const;


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
      clipper::HKL_data<clipper::datatypes::F_sigF<float> > my_fsigf;
      // and this (from the real and imaginary components)
      clipper::HKL_data<clipper::datatypes::F_phi<float> >  my_fphi;
      
      // c.f. get_cell() from a coords file
      clipper::Cell get_cell_for_data(const std::string &file_name) const;
      
      std::pair<bool,clipper::Spacegroup> get_space_group(const std::string &file_name) const;

      clipper::Resolution get_resolution(const clipper::Cell &cell,
					 const std::string &file_name) const;
      void setup_hkls(const std::string &file_name);

      // various ways in which the symmetry can be specified 
      //
      

   public:
      smcif() {};
      smcif(const std::string &file_name) {
	 read_sm_cif(file_name);
      }
      mmdb::Manager *read_sm_cif(const std::string &file_name) const;
      // return success status, true is good
      bool read_data_sm_cif(const std::string &file_name);
      // return an empty map if not possible
      clipper::Xmap<float> map() const;
      bool check_for_f_phis() const; // use sigmaa_maps if we have phis, use
                                     // sigmaa_maps_by_calc_sfs if we don't.
      // calculate maps using fcalc and phi calc in the .cif files
      // return an empty map in first if not possible, 
      std::pair<clipper::Xmap<float>, clipper::Xmap<float> > sigmaa_maps();
      // return an empty map in first if not possible
      std::pair<clipper::Xmap<float>, clipper::Xmap<float> > sigmaa_maps_by_calc_sfs(mmdb::Atom **atom_selection, int n_selected_atoms);
   };

   class simple_sm_u {
   public:
      std::string label; // atom name
      mmdb::realtype u11, u22, u33, u12, u13, u23;
      simple_sm_u() {
	 u11 = 0;
	 u22 = 0;
	 u33 = 0;
	 u12 = 0;
	 u13 = 0;
	 u23 = 0;
      }
      simple_sm_u(const std::string label_in,
		  mmdb::realtype u11_in, mmdb::realtype u22_in, mmdb::realtype u33_in,
		  mmdb::realtype u12_in, mmdb::realtype u13_in, mmdb::realtype u23_in) {
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
