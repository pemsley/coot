
#include "zo-rama.hh"

// this can throw a runtime_error.
//
void
zo::rama_table::read(const std::string &file_name) {
   std::ifstream f(file_name.c_str());
   std::string line;
   // std::cout << "INFO:: reading file " << file_name << std::endl;

   if (f) {
      while(std::getline(f, line)) {
	 std::vector<std::string> bits = coot::util::split_string_no_blanks(line);
	 if (bits.size() == 7) {
	    std::cout << "line: " << line << std::endl;
	    int idx_1 = coot::util::string_to_int(bits[0]);
	    int idx_2 = coot::util::string_to_int(bits[1]);
	    double A_cc = coot::util::string_to_double(bits[3]);
	    double A_cs = coot::util::string_to_double(bits[4]);
	    double A_sc = coot::util::string_to_double(bits[5]);
	    double A_ss = coot::util::string_to_double(bits[6]);
	    zo::rama_coeffs v(idx_1, idx_2, A_cc, A_cs, A_sc, A_ss);
	    rama_vec.push_back(v);
	 }
      }
   } else {
      std::cout << "Warning:: file not found: " << file_name << std::endl;
      throw std::runtime_error("Can't init zo-rama");
   }
}


void
zo::rama_table::test_analytical_derivs() const {

   int width  = 10;
   int height = 10;
   for (int j=0; j<width; j++) {
      for (int i=0; i<height; i++) {
	 zo::realtype phi = M_PI *  2.0 * zo::realtype(i-0.5*height)/zo::realtype(height);  // -pi to +pi
	 zo::realtype psi = M_PI * -2.0 * zo::realtype(j-0.5*width )/zo::realtype(width);
	 zo::realtype sum = 0;
	 std::pair<zo::realtype, zo::realtype> df_n = df_numerical(phi, psi);
	 std::pair<zo::realtype, zo::realtype> df_a = df          (phi, psi);
      }
   }
}

void
zo::rama_table::make_a_png(int width, const std::string &file_name) const {

   
   int height = width; // this is so for a ramachandran plot

   png_bytep *row_pointers = static_cast<png_bytep *>(malloc(sizeof(png_bytep) * height));
   for (int i=0; i<height; i++)
      row_pointers[i] = static_cast<png_byte*>(malloc(width));
   std::vector<std::vector<zo::realtype> > v(height);
   for (int i=0; i<height; i++)
      v[i] = std::vector<zo::realtype>(width, 0);

   for (int j=0; j<width; j++) {
      for (int i=0; i<height; i++) {
	 zo::realtype phi = M_PI *  2.0 * double(i-0.5*height)/double(height);
	 zo::realtype psi = M_PI * -2.0 * double(j-0.5*width )/double(width); // -pi to +pi
	 zo::realtype sum = value(phi,psi);
	 zo::realtype d = expf(sum);
	 v[j][i] = d;
      }
   }

   // replace by gradient for testing
   if (false) {
      for (int j=0; j<width; j++) {
	 for (int i=0; i<height; i++) {
	    zo::realtype phi = M_PI *  2.0 * double(i-0.5*height)/double(height);
	    zo::realtype psi = M_PI * -2.0 * double(j-0.5*width )/double(width); // -pi to +pi
	    std::pair<zo::realtype, zo::realtype> grads = df_numerical(phi,psi);
	    // zo::realtype d = expf(sqrt(grads.first*grads.first + grads.second*grads.second));
	    zo::realtype d = sqrt(grads.first*grads.first + grads.second*grads.second);
	    v[j][i] = d;
	 }
      }
   }

   zo::realtype min_v = 9.99e12;
   zo::realtype max_v = 0.0;
   zo::realtype sum = 0;
   int n = 0;
   for (int j=0; j<width; j++) {
      for (int i=0; i<height; i++) {
	 sum += v[j][i];
	 n++;
	 if (v[j][i] < min_v) min_v = v[j][i];
	 if (v[j][i] > max_v) max_v = v[j][i];
      }
   }
   zo::realtype mean = sum/zo::realtype(n);
   zo::realtype sf = 0.1 / mean;

   for (int j=0; j<width; j++) {
      for (int i=0; i<height; i++) {
	 zo::realtype vv = v[j][i] * sf;
	 vv = 255 * (1 - vv);
	 if (vv < 0) vv = 0;
	 // int pixel_value = std::lround(vv);
	 int pixel_value = int(vv+0.5);
	 row_pointers[j][i] = pixel_value;
      }
   }

   write_png_file(width, height, row_pointers, file_name.c_str());

   // now clean up
   for (int y=0; y<height; y++)
      free(row_pointers[y]);
   free(row_pointers);

}


// fill the table_map with tables
//
void
zo::rama_table_set::init() {

   // Old style - only one table in a file
//    try {
//       table_map[ALL_NON_PRE_PRO] = rama_table("all-non-pre-pro.tab");
//    }
//    catch (const std::runtime_error &rte) {
//       std::cout << "ERROR:: " << rte.what() << std::endl;
//    }

   // read interleaved coeffsolution.dat

   std::string dir1 = coot::package_data_dir();
   std::string dir2 = coot::util::append_dir_dir(dir1, "data");
   std::string dir3 = coot::util::append_dir_dir(dir2, "rama");
   std::string dir4 = coot::util::append_dir_dir(dir3, "zo-tables");

   std::string full = coot::util::append_dir_file(dir4, "coeffsolution.dat");

   std::ifstream f(full.c_str());
   std::string line;

   // std::cout << "INFO:: reading file " << full << std::endl;

   std::vector<std::string> table_type(16);
   table_type[ 0] = "ALL!P";  // all not next Pro
   table_type[ 1] = "ALLP";   // all with next is Pro
   table_type[ 2] = "GLY!P";
   table_type[ 3] = "GLYP";
   table_type[ 4] = "PRO!P";
   table_type[ 5] = "PROP";
   table_type[ 6] = "VI!P";
   table_type[ 7] = "VIP";
   table_type[ 8] = "DN!P";
   table_type[ 9] = "DNP";
   table_type[10] = "ST!P";
   table_type[11] = "STP";
   table_type[12] = "EQ!P";
   table_type[13] = "EQP";
   table_type[14] = "LA!P";
   table_type[15] = "LAP";

   if (f) {
      while(std::getline(f, line)) {
	 std::vector<std::string> bits = coot::util::split_string_no_blanks(line);
	 if (bits.size() == 7) {

	    if (bits[0] != "HEADER") {

	       // std::cout << "line: " << line << std::endl;
	       int idx_1 = coot::util::string_to_int(bits[0]);
	       int idx_2 = coot::util::string_to_int(bits[1]);
	       int idx_res_type = coot::util::string_to_int(bits[2]);
	       double A_cc = coot::util::string_to_double(bits[3]);
	       double A_cs = coot::util::string_to_double(bits[4]);
	       double A_sc = coot::util::string_to_double(bits[5]);
	       double A_ss = coot::util::string_to_double(bits[6]);
	       zo::rama_coeffs v(idx_1, idx_2, A_cc, A_cs, A_sc, A_ss);
	       if (idx_res_type >=0 && idx_res_type <= 15) {
		  std::string tt = table_type[idx_res_type];
		  table_map[tt].rama_vec.push_back(v);
	       }
	    }
	 }
      }
   }
}

std::pair<zo::realtype, zo::realtype>
zo::rama_table_set::df(const std::string &residue_type,
		       const zo::realtype &phi, const zo::realtype &psi) const {
   // residue type is e.g. "ALL!nP"
   std::map<std::string, rama_table>::const_iterator it = table_map.find(residue_type);
   return it->second.df(phi,psi);
}


zo::realtype
zo::rama_table_set::value(const std::string &residue_type,
			  const zo::realtype &phi, const zo::realtype &psi) const {

   // residue type is e.g. "ALL!nP"

   // std::cout << "debug:: in rama_table_set::value residue_type was " << residue_type << std::endl;
   std::map<std::string, rama_table>::const_iterator it = table_map.find(residue_type);
   if (it != table_map.end()) {
      return it->second.value(phi,psi);
   } else {
      std::cout << "ERROR:: unknown residue/table type \"" << residue_type << "\"" << std::endl;
      return 0.0;
   }
}

std::string
zo::rama_table_set::get_residue_type(const std::string &this_residue_type,
				     const std::string &next_residue_type) const {

   std::string r;
   if (next_residue_type == "PRO") {
      r = "ALLP";
      if (this_residue_type == "GLY") r = "GLYP";
      if (this_residue_type == "PRO") r = "PROP";
      if (this_residue_type == "VAL") r = "VIP";
      if (this_residue_type == "ILE") r = "VIP";
      if (this_residue_type == "ASP") r = "DNP";
      if (this_residue_type == "ASN") r = "DNP";
      if (this_residue_type == "SER") r = "STP";
      if (this_residue_type == "THR") r = "STP";
      if (this_residue_type == "GLU") r = "EQP";
      if (this_residue_type == "GLN") r = "EQP";
      if (this_residue_type == "LEU") r = "LAP";
      if (this_residue_type == "ALA") r = "LAP";
   } else {
      r = "ALL!P";
      if (this_residue_type == "GLY") r = "GLY!P";
      if (this_residue_type == "PRO") r = "PRO!P";
      if (this_residue_type == "VAL") r = "VI!P";
      if (this_residue_type == "ILE") r = "VI!P";
      if (this_residue_type == "ASP") r = "DN!P";
      if (this_residue_type == "ASN") r = "DN!P";
      if (this_residue_type == "SER") r = "ST!P";
      if (this_residue_type == "THR") r = "ST!P";
      if (this_residue_type == "GLU") r = "EQ!P";
      if (this_residue_type == "GLN") r = "EQ!P";
      if (this_residue_type == "LEU") r = "LA!P";
      if (this_residue_type == "ALA") r = "LA!P";
   }
   return r;
}
