
#ifndef UPDATING_COORDINATES_MOLECULE_PARAMETERS_T
#define UPDATING_COORDINATES_MOLECULE_PARAMETERS_T

// This is for reading the output of Refmac

class updating_coordinates_molecule_parameters_t {
public:
   int imol;
   std::string pdb_file_name;
   timespec ctime;
   updating_coordinates_molecule_parameters_t() {
      imol = -1;
      ctime.tv_sec = 0;
      ctime.tv_nsec = 0;
   }
   updating_coordinates_molecule_parameters_t(const std::string &file_name) {
      pdb_file_name = file_name;
      ctime.tv_sec = 0;
      ctime.tv_nsec = 0;
   }
};

// This is for updating the difference map when the model has changed.
class updating_model_molecule_parameters_t {
public:
   int imol_coords;
   int imol_map_with_data_attached;
   int imol_map;
   updating_model_molecule_parameters_t() {
      imol_coords = -1;
      imol_map_with_data_attached = -1;
      imol_map = -1;
   }
   updating_model_molecule_parameters_t(int imol_coords_in, int imol_d, int imol_map_in) :
      imol_coords(imol_coords_in), imol_map_with_data_attached(imol_d), imol_map(imol_map_in) {}
};

#endif
