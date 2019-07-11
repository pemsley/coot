
#ifndef UPDATING_COORDINATES_MOLECULE_PARAMETERS_T
#define UPDATING_COORDINATES_MOLECULE_PARAMETERS_T

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

#endif
