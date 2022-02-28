
#ifndef UPDATING_COORDINATES_MOLECULE_PARAMETERS_T
#define UPDATING_COORDINATES_MOLECULE_PARAMETERS_T

#include<string>
#include <sys/types.h>// stat
#include <sys/stat.h>

// This class is for reading the output of Refmac
//
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
   explicit updating_coordinates_molecule_parameters_t(const std::string &file_name) : pdb_file_name(file_name) {
      ctime.tv_sec = 0;
      ctime.tv_nsec = 0;
      imol = -1; // is this used?
   }
   void update_from_stat_info(const struct stat &s) {

#ifdef WINDOWS_MINGW
      ctime.tv_sec = s.st_ctime;
      ctime.tv_nsec = 0.; // not available!? Lets hope not necessary
#else
#ifndef _POSIX_SOURCE
      ctime = s.st_ctimespec; // Mac OS X?
#else
      ctime = s.st_ctim;
#endif
#endif
   }
};

// This class is for updating the difference map when the model has changed.
//
class updating_model_molecule_parameters_t {
public:
   int imol_coords;
   int imol_map_with_data_attached;
   int imol_2fofc_map;  // sigmaA weighted, that is
   int imol_fofc_map; // ditto
   updating_model_molecule_parameters_t() {
      imol_coords = -1;
      imol_map_with_data_attached = -1;
      imol_2fofc_map = -1;
      imol_fofc_map = -1;
   }
   updating_model_molecule_parameters_t(int imol_coords_in, int imol_d, int imol_map_in) :
      imol_coords(imol_coords_in), imol_map_with_data_attached(imol_d), imol_fofc_map(imol_map_in) { imol_2fofc_map = -1;
   }
   updating_model_molecule_parameters_t(int imol_coords_in, int imol_data, int imol_map_2fofc_in, int imol_map_fofc_in) :
      imol_coords(imol_coords_in), imol_map_with_data_attached(imol_data), imol_2fofc_map(imol_map_2fofc_in), imol_fofc_map(imol_map_fofc_in) {}
};

#endif
