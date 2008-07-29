
#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>

#include "mmdb-extras.h"
#include "mmdb.h"

#define BUILT_IN_TESTING

std::string stringify(double x);
std::string stringify(int i);
std::string stringify(unsigned int i);
std::string greg_test(const std::string &file_name);

int test_internal();

int test_alt_conf_rotamers();
int test_wiggly_ligands();
int test_torsion_derivs();
int test_ramachandran_probabilities();
int kdc_torsion_test();

class residue_selection_t {
public:
   CMMDBManager *mol;
   int nSelResidues;
   PCResidue *SelResidues;
   int SelectionHandle;
   void clear_up() {
      mol->DeleteSelection(SelectionHandle);
      delete mol;
      mol = 0;
   } 
};


