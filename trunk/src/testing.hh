
#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>

#define BUILT_IN_TESTING

std::string stringify(double x);
std::string stringify(int i);
std::string stringify(unsigned int i);
std::string greg_test(const std::string &file_name);

int test_internal();

int test_alt_conf_rotamers();
int test_wiggly_ligands();
