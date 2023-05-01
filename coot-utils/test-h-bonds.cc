//
#include <string>
#include <iostream>
#include <iomanip>

#include "utils/coot-utils.hh"
#include "coot-h-bonds.hh"

std::string stuartify(const std::string &s) {

   return coot::util::remove_leading_spaces(s);
}

int main(int argc, char **argv) {

   if (argc > 1) {
      std::string file_name = argv[1];
      int  error_count;
      char error_buf[500];

      mmdb::InitMatType();

      mmdb::Manager *m = new mmdb::Manager;
      mmdb::ERROR_CODE err = m->ReadCoorFile(file_name.c_str());
      if (err) {
         std::cout << "There was an error reading " << file_name.c_str() << ".\n";
         std::cout << "ERROR " << err << " READ: "
                   << mmdb::GetErrorDescription(err) << std::endl;
         m->GetInputBuffer(error_buf, error_count);
         if (error_count >= 0) {
            std::cout << " LINE #" << error_count << "\n " << error_buf << std::endl;
         }
      } else {

         int SelHnd = m->NewSelection();
         m->SelectAtoms(SelHnd, 0, "*",
                        mmdb::ANY_RES, "*",
                        mmdb::ANY_RES, "*",
                        "*", "*", "*", "*");

         // energy lib:
         coot::protein_geometry geom;
         std::string energy_cif_file_name = "ener_lib.cif";

         geom.init_standard();
         geom.read_energy_lib(energy_cif_file_name);

         coot::h_bonds hbs;
         std::vector<coot::h_bond> v = hbs.get(SelHnd, SelHnd, m, geom);

         for (unsigned int i=0; i<v.size(); i++) {
            std::cout << v[i].donor_neigh->GetChainID() << "/"
                      << v[i].donor_neigh->GetSeqNum() << "("
                      << v[i].donor_neigh->GetResName() << ")/"
                      << stuartify(v[i].donor_neigh->name)
                      << "        "
                      << v[i].donor->GetChainID() << "/"
                      << v[i].donor->GetSeqNum() << "("
                      << v[i].donor->GetResName() << ")/"
                      << stuartify(v[i].donor->name)
                      << "        "
                      << v[i].acceptor->GetChainID() << "/"
                      << v[i].acceptor->GetSeqNum() << "("
                      << v[i].acceptor->GetResName() << ")/"
                      << stuartify(v[i].acceptor->name)
                      << "       "
                      << v[i].acceptor_neigh->GetChainID() << "/"
                      << v[i].acceptor_neigh->GetSeqNum() << "("
                      << v[i].acceptor_neigh->GetResName() << ")/"
                      << stuartify(v[i].acceptor_neigh->name)
                      << "\n         "
                      << std::fixed << std::setprecision(2)
                      << v[i].angle_1 << "              "
                      << std::fixed << std::setprecision(2)
                      << v[i].dist    << "                "
                      << std::fixed << std::setprecision(2)
                      << v[i].angle_2 << "\n";
            m->DeleteSelection(SelHnd);
         }
      }
   }
   return 0;
}
