
#include <fstream>
#include "utils/coot-utils.hh"
#include "strand-fragments.hh"

int main(int argc, char **argv) {

   int status = 0;

   std::vector<std::string> pdb_files;
   if (argc > 1) {
      for (int i=1; i<argc; i++) {
         std::string a(argv[i]);
         if (coot::file_exists(a))
            pdb_files.push_back(a);
     }

     std::vector<std::vector<clipper::Coord_orth> > ca_accum;
     std::vector<std::vector<clipper::Coord_orth> > mc_accum;

     if (! pdb_files.empty()) {
         for (std::size_t i=0; i<pdb_files.size(); i++) {
            const std::string &pdb_file = pdb_files[i];
            mmdb::Manager *mol = new mmdb::Manager;
            mol->ReadPDBASCII(pdb_file.c_str());

            std::pair<std::vector<std::vector<clipper::Coord_orth> >, std::vector<std::vector<clipper::Coord_orth> > > v = coot::mol_to_5_residue_strand_fragments(mol);
            std::cout << "INFO:: " << v.first.size() << " fragments in " << pdb_file << std::endl;
            std::move(v.first.begin(),  v.first.end(),  std::back_inserter(ca_accum));
            std::move(v.second.begin(), v.second.end(), std::back_inserter(mc_accum));

            delete mol;
         }
      }

      if (ca_accum.size() > 1) {

         std::cout << "Accumulated " << ca_accum.size() << std::endl;

         std::ofstream f("dist-O-O-1-2.table");

         for (std::size_t i=0; i<ca_accum.size(); i++) {
            for (std::size_t j=0; j<ca_accum.size(); j++) {
               if (i != j) {

                  clipper::RTop_orth rtop(ca_accum[i], ca_accum[j]);
                  double sum = 0.0;
                  for (unsigned int k=0; k<5; k++) {
                      double d_sqrd = (ca_accum[j][k] - ca_accum[i][k].transform(rtop)).lengthsq();
                      double d = sqrt(d_sqrd);
                      sum += d;
                  }
                  std::cout << " " << i << " " << j << " " << 0.2 * sum << std::endl;

                  // file
                  unsigned int ires_1 = 1; // in fragment 0-based.
                  unsigned int ires_2 = 2;
                  unsigned int idx_1 = ires_1 * 4 + 3;  // + 3 indexes O
                  unsigned int idx_2 = ires_2 * 4 + 3;
                  double d_1_sqrd = (mc_accum[i][idx_2] - mc_accum[i][idx_1]).lengthsq();
                  double d_2_sqrd = (mc_accum[j][idx_2] - mc_accum[j][idx_1]).lengthsq();
                  double d_1 = sqrt(d_1_sqrd);
                  double d_2 = sqrt(d_2_sqrd);
                  double d = d_2 - d_1;
                  // f << " " << i << " " << j << " " << d << std::endl;
                  f << " " << i << " " << j << " " << fabs(d) << " for CA-distance " << 0.2 * sum << std::endl;
               }
            }
         }
         f.close();
     }
  }

   return status;

}
