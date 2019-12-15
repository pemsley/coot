
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

      std::vector<std::vector<clipper::Coord_orth> > accum;
      if (! pdb_files.empty()) {
	 for (std::size_t i=0; i<pdb_files.size(); i++) {
	    const std::string &pdb_file = pdb_files[i];
	    mmdb::Manager *mol = new mmdb::Manager;
	    mol->ReadPDBASCII(pdb_file.c_str());

	    std::vector<std::vector<clipper::Coord_orth> > v = coot::mol_to_5_residue_strand_fragments(mol);
	    std::cout << "INFO:: " << v.size() << " fragments in " << pdb_file << std::endl;
	    std::move(v.begin(), v.end(), std::back_inserter(accum));

	    delete mol;
	 }
      }

      if (accum.size() > 1) {

         std::cout << "Accumulated " << accum.size() << std::endl;
         for (std::size_t i=0; i<accum.size(); i++) {
            for (std::size_t j=0; j<accum.size(); j++) {
               clipper::RTop_orth rtop(accum[i], accum[j]);
               double sum = 0.0;
               for (unsigned int k=0; k<5; k++) {
                   double d_sqrd = (accum[j][k] - accum[i][k].transform(rtop)).lengthsq();
                   double d = sqrt(d_sqrd);
                   sum += d;
               }
               std::cout << " " << i << " " << j << " " << 0.2 * sum << std::endl;
            }
         }
     }
  }

   return status;

}
