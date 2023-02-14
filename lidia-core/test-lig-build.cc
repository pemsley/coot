//
#include <fstream>
#include "svg-molecule.hh"
#include "geometry/protein-geometry.hh"
#include "rdkit-interface.hh"

#include "utils/coot-utils.hh"

void
tlc_to_svg_file(const std::string &tlc, coot::protein_geometry &geom) {

   if (tlc == "A4R") return;
   if (tlc == "A85") return;
   if (tlc == "ABZ") return;
   if (tlc == "AC" ) return;
   if (tlc == "AC9") return;
   if (tlc == "AG" ) return;
   if (tlc == "AG1") return;
   if (tlc == "AL" ) return;
   if (tlc == "ALB") return;
   if (tlc == "APW") return;
   if (tlc == "AQZ") return;
   if (tlc == "AR" ) return;
   if (tlc == "ARS") return;
   if (tlc == "ASX") return;
   if (tlc == "ASZ") return;
   if (tlc == "AT" ) return;
   if (tlc == "ATZ") return;
   if (tlc == "AU" ) return;
   if (tlc == "AU3") return;
   if (tlc == "AZZ") return;

   int imol = coot::protein_geometry::IMOL_ENC_ANY;
   mmdb::Residue *r = geom.get_residue(tlc, imol, true, true, 30.0);
   if (r) {

      try {
         bool undelocalize_flag = true;
         // used undelocalize_flag in RDKit::RWMol mol_rw = coot::rdkit_mol(r, rest, "", undelocalize_flag);
         RDKit::RWMol mol_rw = coot::rdkit_mol_sanitized(r, imol, geom);

         // coot::remove_non_polar_Hs(&rdkm); either is good.
         RDKit::MolOps::removeHs(mol_rw);
         RDKit::MolOps::Kekulize(mol_rw);
         int iconf = RDDepict::compute2DCoords(mol_rw, NULL, true);
         RDKit::Conformer &conf = mol_rw.getConformer(iconf);
         RDKit::WedgeMolBonds(mol_rw, &conf);
         svg_molecule_t svg;
         svg.import_rdkit_mol(&mol_rw, iconf);
         std::string s = svg.render_to_svg_string();
         if (true) {
            std::string fn = tlc + std::string(".svg");
            std::ofstream f(fn);
            if (f)
               f << s;
         }
      }
      catch (const Invar::Invariant &e) {
         std::cout << "error " << e.what() << std::endl;
      }
   } else {
      std::cout << "Failed to get residue from dictionary" << std::endl;
   }
}

int main(int argc, char **argv) {

   int status = 0;
   std::string tlc = "CFF";
   coot::protein_geometry geom;
   geom.set_verbose(false);

   if (argc > 1)
      tlc = argv[1];
   if (tlc == "dir") {

      // 20221201-PE 3.8 ms/cif

      const char *e = getenv("CLIBD_MON");
      if (e) {
         std::string dir = coot::util::append_dir_dir(e, "a");
         std::string glob_pattern("*.cif");
         std::vector<std::string> file_names = coot::util::glob_files(dir, glob_pattern);
         unsigned int n_cifs = 0;
         for (unsigned int i=0; i<file_names.size(); i++) {
            const std::string &fn = file_names[i];
            std::string fn_l = coot::util::file_name_non_directory(fn);
            auto parts = coot::util::split_string_no_blanks(fn_l, ".");
            std::string tlc = parts[0];
            // std::cout << "   " << fn << " " << fn_l << " " << tlc << std::endl;
            tlc_to_svg_file(tlc, geom);
            n_cifs++;
            // std::cout << "n_cifs: " << n_cifs << std::endl;
         }
      }

   } else {
      tlc_to_svg_file(tlc, geom);
   }
   return status;
}
