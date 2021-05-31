//
#include <clipper/clipper-ccp4.h>

#include "utils/coot-utils.hh"
#include "coot-map-utils.hh"
#include "coot-map-heavy.hh"

clipper::Coord_orth
get_coords_centre(const std::vector<coot::residue_spec_t> &residues,
                  mmdb::Manager *mol) {
   clipper::Coord_orth sum(0,0,0);
   int n_pt=0;
   for (unsigned int ires=0; ires<residues.size(); ires++) {
      mmdb::Residue *residue_p = coot::util::get_residue(residues[ires], mol);
      if (residue_p) {
         mmdb::PPAtom residue_atoms = 0;
         int n_residue_atoms;
         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            mmdb::Atom *at = residue_atoms[iat];
            clipper::Coord_orth pt(at->x, at->y, at->z);
            sum += pt;
            n_pt++;
         }
      }
   }
   double frac = 1/double(n_pt);
   clipper::Coord_orth c(sum.x()*frac, sum.y()*frac, sum.z()*frac);
   return c;
}

int main(int argc, char **argv) {

   int status = 0;

   if (argc > 4) {
      std::string file_name = argv[1];
      int  error_count;
      char error_buf[500];

      mmdb::InitMatType();

      mmdb::Manager *mol = new mmdb::Manager;
      mmdb::ERROR_CODE err = mol->ReadCoorFile(file_name.c_str());
      if (err) {
         std::cout << "There was an error reading " << file_name.c_str() << ".\n";
         std::cout << "ERROR " << err << " READ: "
                   << mmdb::GetErrorDescription(err) << std::endl;
         mol->GetInputBuffer(error_buf, error_count);
         if (error_count >= 0) {
            std::cout << " LINE #" << error_count << "\n " << error_buf << std::endl;
         }
      } else {

         std::string mtz_file_name = argv[2];
         std::string   f_col = argv[3];
         std::string phi_col = argv[4];
         clipper::Xmap<float> xmap;
         bool r = coot::util::map_fill_from_mtz(&xmap, mtz_file_name, f_col, phi_col, "", 0, 0);
         float radius  = 8.0;

         int imod = 1;
         mmdb::Model *model_p = mol->GetModel(imod);
         mmdb::Chain *chain_p;
         int n_chains = model_p->GetNumberOfChains();

         // Just test a few for now.
         int n_residue_trials = 0;
         int n_max_residue_trials = 60;
         for (int ichain=0; ichain<n_chains; ichain++) {
            chain_p = model_p->GetChain(ichain);
            int nres = chain_p->GetNumberOfResidues();
            mmdb::Residue *residue_p;
            for (int ires=0; ires<nres; ires++) {

               residue_p = chain_p->GetResidue(ires);

               // hack for a few test residues
               n_residue_trials++;
               if (n_residue_trials > n_max_residue_trials)
                  break;
               if (residue_p->GetSeqNum() > 360)
                  break;

               int SelHnd = mol->NewSelection();

               std::vector<coot::residue_spec_t> local_residues =
                  coot::residues_near_residue(coot::residue_spec_t(residue_p),
                                              mol, radius);
               local_residues.push_back(coot::residue_spec_t(residue_p));

               for (unsigned int ilocal=0; ilocal<local_residues.size(); ilocal++) {

                  mol->SelectAtoms(SelHnd, 1,
                                   local_residues[ilocal].chain_id.c_str(),
                                   local_residues[ilocal].res_no,
                                   local_residues[ilocal].ins_code.c_str(),
                                   local_residues[ilocal].res_no,
                                   local_residues[ilocal].ins_code.c_str(),
                                   "*", // any residue name
                                   "*", // atom name
                                   "*", // elements
                                   "*",  // alt loc.
                                   mmdb::SKEY_OR
                                   );
               }
               mmdb::PPAtom atom_sel=NULL;
               int n_sel;
               mol->GetSelIndex(SelHnd, atom_sel, n_sel);

               std::cout << "selected " << n_sel << " atoms " << std::endl;

               std::cout << "------- fffearing -------" << std::endl;
               coot::util::fffear_search f(mol, SelHnd, xmap, radius, true);
               std::vector<std::pair<float, clipper::RTop_orth> > p = f.scored_orientations();

               clipper::Coord_orth coord_centre = get_coords_centre(local_residues, mol);
               std::cout << "coordinates centre: " << coord_centre.format() << std::endl;
               for (unsigned int ires=0; ires<p.size(); ires++) {
                  std::cout << "   " << p[ires].first << " " << p[ires].second.trn().format()
                            << std::endl;
               }

               std::string map_filename = "search-";
               map_filename += coot::util::int_to_string(residue_p->GetSeqNum());
               map_filename += ".map";
               clipper::CCP4MAPfile mapout;
               mapout.open_write(std::string(map_filename));
               mapout.export_xmap(f.get_results_map());
               mapout.close_write();

               mol->DeleteSelection(SelHnd);
            }
         }
      }
   }
   return status;
}
