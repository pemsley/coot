

#include <iostream>
#include <iomanip>

#include "clipper/ccp4/ccp4_map_io.h"
#include "analysis/stats.hh"
#include "atom-selection-container.hh"
#include "coot-map-utils.hh"

class residue_tube_density_container_t {
public:
   std::vector<std::vector<float> > ring_0;
   std::vector<std::vector<float> > ring_1;
   std::vector<std::vector<float> > ring_2;
   coot::residue_spec_t res_spec;
   residue_tube_density_container_t() {}
   residue_tube_density_container_t(const coot::residue_spec_t &rs) : res_spec(rs) {}
   bool empty() const { return ring_0.empty(); }
   void normalize() {

      coot::stats::single s;
      std::vector<std::vector<float> > &ring_ref = ring_1;
      for (unsigned int i=0; i<ring_ref.size(); i++) {
         const std::vector<float> &v_in = ring_ref[i];
         for (unsigned int j=0; j<v_in.size(); j++) {
            s.add(v_in[j]);
         }
      }
      const double mean = s.mean();
      std::cout << "normalize mean " << res_spec << " " << mean << std::endl;

      std::vector<std::reference_wrapper<std::vector<std::vector<float> > > > refs = { ring_0, ring_1, ring_2};
      for (auto &ring : refs) {
         for (unsigned int i=0; i<ring.get().size(); i++) {
            std::vector<float> &v_in = ring.get()[i];
            for (unsigned int j=0; j<v_in.size(); j++) {
               v_in[j] /= mean;
            }
         }
      }

      if (true) { // remove this debugging
         std::string chain_id = res_spec.chain_id;
         if (chain_id == "") chain_id = "---";
         int res_no = res_spec.res_no;
         for (unsigned int i=0; i<ring_0.size(); i++) {
            std::vector<float> &v_in = ring_0[i];
            for (unsigned int j=0; j<v_in.size(); j++) {
               const float &d = v_in[j];
               std::cout << "tube-density-debug ring_0 " << chain_id << " " << res_no
                         << " " << i << " " << j << " density " << d << "\n";
            }
         }
         for (unsigned int i=0; i<ring_1.size(); i++) {
            std::vector<float> &v_in = ring_1[i];
            for (unsigned int j=0; j<v_in.size(); j++) {
               const float &d = v_in[j];
               std::cout << "tube-density-debug ring_1 " << chain_id << " " << res_no
                         << " " << i << " " << j << " density " << d << "\n";
            }
         }
         for (unsigned int i=0; i<ring_2.size(); i++) {
            std::vector<float> &v_in = ring_2[i];
            for (unsigned int j=0; j<v_in.size(); j++) {
               const float &d = v_in[j];
               std::cout << "tube-density-debug ring_2 " << chain_id << " " << res_no
                         << " " << i << " " << j << " density " << d << "\n";
            }
         }
      }
   }
};

std::optional<clipper::Coord_orth>
get_atom_pos(mmdb::Manager *mol, coot::atom_spec_t &atom_spec) {

   mmdb::Atom *at = coot::util::get_atom(atom_spec, mol);
   clipper::Coord_orth co(0,0,0);
   if (at) {
      return clipper::Coord_orth(at->x, at->y, at->z);
   } else {
      return {};
   }
}

std::optional<clipper::Coord_orth>
get_atom_pos(mmdb::Residue *residue_p, const std::string &atom_name) {

   mmdb::Atom **residue_atoms = 0;
   int n_residue_atoms = 0;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      mmdb::Atom *at = residue_atoms[iat];
      if (! at->isTer()) {
         std::string n = at->GetAtomName();
         if (n == atom_name) {
            clipper::Coord_orth co(at->x, at->y, at->z);
            return co;
         }
      }
   }
   return {};
}


residue_tube_density_container_t
get_tube(mmdb::Residue *residue_p, const clipper::Xmap<float> &xmap, const std::vector<double> &radii,
         unsigned int n_length, unsigned int n_ring) {

   // coot::atom_spec_t as_1("A", 12, "", " C4'", "");
   // coot::atom_spec_t as_2("A", 12, "", " C5'", "");
   // coot::atom_spec_t as_r("A", 12, "", " O5'", "");

   residue_tube_density_container_t rtdc;

   std::string at_name_1   = " C5'";
   std::string at_name_2   = " C4'";
   std::string at_name_ref = " O5'";

   std::optional<clipper::Coord_orth> pt_1   = get_atom_pos(residue_p, at_name_1);
   std::optional<clipper::Coord_orth> pt_2   = get_atom_pos(residue_p, at_name_2);
   std::optional<clipper::Coord_orth> pt_ref = get_atom_pos(residue_p, at_name_ref);
   if (pt_1 && pt_2 && pt_ref) {
      rtdc.res_spec = coot::residue_spec_t(residue_p);
      for (unsigned int i=0; i<radii.size(); i++) {
         double radius = radii[i];
         std::vector<std::vector<float> > v =
            coot::util::get_density_on_cylinder(pt_1.value(), pt_2.value(), pt_ref.value(), xmap, radius, n_length, n_ring);
         // bleugh
         if (i == 0) rtdc.ring_0 = v;
         if (i == 1) rtdc.ring_1 = v;
         if (i == 2) rtdc.ring_2 = v;
      }
   } else {
      std::string rn = residue_p->GetResName();
      int res_no = residue_p->GetSeqNum();
      std::string chain_id = residue_p->GetChainID();
      if (false)
         std::cout << "couldn't find atoms " << at_name_1 << " " << at_name_2 << " " << at_name_ref
                   << " in residue " << chain_id << " " << res_no << " with name " << rn << std::endl;
   }
   return rtdc;
}

void tube_analysis(const std::vector<residue_tube_density_container_t> &tubes) {

   // key is based on tube-number, ring index and height index
   //
   std::map<unsigned int, coot::stats::single> s_map;
   std::map<unsigned int, coot::residue_spec_t> spec_map;
   for (const auto &tube_r : tubes) {
      auto tube = tube_r;
      std::vector<std::reference_wrapper<std::vector<std::vector<float> > > > refs =
         { tube.ring_0, tube.ring_1, tube.ring_2};
      for (unsigned int iring=0; iring<refs.size(); iring++) {
         const auto &ring = refs[iring].get();
         for (unsigned int i=0; i<ring.size(); i++) {
            const std::vector<float> &v_in = ring[i];
            for (unsigned int j=0; j<v_in.size(); j++) {
               unsigned int idx = iring * 10000 + i * 100 + j;
               std::cout << "debug idx " << idx << std::endl;
               spec_map[idx] = tube_r.res_spec;
               const float &f = v_in[j];
               s_map[idx].add(f);
            }
         }
      }
   }

   std::map<unsigned int, std::pair<double, double> > mean_sd_map;
   std::map<unsigned int, coot::stats::single>::const_iterator it;
   for (it=s_map.begin(); it!=s_map.end(); it++) {
      unsigned int idx = it->first;
      const auto &s(*it);
      double mean = s.second.mean();
      double var  = s.second.variance();
      double sd   = std::sqrt(var);
      mean_sd_map[idx] = std::make_pair(mean, sd);
      const auto &data = s.second.v;
      double ks = coot::stats::get_kolmogorov_smirnov_vs_normal(data, mean, sd);
      std::cout << "distribution-analysis idx " << idx << " mean " << std::setw(9) << std::left << mean
                << " sd " << std::setw(9) << sd << " n-data " << std::right << std::setw(5) << data.size()
                << " ks " << ks << std::endl;
      if (idx == 10101) {
         for (unsigned int i=0; i<data.size(); i++) {
            std::cout << "idx " << idx << " i " << i << " " << data[i] << std::endl;
         }
      }
   }

   for (const auto &tube_r : tubes) {
      const auto &ring = tube_r.ring_1;
      unsigned int iring = 1;
      double sum_nzz = 0.0;
      unsigned int count = 0;
      for (unsigned int i=0; i<ring.size(); i++) {
         const std::vector<float> &v_in = ring[i];
         for (unsigned int j=0; j<v_in.size(); j++) {
            unsigned int idx = iring * 10000 + i * 100 + j;
            std::map<unsigned int, std::pair<double, double> >::const_iterator it;
            it = mean_sd_map.find(idx);
            if (it == mean_sd_map.end()) {
               std::cout << "idx-find fail " << idx << std::endl;
            } else {
               const double &ref_mean = mean_sd_map.at(idx).first;
               const double &ref_sd   = mean_sd_map.at(idx).second;
               double d = v_in[j];
               double n_z = (d-ref_mean)/ref_sd;
               if (idx   == 1001)
                  std::cout << "debug-nz i " << i << " j " << j << " n_z " << n_z << std::endl;
               sum_nzz += n_z * n_z;
               count++;
            }
         }
      }
      sum_nzz -= static_cast<double>(count);
      std::cout << "result-for-res " << tube_r.res_spec.chain_id << " " << tube_r.res_spec.res_no
                << " " << sum_nzz << std::endl;
   }
}


int main(int argc, char **argv) {

   int status = 0;

   unsigned int n_ring   = 10;
   unsigned int n_length = 6;
   std::vector<double> radii = {0.4, 0.8, 1.2};

   if (argc > 2) {
      std::string map_file_name(argv[1]);
      std::string pdb_file_name(argv[2]);
      try {
         clipper::CCP4MAPfile file;
         file.open_read(map_file_name);
         clipper::Xmap<float> xmap;
         file.import_xmap(xmap);

         // now get the atoms
         atom_selection_container_t asc = get_atom_selection(pdb_file_name, false);

         if (asc.read_success) {

            mmdb::Manager *mol = asc.mol;

            int imod = 1;
            mmdb::Model *model_p = mol->GetModel(imod);
            if (model_p) {

               std::vector<residue_tube_density_container_t> all_the_tubes;
               int n_chains = model_p->GetNumberOfChains();
               for (int ichain=0; ichain<n_chains; ichain++) {

                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                  std::pair<int, int> p = coot::util::min_and_max_residues(chain_p);
                  if (p.first == 9999) {
                     std::string chain_id = chain_p->GetChainID();
                     std::cout << "fail min/max for chain-ID " << chain_id << std::endl;
                  } else {
                     int n_res = chain_p->GetNumberOfResidues();
                     for (int ires=0; ires<n_res; ires++) {
                        mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                        if (residue_p) {
                           residue_tube_density_container_t rtdc = get_tube(residue_p, xmap, radii, n_length, n_ring);
                           if (! rtdc.empty()) {
                              rtdc.normalize();
                              all_the_tubes.push_back(rtdc);
                           }
                        }
                     }
                  }
               }
               tube_analysis(all_the_tubes);
            }
         }
      }
      catch (const clipper::Message_base &exc) {
         std::cout << "WARNING:: failed to open " << map_file_name << std::endl;
      }
   }
   return status;
}
