
#include <vector>
#include <string>
#include <fstream>

#include "clipper/core/xmap.h"
#include "clipper/ccp4/ccp4_map_io.h"
#include "clipper/core/hkl_compute.h"

#include "analysis/stats.hh"
#include "utils/coot-utils.hh"
#include "cablam-markup.hh"
#include "coot-map-utils.hh" // for max_gridding()

class table_info_t {
public:
   std::string ent_file_name;
   std::string chain_id;
   int res_no;
   coot::residue_spec_t spec;
   double dp_prev_to_mid;
   double dp_next_to_mid;
   double dist_proj_point_prev_to_next;
   table_info_t(const std::string &fn, const std::string &chid, int res_no, const coot::residue_spec_t &spec,
                double dd1, double dd2, double dd3) : ent_file_name(fn), chain_id(chid), res_no(res_no),
                                                      spec(spec), dp_prev_to_mid(dd1), dp_next_to_mid(dd2),
                                                      dist_proj_point_prev_to_next(dd3) {}
};

void
make_a_map(const std::vector<table_info_t> &table) {

   std::cout << "Making a map" << std::endl;

   std::string map_file_name = "table_as_map.map";
   double d_1_min =  9999999999.9;
   double d_1_max = -9999999999.9;
   double d_2_min =  9999999999.9;
   double d_2_max = -9999999999.9;
   double d_3_min =  9999999999.9;
   double d_3_max = -9999999999.9;

   coot::stats::single data_1;
   coot::stats::single data_2;
   coot::stats::single data_3;

   for (const auto &entry : table) {
      if (entry.dp_prev_to_mid > d_1_max) d_1_max = entry.dp_prev_to_mid;
      if (entry.dp_prev_to_mid < d_1_min) d_1_min = entry.dp_prev_to_mid;
      if (entry.dp_next_to_mid > d_2_max) d_2_max = entry.dp_next_to_mid;
      if (entry.dp_next_to_mid < d_2_min) d_2_min = entry.dp_next_to_mid;
      if (entry.dist_proj_point_prev_to_next > d_3_max) d_3_max = entry.dist_proj_point_prev_to_next;
      if (entry.dist_proj_point_prev_to_next < d_3_min) d_3_min = entry.dist_proj_point_prev_to_next;
      data_1.add(entry.dp_prev_to_mid);
      data_2.add(entry.dp_next_to_mid);
      data_3.add(entry.dist_proj_point_prev_to_next);
   }

   std::cout << "lims1: " << d_1_min << " " << d_1_max << std::endl;
   std::cout << "lims2: " << d_2_min << " " << d_2_max << std::endl;
   std::cout << "lims3: " << d_3_min << " " << d_3_max << std::endl;

   auto mi_1 = data_1.median_and_iqr();
   auto mi_2 = data_2.median_and_iqr();
   auto mi_3 = data_3.median_and_iqr();

   std::cout << "data_1 " << mi_1.first << " " << mi_1.second << std::endl;
   std::cout << "data_2 " << mi_2.first << " " << mi_2.second << std::endl;
   std::cout << "data_3 " << mi_3.first << " " << mi_3.second << std::endl;

   // limits -5.6 to 5.2
   // limits -5.6 to 5.2
   // limits 2.7 8.0

   clipper::Spacegroup spacegroup(clipper::Spgr_descr("P1"));
   clipper::Cell cell(clipper::Cell_descr(100, 100, 100));
   clipper::Grid_sampling grid_sampling(100, 100, 100);
   clipper::Xmap<float> t(spacegroup, cell, grid_sampling);

   clipper::Xmap_base::Map_reference_index ix;
   for (ix = t.first(); !ix.last(); ix.next() )  { // iterator index.
      t[ix] = 0.0;
   }

   float r1 = 5.2 - -5.6;
   float r2 = 5.2 - -5.6;
   float r3 = 8.0 - 2.7;
   float min_x =  -5.6;
   float min_y =  -5.6;
   float min_z =  -2.7;

   for (const auto &entry : table) {
      int bin_x = static_cast<int>(100 * (entry.dp_prev_to_mid - min_x)/r1);
      int bin_y = static_cast<int>(100 * (entry.dp_next_to_mid - min_y)/r2);
      int bin_z = static_cast<int>(100 * (entry.dist_proj_point_prev_to_next - min_z)/r3);
      clipper::Coord_grid cg(bin_x, bin_y, bin_z);
      int inc = t.get_data(cg) + 1;
      t.set_data(cg, inc);
   }
   clipper::CCP4MAPfile mapout;
   mapout.open_write(map_file_name);
   mapout.export_xmap(t);
   mapout.close_write();

   float mg = coot::util::max_gridding(t);
   std::cout << "INFO:: Map max gridding " << mg << " A/grid-point" << std::endl;
   clipper::Resolution reso(2.0 * mg); // for map -> data, the resolution is half the gridding
   clipper::HKL_info myhkl(spacegroup, cell, reso, true);
   
   clipper::HKL_data< clipper::datatypes::F_phi<float> > fphis(myhkl);
   t.fft_to(fphis);

   std::cout << "fphis num_obs " << fphis.num_obs() << std::endl;

   clipper::HKL_info::HKL_reference_index hri;
   for (hri = fphis.first(); !hri.last(); hri.next()) {
      // float irs = hri.invresolsq();
      // std::cout << hri.hkl().format() << " has reso " << irs <<  " f " << fphis[hri].f() << std::endl;
   }
}

void write_a_big_table(const std::vector<table_info_t> &table) {

   std::cout << "Making a table " << std::endl;

   clipper::Spacegroup spacegroup(clipper::Spgr_descr("P1"));
   clipper::Cell cell(clipper::Cell_descr(100, 100, 100));
   clipper::Grid_sampling grid_sampling(100, 100, 100);
   clipper::Xmap<float> t(spacegroup, cell, grid_sampling);

   clipper::Xmap_base::Map_reference_index ix;
   for (ix = t.first(); !ix.last(); ix.next() )  { // iterator index.
      t[ix] = 0.0;
   }

   float r1 = 5.2 - -5.6;
   float r2 = 5.2 - -5.6;
   float r3 = 8.0 - 2.7;
   float min_x =  -5.6;
   float min_y =  -5.6;
   float min_z =  -2.7;

   for (const auto &entry : table) {
      int bin_x = static_cast<int>(100 * (entry.dp_prev_to_mid - min_x)/r1);
      int bin_y = static_cast<int>(100 * (entry.dp_next_to_mid - min_y)/r2);
      int bin_z = static_cast<int>(100 * (entry.dist_proj_point_prev_to_next - min_z)/r3);
      clipper::Coord_grid cg(bin_x, bin_y, bin_z);
      int inc = t.get_data(cg) + 1;
      t.set_data(cg, inc);
   }

   // reverse of that indexing:
   //
   // p_1 = (bin_x + 0.5) * r1 + min_x
   // p_2 = (bin_y + 0.5) * r2 + min_y
   // p_3 = (bin_z + 0.5) * r3 + min_z

   std::string table_file_name("big.table");
   std::ofstream f(table_file_name.c_str());
   if (f) {
      for (ix = t.first(); !ix.last(); ix.next() )  { // iterator index.
         clipper::Coord_grid cg = ix.coord();
         f << cg.u() << " " << cg.v()  << " " << cg.w() << " " << t[ix] << "\n";
      }
      f.close();
   }
}

void proc_results(const std::string &dir) {

   std::vector<table_info_t> table;

   std::vector<std::string> files = coot::util::glob_files(dir, "*.table");
   for (const auto &ff : files) {
      std::ifstream f(ff.c_str());
      if (f) {
         std::string line;
         while (std::getline(f, line)) {
            std::vector<std::string> parts = coot::util::split_string_no_blanks(line);
            if (parts.size() == 6) {
               try {
                  std::string ch_id = parts[1];
                  int rn = coot::util::string_to_int(parts[2]);
                  coot::residue_spec_t spec(ch_id, rn, "");
                  double dp_prev_to_mid = coot::util::string_to_double(parts[3]);
                  double dp_next_to_mid = coot::util::string_to_double(parts[4]);
                  double dist_proj_point_prev_to_next = coot::util::string_to_double(parts[5]);
                  table_info_t ti(ff, ch_id, rn, spec, dp_prev_to_mid, dp_next_to_mid, dist_proj_point_prev_to_next);
                  table.push_back(ti);
               }
               catch (...) {}
            }
	 }
      }
   }

   std::cout << "table has " << table.size() << " entries" << std::endl;

   make_a_map(table);
   write_a_big_table(table);

}

int main(int argc, char **argv) {

   if (argc < 2) {
      std::cout << "Usage: make-cablam-like-stats <filename>" << std::endl;
      std::cout << "or " << std::endl;
      std::cout << "Usage: make-cablam-like-stats proc-results <tables-dir-name>" << std::endl;
      exit(0);
   } else {

      std::string argv_1 = argv[1];

      if (argv_1 == "proc-results") {
         if (argc == 3) {
            std::string dir = argv[2];
            proc_results(dir);
         }

      } else {

         std::string file_name = argv_1;

         mmdb::Manager *mol = new mmdb::Manager;
         int read_status = mol->ReadCoorFile(file_name.c_str());
         if (read_status == mmdb::Error_NoError) {
            std::vector<coot::cablam_like_geometry_stats_t> v = coot::get_cablam_like_geometry_stats(mol);
            std::string fn_log_file = coot::util::name_sans_extension(file_name) + "-geom-stats.table";
            std::ofstream f(fn_log_file.c_str());
            if (f) {
               for (const auto &item : v) {
                  if (item.residue) {
                     coot::residue_spec_t spec(item.residue);
                     if (! spec.chain_id.empty()) {
                        f << file_name << " " << spec.chain_id << " " << spec.res_no << " "
                          << item.dp_prev_to_mid << " " << item.dp_next_to_mid << " "
                          << item.dist_proj_point_prev_to_next
                          << "\n";
                     }
                  }
               }
               f.close();
            }
         }
      }
   }

   return 0;
}
