/* coot-utils/slurp-map.hh
 *
 * Copyright 2019 by Medical Research Council
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include <iostream>
#include <system_error>
#include <chrono>
#include <zlib.h>

#include "utils/coot-utils.hh"
#include "slurp-map.hh"

// Test on 5778, 10289, 6338

#include <fstream>
#include <iostream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "utils/logging.hh"
extern logging logger;

std::string
coot::util::to_string(slurp_map_result_t r) {
   std::string s;
   if (r == slurp_map_result_t::OK)                      s = "OK";
   if (r == slurp_map_result_t::IS_SLURPABLE_EM_MAP)     s = "IS_SLURPABLE_EM_MAP";
   if (r == slurp_map_result_t::IS_NON_SLURPABLE_EM_MAP) s = "IS_NON_SLURPABLE_EM_MAP";
   if (r == slurp_map_result_t::NOT_AN_EM_MAP)           s = "NOT_AN_EM_MAP";
   if (r == slurp_map_result_t::FAIL)                    s = "FAIL";
   if (r == slurp_map_result_t::FILE_NOT_FOUND)          s = "FILE_NOT_FOUND";
   if (r == slurp_map_result_t::UNRESOLVED)              s = "UNRESOLVED";
   return s;
}

std::ostream&
coot::util::operator<<(std::ostream &s, slurp_map_result_t smr) {
   s << to_string(smr);
   return s;
}


coot::util::slurp_map_result_t
coot::util::is_basic_em_map_file(const std::string &file_name) {

   clipper::Xmap<float> xmap;
   return slurp_fill_xmap_from_map_file(file_name, &xmap, true); // check-only mode
}

#include "voidp-buffer.hh"

coot::util::slurp_map_result_t
coot::util::slurp_fill_xmap_from_map_file(const std::string &file_name,
                                          clipper::Xmap<float> *xmap_p,
                                          bool check_only) { // default arg, false

   auto test_do_something = [] (std::basic_streambuf<char> *data) {
      bool debug = true;
      // std::cout << "data: " << data << std::endl;
      int n_rows = -1;
      int n_cols = -1;
      int n_secs = -1;
      n_cols = *reinterpret_cast<int *>(data);
      n_rows = *reinterpret_cast<int *>(data+4);
      n_secs = *reinterpret_cast<int *>(data+8);
      if (debug) {
         std::cout << "debug:: n_cols " << n_cols << std::endl;
         std::cout << "debug:: n_rows " << n_rows << std::endl;
         std::cout << "debug:: n_sections " << n_secs << std::endl;
      }
   };

   auto slurp_fill_xmap_from_gz_map_file = [] (const std::string &file_name,
                                               clipper::Xmap<float> *xmap_p,
                                               bool check_only) {

      slurp_map_result_t status = slurp_map_result_t::UNRESOLVED;
      struct stat s;
      int fstat = stat(file_name.c_str(), &s);
      if (fstat == 0) {
         gzFile file = gzopen(file_name.c_str(), "rb");
         int z_status = Z_OK;
         voidp_buffer_t buff(4);
         size_t read_pos = 0;
         while (! gzeof(file)) {
            size_t space_remaining = buff.size() - read_pos;
            int bytes_read = gzread(file, (char *)buff.get() + read_pos, space_remaining);
            const char *error_message = gzerror(file, &z_status);
            if ((bytes_read == -1) || z_status != Z_OK) {
               std::cout << "WARNING:: gz read error for " << file_name << " "
                         << error_message << std::endl;
               break;
            }
            read_pos += bytes_read;
            if (buff.size() == read_pos) {
               buff.resize(buff.size() * 2);
            }
         }
         z_status = gzclose_r(file);
         if (z_status != Z_OK) {
            std::cout << "WARNING:: gz close read error for " << file_name << std::endl;
         }
         if (read_pos >= buff.size()) {
            buff.resize(buff.size() + 1);
         }
         *((char *)buff.get() + read_pos) = 0;
         char *data = reinterpret_cast< char *>(buff.get());
         status = slurp_parse_xmap_data(data, xmap_p, check_only); // fill xmap
         std::cout << "DEBUG:: slurp_parse_xmap_data() returned with status "
                   << int(status) << std::endl;
      }
      return status;
   };

   slurp_map_result_t status = slurp_map_result_t::UNRESOLVED;
   if (file_exists(file_name)) {

      bool is_gzip = false;
      std::string ext = file_name_extension(file_name);
      if (ext == ".gz") is_gzip = true;

      if (is_gzip) {
         // This can fail (at the moment) if the axes are not in X,Y,Z order.
         // Also, this returns false when there is PANDDA:: in the header.
         status = slurp_fill_xmap_from_gz_map_file(file_name, xmap_p, check_only);
      } else {
         struct stat s;
         int fstat = stat(file_name.c_str(), &s);
         if (fstat == 0) {
            FILE *fptr = fopen(file_name.c_str(), "rb");
            off_t st_size = s.st_size;
            // std::cout << "st_size: " << st_size << std::endl;
            try {
               // 20231006-PE as it used to be.
               char *space = new char[st_size+1];
               // Happy Path
               size_t st_size_2 = fread(space, st_size, 1, fptr);
               char *data = static_cast<char *>(space);
               fclose(fptr);
               if (st_size_2 == 1) {
                  // Happy Path
                  if (st_size > 1024) {
                     status = slurp_parse_xmap_data(data, xmap_p, check_only); // fill xmap
                  } else {
                     std::cout << "WARNING:: bad read " << file_name << std::endl;
                  }
               } else {
                  std::cout << "WARNING:: bad read " << file_name << std::endl;
               }
               delete [] space;
            }
            catch (const std::bad_alloc &e) {
               std::cout << "WARNING:: out-of-memory " << st_size+1 << " " << e.what() << std::endl;
            }
         }
      }
   } else {
      status = slurp_map_result_t::FILE_NOT_FOUND;
      std::cout << "WARNING:: file does not exist " << file_name << std::endl;
   }

   logger.log(log_t::DEBUG, "slurp_fill_xmap_from_map_file() with check-mode: ",
              check_only, "returning", to_string(status));

   return status;
}

#include "utils/split-indices.hh"

// return value (status) means "is_basic_EM_map" (that's slurpable)
coot::util::slurp_map_result_t
coot::util::slurp_parse_xmap_data(char *data,
                                  clipper::Xmap<float> *xmap_p,
                                  bool check_only) {

   bool debug = false;
   bool print_labels = true;

   slurp_map_result_t status = slurp_map_result_t::UNRESOLVED;
   int n_rows = -1;
   int n_cols = -1;
   int n_secs = -1;
   n_cols = *reinterpret_cast<int *>(data);
   n_rows = *reinterpret_cast<int *>(data+4);
   n_secs = *reinterpret_cast<int *>(data+8);
   if (debug) {
      std::cout << "debug:: n_cols " << n_cols << std::endl;
      std::cout << "debug:: n_rows " << n_rows << std::endl;
      std::cout << "debug:: n_sections " << n_secs << std::endl;
   }

   int mode = -1;
   int data_size = 4; // 4 bytes for a float
   mode = *reinterpret_cast<int *>(data+12);
   if (mode == 0) data_size = 1;
   if (mode == 1) data_size = 2;
   if (mode == 6) data_size = 2; // two-byte integer. Needs fixing in the casting.
   if (debug)
      std::cout << "debug:: slurp_map() mode: " << mode << std::endl;

   int nx_start = -1, ny_start = -1, nz_start = -1;
   int mx = -1, my = -1, mz = -1;
   nx_start = *reinterpret_cast<int *>(data+16);
   ny_start = *reinterpret_cast<int *>(data+20);
   nz_start = *reinterpret_cast<int *>(data+24);
   mx = *reinterpret_cast<int *>(data+28);
   my = *reinterpret_cast<int *>(data+32);
   mz = *reinterpret_cast<int *>(data+36);

   if (debug)
      std::cout << "debug:: slurp_parse_xmap_data(): start and range "
                << nx_start << " " << ny_start << " " << nz_start
                << " range " << mx << " " << my << " " << mz << std::endl;

   float cell_a = 0, cell_b = 0, cell_c = 0;
   float cell_al = 0, cell_be = 0, cell_ga = 0;
   int map_row = -1, map_col = -1, map_sec = -1;

   cell_a  = *reinterpret_cast<float *>(data+40);
   cell_b  = *reinterpret_cast<float *>(data+44);
   cell_c  = *reinterpret_cast<float *>(data+48);
   cell_al = *reinterpret_cast<float *>(data+52);
   cell_be = *reinterpret_cast<float *>(data+56);
   cell_ga = *reinterpret_cast<float *>(data+60);

   if (debug)
      std::cout << "DEBUG:: slurp_parse_xmap_data(): cell " << cell_a << " " << cell_b << " " << cell_c << " "
                << cell_al << " " << cell_be << " " << cell_ga << std::endl;

   // axis order
   map_row = *reinterpret_cast<int *>(data+64);
   map_col = *reinterpret_cast<int *>(data+68);
   map_sec = *reinterpret_cast<int *>(data+72);
   int axis_order_xyz[3];
   axis_order_xyz[map_row-1] = 0;
   axis_order_xyz[map_col-1] = 1;
   axis_order_xyz[map_sec-1] = 2;

   if (debug)
      std::cout << "DEBUG:: slurp_parse_xmap_data(): axis order "
                << map_row << " " << map_col << " " << map_sec << std::endl;

   // At the moment this function only works with simple X Y Z map ordering.
   // So escape with fail status if that is not the case
   bool is_simple_x_y_z_order = false; // initially
   if (map_row == 1)
      if (map_col == 2)
         if (map_sec == 3)
            is_simple_x_y_z_order = true;

   float dmin = 0.0, dmax = 0.0, dmean = 0.0;
   dmax  = *reinterpret_cast<float *>(data+76);
   dmin  = *reinterpret_cast<float *>(data+80);
   dmean = *reinterpret_cast<float *>(data+84);
   int space_group_number = -1;
   space_group_number = *reinterpret_cast<int *>(data+88);
   if (space_group_number == 0) // EM, maybe chimera maps
      space_group_number = 1;

   if (! is_simple_x_y_z_order) {
      if (debug)
         std::cout << "DEBUG:: Not simple X,Y,Z axis order - returning from slurp_parse_xmap_data() now "
                   << std::endl;
      status = slurp_map_result_t::UNRESOLVED;
      if (check_only) {
         auto is_small = [] (float v) { return (fabs(v) < 0.0001); };
         if (mx > 0 && my > 0 && mz > 0) {
            if (space_group_number == 1)
               if(is_small(cell_al-90.0f) && is_small(cell_be-90.0f) && is_small(cell_ga-90.0f))
                  status = slurp_map_result_t::IS_NON_SLURPABLE_EM_MAP;
         }
      }
      return status;
   }

   int size_extended_header = 0;
   size_extended_header = *reinterpret_cast<int *>(data+92);
   char *extra    = data+ 96; // 100 bytes max
   char *ext_type = data+104; // 1 byte
   char *version  = data+108; // 1 byte

   float origin_a = -1, origin_b = -1, origin_c = -1;
   origin_a = *reinterpret_cast<float *>(data+196);
   origin_b = *reinterpret_cast<float *>(data+200);
   origin_c = *reinterpret_cast<float *>(data+204);

   if (check_only) {
      auto is_small = [] (float v) { return (fabs(v) < 0.0001); };
      if (mx > 0 && my > 0 && mz > 0) {
         if (space_group_number == 1)
            if(is_small(cell_al-90.0f) && is_small(cell_be-90.0f) && is_small(cell_ga-90.0f))
               status = slurp_map_result_t::IS_SLURPABLE_EM_MAP;
      }
      return status;
   }

   if (debug)
      std::cout << "debug:: slurp:: origin " << origin_a << " " << origin_b << " " << origin_c << std::endl;

   char *type = data+208; // 4 bytes "MAP "
   int machine_stamp = -1;
   machine_stamp = *reinterpret_cast<int *>(data+212);
   float rmsd;
   rmsd = *reinterpret_cast<float *>(data+216);
   int number_of_labels;
   number_of_labels = *reinterpret_cast<int *>(data+220);
   if (number_of_labels > 10) number_of_labels = 10;
   char *labels = data+224;

   if (print_labels) {
      for(int i=0; i<number_of_labels; i++) {
         char *label = labels + i * 80;
         std::string s(label, 0, 80);
         std::cout << "   :" << s << ":" << std::endl;
      }
   }

   for(int i=0; i<number_of_labels; i++) {
      char *label = labels + i * 80;
      std::string s(label, 0, 80);
      // std::cout << "label length " << s.length() << std::endl;
      if (s.length() >= 8) {
         if (s.substr(0,8) == "PANDDA::") {
            std::cout << "debug:: slurp_parse_xmap_data() found a PANDDA::!" << std::endl;
            return slurp_map_result_t::NOT_AN_EM_MAP;
         } else {
            std::cout << "debug:: " << s << " does not start with PANDDA::" << std::endl;
            if (s.length() >= 16) {
               if (s.substr(0,16) == "APPLY-SYMMETRY::") {
                  return slurp_map_result_t::NOT_AN_EM_MAP;
               }
            }
         }
      }
   }

   char *map_data = data + size_extended_header + 1024; // points to the start of the grid (of 4 byte floats)

   int index_axis_order[3];
   index_axis_order[axis_order_xyz[0]] = mx;
   index_axis_order[axis_order_xyz[1]] = my;
   index_axis_order[axis_order_xyz[2]] = mz;

   clipper::Cell_descr cell_descr(cell_a, cell_b, cell_c, cell_al, cell_be, cell_ga);
   clipper::Cell cell(cell_descr);
   clipper::Spgr_descr sgd(space_group_number);
   clipper::Spacegroup space_group(sgd);
   clipper::Grid_sampling grid_sampling(index_axis_order[0], index_axis_order[1], index_axis_order[2]);
   clipper::Grid grid(mx, my, mz);
   std::cout << "DEBUG:: init xmap with " << space_group.symbol_hm() << " " << cell.format() << " " << grid_sampling.format()
             << std::endl;
   xmap_p->init(space_group, cell, grid_sampling);
   clipper::Coord_grid coord_grid_min(nx_start,ny_start,nz_start);
   clipper::Coord_grid coord_grid_max(index_axis_order[0]-1, index_axis_order[1]-1, index_axis_order[2]-1);

   clipper::Grid_range gr(coord_grid_min, coord_grid_max);

   clipper::Xmap<float>::Map_reference_coord i0 = clipper::Xmap<float>::Map_reference_coord(*xmap_p, gr.min());
   clipper::Coord_grid gr_max = gr.max();
   clipper::Xmap<float> &xmap = *xmap_p;

   if (debug) {
      std::cout << "debug:: slurp_map() grid " << grid.format() << std::endl;
      std::cout << "debug:: slurp_map() gr "   << gr.format() << std::endl;
      std::cout << "debug:: slurp_map() i0 "   << i0.coord().format() << " gr_max " << gr_max.format() << std::endl;
   }

   // mrc.set_coord() is slow - it can be multi-threaded.
   // much faster than using conventional map loading though.

   // 20240421-PE this crashes with tomogram emd_43330

   auto fill_map_sections_by_grid = [data_size] (std::pair<unsigned int, unsigned int> start_stop_section_index,
                                clipper::Xmap<float> *xmap,
                                int n_secs, int n_rows, int n_cols,
                                int nx_start, int ny_start, int nz_start,
                                int *axis_order_xyz,
                                const char *map_data,
                                std::atomic<bool> &print_lock) {

      // 20240421-PE Note to self: how about trying to get something like this to work? (c.f. mini-texture.cc)

      unsigned int section_start = start_stop_section_index.first;
      unsigned int section_stop  = start_stop_section_index.second;
      clipper::Grid_sampling gs = xmap->grid_sampling();
      clipper::Coord_grid cg_0(0,0, section_start);
      clipper::Coord_grid cg_1(gs.nu()-1, gs.nv()-1, section_stop-1);
      clipper::Grid_map grid(cg_0, cg_1);
      clipper::Xmap_base::Map_reference_coord ix(*xmap, grid.min()), iu, iv, iw;
      unsigned int offset = section_start * n_rows * n_cols;
      unsigned int c_u = 0;
      for ( iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) {
         unsigned int c_v = 0;
         for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
            unsigned int c_w = 0;
            for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
               unsigned int grid_offset = c_u * gs.nu() * gs.nv() + c_v * gs.nv() + c_w;
               const float f = *reinterpret_cast<const float *>(map_data + 4 * offset + 4 * grid_offset);

               if (true) { // debugging
                  bool unlocked = false;
                  while (! print_lock.compare_exchange_weak(unlocked, true)) {
                     std::this_thread::sleep_for(std::chrono::microseconds(1));
                     unlocked = false;
                  }
                  std::cout << " " << iw.coord().format() << " " << f << std::endl;
                  print_lock = false;
               }

               (*xmap)[iw] = f;
               c_w++;
            }
            c_v++;
         }
         c_u++;
      }
   };


   auto fill_map_sections = [data_size] (std::pair<unsigned int, unsigned int> start_stop_section_index,
                                clipper::Xmap<float> *xmap,
                                int n_secs, int n_rows, int n_cols,
                                int nx_start, int ny_start, int nz_start,
                                int *axis_order_xyz,
                                const char *map_data,
                                std::atomic<bool> &print_lock) {


                               int offset = start_stop_section_index.first * n_rows * n_cols;
                               int crs[3];  // col,row,sec coordinate
                               clipper::Xmap<float>::Map_reference_coord mrc(*xmap);
                               // this is unsigned int because split_indices returns unsigned int pairs.

                               if (false) { // debugging
                                  bool unlocked = false;
                                  while (! print_lock.compare_exchange_weak(unlocked, true)) {
                                     std::this_thread::sleep_for(std::chrono::microseconds(1));
                                     unlocked = false;
                                  }
                                  std::cout << "start stop sections: "
                                            << start_stop_section_index.first << " "
                                            << start_stop_section_index.second << " "
                                            << std::endl;
                                  print_lock = false;
                               }

                               for (unsigned int isec = start_stop_section_index.first; isec < start_stop_section_index.second; isec++ ) {

                                  crs[2] = isec + nz_start;
                                  for ( int irow = 0; irow < n_rows; irow++ ) {
                                     crs[1] = irow + ny_start;
                                     for ( int icol = 0; icol < n_cols; icol++ ) {
                                        crs[0] = icol + nx_start;
                                        mrc.set_coord(clipper::Coord_grid(crs[axis_order_xyz[0]],
                                                                          crs[axis_order_xyz[1]],
                                                                          crs[axis_order_xyz[2]]));
                                        if (data_size == 4) { // typical floating point data
                                           float f = *reinterpret_cast<const float *>(map_data + 4 * offset);
                                           if (f >  10000) f = 0;
                                           if (f < -10000) f = 0;
                                           (*xmap)[mrc] = f;
                                        }
                                        if (data_size == 1) {
                                           const char *c = reinterpret_cast<const char *>(map_data + offset);
                                           int ii = *c;
                                           (*xmap)[mrc] = static_cast<float>(ii);
                                        }
                                        offset++;
                                     }
                                  }
                               }

                               if (false) { // debugging
                                  bool unlocked = false;
                                  while (! print_lock.compare_exchange_weak(unlocked, true)) {
                                     std::this_thread::sleep_for(std::chrono::microseconds(1));
                                     unlocked = false;
                                  }
                                  std::cout << "DEBUG:: slurping: done " << start_stop_section_index.first << " "
                                            << start_stop_section_index.second << "\n";
                                  print_lock = false;
                               }
                            };


   bool single_thread_slurp = false;
   unsigned int n_threads = coot::get_max_number_of_threads();

   if (n_threads < 2) single_thread_slurp = true;

   if (single_thread_slurp) {
      int offset = 0;
      int crs[3];  // col,row,sec coordinate
      clipper::Xmap<float>::Map_reference_coord mrc(xmap);
      std::cout << "info:: n_sections: " << n_secs  << std::endl;
      for ( int isec = 0; isec < n_secs; isec++ ) {
         crs[2] = isec + nz_start;
         for ( int irow = 0; irow < n_rows; irow++ ) {
            crs[1] = irow + ny_start;
            for ( int icol = 0; icol < n_cols; icol++ ) {
               crs[0] = icol + nx_start;
               mrc.set_coord( clipper::Coord_grid( crs[axis_order_xyz[0]], crs[axis_order_xyz[1]], crs[axis_order_xyz[2]] ) );
               float f = *reinterpret_cast<float *>(map_data + 4 * offset);
               if (f >  10000) f = 0;
               if (f < -10000) f = 0;
               xmap[mrc] = f;
               offset++;
            }
         }
      }
      status = slurp_map_result_t::OK;
   } else {
      std::vector<std::pair<unsigned int, unsigned int> > airs = atom_index_ranges(n_secs, n_threads);
      std::vector<std::thread> threads;
      std::atomic<bool> print_lock(false);
      try {
         for (auto air : airs) {
            if (debug)
               std::cout << "DEBUG:: thread fill sections " << air.first << " to " << air.second << std::endl;
            threads.push_back(std::thread(fill_map_sections, air, &xmap, n_secs, n_rows, n_cols,
                                          nx_start, ny_start, nz_start,
                                          axis_order_xyz, map_data, std::ref(print_lock)));
         }
         for (std::size_t i=0; i<airs.size(); i++)
            threads[i].join();
         if (debug)
            std::cout << "DEBUG:: thread join finished" << std::endl;
         status = slurp_map_result_t::OK;
      }
      catch (const std::system_error &e) {
         std::cout << "ERROR:: std::system_error: " << e.what() << std::endl;
      }
   }

   if (debug)
      std::cout << "DEBUG:: coot::util::slurp_parse_xmap_data(): returning status " << int(status) << std::endl;

   return status;
}

// PANDDA::, that is.
bool
coot::util::map_labels_contain_PANDDA(const std::string &file_name) {

   bool status = false;
   std::vector<std::string> v = get_map_labels(file_name);
   for (const auto &s : v) {
      if (s.length() >= 8) {
         if (s.substr(0,8) == "PANDDA::") {
            status = true;
            break;
         }
      }
   }
   return status;

}

std::vector<std::string>
coot::util::get_map_labels(const std::string &file_name) {

   std::vector<std::string> v;

   auto slurp_parse_get_labels = [] (char *data) {

      std::vector<std::string> v;
      bool debug = false;

      int n_rows = -1;
      int n_cols = -1;
      int n_secs = -1;
      n_cols = *reinterpret_cast<int *>(data);
      n_rows = *reinterpret_cast<int *>(data+4);
      n_secs = *reinterpret_cast<int *>(data+8);
      if (debug) {
         std::cout << "debug:: n_cols " << n_cols << std::endl;
         std::cout << "debug:: n_rows " << n_rows << std::endl;
         std::cout << "debug:: n_sections " << n_secs << std::endl;
      }

      int mode = -1;
      int data_size = 4; // 4 bytes for a float
      mode = *reinterpret_cast<int *>(data+12);
      if (mode == 0) data_size = 1;
      if (mode == 1) data_size = 2;
      if (mode == 6) data_size = 2; // two-byte integer. Needs fixing in the casting.
      if (debug)
         std::cout << "debug:: slurp_parse_get_labels() mode: " << mode << std::endl;

      int nx_start = -1, ny_start = -1, nz_start = -1;
      int mx = -1, my = -1, mz = -1;
      nx_start = *reinterpret_cast<int *>(data+16);
      ny_start = *reinterpret_cast<int *>(data+20);
      nz_start = *reinterpret_cast<int *>(data+24);
      mx = *reinterpret_cast<int *>(data+28);
      my = *reinterpret_cast<int *>(data+32);
      mz = *reinterpret_cast<int *>(data+36);

      if (debug)
         std::cout << "debug:: start and range " << nx_start << " " << ny_start << " " << nz_start
                   << " range " << mx << " " << my << " " << mz << std::endl;
      float cell_a = 0, cell_b = 0, cell_c = 0;
      float cell_al = 0, cell_be = 0, cell_ga = 0;
      int map_row = -1, map_col = -1, map_sec = -1;

      cell_a  = *reinterpret_cast<float *>(data+40);
      cell_b  = *reinterpret_cast<float *>(data+44);
      cell_c  = *reinterpret_cast<float *>(data+48);
      cell_al = *reinterpret_cast<float *>(data+52);
      cell_be = *reinterpret_cast<float *>(data+56);
      cell_ga = *reinterpret_cast<float *>(data+60);

      if (debug)
         std::cout << "debug:: cell " << cell_a << " " << cell_b << " " << cell_c << " "
                   << cell_al << " " << cell_be << " " << cell_ga << std::endl;

      // axis order
      map_row = *reinterpret_cast<int *>(data+64);
      map_col = *reinterpret_cast<int *>(data+68);
      map_sec = *reinterpret_cast<int *>(data+72);
      int axis_order_xyz[3];
      axis_order_xyz[map_row-1] = 0;
      axis_order_xyz[map_col-1] = 1;
      axis_order_xyz[map_sec-1] = 2;

      if (debug)
         std::cout << "debug:: axis order " << map_row << " " << map_col << " " << map_sec << std::endl;

      // At the moment this function only works with simple X Y Z map ordering.
      // So escape with fail status if that is not the case
      bool is_simple_x_y_z_order = false; // initially
      if (map_row == 1)
         if (map_col == 2)
            if (map_sec == 3)
               is_simple_x_y_z_order = true;

      float dmin = 0.0, dmax = 0.0, dmean = 0.0;
      dmax  = *reinterpret_cast<float *>(data+76);
      dmin  = *reinterpret_cast<float *>(data+80);
      dmean = *reinterpret_cast<float *>(data+84);
      int space_group_number = -1;
      space_group_number = *reinterpret_cast<int *>(data+88);
      if (space_group_number == 0) // EM, maybe chimera maps
         space_group_number = 1;
      int size_extended_header = 0;
      size_extended_header = *reinterpret_cast<int *>(data+92);
      char *extra = data+96; // 100 bytes max
      char *ext_type = data+104; // 1 byte
      char *version = data+108;  // 1 byte

      float origin_a = -1, origin_b = -1, origin_c = -1;
      origin_a = *reinterpret_cast<float *>(data+196);
      origin_b = *reinterpret_cast<float *>(data+200);
      origin_c = *reinterpret_cast<float *>(data+204);

      if (debug)
         std::cout << "origin " << origin_a << " " << origin_b << " " << origin_c << std::endl;

      char *type = data+208; // 4 bytes "MAP "
      int machine_stamp = *reinterpret_cast<int *>(data+212);
      float rmsd = *reinterpret_cast<float *>(data+216);
      int number_of_labels = *reinterpret_cast<int *>(data+220);
      if (number_of_labels > 10) number_of_labels = 10;
      char *labels = data+224;

      for(int i=0; i<number_of_labels; i++) {
         char *label = labels + i * 80;
         std::string s(label, 0, 80);
         v.push_back(s);
      }
      return v;
   };

   struct stat s;
   int fstat = stat(file_name.c_str(), &s);
   if (fstat == 0) {
      FILE *fptr = fopen(file_name.c_str(), "rb");
      off_t st_size = s.st_size;
      std::cout << "get_map_labels(): st_size: " << st_size << std::endl;
      try {
         // 20231006-PE as it used to be.
         char *space = new char[st_size+1];
         // Happy Path
         size_t st_size_2 = fread(space, st_size, 1, fptr);
         char *data = static_cast<char *>(space);
         fclose(fptr);
         if (st_size_2 == 1) {
            // Happy Path
            if (st_size > 1024) {
               v = slurp_parse_get_labels(data);
            }
         }
      }
      catch (const std::bad_alloc &e) {
         std::cout << "WARNING:: out-of-memory " << st_size+1 << " " << e.what() << std::endl;
      }
   }
   return v;
}



#if 0

int main(int argc, char **argv) {

   int status = 0;

   if (argc > 1) {

      std::string file_name(argv[1]);

      try {
         auto tp_0 = std::chrono::high_resolution_clock::now();
         clipper::CCP4MAPfile file;
         clipper::Xmap<float> xmap;
         file.open_read(file_name);
         file.import_xmap(xmap);
         file.close_read();

         auto tp_1 = std::chrono::high_resolution_clock::now();
         auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
         std::cout << "Clipper read of map: " << d10 << " milliseconds" << std::endl;

         struct stat s;
         int fstat = stat(file_name.c_str(), &s);
         if (fstat == 0) {
            // Happy Path
            int st_size = s.st_size;
            std::cout << "size " << st_size << " bytes" << std::endl;
            auto tp_2 = std::chrono::high_resolution_clock::now();
            FILE *fptr = fopen(file_name.c_str(), "rb");
            void *space = new void [st_size+1];
            auto tp_3 = std::chrono::high_resolution_clock::now();
            size_t st_size_2 = fread(space, st_size, 1, fptr);
            char *data = static_cast<char *>(space);
            fclose(fptr);
            auto tp_4 = std::chrono::high_resolution_clock::now();
            // now act on st_size bytes of space
            // std::cout << "st_size " << st_size << " st_size_2 " << st_size_2 << std::endl;
            if (st_size_2 == 1) {
               // Happy Path
               if (st_size > 1024) {
                  clipper::Xmap<float> xmap;
                  slurp_parse_xmap_data(data, &xmap); // fill xmap
                  auto tp_5 = std::chrono::high_resolution_clock::now();
                  clipper::CCP4MAPfile out_file;
                  out_file.open_write("slurp-parsed.map");
                  out_file.export_xmap(xmap);
                  out_file.close_write();
                  auto tp_6 = std::chrono::high_resolution_clock::now();
                  auto d32 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_3 - tp_2).count();
                  auto d43 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_4 - tp_3).count();
                  auto d54 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_5 - tp_4).count();
                  auto d65 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_6 - tp_5).count();
                  std::cout << "malloc() for map:  " << d32 << " milliseconds" << std::endl;
                  std::cout << "Slurp read of map: " << d43 << " milliseconds" << std::endl;
                  std::cout << "Xmap construction: " << d54 << " milliseconds" << std::endl;
                  std::cout << "Xmap write: " << d65 << " milliseconds" << std::endl;
               } else {
                  std::cout << "File " << file_name << " was not a map" << std::endl;
               }
            }
            delete [] space;
         } else {
            std::cout << "stat() failed for " << file_name << std::endl;
         }

      }
      catch (const clipper::Message_base &exc) {
         std::cout << "WARNING:: failed to open " << file_name << std::endl;
      }

   }

   return status;

}

#endif
