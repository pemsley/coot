/*
 * coot-utils/make-acedrg-sqlite.cc
 *
 * Copyright 2026 by Global Phasing Ltd.
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "acedrg-sqlite-tables.hh"

int main(int argc, char **argv) {

   bool if_missing = false;
   std::vector<std::string> args;
   for (int i=1; i<argc; i++) {
      std::string a(argv[i]);
      if (a == "--if-missing") if_missing = true;
      else args.push_back(a);
   }
   if (args.empty()) {
      std::cout << "Usage: coot-make-acedrg-sqlite [--if-missing] "
                << "<acedrg-ascii-tables-dir> [<output-data-dir>]\n"
                << "   default output: " << coot::acedrg_sqlite_tables::default_data_dir()
                << std::endl;
      return if_missing ? 0 : 2;
   }
   std::string tables_dir = args[0];
   std::string out_dir = (args.size() > 1) ? args[1] : coot::acedrg_sqlite_tables::default_data_dir();

   if (if_missing) {
      std::filesystem::path db = std::filesystem::path(out_dir) / "acedrg.sqlite";
      if (std::filesystem::exists(db)) {
         std::cout << "coot-make-acedrg-sqlite: " << db.string() << " exists - nothing to do" << std::endl;
         return 0;
      }
   }
   bool status = false;
   try {
      status = coot::acedrg_sqlite_tables::build(tables_dir, out_dir);
   }
   catch (const std::exception &e) {
      std::cout << "coot-make-acedrg-sqlite: " << e.what() << std::endl;
   }
   if (status)
      std::cout << "coot-make-acedrg-sqlite: wrote " << out_dir << "/acedrg.sqlite" << std::endl;
   // with --if-missing a failure is tolerated (build machines without the
   // tables, or without write access): exit 0 so "make" carries on.
   return status ? 0 : (if_missing ? 0 : 1);
}
