/* lidia-core/bond-record-container-t.hh
 * 
 * Copyright 2016 by Medical Research Council
 * Author: Paul Emsley
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

#ifdef USE_SQLITE3
#include <sqlite3.h>
#endif // USE_SQLITE3

#include "cod-atom-type-t.hh"
#include "bond-table-record-t.hh"
// for validation
#include "geometry/protein-geometry.hh"

namespace cod {

      // --------------- record container ---------------------

   class bond_record_container_t {
      std::string::size_type get_max_atom_type_width() const;
      std::map<std::string, unsigned int> level_4_atom_types_map; // convert from acedrg-table 
                                                                  // COD string to its index
      bool write_atom_type_indices(const std::string &file_name) const;

      void sort() { std::sort(bonds.begin(), bonds.end()); }
      void fill_level_4_atom_map();
      void fill_bonds_map();

      void fill_cod_atom_type_map();
      std::set<atom_type_t> cod_atom_type_set;
      std::map<atom_type_t, int> cod_atom_type_map; // for extraction of atom type indices

      std::vector<std::string>
      read_atom_type_indices(const std::string &atom_type_indices_file_name) const;
      
      bool read_bonds(const std::string &bonds_file_name,
		      const std::vector<std::string> &types);

      // can throw std::runtime_error
      unsigned int get_atom_index(const std::string &at_name_1,
				  const coot::dictionary_residue_restraints_t &rest) const;
      unsigned int get_atom_index(const std::string &at_name_1,
				  const RDKit::RWMol &mol) const;

      std::vector<bool> get_is_hydrogen_flags(const RDKit::RWMol &rdkm) const;
      
      double get_bond_distance_from_model(const std::string &at_name_1,
					  const std::string &at_name_2,
					  mmdb::Residue *res) const;

      bond_table_record_t
      get_cod_bond_from_table(const atom_type_t &cod_type_1,
			      const atom_type_t &cod_type_2) const;
      
      // which uses
      bond_table_record_t
      make_bond_from_level_3_vector(const atom_type_t &cod_type_1,
				    const atom_type_t &cod_type_2,
				    const std::vector<bond_table_record_t> &v,
				    bond_table_record_t::approximation_level_t al) const;

      // generalization
      //
      bond_table_record_t
      make_bond_from_level_2_map(const atom_type_t &cod_type_1,
				 const atom_type_t &cod_type_2,
				 const std::map<std::string, std::map<std::string, std::vector<bond_table_record_t> > > &l3_map,
				 bond_table_record_t::approximation_level_t approx_level) const;
      

      bond_table_record_t
      consolidate_bonds(const atom_type_t &cod_type_1,
			const atom_type_t &cod_type_2,
			const std::vector<bond_table_record_t> &lb,
			bond_table_record_t::approximation_level_t approx_level) const;

      void t3_miss_diagnose(const atom_type_t &cod_type_1,
			    const atom_type_t &cod_type_2) const;
      
#ifdef USE_SQLITE3
      static int db_callback(void *NotUsed, int argc, char **argv, char **azColName);
      sqlite3 *make_sqlite_db(const std::string &file_name);
      bool db_add_level_4_types(sqlite3 *db);
#endif // USE_SQLITE3

      clipper::Coord_orth co(mmdb::Atom *at) const { return clipper::Coord_orth(at->x, at->y, at->z); }

   public:
      bond_record_container_t() {}
      bond_record_container_t(const std::string &table_dir_name) {
	 read_acedrg_table_dir(table_dir_name);
      }
      bool read_acedrg_table(const std::string &file_name);

      // fill the bonds table
      // 
      void read_acedrg_table_dir(const std::string &table_dir_name);

      std::vector<bond_table_record_t> bonds;

      // Level-4 atom type is the (outer) key and (outer) value
      // is a map for which the key is the second level-4 atom type,
      // returning a bond record.
      // std::map<std::string, std::map<std::string, bond_table_record_t> > bonds_map;

      // the key is a level-3 atom type, the value is
      // a map for which the key is a level-3 atom type and the value of that
      // map is a vector of bond_table_record_t.
      // 
      // std::map<std::string, std::map<std::string, std::vector<bond_table_record_t> > >
      // bonds_map;

      // same again, but with a level-2 type to index first
      //
      // bonds_map[l2][l2][l3][l3]
      // 
      std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, std::vector<bond_table_record_t> > > > > bonds_map;

      void add(const bond_table_record_t &rec) {
	 bonds.push_back(rec);
      }

      void add_table(const bond_record_container_t &brc);

      unsigned int size() { return bonds.size(); }
      bool write(const std::string &atom_indices_file_name,
		 const std::string &bonds_file_name) const;

      bool read(const std::string &atom_type_indices_file_name,
		const std::string &bonds_file_name);
      void check() const;
      void validate(mmdb::Residue *res, const coot::dictionary_residue_restraints_t &rest) const;
      void make_db(const std::string &file_name); // probably not const
   };
}
