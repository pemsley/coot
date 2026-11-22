#ifndef COOT_UTILS_ACEDRG_SQLITE_TABLES_HH
#define COOT_UTILS_ACEDRG_SQLITE_TABLES_HH

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "geometry/protein-geometry.hh"

namespace gemmi { struct ChemComp; }

namespace coot {

   //! \brief AceDRG bond/angle ideal-value lookups from a SQLite conversion
   //! of the AceDRG tables, using gemmi's AcedrgTables as the engine.
   //!
   //! See 2026-08-14-acedrg-sqlite-tables-in-coot-design.md. Usable only
   //! when the data directory (default: XDG cache dir Coot/acedrg-tables)
   //! holds the cheap ASCII tables and acedrg.sqlite, as produced by
   //! build() / the coot-make-acedrg-sqlite tool.
   class acedrg_sqlite_tables {
   public:
      acedrg_sqlite_tables();
      ~acedrg_sqlite_tables();
      acedrg_sqlite_tables(const acedrg_sqlite_tables &) = delete;
      acedrg_sqlite_tables& operator=(const acedrg_sqlite_tables &) = delete;

      //! xdg_t::get_cache_home()/acedrg-tables
      static std::string default_data_dir();

      //! @return false if this build has only the stubs (no SQLite3, or
      //!         a gemmi without AcedrgTables::insert_bond_row())
      static bool compiled_with_support();

      //! open the tables in the given data directory. Returns success.
      bool init(const std::string &acedrg_data_dir);
      //! init() on default_data_dir()
      bool init();
      bool is_usable() const;

      //! prefetch this molecule's rows from SQLite, then run gemmi's
      //! fill_restraints() on cc. Returns false if not usable.
      //!
      //! @param acedrg_atom_types if not null, set to the AceDRG atom type (the
      //!        COD class, e.g. "C[5](N[5]C[5]C[6])(C[5]C[5]HH)(O){1|C<3>,2|C<4>,4|H<1>}")
      //!        of each atom of cc, in the order of cc.atoms
      bool fill_chemcomp(gemmi::ChemComp &cc, std::vector<std::string> *acedrg_atom_types = nullptr);

      //! Return (status, restraints) where the restraints contain only the
      //! AceDRG-derived bond and angle restraints for the molecule described
      //! by restraints_in (atoms, bond connectivity and bond-type strings
      //! are read from restraints_in; values are ignored). Suitable for
      //! dictionary_residue_restraints_t::conservatively_replace_with().
      //! The atom_info of the restraints is a copy of that of restraints_in,
      //! with the AceDRG atom types (acedrg_atom_type) set.
      std::pair<bool, dictionary_residue_restraints_t>
      make_bond_and_angle_restraints(const dictionary_residue_restraints_t &restraints_in);

      //! Convert the AceDRG ASCII tables: copy the cheap table files into
      //! output_data_dir and build output_data_dir/acedrg.sqlite from the
      //! heavy bond/angle tables. Returns success.
      static bool build(const std::string &acedrg_ascii_tables_dir,
                        const std::string &output_data_dir);

   private:
      class impl;
      std::unique_ptr<impl> pimpl;
   };
}

#endif // COOT_UTILS_ACEDRG_SQLITE_TABLES_HH
