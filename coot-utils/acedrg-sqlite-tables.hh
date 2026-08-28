#ifndef COOT_UTILS_ACEDRG_SQLITE_TABLES_HH
#define COOT_UTILS_ACEDRG_SQLITE_TABLES_HH

#include <memory>
#include <string>
#include <utility>

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

      //! open the tables in the given data directory. Returns success.
      bool init(const std::string &acedrg_data_dir);
      //! init() on default_data_dir()
      bool init();
      bool is_usable() const;

      //! prefetch this molecule's rows from SQLite, then run gemmi's
      //! fill_restraints() on cc. Returns false if not usable.
      bool fill_chemcomp(gemmi::ChemComp &cc);

      //! Return (status, restraints) where the restraints contain only the
      //! AceDRG-derived bond and angle restraints for the molecule described
      //! by restraints_in (atoms, bond connectivity and bond-type strings
      //! are read from restraints_in; values are ignored). Suitable for
      //! dictionary_residue_restraints_t::conservatively_replace_with().
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
