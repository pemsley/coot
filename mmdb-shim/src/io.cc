// mmdb-shim I/O — Manager read/write via gemmi (architecture B: gemmi is the
// live store, so read = parse + build wrapper tree; write = serialize st).
// Compiled once (links libgemmi_cpp) so the heavy gemmi write/read headers stay
// out of the many Coot TUs that include <mmdb2/mmdb_manager.h>.
#define COOT_USE_MMDB_SHIM 1
#include <mmdb2/mmdb_manager.h>

#include <gemmi/mmread.hpp>   // read_structure_file (auto-detect)
#include <gemmi/pdb.hpp>      // read_pdb_file
#include <gemmi/to_pdb.hpp>   // write_pdb
#include <gemmi/to_mmcif.hpp> // make_mmcif_document
#include <gemmi/to_cif.hpp>   // write_cif_to_stream

#include <fstream>

namespace mmdb {

// Rebuild the wrapper tree from a freshly loaded gemmi::Structure. We do NOT run
// gemmi's setup_entities()/subchain splitting — MMDB parity wants the raw chains
// as they appear in the file.
ERROR_CODE Manager::ReadPDBASCII(cpstr fname) {
  try {
    st = gemmi::read_pdb_file(fname);
  } catch (const std::exception &) {
    return Error_CantOpenFile;
  }
  build_from_gemmi();
  return Error_NoError;
}

ERROR_CODE Manager::ReadCoorFile(cpstr fname) {
  try {
    st = gemmi::read_structure_file(fname);
  } catch (const std::exception &) {
    return Error_CantOpenFile;
  }
  build_from_gemmi();
  return Error_NoError;
}

ERROR_CODE Manager::WritePDBASCII(cpstr fname) {
  std::ofstream os(fname);
  if (!os) return Error_CantOpenFile;
  gemmi::write_pdb(st, os);
  return Error_NoError;
}

ERROR_CODE Manager::WriteCIFASCII(cpstr fname) {
  std::ofstream os(fname);
  if (!os) return Error_CantOpenFile;
  gemmi::cif::Document doc = gemmi::make_mmcif_document(st);
  gemmi::cif::write_cif_to_stream(os, doc);
  return Error_NoError;
}

} // namespace mmdb
