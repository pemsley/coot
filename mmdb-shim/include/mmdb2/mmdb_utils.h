// mmdb-shim public header. When COOT_USE_MMDB_SHIM is defined, resolves to the
// gemmi-backed shim; otherwise falls through to the real MMDB header via
// #include_next (requires this dir to precede real mmdb2 on the include path).
#ifndef COOT_MMDB_SHIM_mmdb_utils_H
#define COOT_MMDB_SHIM_mmdb_utils_H
#  ifdef COOT_USE_MMDB_SHIM
#    include "_shim_impl.hh"
#  else
#    include_next <mmdb2/mmdb_utils.h>
#  endif
#endif
