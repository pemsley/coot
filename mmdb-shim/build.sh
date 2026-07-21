#!/usr/bin/env bash
# Build + run the mmdb-shim tests.
#   macro ON  -> gemmi-backed shim (mmdb-shim/include/mmdb2 + _shim_impl.hh)
#   macro OFF -> falls through to real MMDB via #include_next (coexistence check)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
GEMMI="${GEMMI_INC:-$(brew --prefix gemmi)/include}"
MMDB="${MMDB_INC:-/opt/homebrew/Cellar/mmdb2/2.0.22/include}"
CXX=/usr/bin/clang++
STD="-std=c++17 -O0 -g -Wall"

echo "=== core unit test ==="
$CXX $STD -I"$GEMMI" "$HERE/core/test_core.cc" -o "$HERE/core/test_core" && "$HERE/core/test_core" | tail -1

echo "=== shim test (COOT_USE_MMDB_SHIM) ==="
$CXX $STD -DCOOT_USE_MMDB_SHIM -I"$HERE/include" -I"$GEMMI" \
  "$HERE/test/test_shim.cc" -o "$HERE/test/test_shim" && "$HERE/test/test_shim" | tail -1

echo "=== UDData engine test ==="
$CXX $STD -DCOOT_USE_MMDB_SHIM -I"$HERE/include" -I"$GEMMI" \
  "$HERE/test/test_udd.cc" -o "$HERE/test/test_udd" && "$HERE/test/test_udd" | tail -1

echo "=== selection engine test ==="
$CXX $STD -DCOOT_USE_MMDB_SHIM -I"$HERE/include" -I"$GEMMI" \
  "$HERE/test/test_sel.cc" -o "$HERE/test/test_sel" && "$HERE/test/test_sel" | tail -1

echo "=== contacts test (SelectSphere + SeekContacts via gemmi NeighborSearch) ==="
$CXX $STD -DCOOT_USE_MMDB_SHIM -I"$HERE/include" -I"$GEMMI" \
  "$HERE/src/contacts.cc" "$HERE/test/test_contacts.cc" \
  -L"$(brew --prefix gemmi)/lib" -lgemmi_cpp -lz -o "$HERE/test/test_contacts" \
  && "$HERE/test/test_contacts" | tail -1

echo "=== I/O round-trip test (gemmi read/write; links libgemmi_cpp) ==="
PDB="${PDB:-$HERE/../1hr2_final.pdb}"
$CXX $STD -DCOOT_USE_MMDB_SHIM -I"$HERE/include" -I"$GEMMI" \
  "$HERE/src/io.cc" "$HERE/test/test_io.cc" \
  -L"$(brew --prefix gemmi)/lib" -lgemmi_cpp -lz -o "$HERE/test/test_io" \
  && "$HERE/test/test_io" "$PDB" | tail -1

echo "=== leaf integration test (real Coot mmdb idioms vs gemmi ground truth) ==="
$CXX $STD -DCOOT_USE_MMDB_SHIM -I"$HERE/include" -I"$GEMMI" \
  "$HERE/src/io.cc" "$HERE/test/test_leaf.cc" \
  -L"$(brew --prefix gemmi)/lib" -lgemmi_cpp -lz -o "$HERE/test/test_leaf" \
  && "$HERE/test/test_leaf" | tail -1

if [ -x "$HERE/../mmdb-recon/ast/build2/mmdb_tool" ]; then
  echo "=== field-access rewriter test (synthetic) ==="
  bash "$HERE/test/rewrite/run_rewrite_test.sh" | tail -1
fi

echo "=== coexistence: macro OFF must compile against real MMDB ==="
printf '#include <mmdb2/mmdb_manager.h>\nint main(){mmdb::Manager*m=new mmdb::Manager();int n=m->GetNumberOfModels();delete m;return n;}\n' > /tmp/mmdb_fallback.cc
$CXX $STD -I"$HERE/include" -I"$MMDB" -c /tmp/mmdb_fallback.cc -o /tmp/mmdb_fallback.o \
  && echo "  ok: fallback compiles against real MMDB"
