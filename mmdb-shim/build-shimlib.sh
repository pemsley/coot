#!/usr/bin/env bash
# Build the shim's compiled units (io.cc + contacts.cc) into a static lib that
# the real Coot build links against. Header-only parts come via -Iinclude.
#   -> mmdb-shim/lib/libmmdbshim.a
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
GEMMI="${GEMMI_INC:-$(brew --prefix gemmi)/include}"
CXX="${CXX:-/usr/bin/clang++}"
mkdir -p "$HERE/lib" "$HERE/obj"

for u in io contacts; do
  "$CXX" -std=c++17 -O2 -fPIC -DCOOT_USE_MMDB_SHIM -I"$HERE/include" -I"$GEMMI" \
    -c "$HERE/src/$u.cc" -o "$HERE/obj/$u.o"
done
ar rcs "$HERE/lib/libmmdbshim.a" "$HERE/obj/io.o" "$HERE/obj/contacts.o"
echo "built $HERE/lib/libmmdbshim.a"
echo
echo "To wire into the autotools build (before ./configure), e.g.:"
echo "  export CPPFLAGS=\"-I$HERE/include -DCOOT_USE_MMDB_SHIM \$CPPFLAGS\""
echo "  export LDFLAGS=\"-L$HERE/lib -L\$(brew --prefix gemmi)/lib \$LDFLAGS\""
echo "  export LIBS=\"-lmmdbshim -lgemmi_cpp \$LIBS\""
echo "  # (apply the field rewrite on a branch first; note clipper/ssm still"
echo "  #  see the shim headers in THIS version — expected, per project notes)"
