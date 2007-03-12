
lib_LTLIBRARIES=libccp4mg-utils.la

include_HEADERS = geomutil.h mgtree.h cartesian.h matrix.h quat.h ppmutil.h

EXTRA_DIST = convert.sh get-new-utils.sh

INCLUDES =  \
	$(MMDB_CXXFLAGS)

# for ppmutil.cc compilation
# CXXFLAGS = -DUSE_LIBPNG -I/sw/include


libccp4mg_utils_la_SOURCES = geomutil.cc mgtree.cc cartesian.cc matrix.cc quat.cc \
	ppmutil.cc

