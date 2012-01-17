
EXTRA_DIST = Makefile.simple 
INCLUDES = $(MMDB_CXXFLAGS) $(CLIPPER_CXXFLAGS)
check_PROGRAMS = diff-vec

diff_vec_SOURCES = diff-vec.cc
diff_vec_LDADD =        \
    $(CLIPPER_LIBS) 	\
    $(MMDB_LIBS)
