
EXTRA_DIST = Makefile.simple 


INCLUDES = -I$(top_srcdir)/auxil/sec-str \
	$(MMDB_CXXFLAGS) 	         \
	$(CLIPPER_CXXFLAGS) 		


check_PROGRAMS = strand-analysis

strand_analysis_LDADD = \
    $(CLIPPER_LIBS) 	\
    $(MMDB_LIBS)


strand_analysis_SOURCES = strand-analysis.cc
