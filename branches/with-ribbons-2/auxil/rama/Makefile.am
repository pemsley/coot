
EXTRA_DIST = Makefile.simple fit-rama.scm run.scm test-arr.scm

INCLUDES =  $(CLIPPER_CXXFLAGS) $(MMDB_CXXFLAGS)


check_PROGRAMS = rama-data rama-doughnut 

rama_data_SOURCES = rama-data.cc
rama_data_LDADD =  $(MMDB_LIBS)	$(CLIPPER_LIBS) 

rama_doughnut_SOURCES = rama-doughnut.cc
rama_doughnut_LDADD =  $(MMDB_LIBS) $(CLIPPER_LIBS) 
