

lib_LTLIBRARIES=libcoot-skeleton.la

INCLUDES =-I$(top_srcdir)/coords -I$(top_srcdir)/angles \
	 -I$(top_srcdir)/atom-utils  $(MMDB_CXXFLAGS) $(CLIPPER_CXXFLAGS) 

include_HEADERS = graphical_skel.h dm-utils.h BuildCas.h

libskeleton_a_sources =  \
	dm-utils.cc \
	graphical_skel.cc \
	BuildCas.cc

# No more fc_skeleton
#	fc_skeleton.cc fc_skeleton.h


libcoot_skeleton_la_SOURCES = $(libskeleton_a_sources)


# bin_PROGRAMS = ncssearch test_various test_torsion test_theta_2

# ncssearch_SOURCES = ncs-search.cc ncs-search.h

# test_various_SOURCES = test_various.cc 

# test_torsion_SOURCES = test-torsion.cc

# test_theta_2_SOURCES = test_theta_2.cc

# ncssearch_LDADD = -L. -lcoot-skeleton \
# 	-L../coords -lcoot-coords	\
# 	-L../angles -lcoot-AngleInfo \
# 	$(CLIPPER_LIBS) 	\
# 	$(MMDB_LIBS)

# test_torsion_LDADD = -L. -lcoot-skeleton 	\
# 		-L../coords -lcoot-coords	\
# 		-L../angles -lcoot-AngleInfo \
# 		$(CLIPPER_LIBS) 	\
# 		$(MMDB_LIBS)

# test_theta_2_LDADD = -L../angles -lcoot-AngleInfo

# sktest_SOURCES = sktest.cpp
# sktest_LDADD = -L. -lskeleton $(CLIPPER_LIBS) 

# test_various_LDADD = -L. -lcoot-skeleton \
# 	-L../coords -lcoot-coords	\
# 	-L../angles -lcoot-AngleInfo \
# 	$(CLIPPER_LIBS) 	\
# 	$(MMDB_LIBS)

