
SURF_SOURCES = Vectors.cpp CIsoSurface.cpp

COORD_INCLUDES = -I$(top_srcdir)/coords

INCLUDES = $(COORD_INCLUDES) $(CLIPPER_CXXFLAGS) 


lib_LTLIBRARIES = libcoot-isosurface.la

include_HEADERS = Vectors.h CIsoSurface.h

libcoot_isosurface_la_SOURCES = $(SURF_SOURCES)
