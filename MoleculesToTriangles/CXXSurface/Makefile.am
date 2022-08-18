
ACLOCAL_AMFLAGS="-Im4"
lib_LTLIBRARIES=libMoleculesToTrianglesCXXSurface.la

libMoleculesToTrianglesCXXSurface_la_SOURCES=CXXCircle.cpp CXXCircleNode.cpp \
        CXXCoord.cpp CXXFortranFile.cpp CXXNewHood.cpp CXXSphereElement.cpp \
	CXXSphereNode.cpp CXXSphereTriangle.cpp CXXSphereTriangleEdge.cpp CXXSurface.cpp CXXSurfaceVertex.cpp \
	CXXTorusElement.cpp CXXTorusNode.cpp CXXTriangle.cpp CXXQADSurface.cpp  CXXQADSurfaceStaticData.cpp \
	CXXCreator.cpp CXXSpace.cpp CXXFFTSolventMap.cpp CXXChargeTable.cpp CXXException.cpp \
	CXXUtils.cpp CXXBall.cpp CXXSurfaceMaker.cpp

libMoleculesToTrianglesCXXSurface_la_CPPFLAGS=${CPPFLAGS} -I@MMDB_INCLUDE_DIR@ -I@CLIPPER_INCLUDE_DIR@ $(OPENMP_CFLAGS)
# AC_OPENMP needs to be used when linking programs:
# https://www.gnu.org/software/autoconf/manual/autoconf-2.67/html_node/Generic-Compiler-Characteristics.html
# libMoleculesToTrianglesCXXSurface_la_LIBADD=-lrfftw -lfftw -lmmdb2 -lccp4c -lclipper-core -lclipper-ccp4 -lgomp
libMoleculesToTrianglesCXXSurface_la_LIBADD=$(FFTW2_LIBS) -lmmdb2 -lccp4c -lclipper-core -lclipper-ccp4 $(OPENMP_CFLAGS)

# was this: otherincludedir = $(includedir)/MoleculesToTriangles/CXXSurface
otherincludedir = $(top_srcdir)/MoleculesToTriangles/CXXSurface
otherinclude_HEADERS = CXXSurfaceMaker.h CXXSurface.h CXXCoord.h CXXSurfaceVertex.h \
      TokenIterator.h CXXTriangle.h CXXUtils.h CXXChargeTable.h CXXCreator.h CXXException.h \
      CXXSpace.h CXXFFTSolventMap.h

bin_PROGRAMS=ElectroSurface
ElectroSurface_SOURCES=ElectroSurface.cpp
ElectroSurface_LDADD=libMoleculesToTrianglesCXXSurface.la $(CLIPPER_LIBS) $(MMDB_LIBS)
ElectroSurface_CPPFLAGS=${CPPFLAGS} -std=c++11 $(MMDB_CXXFLAGS) $(CLIPPER_CXXFLAGS) $(OPENMP_CFLAGS)
ElectroSurface_LDFLAGS=${LDFLAGS} $(OPENMP_CXXFLAGS) $(CLIPPER_LIBS) $(MMDB_LIBS)

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = MoleculesToTrianglesCXXSurface.pc

