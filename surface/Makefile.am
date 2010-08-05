
INCLUDES = -I$(top_srcdir)/coot-utils $(MMDB_CXXFLAGS) $(CLIPPER_CXXFLAGS) $(FFTW_CXXFLAGS)

AM_CPPFLAGS = $(OPENMP_CFLAGS)

lib_LTLIBRARIES = libccp4mg-surface.la

EXTRA_DIST = surface.cpp surface.h \
	CXXTorusTriangle.cpp \
	CXXUtils.cpp \
	CXXUtils.h \
	ElectroSurface.cpp


libccp4mg_surface_la_SOURCES = \
	CXXSphereElement.cpp \
	CXXSphereElement.h \
	CXXChargeTable.cpp		CXXSphereNode.cpp \
	CXXChargeTable.h		CXXSphereNode.h \
	CXXCircle.cpp			CXXSphereTriangle.cpp \
	CXXCircle.h			CXXSphereTriangle.h \
	CXXCircleNode.cpp		CXXSphereTriangleEdge.cpp  \
	CXXCircleNode.h			CXXSphereTriangleEdge.h \
	CXXCoord.cpp			CXXSurface.cpp \
	CXXCoord.h			CXXSurface.h \
	CXXCreator.cpp			\
	CXXCreator.h			CXXSurfaceVertex.cpp \
	CXXException.cpp		CXXSurfaceVertex.h \
	CXXException.h			CXXTorusElement.cpp \
	CXXFFTSolventMap.cpp		CXXTorusElement.h \
	CXXFFTSolventMap.h		CXXTorusNode.cpp \
	CXXFortranFile.cpp		CXXTorusNode.h \
	CXXFortranFile.h		 \
	CXXTorusTriangle.h \
	CXXTorusTriangle.h.h \
	CXXNewHood.cpp			CXXTriangle.cpp \
	CXXNewHood.h			CXXTriangle.h \
	CXXNode.cpp\
	CXXNode.h			TokenIterator.h \
	CXXQADSurface.cpp		\
	CXXQADSurface.h	\
	CXXBall.cpp CXXBall.h       \
	CXXQADSurfaceStaticData.cpp \
	CXXSpace.cpp				\
	CXXSpace.h \
	CXXMatrix4x4.h \
	CXXPointyBit.h \
	CXXSphereFlatTriangle.h

libccp4mg_surface_la_LIBADD = $(CLIPPER_LIBS) $(MMDB_LIBS)


# CXXAtom.h CXXAtom.cpp

