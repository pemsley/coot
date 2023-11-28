
ACLOCAL_AMFLAGS="-Im4"
lib_LTLIBRARIES=libMoleculesToTrianglesCXXSurface.la

include_HEADERS = \
ctpl_boost.h	  CXXCoord.h	      CXXNewHood.h	       CXXSphereNode.h		CXXTorusElement.h     surface.h        \
CXXAlloc.h	  CXXCreator.h	      CXXPointyBit.h	       CXXSphereTriangleEdge.h	CXXTorusNode.h	      TokenIterator.h  \
CXXBall.h	  CXXException.h      CXXQADSurface.h	       CXXSphereTriangle.h	CXXTorusTriangle.h.h                   \
CXXChargeTable.h  CXXFFTSolventMap.h  CXXSpace.h	       CXXSurface.h		CXXTriangle.h                          \
CXXCircle.h	  CXXFortranFile.h    CXXSphereElement.h       CXXSurfaceMaker.h	CXXUtils.h                             \
CXXCircleNode.h   CXXMatrix4x4.h      CXXSphereFlatTriangle.h  CXXSurfaceVertex.h	CXXVCN.h

libMoleculesToTrianglesCXXSurface_la_SOURCES=CXXCircle.cpp CXXCircleNode.cpp \
        CXXCoord.cpp CXXFortranFile.cpp CXXNewHood.cpp CXXSphereElement.cpp \
	CXXSphereNode.cpp CXXSphereTriangle.cpp CXXSphereTriangleEdge.cpp CXXSurface.cpp CXXSurfaceVertex.cpp \
	CXXTorusElement.cpp CXXTorusNode.cpp CXXTriangle.cpp CXXQADSurface.cpp  CXXQADSurfaceStaticData.cpp \
	CXXCreator.cpp CXXSpace.cpp CXXFFTSolventMap.cpp CXXChargeTable.cpp CXXException.cpp \
	CXXUtils.cpp CXXBall.cpp CXXSurfaceMaker.cpp

libMoleculesToTrianglesCXXSurface_la_LIBADD = $(CLIPPER_LIBS) $(MMDB_LIBS)

libMoleculesToTrianglesCXXSurface_la_LDFLAGS=$(SHARED_LDFLAGS)

AM_CXXFLAGS=$(OPENMP_CFLAGS) $(COOT_WINDOWS_CFLAGS)

