
ACLOCAL_AMFLAGS="-Im4"
AM_CPPFLAGS = -I$(top_srcdir) $(BOOST_CPPFLAGS) $(MMDB_CFLAGS)
AM_LDFLAGS = $(BOOST_LDFLAGS)

lib_LTLIBRARIES=libMoleculesToTrianglesCXXClasses.la libBoostedMoleculesToTrianglesCXXClasses.la

libMoleculesToTrianglesCXXClasses_la_SOURCES=AtomPropertyRampColorRule.cpp	 \
        DirectedPositionedEntity.cpp						 \
	Representation.cpp BallsPrimitive.cpp DiscreteSegment.cpp		 \
	RepresentationInstance.cpp BondsPrimitive.cpp DisplayPrimitive.cpp	 \
	RotatedTranslatedScaledEntity.cpp BoxSectionPrimitive.cpp Light.cpp	 \
	SceneSetup.cpp CXXMatrix.cpp MolecularRepresentation.cpp		 \
	SecondaryColorScheme.cpp Camera.cpp MyMolecule.cpp	                 \
	SolidColorRule.cpp ColorRule.cpp NRStuff.cpp SpherePrimitive.cpp	 \
	ColorScheme.cpp Polyhedron.cpp SticksPrimitive.cpp CompoundSelection.cpp \
	PositionedEntity.cpp SurfacePrimitive.cpp CylinderPoint.cpp Renderer.cpp \
	VertexColorNormalPrimitive.cpp CylindersPrimitive.cpp RendererGL.cpp	 \
	oglPolyhedron.cpp MolecularRepresentationInstance.cpp Quaternion.cpp	 \
	Vector2.cpp Vector3.cpp RendererGLSL.cpp lfit.c DishyBase.cpp            \
        FlatFanPrimitive.cpp

# I want pkg-config/configure to make MMDB_CPPFLAGS, not MMDB_CFLAGS (but, after some reading,
# I don't think that's how pkg-config works - oh well).
#
# is $(CPPFLAGS) needed here?
libMoleculesToTrianglesCXXClasses_la_CPPFLAGS=$(AM_CPPFLAGS) \
   $(MoleculesToTrianglesCXXSurface_CFLAGS) $(MMDB_CFLAGS) \
   $(OPENMP_CXXFLAGS)

# rather than specifying the libs directly, you should test for clipper in configure
# and use $(CLIPPER_LIBS) here
#
# libMoleculesToTrianglesCXXClasses_la_LIBADD=$(MoleculesToTrianglesCXXSurface_LIBS) \
# -lclipper-ccp4 -lclipper-core  -lrfftw -lfftw -lccp4c -lmmdb2 \
# $(GL_LIBS) $(GLU_LIBS)

libMoleculesToTrianglesCXXClasses_la_LIBADD=$(MoleculesToTrianglesCXXSurface_LIBS) \
   $(CLIPPER_LIBS) $(MMDB_LIBS) \
   $(GL_LIBS) $(GLU_LIBS)

# otherincludedir = $(includedir)/MoleculesToTriangles/CXXClasses
# otherincludedir = $(top_srcdir)/MoleculesToTriangles/CXXClasses
# other
include_HEADERS = Quaternion.hpp Vector.hpp Vector2.hpp \
   AtomPropertyRampColorRule.h 	DirectedPositionedEntity.h		RendererImplementations.h       \
   BallsPrimitive.h		DiscreteSegment.h			Representation.h                \
   BondsPrimitive.h		DisplayPrimitive.h			RepresentationInstance.h        \
   BoxSectionPrimitive.h	GLIndexType.h				RotatedTranslatedScaledEntity.h \
   CXXMatrix.h			Light.h					SceneSetup.h                    \
   Callbacks.h			LinesPrimitive.h			SecondaryColorScheme.h          \
   Camera.h			MolecularRepresentation.h		SolidColorRule.h                \
   CameraPort.h			MolecularRepresentationInstance.h	SpherePrimitive.h               \
   ColorNamesMap.h		MyMolecule.h				SticksPrimitive.h               \
   ColorRule.h			NRStuff.h				SurfacePrimitive.h              \
   ColorScheme.h		Polyhedron.h				VertexColorNormalPrimitive.h    \
   CompoundSelection.h		PositionedEntity.h			lfit.h                          \
   CylinderPoint.h		Renderer.h				oglPolyhedron.h                 \
   CylindersPrimitive.h		RendererGL.h                            RendererGLSL.hpp

# PointLight.hpp: $(top_srcdir)/PointLight.frag $(top_srcdir)/PointLight.vert $(top_srcdir)/ShaderPacker.py
# python $(top_srcdir)/ShaderPacker.py PointLight $(top_srcdir)

libBoostedMoleculesToTrianglesCXXClasses_la_SOURCES=BoostWrappers.cpp
libBoostedMoleculesToTrianglesCXXClasses_la_CPPFLAGS= ${AM_CPPFLAGS} \
     $(MoleculesToTrianglesCXXSurface_CFLAGS)
libBoostedMoleculesToTrianglesCXXClasses_la_LIBADD=-lMoleculesToTrianglesCXXClasses \
     $(MoleculesToTrianglesCXXSurface_LIBS) \
     $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB) \
     $(CLIPPER_LIBS) $(MMDB_LIBS) \
     $(GL_LIBS)

# -lrfftw -lfftw -lclipper-core -lclipper-ccp4 -lccp4c -lmmdb2

check_PROGRAMS=DemoApp

DemoApp_SOURCES=DemoApp.cpp CameraPort.cpp

#DemoApp_LDADD=-lMoleculesToTrianglesCXXClasses \
#   $(MoleculesToTrianglesCXXSurface_LIBS) \
#   -lrfftw -lfftw -lclipper-core -lclipper-ccp4 -lccp4c -lmmdb2 \
#   $(GL_LIBS) $(GLU_LIBS) $(glfw3_LIBS)
#
# We don't need to specify those extra libraries, using _LIBADD above means that the MoleculesToTrianglesCXXClasses
# knows what it needs to link when assembling the executable (perhaps at run-time).
#
DemoApp_LDADD=-lMoleculesToTrianglesCXXClasses $(glfw3_LIBS)

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = MoleculesToTrianglesCXXClasses.pc
