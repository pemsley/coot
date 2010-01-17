

CCP4MGHOME=$HOME/ccp4mg/cvs

# clean up doublicates at some point!!! FIXME

# this is the new ones for ribbons...
# 1.) dir util
cp $CCP4MGHOME/util/cartesian.cc  \
   $CCP4MGHOME/util/connect.cc    \
   $CCP4MGHOME/util/CParamsManager.cc \
   $CCP4MGHOME/util/geomutil.cc    \
   $CCP4MGHOME/util/matrix.cc	\
   $CCP4MGHOME/util/mginterrupt.cc 	\
   $CCP4MGHOME/util/mgtree.cc	\
   $CCP4MGHOME/util/mgutil.cc 	\
   $CCP4MGHOME/util/plane.cc 	\
   $CCP4MGHOME/util/quat.cc	\
   $CCP4MGHOME/util/redirect.cc \
   $CCP4MGHOME/util/volume.cc util/ 

cp $CCP4MGHOME/util/cartesian.h  \
   $CCP4MGHOME/util/connect.h    \
   $CCP4MGHOME/util/CParamsManager.h \
   $CCP4MGHOME/util/geomutil.h    \
   $CCP4MGHOME/util/matrix.h	\
   $CCP4MGHOME/util/mginterrupt.h 	\
   $CCP4MGHOME/util/mgtree.h	\
   $CCP4MGHOME/util/mgutil.h 	\
   $CCP4MGHOME/util/plane.h 	\
   $CCP4MGHOME/util/quat.h	\
   $CCP4MGHOME/util/redirect.h \
   $CCP4MGHOME/util/volume.h util/ 

# 2.) dir mmut
cp $CCP4MGHOME/mmut/mman_base.cc  \
   $CCP4MGHOME/mmut/mman_manager.cc    \
   $CCP4MGHOME/mmut/mmut_basepairs.cc    \
   $CCP4MGHOME/mmut/mmut_bonds.cc    \
   $CCP4MGHOME/mmut/mmut_connectivity.cc    \
   $CCP4MGHOME/mmut/mmut_contact.cc    \
   $CCP4MGHOME/mmut/mmut_hbond.cc    \
   $CCP4MGHOME/mmut/mmut_lipids.cc    \
   $CCP4MGHOME/mmut/mmut_manager.cc    \
   $CCP4MGHOME/mmut/mmut_morph.cc    \
   $CCP4MGHOME/mmut/mmut_nma.cc    \
   $CCP4MGHOME/mmut/mmut_sasarea.cc    \
   $CCP4MGHOME/mmut/mmut_sbase.cc    \
   $CCP4MGHOME/mmut/mmut_secstr.cc    \
   $CCP4MGHOME/mmut/mmut_spline.cc    \
   $CCP4MGHOME/mmut/mmut_util.cc mmut/

cp $CCP4MGHOME/mmut/mman_base.h  \
   $CCP4MGHOME/mmut/mman_manager.h    \
   $CCP4MGHOME/mmut/mmut_basepairs.h    \
   $CCP4MGHOME/mmut/mmut_bonds.h    \
   $CCP4MGHOME/mmut/mmut_connectivity.h    \
   $CCP4MGHOME/mmut/mmut_contact.h    \
   $CCP4MGHOME/mmut/mmut_hbond.h    \
   $CCP4MGHOME/mmut/mmut_lipids.h    \
   $CCP4MGHOME/mmut/mmut_manager.h    \
   $CCP4MGHOME/mmut/mmut_morph.h    \
   $CCP4MGHOME/mmut/mmut_nma.h    \
   $CCP4MGHOME/mmut/mmut_sasarea.h    \
   $CCP4MGHOME/mmut/mmut_sbase.h    \
   $CCP4MGHOME/mmut/mmut_secstr.h    \
   $CCP4MGHOME/mmut/mmut_util.h \
   $CCP4MGHOME/mmut/splineinfo.h mmut/
 
# fix mmut/mmut_sbase.h/cc
cp mmut/mmut_sbase.h  mmut/mmut_sbase.h.orig
cp mmut/mmut_sbase.cc mmut/mmut_sbase.cc.orig
sed -e 's/CLink/MGCLink/g' mmut/mmut_sbase.h.orig  > mmut/mmut_sbase.h
sed -e 's/CLink/MGCLink/g' mmut/mmut_sbase.cc.orig > mmut/mmut_sbase.cc

# 3.) dir pygl
cp $CCP4MGHOME/pygl/atom_util.cc  \
   $CCP4MGHOME/pygl/bezier.cc    \
   $CCP4MGHOME/pygl/billboard.cc    \
   $CCP4MGHOME/pygl/build_tree_primitives.cc    \
   $CCP4MGHOME/pygl/catmull.cc    \
   $CCP4MGHOME/pygl/cdisplayobject.cc    \
   $CCP4MGHOME/pygl/cprimitive.cc    \
   $CCP4MGHOME/pygl/font_info.cc    \
   $CCP4MGHOME/pygl/font_util.cc    \
   $CCP4MGHOME/pygl/freetype_dl.c    \
   $CCP4MGHOME/pygl/freetype_font.cc    \
   $CCP4MGHOME/pygl/help.cc    \
   $CCP4MGHOME/pygl/lincrv.cc    \
   $CCP4MGHOME/pygl/ppmutil.cc    \
   $CCP4MGHOME/pygl/psutil.cc    \
   $CCP4MGHOME/pygl/rgbreps.cc    \
   $CCP4MGHOME/pygl/sphere.cc    \
   $CCP4MGHOME/pygl/splineset.cc    \
   $CCP4MGHOME/pygl/subdivide.cc    \
   $CCP4MGHOME/pygl/symmetry.cc    \
   $CCP4MGHOME/pygl/text.cc    \
   $CCP4MGHOME/pygl/texture.cc    \
   $CCP4MGHOME/pygl/win_font.cc    \
   $CCP4MGHOME/pygl/x11_font.cc    \
   $CCP4MGHOME/pygl/zsortps.cc pygl/ 
# we dont get this as we dont want any QT things and it contains
# only a dummy function currently
#   $CCP4MGHOME/pygl/qt-text.cc    \

cp $CCP4MGHOME/pygl/atom_util.h  \
   $CCP4MGHOME/pygl/cbuild.h    \
   $CCP4MGHOME/pygl/catmull.h    \
   $CCP4MGHOME/pygl/cdisplayobject.h    \
   $CCP4MGHOME/pygl/cprimitive.h    \
   $CCP4MGHOME/pygl/font_info.h    \
   $CCP4MGHOME/pygl/font_util.h    \
   $CCP4MGHOME/pygl/freetype_font.h    \
   $CCP4MGHOME/pygl/freetype_dl.h    \
   $CCP4MGHOME/pygl/help.h    \
   $CCP4MGHOME/pygl/lincrv.h    \
   $CCP4MGHOME/pygl/ppmutil.h    \
   $CCP4MGHOME/pygl/psutil.h    \
   $CCP4MGHOME/pygl/rgbreps.h    \
   $CCP4MGHOME/pygl/sphere.h    \
   $CCP4MGHOME/pygl/sphere_0.h    \
   $CCP4MGHOME/pygl/sphere_1.h    \
   $CCP4MGHOME/pygl/sphere_2.h    \
   $CCP4MGHOME/pygl/sphere_3.h    \
   $CCP4MGHOME/pygl/sphere_4.h    \
   $CCP4MGHOME/pygl/sphere_arrays.h    \
   $CCP4MGHOME/pygl/subdivide.h    \
   $CCP4MGHOME/pygl/symmetry.h    \
   $CCP4MGHOME/pygl/texture.h    \
   $CCP4MGHOME/pygl/win_font.h    \
   $CCP4MGHOME/pygl/x11_font.h pygl/ 

# 4.) dir mgapp
cp $CCP4MGHOME/mgapp/mgapp_base.cc  \
   $CCP4MGHOME/mgapp/mg_colour.cc mgapp/ 

cp $CCP4MGHOME/mgapp/mgapp_base.h  \
   $CCP4MGHOME/mgapp/mg_colour.h mgapp/ 

# insert lines in mg_colour.h
cp mgapp/mg_colour.h mgapp/mg_colour.h.orig
sed -e 's/#include "atom_util.h"/#include "atom_util.h"\n\/\/ BL says:: this is in mmdb_atom.h in MG\n#define SSE_Bridge        7\n#define SSE_Bend          8/' mgapp/mg_colour.h.orig  > mgapp/mg_colour.h
