
install_top_dir=$(HOME)/autobuild/Linux-tiresias-pre-release-gtk2-python

prfx = $(install_top_dir)
CXX = g++
INCLUDES = -I$(prfx)/include/python2.7 -I$(prfx)/include
LIBS=$(prfx)/lib/libcoot-geometry.la -L$(prfx)/lib -lpython2.7

libs: _coot-libs.so

_coot-libs.so: libcoot-libs.la

libcoot-libs.la: wrap-coot-libs.lo restraints.lo
#	libtool --tag=CXX --mode=link $(CXX) $^ -o $@ $(LIBS)
#	ld -shared -fPIC $^ -o $@
	libtool --tag=CXX --mode=link $(CXX)  $^ -o $@ $(LIBS)


wrap-coot-libs.lo: wrap-coot-libs.cc
	libtool --tag=CXX --mode=compile $(CXX) $(INCLUDES) $^ -c -o $@

wrap-coot-libs.cc: wrap-coot-libs.i
	swig -o wrap-coot-libs.cc -python -c++ wrap-coot-libs.i

restraints.lo: restraints.cc
	libtool --tag=CXX --mode=compile $(CXX) $(INCLUDES) $^ -c -o $@

clean:
	rm -f restraints.lo restraints.o example.o _example.so example_wrap.o \
	      wrap-coot-libs.lo wrap-coot-libs.o example-wrap.o

