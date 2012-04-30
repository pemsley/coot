

install_top_dir=$(HOME)/autobuild/build-coot+rdkit-pre-release-gtk2-python


prfx = $(install_top_dir)
CXX = g++
HACK_INCLUDES = -I../coot/analysis
INCLUDES = -I$(prfx)/include/python2.7 -I$(prfx)/include $(HACK_INCLUDES) -I$(prfx)/include/rdkit
LIBS=$(prfx)/lib/libcoot-geometry.la $(prfx)/lib/libcoot-analysis.la -L$(prfx)/lib -lGraphMol -lpython2.7

libs: _coot_libs.so

_coot_libs.so: wrap-coot-libs.o restraints.o
	libtool --mode=link g++ -shared -fPIC $^ -o $@ $(LIBS)

wrap-coot-libs.o: wrap-coot-libs.cc
	libtool --mode=compile $(CXX) -fPIC  $(INCLUDES) -c $^ -o $@

restraints.o: restraints.cc
	libtool --mode=compile $(CXX) -fPIC $(INCLUDES) -c $^ -o $@

wrap-coot-libs.cc: wrap-coot-libs.i restraints.hh
	swig  $(INCLUDES) -o wrap-coot-libs.cc -python -shadow -c++ wrap-coot-libs.i

clean:
	rm -f restraints.lo restraints.o example.o _example.so example_wrap.o \
	      wrap-coot-libs.cc wrap-coot-libs.lo wrap-coot-libs.o example-wrap.o

test:
	# python test-coot-libs.py
	python p2.py

