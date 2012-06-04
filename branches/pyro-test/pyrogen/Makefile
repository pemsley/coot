
# It seems that we need to setup PYTHONPATH correctly to properly link
# _coot_libs.so?


install_top_dir=$(HOME)/autobuild/build-coot+rdkit-pre-release-gtk2-python


prfx = $(install_top_dir)
CXX = g++ -Wreturn-type
HACK_INCLUDES = -I../coot/analysis
INCLUDES = -I$(prfx)/include/python2.7 -I$(prfx)/include $(HACK_INCLUDES) -I$(prfx)/include/rdkit
LIBS=$(prfx)/lib/libcoot-ideal.la \
     $(prfx)/lib/libcoot-coord-utils.la \
     $(prfx)/lib/libcoot-geometry.la \
     $(prfx)/lib/libcoot-analysis.la \
     $(prfx)/lib/libcoot-lidia-core.la \
	 -L$(prfx)/lib -lGraphMol -lboost_python -lpython2.7

libs: _coot_libs.so restraints_boost.so

# restraints-boost.o
_coot_libs.so: wrap-coot-libs.o restraints.o py-restraints.o 
	libtool --mode=link g++ -shared -fPIC $^ -o $@ $(LIBS) 

restraints_boost.so: restraints-boost.o py-restraints.o
	libtool --mode=link g++ -shared -fPIC $^ -o $@ $(LIBS) 

wrap-coot-libs.o: wrap-coot-libs.cc
	libtool --mode=compile $(CXX) -fPIC  $(INCLUDES) -c $^ -o $@

restraints.o: restraints.cc restraints.hh py-restraints.hh
	libtool --mode=compile $(CXX) -fPIC $(INCLUDES) -c $< -o $@

py-restraints.o: py-restraints.cc py-restraints.hh
	libtool --mode=compile $(CXX) -fPIC $(INCLUDES) -c $< -o $@

restraints-boost.o: restraints-boost.cc restraints.hh
	libtool --mode=compile $(CXX) -fPIC $(INCLUDES) -c $< -o $@

wrap-coot-libs.cc: wrap-coot-libs.i restraints.hh
	swig  $(INCLUDES) -o wrap-coot-libs.cc -python -shadow -c++ wrap-coot-libs.i

test-call-python: test-call-python.cc
	libtool --mode=link $(CXX) -fPIC $(INCLUDES) $< -o $@ $(LIBS)

clean:
	rm -f _coot_libs.so coot_libs.py \
              restraints.lo restraints.o example.o _example.so example_wrap.o \
	      wrap-coot-libs.cc wrap-coot-libs.lo wrap-coot-libs.o example-wrap.o

test:
	# python test-coot-libs.py
	python p2.py

boost-example.o: boost-example.cc
	libtool --mode=compile $(CXX) -fPIC $(INCLUDES) -c $< -o $@

boost-example-defs.o: boost-example-defs.cc
	libtool --mode=compile $(CXX) -fPIC $(INCLUDES) -c $< -o $@

mytest.so: boost-example.o boost-example-defs.o
	libtool --mode=link g++ -shared -fPIC $^ -o $@ $(LIBS)
