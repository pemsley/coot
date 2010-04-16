

GOO_CFLAGS=`PKG_CONFIG_PATH=/usr/lib/pkgconfig:$(HOME)/goocanvas/lib/pkgconfig pkg-config --cflags goocanvas`
GOO_LDFLAGS=`PKG_CONFIG_PATH=/usr/lib/pkgconfig:$(HOME)/goocanvas/lib/pkgconfig pkg-config --libs goocanvas`

EXTRA_FLAGS=-Wchar-subscripts -Wcomment -Wformat -Winit-self -Wmissing-braces -Wparentheses -Wreturn-type -Wswitch -Wuninitialized -Wundef -Werror -O2

GTK_CFLAGS=`pkg-config --cflags gtk+-2.0` 
GTK_LDFLAGS=`pkg-config --libs gtk+-2.0`

CAIRO_CFLAGS=`pkg-config --cflags cairo` 
# CAIRO_LDFLAGS=`pkg-config --libs cairo` in GOO_LDFLAGS already

CRYST=/home/emsley/autobuild/f12/Linux-jackal-f12-pre-release-gtk2
CRYST=$(HOME)/autobuild/Linux-scylla-pre-release-gtk2

GSL_CFLAGS=-I$(CRYST)/include
GSL_LDFLAGS=$(CRYST)/lib/libgsl.la  $(CRYST)/lib/libgslcblas.la -lm

MMDB=$(CRYST)
# MMDB=$(HOME)/test/mmdb

MMDB_CXXFLAGS=-I$(MMDB)/include/mmdb
MMDB_LDFLAGS=$(MMDB)/lib/libmmdb.la

CLIPPER=$(CRYST)
CLIPPER_CXXFLAGS=-I$(CLIPPER)/include
CLIPPER_LDFLAGS=$(CLIPPER)/lib/libclipper-core.la \
        $(CLIPPER)/lib/librfftw.la \
        $(CLIPPER)/lib/libfftw.la


CXXFLAGS=-g -O2 $(EXTRA_FLAGS) $(MMDB_CXXFLAGS) $(CLIPPER_CXXFLAGS) 

lbg: lig-build.o lbg.o lbg-callbacks.o lbg-search.o some-coot-utils.o lbg-molfile.o residue-circles.o
	libtool --mode=link $(CXX)  -Wl,--export-dynamic -o $@ $^ \
	$(CLIPPER_LDFLAGS) $(MMDB_LDFLAGS) $(GTK_LDFLAGS) $(GOO_LDFLAGS) $(GSL_LDFLAGS)

lig-build.o: lig-build.cc
	libtool --tag=CXX --mode=compile $(CXX) $(CXXFLAGS) $(GTK_CFLAGS) $(GOO_CFLAGS) -c $^

lbg-search.o: lbg-search.cc
	libtool --tag=CXX --mode=compile $(CXX) $(CXXFLAGS)  $(GTK_CFLAGS) $(GOO_CFLAGS) -c $^

residue-circles.o: residue-circles.cc
	libtool --tag=CXX --mode=compile $(CXX) $(CXXFLAGS)  $(GTK_CFLAGS) $(GOO_CFLAGS) $(GSL_CFLAGS) -c $^

lbg.o: lbg.cc
	libtool --tag=CXX --mode=compile $(CXX) $(CXXFLAGS) $(GTK_CFLAGS) $(GOO_CFLAGS) $(CAIRO_CFLAGS) -c $^

lbg-molfile.o: lbg-molfile.cc
	libtool --tag=CXX --mode=compile $(CXX) $(CXXFLAGS) -c $^

test-read-molfile.o: test-read-molfile.cc
	libtool --tag=CXX --mode=compile $(CXX) $(CXXFLAGS) -c $^

some-coot-utils.o: some-coot-utils.cc
	libtool --tag=CXX --mode=compile $(CXX) $(CXXFLAGS) -c $^

lbg-callbacks.o: lbg-callbacks.cc
	libtool --tag=CXX --mode=compile $(CXX) $(CXXFLAGS) $(GTK_CFLAGS) $(GOO_CFLAGS) -c $^

example: example.c
	$(CXX) $(CXXFLAGS) -export-dynamic -o $@ $^ `pkg-config --cflags --libs gtk+-2.0`

test-lig-build: lig-build.o test-lig-build.o
	libtool --mode=link $(CXX)  -Wl,--export-dynamic -o $@ $^ $(GTK_LDFLAGS) $(GOO_LDFLAGS)

test-lig-build.o: test-lig-build.cc 
	libtool --tag=CXX --mode=compile $(CXX) $(CXXFLAGS) $(GTK_CFLAGS) $(GOO_CFLAGS) -c $^

test-inh: test-inh.cc
	libtool --mode=link $(CXX) $(CXXFLAGS) $(GTK_CFLAGS) $(GOO_CFLAGS) $^ -o $@ 

test-read-molfile: test-read-molfile.o lbg-molfile.o
	libtool --mode=link $(CXX) $(CXXFLAGS) $(GTK_CFLAGS) $(GOO_CFLAGS) $^ -o $@ 

dist:
	tar czf lbg.tar.gz lbg.glade *.cc *.hh *.svg coot-icon.png coot-favicon.png notes Makefile

