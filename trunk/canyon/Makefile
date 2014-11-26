

COOT_INSTALL=$(HOME)/autobuild/build-coot+rdkit-pre-release-gtk2-python

CXXFLAGS=-I$(COOT_INSTALL)/include -I$(COOT_INSTALL)/include/coot -g

LDFLAGS=-L$(COOT_INSTALL)/lib -g -lcoot-coord-utils -lclipper-core -lmmdb2

canyon: canyon.o quaternion.o
	libtool --mode=link $(CXX) -o $@ $^ $(LDFLAGS)

quaternion.o: quaternion.cc quaternion.hh
	libtool --mode=compile $(CXX) $(CXXFLAGS) -c $<

canyon.o: canyon.cc canyon.hh
	libtool --mode=compile $(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f canyon canyon.o canyon.lo

dist:
	d=$$(date +%Y-%m-%d-%H%M) && echo $$d && tar czf canyon-$$d.tar.gz canyon.cc canyon.hh quaternion.cc quaternion.hh Makefile
