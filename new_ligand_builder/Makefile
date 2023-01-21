
# to find clipper.pc, you will need to do something like
# export PKG_CONFIG_PATH=$COOT_DIR/lib/pkgconfig

COOT_DIR=$(HOME)/freshcoot/experimental_systemlibsbuild/prefix

CXXFLAGS= \
   -I$(COOT_DIR)/include \
   -I$(COOT_DIR)/include/coot \
   `pkg-config --cflags gtk4` \
   -ggdb -O2 -march=native -std=c++11
# `pkg-config --cflags epoxy` 
# `pkg-config --cflags clipper` 

LDLIBS=-fuse-ld=mold \
       -L./ligand_editor_canvas/target/release \
       -L$(COOT_DIR)/lib \
       `pkg-config --libs gtk4`  \
        -g \
        -Wl,--rpath=$(COOT_DIR)/lib \
        -lfreetype -lm -ldl -lligand_editor_canvas

	   # -lcoot-ligand -lcoot-ideal -lcoot-coords -lcoot-map-utils -lcoot-coord-utils -lcoot-geometry -lmmdb2 
       #`pkg-config --libs epoxy` 
       #`pkg-config --libs clipper` 

HEADERS = \
    ligand-builder.hpp \
	ligand_editor_canvas.hpp

SOURCES = \
    ligand-builder.cpp \
    main.cpp

OBJS := $(patsubst %.cc,%.o,$(SOURCES))

%.o: %.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

BUILD_RUST:
	pushd ligand_editor_canvas && cargo build -r && popd 

new-ligand-editor: BUILD_RUST $(OBJS) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDLIBS)
