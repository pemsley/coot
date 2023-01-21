
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
    ligand-builder.hpp

SOURCES = \
    ligand-builder.cpp \
    main.cpp

RUST_SOURCES = \
    ligand_editor_canvas/src/ffi.rs \
    ligand_editor_canvas/src/lib.rs \
    ligand_editor_canvas/src/model.rs \
    ligand_editor_canvas/src/widget.rs \
    ligand_editor_canvas/build.rs \
    ligand_editor_canvas/cbindgen.toml \
    ligand_editor_canvas/Cargo.toml 


OBJS := $(patsubst %.cc,%.o,$(SOURCES))

%.o: %.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

new-ligand-editor: ligand_editor_canvas.hpp ligand_editor_canvas/target/release/libligand_editor_canvas.a $(OBJS) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDLIBS)

ligand_editor_canvas.hpp: $(RUST_SOURCES)
	./build_rust_lib.sh