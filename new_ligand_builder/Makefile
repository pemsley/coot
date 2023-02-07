
# to find clipper.pc, you will need to do something like
# export PKG_CONFIG_PATH=$COOT_DIR/lib/pkgconfig

COOT_DIR=$(HOME)/freshcoot/experimental_systemlibsbuild/prefix

CXXFLAGS= \
   -I$(COOT_DIR)/include \
   -I$(COOT_DIR)/include/coot \
   -I/usr/include/rdkit \
   `pkg-config --cflags gtk4` \
   -ggdb -O2 -march=native -std=c++17
# `pkg-config --cflags epoxy` 
# `pkg-config --cflags clipper` 

LDLIBS=-fuse-ld=mold \
       -L$(COOT_DIR)/lib \
       `pkg-config --libs gtk4`  \
        -g \
        -Wl,--rpath=$(COOT_DIR)/lib \
        -lfreetype -lm -ldl -lRDKitSmilesParse -lRDKitMolDraw2D -lRDKitGraphMol

	   # -lcoot-ligand -lcoot-ideal -lcoot-coords -lcoot-map-utils -lcoot-coord-utils -lcoot-geometry -lmmdb2 
       #`pkg-config --libs epoxy` 
       #`pkg-config --libs clipper` 

HEADERS = \
    ligand-builder.hpp \
    ligand_editor_canvas/model.hpp \
    ligand_editor_canvas/tools.hpp \
    ligand_editor_canvas.hpp

SOURCES = \
    ligand-builder.cpp \
    ligand_editor_canvas/model.cpp \
    ligand_editor_canvas/tools.cpp \
    ligand_editor_canvas.cpp \
    main.cpp



OBJS := $(patsubst %.cpp,%.o,$(SOURCES))

%.o: %.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

new-ligand-editor: $(OBJS) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDLIBS)
