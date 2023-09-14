
# to find clipper.pc, you will need to do something like
# export PKG_CONFIG_PATH=$COOT_DIR/lib/pkgconfig

COOT_DIR=$(HOME)/freshcoot/experimental_systemlibsbuild/prefix2

CXXFLAGS= \
   -I$(COOT_DIR)/include \
   -I$(COOT_DIR)/include/coot \
   -I/usr/include/rdkit \
   -I/usr/include/python3.11 \
   -DUSE_PYTHON \
   `pkg-config --cflags gtk4` \
   -ggdb -O2 -march=native -std=c++17
# `pkg-config --cflags epoxy` 
# `pkg-config --cflags clipper` 

LDLIBS=-fuse-ld=mold \
       -L$(COOT_DIR)/lib \
       `pkg-config --libs gtk4`  \
        -g \
        -Wl,--export-dynamic,--rpath=$(COOT_DIR)/lib \
        -lfreetype -lm -ldl \
        -lcootsumo -lcoot-geometry -lcoot-lidia-core \
        -lRDKitDepictor -lRDKitSmilesParse -lRDKitMolDraw2D -lRDKitGraphMol -lRDKitRDGeneral -lRDKitFileParsers -lRDKitSubstructMatch \
        -lboost_iostreams -lboost_serialization -lboost_python311 -lpython3.11

	   # -lcoot-ligand -lcoot-ideal -lcoot-coords -lcoot-map-utils -lcoot-coord-utils -lcoot-geometry -lmmdb2 
       #`pkg-config --libs epoxy` 
       #`pkg-config --libs clipper` 

HEADERS = \
    ligand_builder_state.hpp \
    ligand_builder_ui.hpp \
    ligand_builder_utils.hpp \
    ligand_builder_generators.hpp \
    ligand_editor_canvas/model.hpp \
    ligand_editor_canvas/core.hpp \
    ligand_editor_canvas/tools.hpp \
    ligand_editor_canvas.hpp \

SOURCES = \
    ligand_builder_state.cpp \
    ligand_builder_signals.cpp \
    ligand_builder_ui.cpp \
    ligand_builder_utils.cpp \
    ligand_builder_generators.cpp \
    ligand_editor_canvas/model.cpp \
    ligand_editor_canvas/core.cpp \
    ligand_editor_canvas/tools.cpp \
    ligand_editor_canvas.cpp \
    main.cpp



OBJS := $(patsubst %.cpp,%.o,$(SOURCES))

%.o: %.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

layla: $(OBJS) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDLIBS)
