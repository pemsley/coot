# buccaneer makefile

SHARED = \
buccaneer-build \
buccaneer-correct \
buccaneer-filter \
buccaneer-find \
buccaneer-grow \
buccaneer-known \
buccaneer-lib \
buccaneer-link \
buccaneer-merge \
buccaneer-ncsbuild \
buccaneer-prep \
buccaneer-prot \
buccaneer-prune \
buccaneer-sequence \
buccaneer-tidy \
buccaneer-util \
k2c_activations \
k2c_convolution_layers \
k2c_core_layers \
k2c_embedding_layers \
k2c_helper_functions \
k2c_merge_layers \
k2c_normalization_layers \
k2c_pooling_layers \
k2c_recurrent_layers \
protein_db \
ml-grow \
ml-grow-next \
ml-grow-prev \
simplex-lib \
simulate-lib \

ORIG_ONLY = \
buccaneer-join \

NN_ONLY = \
buccaneer-NN-build-cyc \
buccaneer-NN-build \
buccaneer-NN-features \
buccaneer-NN-join \
buccaneer-NN-model \
buccaneer-NN-Select \

LIBS = \
-l:libccp4c.so.8.0 \
-l:libclipper-ccp4.so.2 \
-l:libclipper-contrib.so.2 \
-l:libclipper-core.so.2 \
-l:libclipper-minimol.so.2 \
-l:libclipper-mmdb.so.2 \
-l:libfftw.so.2 \
-l:libmmdb2.so.0 \
-l:librfftw.so.2 \
-lm \
-lstdc++ \

FLAGS = \
-std=c++11 \
-O2 \
-Wall \
-Wno-sign-compare \
-fPIC \
-ftemplate-depth-50 \
-D_GLIBCXX_USE_CXX11_ABI=0 \
-DBUCCANEER_PROFILE=1 \

CFLAGS = ${FLAGS} -I${CCP4}/include -c
LFLAGS = ${FLAGS} -L${CCP4}/lib ${LIBS}

SHARED_OBJS = $(SHARED:=.o)
ORIG_OBJS = $(ORIG_ONLY:=.o)
NN_OBJS = $(NN_ONLY:=.o)

TARGETS = cbuccaneer cbuccaneerNN csequins

all: ${TARGETS}

cbuccaneer: ${SHARED_OBJS} ${ORIG_OBJS} cbuccaneer.o
	g++ ${SHARED_OBJS} ${ORIG_OBJS} cbuccaneer.o -o $@ ${LFLAGS}

cbuccaneerNN: ${SHARED_OBJS} ${NN_OBJS} cbuccaneerNN.o
	g++ ${SHARED_OBJS} ${NN_OBJS} cbuccaneerNN.o -o $@ ${LFLAGS}

csequins: ${SHARED_OBJS} ${ORIG_OBJS} csequins.o
	g++ ${SHARED_OBJS} ${ORIG_OBJS} csequins.o -o $@ ${LFLAGS}

cbuccaneer.o: cbuccaneer.cpp
	g++ ${CFLAGS} $<

cbuccaneerNN.o: cbuccaneerNN.cpp
	g++ ${CFLAGS} $<

csequins.o: csequins.cpp
	g++ ${CFLAGS} $<

buccaneer-NN-join.o: buccaneer-NN-join.cpp buccaneer-join.h
	g++ ${CFLAGS} $<

%.o: %.cpp %.h
	g++ ${CFLAGS} $<

clean:
	rm *.o ${TARGETS}
