#
AM_CPPFLAGS = \
	-I$(top_srcdir)     \
	$(CLIPPER_CXXFLAGS) \
	$(MMDB_CXXFLAGS)    \
	$(COOT_WINDOWS_CFLAGS)

# what should the library be called? libcoot-buccaneer-ml-growing.la perhaps?
lib_LTLIBRARIES=libcoot-ml-growing.la

libcoot_ml_growing_la_SOURCES = \
   buccaneer-build.cpp    buccaneer-merge.cpp	      buccaneer-prot.cpp	   k2c_convolution_layers.cpp	 ml-grow-next.cpp \
   buccaneer-correct.cpp  buccaneer-ncsbuild.cpp      buccaneer-prune.cpp	   k2c_core_layers.cpp		 ml-grow-prev.cpp \
   buccaneer-filter.cpp   buccaneer-NN-build.cpp      buccaneer-sequence.cpp       k2c_embedding_layers.cpp	 protein_db.cpp   \
   buccaneer-find.cpp     buccaneer-NN-build-cyc.cpp  buccaneer-tidy.cpp	   k2c_helper_functions.cpp	 simplex-lib.cpp  \
   buccaneer-grow.cpp     buccaneer-NN-features.cpp   buccaneer-util.cpp	   k2c_merge_layers.cpp		 simulate-lib.cpp \
                          buccaneer-NN-join.cpp	                 	           k2c_normalization_layers.cpp \
   buccaneer-known.cpp    buccaneer-NN-model.cpp                                   k2c_pooling_layers.cpp       \
   buccaneer-lib.cpp      buccaneer-NN-Select.cpp                		   k2c_recurrent_layers.cpp     \
   buccaneer-link.cpp     buccaneer-prep.cpp	      k2c_activations.cpp	   ml-grow.cpp

libcoot_ml_growing_la_LIBADD = $(CLIPPER_LIBS)

libcoot_ml_growing_la_LDFLAGS = $(SHARED_LDFLAGS)

pkginclude_HEADERS = \
 buccaneer-build.h    buccaneer-known.h	        buccaneer-NN-build.h	buccaneer-prune.h     ml-grow.h              \
 buccaneer-correct.h  buccaneer-lib.h	        buccaneer-NN-features.h	buccaneer-sequence.h  ml-grow-next.h         \
 buccaneer-filter.h   buccaneer-link.h	        buccaneer-NN-model.h	buccaneer-tidy.h      ml-grow-prev.h         \
 buccaneer-find.h     buccaneer-merge.h	        buccaneer-NN-Select.h	buccaneer-util.h      protein_db.h           \
 buccaneer-grow.h     buccaneer-ncsbuild.h      buccaneer-prep.h		k2c_include.h	      simplex-lib.h  \
 buccaneer-join.h     buccaneer-NN-build-cyc.h  buccaneer-prot.h		k2c_tensor_include.h  simulate-lib.h

pkgincludedir = $(includedir)/coot/buccaneer_ml_growing



libexec_PROGRAMS = cbuccaneerNN csequins

cbuccaneerNN_SOURCES = cbuccaneerNN.cpp

cbuccaneerNN_LDADD = ./libcoot-ml-growing.la $(CLIPPER_LIBS)

csequins_SOURCES = csequins.cpp

csequins_LDADD = ./libcoot-ml-growing.la $(CLIPPER_LIBS)
