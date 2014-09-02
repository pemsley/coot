## Process this file with automake to produce Makefile.in

EXTRA_DIST = post-glade post-gtk2-glade.sh

INCLUDES = \
	-DPACKAGE_DATA_DIR=\""$(datadir)"\" \
	-DPACKAGE_LOCALE_DIR=\""$(prefix)/$(DATADIRNAME)/locale"\" \
	@PACKAGE_CFLAGS@

bin_PROGRAMS = restraints-editor example

example_SOURCES = example.c
example_LDADD = @PACKAGE_LIBS@

restraints_editor_SOURCES = \
	main.cc \
	support.c support.h \
	interface.c interface.h \
	callbacks.c callbacks.h

restraints_editor_LDADD = @PACKAGE_LIBS@ $(INTLLIBS)

