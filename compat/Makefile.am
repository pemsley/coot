

# these compat files depend on libintl.h - which is part of glibc, 
# which is not necessarily installed in /usr (e.g. not on a Mac)
# so the compiler doesn't find it.  It *is* in CCP4

lib_LTLIBRARIES = libcoot-compat.la
pkginclude_HEADERS = coot-getopt.h
pkgincludedir = $(includedir)/coot/compat

libcoot_compat_la_SOURCES = getopt.c getopt1.c

EXTRA_DIST = coot-sysdep.h sleep-fixups.h

