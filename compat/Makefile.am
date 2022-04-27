
lib_LTLIBRARIES = libcoot-compat.la
pkginclude_HEADERS = coot-getopt.h coot-sysdep.h coot-sysdep-c.h getopt_int.h coot-win32.h coot-posix.h

pkgincludedir = $(includedir)/coot/compat

libcoot_compat_la_SOURCES = getopt.cc getopt1.cc coot-sysdep-c.cc
if COOT_BUILD_POSIX
libcoot_compat_la_SOURCES += coot-posix.cc
endif
if COOT_BUILD_WINDOWS
libcoot_compat_la_SOURCES += coot-win32.cc
endif

libcoot_compat_la_LDFLAGS = $(SHARED_LDFLAGS)
