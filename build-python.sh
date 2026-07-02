#!/bin/bash

# build-python.sh
# Extracted from build-it-3-3: builds libffi, openssl (if needed), and Python.
#
# Edit the variables below to match your environment.

# ========================= Configuration =========================

# Where everything gets installed
install_top_dir=${install_top_dir:=$HOME/test-python}

# Where source tarballs are downloaded/cached
sources_dir=${sources_dir:=$HOME/autobuild/sources}

# Where build logs go
log_dir=${log_dir:=$HOME/autobuild/logs}

# Versions
python_version=python3.14
pver=3.14.2
libffi_version=3.4.6
gettext_version=0.22.5

# Set to "true" to skip building components that are already installed
skip_if_installed=true

# OpenSSL prefix — use system/Homebrew OpenSSL rather than building it.
# On macOS with Homebrew: /opt/homebrew/opt/openssl@3
# On Linux with system OpenSSL: /usr
openssl_prefix=${openssl_prefix:=/opt/homebrew/opt/openssl@3}

# ========================= Derived variables =========================

OS=$(uname)
MAKE=${MAKE:=make}
WGET="wget -nv -N -P ${sources_dir}"
generic_prefix="--prefix=$install_top_dir"

mkdir -p "$sources_dir" "$log_dir" "$install_top_dir"

export PATH=$install_top_dir/bin:$PATH
export PKG_CONFIG_PATH=$install_top_dir/lib/pkgconfig:$install_top_dir/lib64/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$install_top_dir/lib:$install_top_dir/lib64:$LD_LIBRARY_PATH

# Detect distro for platform-specific logic
dist_name=unknown
if [ $OS = Linux ] ; then
   if [ -e /etc/os-release ] ; then
      dist_name=$(grep ^ID= /etc/os-release | sed -e 's/ID=//' -e 's/"//g')
   fi
fi

# ========================= libffi =========================

build_libffi=true
if [ "$skip_if_installed" = true ] ; then
   if [ -e $install_top_dir/lib/pkgconfig/libffi.pc ] || \
      [ -e $install_top_dir/include/ffi.h ] ; then
      build_libffi=false
      echo "libffi already installed, skipping"
   fi
fi

if [ "$build_libffi" = true ] ; then
   echo "========================= building libffi $libffi_version ========================="
   (
   unset LD_LIBRARY_PATH
   ${WGET} https://github.com/libffi/libffi/releases/download/v$libffi_version/libffi-$libffi_version.tar.gz
   if [ -e ${sources_dir}/libffi-$libffi_version.tar.gz ] ; then
      tar xf ${sources_dir}/libffi-$libffi_version.tar.gz
   else
      echo "ERROR: failed to get libffi tarball"
      exit 1
   fi
   cd libffi-$libffi_version
   if ./configure --prefix=$install_top_dir --enable-portable-binary --enable-shared --disable-static ; then
      if $MAKE ; then
         make install
      fi
   fi
   cd -
   ) > $log_dir/libffi.txt 2>&1
   echo "done libffi"
fi

# ========================= openssl =========================

echo "Using OpenSSL from $openssl_prefix"
if [ ! -d "$openssl_prefix/include/openssl" ] ; then
   echo "WARNING: $openssl_prefix/include/openssl not found — Python may build without ssl support"
fi

# ========================= gettext (macOS only, for libintl) =========================

build_gettext=false
if [ $OS = Darwin ] ; then
   if [ ! -e $install_top_dir/include/libintl.h ] ; then
      build_gettext=true
   fi
fi

if [ "$build_gettext" = true ] ; then
   echo "========================= building gettext $gettext_version ========================="
   (
   unset LD_LIBRARY_PATH
   ${WGET} http://www.mirrorservice.org/sites/ftp.gnu.org/gnu/gettext/gettext-$gettext_version.tar.gz
   tar xf $sources_dir/gettext-$gettext_version.tar.gz
   cd gettext-$gettext_version
   ./configure --prefix=$install_top_dir --disable-java
   $MAKE && make install
   cd -
   ) > $log_dir/gettext.txt 2>&1
   echo "done gettext"
else
   echo "not building gettext"
fi

# ========================= Python =========================

build_python=true
if [ "$skip_if_installed" = true ] ; then
   if [ -e $install_top_dir/bin/$python_version ] ; then
      installed_ver=$($install_top_dir/bin/$python_version -V 2>&1)
      if [ "$installed_ver" = "Python $pver" ] ; then
         if $install_top_dir/bin/$python_version -c 'import ctypes' 2>/dev/null ; then
            echo "$python_version ($pver) already installed with working ctypes, skipping"
            build_python=false
         else
            echo "$python_version ctypes failure, rebuilding"
         fi
      else
         echo "wrong Python version ($installed_ver), rebuilding"
      fi
   fi
fi

if [ "$build_python" = true ] ; then
   echo "========================= building Python $pver ========================="
   (
   date
   unset LD_LIBRARY_PATH

   ${WGET} https://www.python.org/ftp/python/$pver/Python-$pver.tgz
   if [ ! -e ${sources_dir}/Python-$pver.tgz ] ; then
      echo "ERROR: failed to get Python tarball"
      exit 1
   fi

   tar xf ${sources_dir}/Python-$pver.tgz
   cd Python-$pver

   python_ldflags=""
   python_libs=""
   python_spec_flags=""
   if [ "$dist_name" = fedora ] ; then
      python_ldflags="LDFLAGS=-Wl,--rpath=$install_top_dir/lib"
      python_libs="LIBS=-L$install_top_dir/lib64"
   elif [ $OS = Linux ] ; then
      python_ldflags="LDFLAGS=-Wl,--rpath=$install_top_dir/lib"
      python_libs="LIBS=-L$install_top_dir/lib64"
   fi

   echo "INFO:: ./configure $generic_prefix --with-ensurepip=install --enable-shared" \
        "--with-openssl=$openssl_prefix CPPFLAGS=-I$install_top_dir/include" \
        "$python_libs $python_ldflags $python_spec_flags"

   if ./configure $generic_prefix --with-ensurepip=install --enable-shared \
                  --with-openssl=$openssl_prefix \
                  CPPFLAGS=-I$install_top_dir/include \
                  $python_libs $python_ldflags $python_spec_flags ; then
      echo "INFO:: make for python"
      if $MAKE ; then
         echo "INFO:: make install for python"
         $MAKE install
      fi
   fi
   cd -
   ) > $log_dir/python.txt 2>&1

   # symlink python3 and python
   if [ -e $install_top_dir/bin/$python_version ] ; then
      ln -sf $install_top_dir/bin/$python_version $install_top_dir/bin/python3
      ln -sf $install_top_dir/bin/python3 $install_top_dir/bin/python
   fi

   echo "done Python $pver"

   # upgrade pip and install setuptools
   echo "upgrading pip..."
   (
   $install_top_dir/bin/$python_version -m pip install --upgrade pip
   $install_top_dir/bin/$python_version -m pip install setuptools
   ) > $log_dir/pip-upgrade.txt 2>&1
   echo "done pip upgrade"
else
   echo "not building Python"
fi

# ========================= test =========================

echo "testing $install_top_dir/bin/$python_version ..."
if $install_top_dir/bin/$python_version -c 'import ctypes; print("ctypes OK")' ; then
   echo "Python is working"
else
   echo "ERROR: Python not working"
   exit 2
fi
