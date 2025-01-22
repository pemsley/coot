
;;; Code:

;;; how do I add "-Wno-unknown-pragmas" here?

;; remove USE_GUILE=1 for now
;; "USE_MOLECULES_TO_TRIANGLES=1"
;; "HAVE_GOOCANVAS=1" 
;;
(setq flycheck-gcc-definitions   (list "HAVE_CXX_THREAD=1" "HAVE_GSL=1" "VERSION=\"0.9.9-pre\"" "HAVE_BOOST=1" "HAVE_BOOST_BASED_THREAD_POOL_LIBRARY=1" "PKGDATADIR=11111" "USE_PYTHON=1" "PYTHONDIR=/asdf" "PKGPYTHONDIR=asdf" "USE_LIBCURL=1" "BUILD_CURLEW=1" "COOT_BUILD_INFO_STRING=x" "_OPENMP=1" "MAKE_ENHANCED_LIGAND_TOOLS=1" "WITH_SOUND=1" "USE_GEMMI=1"))
(setq flycheck-clang-definitions (list "HAVE_CXX_THREAD=1" "HAVE_GSL=1" "VERSION=\"0.9.9-pre\"" "HAVE_BOOST=1" "HAVE_BOOST_BASED_THREAD_POOL_LIBRARY=1" "PKGDATADIR=11111" "USE_PYTHON=1" "PYTHONDIR=/asdf" "PKGPYTHONDIR=asdf" "USE_LIBCURL=1" "BUILD_CURLEW=1" "COOT_BUILD_INFO_STRING=x" "_OPENMP=1" "MAKE_ENHANCED_LIGAND_TOOLS=1" "WITH_SOUND=1" "USE_GEMMI=1"))

(setq flycheck-gcc-args   '("-Wno-unknown-pragmas" "-std=c++17"))
(setq flycheck-clang-args '("-Wno-unknown-pragmas" "-std=c++17"))
(setq flycheck-gcc-language-standard   "c++17")
(setq flycheck-clang-language-standard "c++17")

(setq build-path-list
  (list

   "." ".." "../.."
   "/usr/lib/x86_64-linux-gnu/glib-2.0/include"
   "/usr/local/include"
   "/usr/local/include/coot"
   "/usr/local/include/rdkit"
   "/usr/local/include/boost"
   "/usr/local/include/python3.11"
   "/usr/local/include/MoleculesToTriangles"
   "/usr/local/include/gsl"
   "/Users/pemsley/python3/include/python3.8"
   "/Users/pemsley/glm/include"
   "/Users/pemsley/ogg-vorbis/include"
   "/Users/pemsley/assimp/include"

   "/Users/pemsley/python3/include/python3.9"
   "/Users/pemsley/autobuild/build-for-chapi-gtk4/lib/python3.10/site-packages/nanobind/include"
   "/Users/pemsley/autobuild/Linux-penelope-gtk4-python/include"

   ;; from https://github.com/Wilfred/flycheck-pkg-config/issues/2
   "/Users/pemsley/gtk4/include/gtk-4.0"
   ;; "/usr/include/gtk-4.0"
   "/usr/include/graphene-1.0"
   "/Users/pemsley/gtk/lib/x86_64-linux-gnu/graphene-1.0/include" ;; for graphene-config.h
   "/usr/local/include/libpng16" "/usr/local/include/gdk-pixbuf-2.0"
   "/usr/local/include/libdrm"
   "/usr/local/include/harfbuzz" "/usr/local/include/freetype2"
   "/usr/lib/glib-2.0/include" "/usr/local/include/glib-2.0"
   "/usr/local/Cellar/glib/2.82.4/lib/glib-2.0/include"
   "/usr/local/Cellar/graphene/1.10.8/include/graphene-1.0"
   "/usr/local/Cellar/graphene/1.10.8/lib/graphene-1.0/include"
   "/usr/local/include/libpng16" "/usr/local/include/pixman-1"
   "/usr/local/include/cairo" "/usr/local/include/atk-1.0" "/usr/local/include/pango-1.0"
   "/usr/local/include/gio-unix-2.0"
   "/usr/lib/dbus-1.0/include" "/usr/local/include/dbus-1.0"
   "/usr/local/include/at-spi-2.0"
   "/usr/local/include/at-spi2-atk/2.0"
   "/usr/local/Cellar/python@3.12/3.12.8//Frameworks/Python.framework/Versions/3.12/include/python3.12"
   ))


(setq flycheck-gcc-include-path build-path-list)

(setq flycheck-clang-include-path build-path-list)

;;; setup-flycheck.el ends here
