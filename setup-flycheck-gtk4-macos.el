
;;; Code:

;;; how do I add "-Wno-unknown-pragmas" here?

;; remove USE_GUILE=1 for now
;; "USE_MOLECULES_TO_TRIANGLES=1"
;; "HAVE_GOOCANVAS=1"
;;
(setq flycheck-clang-definitions (list "HAVE_CXX_THREAD=1" "HAVE_GSL=1" "VERSION=\"0.9.9-pre\"" "HAVE_BOOST=1" "HAVE_BOOST_BASED_THREAD_POOL_LIBRARY=1" "PKGDATADIR=11111" "USE_PYTHON=1" "PYTHONDIR=/asdf" "PKGPYTHONDIR=asdf" "USE_LIBCURL=1" "BUILD_CURLEW=1" "COOT_BUILD_INFO_STRING=x" "_OPENMP=1" "MAKE_ENHANCED_LIGAND_TOOLS=1" "WITH_SOUND=1" "USE_GEMMI=1" "USE_MOLECULES_TO_TRIANGLES=1"))

(setq flycheck-clang-args '("-Wno-unknown-pragmas" "-std=c++17"))
(setq flycheck-clang-language-standard "c++17")

(setq build-path-list
  (list

   "." ".." "../.."

   "/opt/homebrew/include"
   "/opt/homebrew/include/coot"
   "/opt/homebrew/include/rdkit"
   "/opt/homebrew/include/boost"
   "/opt/homebrew/include/python3.14"
   "/opt/homebrew/include/MoleculesToTriangles"
   "/opt/homebrew/include/gsl"
   "/opt/homebrew/include/gtk-4.0"
   "/opt/homebrew/Cellar/python@3.13/3.13.5/Frameworks/Python.framework/Versions/3.13/include/python3.13"
   "/opt/homebrew/Cellar/nanobind//2.7.0/share/nanobind/include"
   "/opt/homebrew/Cellar/glm/1.0.1/include"
   "/opt/homebrew/Cellar/gtk4/4.18.6/include/gtk-4.0"
   "/opt/homebrew/Cellar/glib/2.84.3/include/glib-2.0"
   "/opt/homebrew/Cellar/glib/2.84.3/lib/glib-2.0/include"
   "/opt/homebrew/Cellar/cairo/1.18.4/include/cairo"
   "/opt/homebrew/Cellar/pango/1.56.3/include/pango-1.0"
   "/opt/homebrew/Cellar/harfbuzz/11.2.1/include/harfbuzz"
   "/opt/homebrew/Cellar/gdk-pixbuf/2.42.12_1/include/gdk-pixbuf-2.0"
   "/opt/homebrew/Cellar/graphene/1.10.8/include/graphene-1.0"
   "/opt/homebrew/Cellar/graphene/1.10.8/lib/graphene-1.0/include"
   "/opt/homebrew/Cellar/freetype/2.13.3/include/freetype2"

   "/Users/pemsley/glm/include"
   "/Users/pemsley/ogg-vorbis/include"
   "/Users/pemsley/assimp/include"

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
   "/usr/local/include/glib-2.0"
   "/usr/local/Cellar/glib/2.82.4/lib/glib-2.0/include"
   "/usr/local/Cellar/graphene/1.10.8/include/graphene-1.0"
   "/usr/local/Cellar/graphene/1.10.8/lib/graphene-1.0/include"
   "/usr/local/include/libpng16" "/usr/local/include/pixman-1"
   "/usr/local/include/cairo" "/usr/local/include/atk-1.0" "/usr/local/include/pango-1.0"
   "/usr/local/include/gio-unix-2.0"
   "/usr/local/include/dbus-1.0"
   "/usr/local/include/at-spi-2.0"
   "/usr/local/include/at-spi2-atk/2.0"
   "/usr/local/Cellar/python@3.12/3.12.8//Frameworks/Python.framework/Versions/3.12/include/python3.12"
   ))


(setq flycheck-gcc-include-path build-path-list)

(setq flycheck-clang-include-path build-path-list)

;;; setup-flycheck.el ends here
