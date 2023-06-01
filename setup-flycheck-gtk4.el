
;;; Code:

;;; how do I add "-Wno-unknown-pragmas" here?

;; remove USE_GUILE=1 for now
;; "USE_MOLECULES_TO_TRIANGLES=1"
;; "HAVE_GOOCANVAS=1" 
;;
(setq flycheck-gcc-definitions   (list "HAVE_CXX_THREAD=1" "HAVE_GSL=1" "VERSION=\"0.9.9-pre\"" "HAVE_BOOST=1" "PKGDATADIR=11111" "USE_PYTHON=1" "PYTHONDIR=/asdf" "PKGPYTHONDIR=asdf" "USE_LIBCURL=1" "BUILD_CURLEW=1" "COOT_BUILD_INFO_STRING=x"  "_OPENMP=1"))
(setq flycheck-clang-definitions (list "HAVE_CXX_THREAD=1" "HAVE_GSL=1" "VERSION=\"0.9.9-pre\"" "HAVE_BOOST=1" "PKGDATADIR=11111" "USE_PYTHON=1" "PYTHONDIR=/asdf" "PKGPYTHONDIR=asdf" "USE_LIBCURL=1" "BUILD_CURLEW=1" "COOT_BUILD_INFO_STRING=x" "_OPENMP=1"))

(setq flycheck-gcc-args   "-Wno-unknown-pragmas")
(setq flycheck-clang-args "-Wno-unknown-pragmas")

(setq build-path-list
  (list

   "." ".."
   "/usr/lib/x86_64-linux-gnu/glib-2.0/include"
   "/home/paule/autobuild/Linux-penelope-gtk4/include"
   "/home/paule/autobuild/Linux-penelope-gtk4/include/coot"
   "/home/paule/autobuild/Linux-penelope-gtk4/include/rdkit"
   "/home/paule/autobuild/Linux-penelope-gtk4/include/boost"
   "/home/paule/autobuild/Linux-penelope-gtk4/include/python3.9"
   "/home/paule/autobuild/Linux-penelope-gtk4/include/MoleculesToTriangles"
   "/home/paule/python3/include/python3.8"
   "/home/paule/glm/include"
   "/home/paule/goocanvas/include/goocanvas-2.0"
   "/home/paule/ogg-vorbis/include"
   "/home/paule/assimp/include"

   "/home/paule/python3/include/python3.9"
   "/home/paule/autobuild/Linux-penelope-gtk4-python/include"

   ;; from https://github.com/Wilfred/flycheck-pkg-config/issues/2
   "/usr/include/gtk-4.0"
   "/usr/include/graphene-1.0"
   "/home/paule/gtk/lib/x86_64-linux-gnu/graphene-1.0/include" ;; for graphene-config.h
   "/usr/include/libpng16" "/usr/include/gdk-pixbuf-2.0"
   "/usr/include/libdrm"
   "/usr/include/harfbuzz" "/usr/include/freetype2"
   "/usr/lib/glib-2.0/include" "/usr/include/glib-2.0"
   "/usr/include/libpng16" "/usr/include/pixman-1"
   "/usr/include/cairo" "/usr/include/atk-1.0" "/usr/include/pango-1.0"
   "/usr/include/gio-unix-2.0"
   "/usr/lib/dbus-1.0/include" "/usr/include/dbus-1.0"
   "/usr/include/at-spi-2.0"
   "/usr/include/at-spi2-atk/2.0"
   ))


(setq flycheck-gcc-include-path build-path-list)

(setq flycheck-clang-include-path build-path-list)

;;; setup-flycheck.el ends here
