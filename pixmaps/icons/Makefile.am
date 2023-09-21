
EXTRA_DIST = hicolor_apps_48x48_coot.png   \
	     hicolor_apps_64x64_coot.png   \
	     hicolor_apps_96x96_coot.png   \
	     hicolor_apps_256x256_coot.png \
	     hicolor_apps_scalable_coot.svg \
		 hicolor_apps_scalable_coot-layla.svg \
		 hicolor_actions_scalable_Stock-dialog-warning.svg \
		 hicolor_actions_scalable_Stock-dialog-information.svg

# This needs to be installed in Coot.app/Contents/Resources, only on macOS though.
# Currently this provides a blank icon (i.e. it doesn't work) when viewing the Coot.app
# folder  with the Finder.
# To make the .icns file:
# https://gist.github.com/jamieweavis/b4c394607641e1280d447deed5fc85fc
# but am I doing something wrong with the installation?
MACOS_ICON = Coot.icns

# is this install-data-local thing needed now? Check with other software that install
# desktop icons.
#
install-data-local:
	for icon in $(EXTRA_DIST); do \
		THEME=`   echo $$icon | cut -d_ -f1`; \
		CONTEXT=` echo $$icon | cut -d_ -f2`; \
		SIZE=`    echo $$icon | cut -d_ -f3`; \
		ICONFILE=`echo $$icon | cut -d_ -f4`; \
		mkdir -p $(DESTDIR)$(datadir)/icons/$$THEME/$$SIZE/$$CONTEXT; \
		$(INSTALL_DATA) $(srcdir)/$$icon $(DESTDIR)$(datadir)/icons/$$THEME/$$SIZE/$$CONTEXT/$$ICONFILE; \
	done;
