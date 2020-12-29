

if True:
    if coot_gui_api.main_menubar():
        menu = coot_menubar_menu("PISA")
        submenu_pisa = gtk.Menu()
        menuitem_pisa = gtk.MenuItem("PISA Assemblies...")

        add_simple_coot_menu_menuitem(
            submenu_pisa, "PISA assemblies...",
            lambda func:
            molecule_chooser_gui("Choose molecule for PISA assembly analysis",
                                 lambda imol:
                                 pisa_assemblies(imol)))

        add_simple_coot_menu_menuitem(
            submenu_pisa, "PISA interfaces...",
            lambda func:
            molecule_chooser_gui("Choose molecule for PISA interface analysis",
                                 lambda imol:
                                 pisa_interfaces(imol)))
