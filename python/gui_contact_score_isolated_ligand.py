

add_simple_coot_menu_menuitem(coot_menubar_menu("Ligand"),
			      "Isolated dots for this ligand",
                              lambda func: contact_score_ligand_func())

add_simple_coot_menu_menuitem(coot_menubar_menu("Ligand"),
                              "Coot Isolated Ligand Dots",
                              lambda func: coot_contact_dots_ligand_func())
