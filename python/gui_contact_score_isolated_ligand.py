

add_simple_coot_menu_menuitem(coot_menubar_menu("Ligand"),
			      "Isolated Molprobity Dots for this Ligand",
                              lambda func: contact_score_ligand_func())

add_simple_coot_menu_menuitem(coot_menubar_menu("Ligand"),
                              "Isolated Coot Ligand Dots for this Ligand",
                              lambda func: coot_contact_dots_ligand_func())
