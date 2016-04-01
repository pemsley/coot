

add_simple_coot_menu_menuitem(coot_menubar_menu("Ligand"),
			      "Isolated dots for this ligand",
                              lambda func: contact_score_ligand_func())
