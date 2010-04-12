      int j = i+1;
      double r = 30;
      double theta_1 = i * angle_step * DEG_TO_RAD;
      double theta_2 = j * angle_step * DEG_TO_RAD;
      double pt_1_x = r *sin(theta_1) + dx_cen;
      double pt_1_y = r *cos(theta_1) + dy_cen;
      double pt_2_x = r *sin(theta_2) + dx_cen;
      double pt_2_y = r *cos(theta_2) + dy_cen;
      goo_canvas_polyline_new_line(root,
				   pt_1_x, pt_1_y,
				   pt_2_x, pt_2_y,
				   NULL);
      int atom_1 = mol.add_atom(pt_1_x, pt_1_y, " C", 0);
      int atom_2 = mol.add_atom(pt_2_x, pt_2_y, " C", 0);
      mol.add_bond(atom_1, atom_2, cen, lig_build::bond_t::SINGLE_BOND);
