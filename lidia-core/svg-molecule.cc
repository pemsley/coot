
#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#include "svg-molecule.hh"
#include "rdkit-interface.hh"

// put this into (new file) svg-molecule.cc when you're happy with it
void
svg_molecule_t::import_rdkit_mol(RDKit::ROMol *rdkm, int iconf) {

   // we don't fiddle with the hydrogens here, we just draw what we are given.
   //
   // So typically, user has called Chem.RemoveHs() or remove_non_polar_Hs()
   // and MolOps::Kekulize() and WedgeMolBonds() before calling this function.
   
   int n_conf  = rdkm->getNumConformers();
   if (iconf < n_conf) {
      const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();

      RDKit::Conformer &conf = rdkm->getConformer(iconf);
      unsigned int n_mol_atoms = rdkm->getNumAtoms();

      // determine the centre correction
      double sum_x = 0;
      double sum_y = 0;
      double min_y = 9e9;
      for (unsigned int iat=0; iat<n_mol_atoms; iat++) {
	 const RDKit::Atom *at_p = (*rdkm)[iat];
	 RDGeom::Point3D &r_pos = conf.getAtomPos(iat);
	 sum_x += r_pos.x;
	 sum_y += r_pos.y;
	 if (r_pos.y < min_y)
	    min_y = r_pos.y;
      }

      // set the scale correction
      // 
      std::vector<double> bond_lengths;
      for (unsigned int i=0; i<rdkm->getNumBonds(); i++) {
	 const RDKit::Bond *bond_p = rdkm->getBondWithIdx(i);
	 int idx_1 = bond_p->getBeginAtomIdx();
	 int idx_2 = bond_p->getEndAtomIdx();
	 if ( (*rdkm)[idx_1]->getAtomicNum() != 1) {
	    if ( (*rdkm)[idx_2]->getAtomicNum() != 1) {
	       RDGeom::Point3D &r_pos_1 = conf.getAtomPos(idx_1);
	       RDGeom::Point3D &r_pos_2 = conf.getAtomPos(idx_2);
	       clipper::Coord_orth p1(r_pos_1.x, r_pos_1.y, r_pos_1.z);
	       clipper::Coord_orth p2(r_pos_2.x, r_pos_2.y, r_pos_2.z);
	       double l = clipper::Coord_orth::length(p1, p2);
	       bond_lengths.push_back(l);
	    }
	 }
      }
      if (bond_lengths.size() > 0) {
	 std::sort(bond_lengths.begin(), bond_lengths.end());
	 int index = bond_lengths.size()/2;
	 double bll = bond_lengths[index];
	 double scale = 1.0/bll;
	 // scale_correction.first = 1;
	 // scale_correction.second = scale;
      }

      for (unsigned int iat=0; iat<n_mol_atoms; iat++) {
	 const RDKit::Atom *at_p = (*rdkm)[iat];
	 RDGeom::Point3D &r_pos = conf.getAtomPos(iat);
	 std::string name = "";
	 try {
	    at_p->getProp("name", name);
	 }
	 catch (const KeyErrorException &kee) {
	    // we don't need to see these.  We get them when reading an mdl file
	    // (for example).
	    // std::cout << "caught no-name for atom exception in import_rdkit_mol(): "
	    // <<  kee.what() << std::endl;
	 }
	 try {
	    at_p->getProp("molFileAlias", name);
	 }
	 catch (const KeyErrorException &kee) { }
	 // lig_build::pos_t pos = input_coords_to_canvas_coords(cp);
	 lig_build::pos_t pos(r_pos.x, r_pos.y);
	 int n = at_p->getAtomicNum();
	 std::string element = tbl->getElementSymbol(n);
	 int charge = at_p->getFormalCharge();
	 if (false)
	    std::cout << "found atom " << iat <<  " with element \"" << element
		      << "\" with formal charge " << charge << std::endl;
	 svg_atom_t mol_at(pos, element, charge);

	 mol_at.charge = charge;
	 if (! name.empty())
	    mol_at.atom_name = name;
	 std::pair<bool,int> added = add_atom(mol_at);
      }

      unsigned int n_bonds = rdkm->getNumBonds();
      // std::cout << "considering " << n_bonds << " bonds" << std::endl;
      for (unsigned int ib=0; ib<n_bonds; ib++) {
	 const RDKit::Bond *bond_p = rdkm->getBondWithIdx(ib);
	 int idx_1 = bond_p->getBeginAtomIdx();
	 int idx_2 = bond_p->getEndAtomIdx();
	 lig_build::bond_t::bond_type_t bt = coot::convert_bond_type(bond_p->getBondType());
	 if (false)
	    std::cout << "on import got bt " << bt << " from " << bond_p->getBondType()
		      << std::endl;

	 try { 
	    const svg_atom_t &cat1 = atoms[idx_1];
	    const svg_atom_t &cat2 = atoms[idx_2];
	    bool shorten_first  = false;
	    bool shorten_second = false;
	    if ( (*rdkm)[idx_1]->getAtomicNum() != 6) {
	       shorten_first = true;
	    } else {
	       // how many bonds does this atom have?
	       // if (mol[idx_1].getNumBonds() == 1)
	       //  shorten_first = true;
	    }
	    if ( (*rdkm)[idx_2]->getAtomicNum() != 6) {
	       shorten_second = true;
	    } else {
	    }

	    // do we need to pass idx_1 and idx_2, neighbour info?
            // 	    cairo_bond_t bond(idx_1, idx_2, cat1, cat2, shorten_first, shorten_second,
            // 			      bt, empty, empty);
	    svg_bond_t bond(idx_1, idx_2, bt);
	    RDKit::Bond::BondDir bond_dir = bond_p->getBondDir();
	    if (false)
	       std::cout << "bond " << ib << ":  type " << bt
			 << " between " << idx_1 << " at "
			 << conf.getAtomPos(idx_1)
			 << " and " << idx_2 << " at "
			 << conf.getAtomPos(idx_2)
			 << " dir " << bond_dir << std::endl;
	    if (bond_dir != RDKit::Bond::NONE) {
	       if (bond_dir == RDKit::Bond::BEGINWEDGE) {
		  bond.set_bond_type(lig_build::bond_t::OUT_BOND);
	       }
	       if (bond_dir == RDKit::Bond::BEGINDASH)
		  bond.set_bond_type(lig_build::bond_t::IN_BOND);
	    }
	    add_bond(bond);
	 }
	 catch (...) {
	    std::cout << "WARNING:: problem. scrambled input molecule? numbers of atoms: ";
	    std::cout << rdkm->getNumAtoms() << " vs "
		      << get_number_of_atoms_including_hydrogens() << std::endl;
	 }
      }
      assign_ring_centres();
   }

   median_bond_length_ = median_bond_length();
}


// This is for a double bond in a ring
//
// pos_1 and pos_2 are the coordinates of the atoms, not the pre-shorted
// coordinates
//
std::string
svg_bond_t::draw_double_in_ring_bond(const lig_build::pos_t &pos_1_in,
                                     const lig_build::pos_t &pos_2_in,
                                     bool shorten_first,
                                     bool shorten_second,
                                     const lig_build::pos_t &centre,
                                     double scale, bool dashed_inner) {

   std::string s;

   lig_build::pos_t pos_1 = pos_1_in;
   lig_build::pos_t pos_2 = pos_2_in;

   double shorten_fraction = 0.74;
   if (shorten_first)
      pos_1 = lig_build::pos_t::fraction_point(pos_2_in, pos_1_in, shorten_fraction);
   if (shorten_second)
      pos_2 = lig_build::pos_t::fraction_point(pos_1_in, pos_2_in, shorten_fraction);

   // outside "normal" bond
   lig_build::pos_t p1 = svg_molecule_t::mol_coords_to_svg_coords(pos_1, centre, scale);
   lig_build::pos_t p2 = svg_molecule_t::mol_coords_to_svg_coords(pos_2, centre, scale);

   // for inner bond
   std::pair<lig_build::pos_t, lig_build::pos_t> p = 
      make_double_aromatic_short_stick(pos_1_in, pos_2_in, shorten_first, shorten_second);

   // cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
   // cairo_move_to(cr, p1.x, p1.y);
   // cairo_line_to(cr, p2.x, p2.y);
   // cairo_stroke(cr);

   s += make_bond_line_string(p1, p2);

   if (dashed_inner) {
      double dashlength = 0.015; // 0.01 is also fine
      // cairo_set_dash(cr, &dashlength, 1, 0);
      std::cout << "draw_double_in_ring_bond(): set dash inner here " << std::endl; // 20221130-PE 
   }
   p1 = svg_molecule_t::mol_coords_to_svg_coords(p.first,  centre, scale);
   p2 = svg_molecule_t::mol_coords_to_svg_coords(p.second, centre, scale);

   // cairo_move_to(cr, p1.x, p1.y);
   // cairo_line_to(cr, p2.x, p2.y);
   // cairo_stroke(cr);

   s += make_bond_line_string(p1, p2);

   if (dashed_inner) {
      //cairo_set_dash(cr, NULL, 0, 0); // restore
   }

   // std::cout << "draw_double_in_ring_bond return s " << s << std::endl;
   return s;
}

// not in-ring
//
// shorten_first and shorten_second are set depending on the element of the atoms and number of bonds.
std::string
svg_bond_t::draw_double_bond(const lig_build::atom_t &at_1,
                             const lig_build::atom_t &at_2,
                             bool shorten_first, bool shorten_second,
                             const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_first_atom,
                             const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
                             const lig_build::pos_t &centre, double scale) {

   std::string s;

   lig_build::pos_t pos_1 = at_1.atom_position;
   lig_build::pos_t pos_2 = at_2.atom_position;
   double shorten_fraction = 0.74;
   if (shorten_first)
      pos_1 = lig_build::pos_t::fraction_point(at_2.atom_position, at_1.atom_position, shorten_fraction);
   if (shorten_second)
      pos_2 = lig_build::pos_t::fraction_point(at_1.atom_position, at_2.atom_position, shorten_fraction);

   // the boring case (symmetrically, same-length both displaced) (e.g. N=N, C=C) occurs when:
   // not
   // at least one C atom with at least one bond and not one C atom with more than one bond.
   //
   bool do_boring_case = false;

   if ((other_connections_to_second_atom.size() == 0) &&
       (other_connections_to_first_atom.size()  == 0)) {
      do_boring_case = true;
   }

   unsigned int n_C = 0;
   if (at_1.element == "C") n_C++;
   if (at_2.element == "C") n_C++;

   if (n_C == 0) {
      do_boring_case = true;
   } else {
      // more than 1 Carbon
      if (n_C == 2) {
	 if (other_connections_to_first_atom.size() == 0)
	    if (other_connections_to_second_atom.size() == 0)
	       do_boring_case = true;
      }
      if (n_C == 1) {
	 if (at_1.element == "C")
	    if (other_connections_to_first_atom.size() > 1)
	       do_boring_case = true;
	 if (at_2.element == "C")
	    if (other_connections_to_second_atom.size() > 1)
	       do_boring_case = true;
      }
   }


   if (do_boring_case) {

      // std::cout << "boring case" << std::endl;

      // boring case :-)

      std::pair<std::pair<lig_build::pos_t, lig_build::pos_t>, std::pair<lig_build::pos_t, lig_build::pos_t> > p =
	 make_double_bond(at_1.atom_position, at_2.atom_position, shorten_first, shorten_second);

      lig_build::pos_t p1 = svg_molecule_t::mol_coords_to_svg_coords(p.first.first,  centre, scale);
      lig_build::pos_t p2 = svg_molecule_t::mol_coords_to_svg_coords(p.first.second, centre, scale);

      // cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
      // cairo_move_to(cr, p1.x, p1.y);
      // cairo_line_to(cr, p2.x, p2.y);
      // cairo_stroke(cr);

      s += make_bond_line_string(p1, p2);

      p1 = svg_molecule_t::mol_coords_to_svg_coords(p.second.first,  centre, scale);
      p2 = svg_molecule_t::mol_coords_to_svg_coords(p.second.second, centre, scale);

      // cairo_move_to(cr, p1.x, p1.y);
      // cairo_line_to(cr, p2.x, p2.y);
      // cairo_stroke(cr);

      s += make_bond_line_string(p1, p2);

   } else {

      // std::cout << "dobule bond elegant case" << std::endl;
      // elegant case, offset, shorten a bond
      std::pair<std::pair<lig_build::pos_t, lig_build::pos_t>, std::pair<lig_build::pos_t, lig_build::pos_t> > bonds =
	 make_double_bond(at_1.atom_position, at_2.atom_position, shorten_first, shorten_second,
			  other_connections_to_first_atom, other_connections_to_second_atom);

      lig_build::pos_t p1 = svg_molecule_t::mol_coords_to_svg_coords(bonds.first.first,  centre, scale);
      lig_build::pos_t p2 = svg_molecule_t::mol_coords_to_svg_coords(bonds.first.second, centre, scale);
      // cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
      // cairo_move_to(cr, p1.x, p1.y);
      // cairo_line_to(cr, p2.x, p2.y);

      s += make_bond_line_string(p1, p2);

      p1 = svg_molecule_t::mol_coords_to_svg_coords(bonds.second.first,  centre, scale);
      p2 = svg_molecule_t::mol_coords_to_svg_coords(bonds.second.second, centre, scale);
      // cairo_move_to(cr, p1.x, p1.y);
      // cairo_line_to(cr, p2.x, p2.y);
      // cairo_stroke(cr);

      s += make_bond_line_string(p1, p2);
   }
   return s;
}



// static
lig_build::pos_t
svg_molecule_t::svg_molecule_t::mol_coords_to_svg_coords(const lig_build::pos_t &pos_1,
                                                         const lig_build::pos_t &centre,
                                                         double scale) {

   lig_build::pos_t p1 = (pos_1 - centre) * scale;
   p1.y = -p1.y; // canvas is upside down c.f. normal/real-world/molecule coordinates
   p1 += lig_build::pos_t(0.5,0.5);
   return p1;

}

std::string
svg_bond_t::make_bond_line_string(const lig_build::pos_t &p1, const lig_build::pos_t &p2) const {

   double sf = 400.0; // scale factor
   std::string s;
   s += "   <line x1=\"";
   s += std::to_string(sf * p1.x);
   s += "\" y1=\"";
   s += std::to_string(sf * p1.y);
   s += "\" x2=\"";
   s += std::to_string(sf * p2.x);
   s += "\" y2=\"";
   s += std::to_string(sf * p2.y);
   s += "\"";
   s += " style=\"stroke:#202020; stroke-width:2; fill:none; stroke-linecap:round;\" />\n";
   return s;
}


std::string
svg_bond_t::draw_bond(const svg_atom_t &at_1, const svg_atom_t &at_2,
                      bool at_1_in_ring_flag, bool at_2_in_ring_flag,
                      lig_build::bond_t::bond_type_t bt,
                      bool shorten_first, bool shorten_second,
                      const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_first_atom,
                      const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
                      const lig_build::pos_t &centre, double scale) {


   std::string s;

   lig_build::pos_t pos_1_in = at_1.atom_position;
   lig_build::pos_t pos_2_in = at_2.atom_position;

   lig_build::pos_t pos_1 = pos_1_in;
   lig_build::pos_t pos_2 = pos_2_in;

   

   // fraction_point() returns a point that is (say) 0.8 of the way
   // from p1 (first arg) to p2 (second arg).
   //
   double shorten_fraction = 0.74;
   if (shorten_first)
      pos_1 = lig_build::pos_t::fraction_point(pos_2_in, pos_1_in, shorten_fraction);
   if (shorten_second)
      pos_2 = lig_build::pos_t::fraction_point(pos_1_in, pos_2_in, shorten_fraction);

   // std::cout << "draw_bond for bt " << bt << " between " << at_1 << " and " << at_2 << std::endl;

   switch (bt) {
   case lig_build::bond_t::SINGLE_BOND:
   case lig_build::bond_t::SINGLE_OR_DOUBLE:
   case lig_build::bond_t::SINGLE_OR_AROMATIC:
   case lig_build::bond_t::DELOC_ONE_AND_A_HALF:
   case lig_build::bond_t::BOND_ANY:
      {
         //  std::cout << "draw_bond single" << std::endl;
	 // push down the decision making for the shortening
	 bool at_1_is_singleton = false;
	 bool at_2_is_singleton = false;
	 if (other_connections_to_first_atom.size()  == 0) at_1_is_singleton = true;
	 if (other_connections_to_second_atom.size() == 0) at_2_is_singleton = true;
	 std::pair<lig_build::pos_t, lig_build::pos_t> bp =
	    coords_for_single_bond(at_1, at_2, at_1_is_singleton, at_2_is_singleton);
	 pos_1 = bp.first;
	 pos_2 = bp.second;
         lig_build::pos_t p1 = svg_molecule_t::mol_coords_to_svg_coords(pos_1, centre, scale);
         lig_build::pos_t p2 = svg_molecule_t::mol_coords_to_svg_coords(pos_2, centre, scale);
	 std::cout << "pos_1 " << pos_1 << " pos_2 " << pos_2
		   << " bp.first " << bp.first << " bp.second " << bp.second
		   << " p1 " << p1 << " p2 " << p2 << std::endl;

	 // cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
	 // cairo_move_to(cr, p1.x, p1.y);
	 // cairo_line_to(cr, p2.x, p2.y);
	 // cairo_stroke(cr);

         std::string bond_string = make_bond_line_string(p1, p2);

         s += bond_string;
      }
      break;

   case lig_build::bond_t::DOUBLE_BOND:

   case lig_build::bond_t::DOUBLE_OR_AROMATIC:
      {
	 if (have_centre_pos()) {
            // std::cout << "draw_double in ring bond " << std::endl;
	    s += draw_double_in_ring_bond(pos_1_in, pos_2_in, shorten_first, shorten_second,
                                          centre, scale);
	 } else {
            // std::cout << "draw_double bond " << std::endl;
	    s += draw_double_bond(at_1, at_2,
			     shorten_first, shorten_second,
			     other_connections_to_first_atom,
			     other_connections_to_second_atom,
			     centre, scale);
	 }
      }
      break;

   // other bond types

   case lig_build::bond_t::AROMATIC_BOND:
      if (have_centre_pos()) {
         std::cout << "draw aromatic bond " << std::endl;
	 bool dashed_inner = true;
	 s += draw_double_in_ring_bond(pos_1_in, pos_2_in, shorten_first, shorten_second,
                                       centre, scale, dashed_inner);
      }

      break;

   case lig_build::bond_t::TRIPLE_BOND:
      {
	 lig_build::pos_t buv = (pos_2-pos_1).unit_vector();
	 lig_build::pos_t buv_90 = buv.rotate(90);
	 double small = 0.0125/scale;
	 lig_build::pos_t p1 = pos_1 + buv_90 * small;
	 lig_build::pos_t p2 = pos_2 + buv_90 * small;
	 lig_build::pos_t p3 = pos_1;
	 lig_build::pos_t p4 = pos_2;
	 lig_build::pos_t p5 = pos_1 - buv_90 * small;
	 lig_build::pos_t p6 = pos_2 - buv_90 * small;

	 lig_build::pos_t sc_p1 = svg_molecule_t::mol_coords_to_svg_coords(p1, centre, scale);
	 lig_build::pos_t sc_p2 = svg_molecule_t::mol_coords_to_svg_coords(p2, centre, scale);
	 lig_build::pos_t sc_p3 = svg_molecule_t::mol_coords_to_svg_coords(p3, centre, scale);
	 lig_build::pos_t sc_p4 = svg_molecule_t::mol_coords_to_svg_coords(p4, centre, scale);
	 lig_build::pos_t sc_p5 = svg_molecule_t::mol_coords_to_svg_coords(p5, centre, scale);
	 lig_build::pos_t sc_p6 = svg_molecule_t::mol_coords_to_svg_coords(p6, centre, scale);

	 // cairo_move_to(cr, sc_p1.x, sc_p1.y);
	 // cairo_line_to(cr, sc_p2.x, sc_p2.y);
	 // cairo_stroke(cr);
	 // cairo_move_to(cr, sc_p3.x, sc_p3.y);
	 // cairo_line_to(cr, sc_p4.x, sc_p4.y);
	 // cairo_stroke(cr);
	 // cairo_move_to(cr, sc_p5.x, sc_p5.y);
	 // cairo_line_to(cr, sc_p6.x, sc_p6.y);
	 // cairo_stroke(cr);
      }
      break;

   case IN_BOND:
      {
	 // set of lines
	 std::vector<std::pair<lig_build::pos_t, lig_build::pos_t> > vp =
	    lig_build::pos_t::make_wedge_in_bond(pos_1, pos_2);
	 if (vp.size()) {
	    // cairo_set_source_rgb(0.1, 0.1, 0.1);
	    for (unsigned int i=0; i<vp.size(); i++) { 

	       lig_build::pos_t p1 = svg_molecule_t::mol_coords_to_svg_coords(vp[i].first,  centre, scale);
	       lig_build::pos_t p2 = svg_molecule_t::mol_coords_to_svg_coords(vp[i].second, centre, scale);

	       // cairo_move_to(p1.x, p1.y);
	       // cairo_line_to(p2.x, p2.y);

               std::string bond_string = make_bond_line_string(p1, p2);
               s += bond_string;
	    }
            // cairo_stroke(cr);
	 }
      }
      break;

   case OUT_BOND:
      {

	 bool done_darted = false;
	 if (other_connections_to_second_atom.size() > 0) {
	    // don't make sheared or darted wedge bonds to atoms that are not Carbon
	    // don't make sheared or darted wedge bond to a C that is tetrahedral
	    bool draw_dart_or_wedge = false;
	    if (at_2.element == "C")
	       if (other_connections_to_second_atom.size() <= 2)
		  draw_dart_or_wedge = true;
	    if (draw_dart_or_wedge) {
               std::string bond_string = draw_sheared_or_darted_wedge_bond(pos_1, pos_2, other_connections_to_second_atom,
                                                                           centre, scale);
               s += bond_string;
	       done_darted = true;
	    }
	 }

	 if (! done_darted) {
	    // filled shape (normal wedge)
	    std::vector<lig_build::pos_t> v = lig_build::pos_t::make_wedge_out_bond(pos_1, pos_2);

	    // lig_build::pos_t p = svg_molecule_t::mol_coords_to_svg_coords(v[0], centre, scale);
	    // cairo_move_to(cr, p.x, p.y);
	    // cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
	    for (unsigned int i=1; i<v.size(); i++) {
	       // lig_build::pos_t p_i = svg_molecule_t::mol_coords_to_svg_coords(v[i], centre, scale);
	       // cairo_line_to(cr, p_i.x, p_i.y);
               std::string bond_string = "   <polygon points=\"";
               for (unsigned int i=0; i<v.size(); i++) {
                  double sf = 400.0; // scale_factor
                  lig_build::pos_t p_i = svg_molecule_t::mol_coords_to_svg_coords(v[i], centre, scale);         
                  bond_string += std::to_string(sf * p_i.x);
                  bond_string += ",";
                  bond_string += std::to_string(sf * p_i.y);
                  bond_string += " ";
               }
               bond_string += "\" style=\"fill:#2020202;\" />\n";
               s += bond_string;
	    }
	    // cairo_close_path(cr);
	    // cairo_fill(cr);
	    // cairo_stroke(cr);
	 }
      }
      break;
   case BOND_UNDEFINED:
      break;
   }
   return s;

}

std::string
svg_bond_t::draw_sheared_or_darted_wedge_bond(const lig_build::pos_t &pos_1,
                                              const lig_build::pos_t &pos_2,
                                              const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
                                              const lig_build::pos_t &centre,
                                              double scale) const {

   std::string s;

   std::vector<lig_build::pos_t> v = coords_for_sheared_or_darted_wedge_bond(pos_1, pos_2, other_connections_to_second_atom);

   // if there is a triple bond connected to the second atom, we want to draw an ordinary
   // bond between the atom postions as well as the wedge (because the wedge is shortened
   // for aesthetic reasons).
   //
   if (other_connections_to_second_atom.size() == 1) {
      const lig_build::pos_t  &third_atom_pos = other_connections_to_second_atom[0].first.atom_position;
      const lig_build::bond_t &third_bond     = other_connections_to_second_atom[0].second;
      if (third_bond.get_bond_type() == lig_build::bond_t::TRIPLE_BOND) {
	 lig_build::pos_t p1 = svg_molecule_t::mol_coords_to_svg_coords(pos_1, centre, scale);
	 lig_build::pos_t p2 = svg_molecule_t::mol_coords_to_svg_coords(pos_2, centre, scale);

	 // cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
	 // cairo_move_to(cr, p1.x, p1.y);
	 // cairo_line_to(cr, p2.x, p2.y);
	 // cairo_stroke(cr);

         s += make_bond_line_string(p1, p2);

      }
   }

   // std::cout << "in draw_sheared_or_darted_wedge_bond() v.size " << v.size() << std::endl;

   // // darts have 5 points, sheared bonds have 4 points. We don't want to stroke the
   // // darts (I think).
   // //
   // if (v.size() == 4) {
   //    // stroked shape
   //    //
   //    lig_build::pos_t p = svg_molecule_t::mol_coords_to_svg_coords(v[0], centre, scale);
   //    cairo_move_to(cr, p.x, p.y);
   //    for (unsigned int i=1; i<v.size(); i++) {
   //       lig_build::pos_t p_i = svg_molecule_t::mol_coords_to_svg_coords(v[i], centre, scale);
   //       cairo_line_to(cr, p_i.x, p_i.y);
   //    }
   //    cairo_close_path(cr);
   //    cairo_stroke(cr);
   // }

   // // filled shape
   // //
   // lig_build::pos_t p = svg_molecule_t::mol_coords_to_svg_coords(v[0], centre, scale);
   // cairo_move_to(cr, p.x, p.y);
   // for (unsigned int i=1; i<v.size(); i++) {
   //    lig_build::pos_t p_i = svg_molecule_t::mol_coords_to_svg_coords(v[i], centre, scale);
   //    cairo_line_to(cr, p_i.x, p_i.y);
   // }
   // cairo_close_path(cr);
   // cairo_fill(cr);
   // cairo_stroke(cr);

   double sf = 400.0; // scale_factor

   if (v.size() > 4) {
      s += "<polygon points=\"";
      for (unsigned int i=0; i<v.size(); i++) {
         lig_build::pos_t p_i = svg_molecule_t::mol_coords_to_svg_coords(v[i], centre, scale);         
         s += std::to_string(sf * p_i.x);
         s += ",";
         s += std::to_string(sf * p_i.y);
         s += " ";
      }
      s += "\" style=\"fill:#2020202;\" />\n";
   }

   return s;
}



void
svg_atom_t::set_colour() {

   colour = "grey";
   if (element == "C") colour = "black";
   if (element == "O") colour = "red";
   if (element == "N") colour = "blue";
   if (element == "S") colour = "#bbbb00";
   if (element == "F") colour = "green";
   if (element == "Cl") colour = "green";
   if (element == "Br") colour = "brown";
   if (element == "I") colour = "purple";
   if (element == "P") colour = "orange";
   if (element == "Fe") colour = "brown";
   if (element == "H") colour = "lightgrey";

}

std::string
svg_atom_t::make_text_item(const lig_build::atom_id_info_t &atom_id_info,
                           const lig_build::pos_t &centre, double scale,
                           double median_bond_length) const {

   double sf = 400.0; // scale factor

   std::string s;

   for (unsigned int i=0; i<atom_id_info.n_offsets(); i++) {

      // cairo_set_font_size(cr, 0.44 * scale * median_bond_length);
      lig_build::pos_t p = svg_molecule_t::mol_coords_to_svg_coords(atom_position, centre, scale);
      p += atom_id_info.offsets[i].tweak * scale * 0.030 * median_bond_length;

      // should these positions depend on the median_bond_length_?

      if (atom_id_info[i].text_pos_offset == lig_build::offset_text_t::UP)
	 p.y -= 0.36 * scale * median_bond_length;
      if (atom_id_info[i].text_pos_offset == lig_build::offset_text_t::DOWN)
	 p.y += 0.36 * scale * median_bond_length;

      if (atom_id_info.size_hint == -1) {
	 // cairo_set_font_size(cr, 0.44 * scale * 0.7 * median_bond_length);
      }

      if (atom_id_info.offsets[i].subscript) {
	 p.y += 0.2 * scale * median_bond_length;
	 // cairo_set_font_size(cr, 0.66 * scale * 0.533 * median_bond_length);
      }
      if (atom_id_info.offsets[i].superscript) {
	 // cairo_set_font_size(cr, 0.66 * scale * 0.533 * median_bond_length);
	 p.y -= 0.2 * scale * median_bond_length;
      }

      if (false)
	 std::cout << "Rendering tweak " << i << " :" << atom_id_info[i].text
		   << ": with tweak " << atom_id_info[i].tweak
		   << ": with size_hint " << atom_id_info.size_hint
		   << " with text_pos_offset " << atom_id_info[i].text_pos_offset
		   << " at pos " << p << std::endl;

      std::string txt = atom_id_info.offsets[i].text;
      if (txt == std::string("−")) {
	 // std::cout << "--- convert unicode! " << std::endl;
	 txt = "-";
      }
      if (true) {
	 if (txt.size() > 0) {

	    // cairo_text_extents_t te;
	    // // we need to "centre" the first letter of the text, ie, the N, not the NH.
	    // std::string t0(txt.substr(0,1));
	    // cairo_text_extents(cr, t0.c_str(), &te);
	    // cairo_move_to(cr, p.x, p.y);
	    // cairo_rel_move_to(cr,
	    //     	      - te.x_bearing - te.width / 2,
	    //     	      - te.y_bearing - te.height / 2);
	    // if (false)
	    //    std::cout << "show text \"" << txt << "\" at " << p << " for t0 "
	    //     	 << t0 << " with move-rel "
	    //     	 << te.x_bearing << " - " << te.width  << "/2 "
	    //     	 << te.y_bearing << " - " << te.height << "/2 "
	    //     	 << std::endl;
	    // cairo_show_text(cr, txt.c_str());
	    // cairo_stroke(cr);

            // 20230215-PE upated
            double x_fudge = -sf * 0.50 / 50.0;
            double y_fudge =  sf * 0.65  / 50.0;

            std::string default_font_size = "\"0.8em\"";
            std::string font_size = default_font_size;
            if (atom_id_info.offsets[i].superscript) font_size = "\"0.6em\"";
            if (atom_id_info.offsets[i].subscript)   font_size = "\"0.6em\"";
            if (atom_id_info.offsets[i].subscript) x_fudge += sf * 0.002;
            if (txt == "-") font_size = "\"1.0em\"";
            if (txt == "-") x_fudge += sf * 0.005;

            std::string atom_string;
            atom_string += "   <text x=\"";
            atom_string += std::to_string(sf * p.x + x_fudge);
            atom_string += "\" y=\"";
            atom_string += std::to_string(sf * p.y + y_fudge);
            atom_string += "\" font-family=\"Helvetica, sans-serif\" font-size=" + font_size + " fill=\"";
            atom_string += colour;
            atom_string += "\">";
            atom_string += txt;
            atom_string += "</text>\n";

            // std::cout << "atom_string: " << atom_string;

            s += atom_string;
	 } else {
	    std::cout << "oops empty text!" << std::endl;
	 }
      }
   }
   return s;
}


std::string
svg_molecule_t::render_to_svg_string() {

   auto make_bond_comment = [] (unsigned int bond_idx, const svg_bond_t &bond) {
      std::string s("<!-- ");
      s += "Bond ";
      s += std::to_string(bond_idx);
      s += " atom ";
      s += std::to_string(bond.get_atom_1_index());
      s += " to atom ";
      s += std::to_string(bond.get_atom_2_index());
      s += " ";
      if (bond.get_bond_type() == lig_build::bond_t::SINGLE_BOND)          s += std::string("single");
      if (bond.get_bond_type() == lig_build::bond_t::DOUBLE_BOND)          s += std::string("double");
      if (bond.get_bond_type() == lig_build::bond_t::TRIPLE_BOND)          s += std::string("triple");
      if (bond.get_bond_type() == lig_build::bond_t::SINGLE_OR_DOUBLE)     s += std::string("single-or-double");
      if (bond.get_bond_type() == lig_build::bond_t::SINGLE_OR_AROMATIC)   s += std::string("single-or-aromatic");
      if (bond.get_bond_type() == lig_build::bond_t::DELOC_ONE_AND_A_HALF) s += std::string("deloc-one-and-a-half");
      if (bond.get_bond_type() == lig_build::bond_t::IN_BOND)              s += std::string("in-bond");
      if (bond.get_bond_type() == lig_build::bond_t::OUT_BOND)             s += std::string("out-bond");
      if (bond.get_bond_type() == lig_build::bond_t::BOND_ANY)             s += std::string("bond-any");
      s += std::string(" -->\n");
      return s;
   };

   double sf = 400.0; // scale factor
   std::string s;
   std::string svg_header_1 = "<svg xmlns=\"http://www.w3.org/2000/svg\"\n    xmlns:xlink=\"http://www.w3.org/1999/xlink\" ";
   std::string svg_header_2 = ">\n";
   std::string svg_footer = "</svg>\n";

   // determine the viewBox

   std::string viewBox_string;
   if (! atoms.empty()) {
      lig_build::pos_t centre = get_ligand_centre();
      double scale = get_scale();
      float min_x =  100000.0;
      float min_y =  100000.0;
      float max_x = -100000.0;
      float max_y = -100000.0;
      for (unsigned int iat=0; iat<atoms.size(); iat++) {
         const auto &atom_pos = atoms[iat].atom_position;
         lig_build::pos_t pos = mol_coords_to_svg_coords(atom_pos, centre, scale) * sf;
         if (pos.x < min_x) min_x = pos.x;
         if (pos.y < min_y) min_y = pos.y;
         if (pos.x > max_x) max_x = pos.x;
         if (pos.y > max_y) max_y = pos.y;
      }
      // now adjust so that the labels can fit
      min_x -= 0.0;
      min_y -= 0.0;

      float width  = max_x - min_x;
      float height = max_y - min_y;
      viewBox_string = "viewBox=" + std::string("\"") +
         std::to_string(min_x) + std::string(" ") +
         std::to_string(min_y) + std::string(" ") +
         std::to_string(max_x) + std::string(" ") +
         std::to_string(max_y) + std::string("\"");
   }

   s += svg_header_1;
   s += viewBox_string;
   s += svg_header_2;

   std::cout << "viewBox: " << viewBox_string << std::endl;

   // just testing that I can see something. No longer needed because I can
   // s += "   <rect x=\"10\" y=\"10\" width=\"10\" height=\"10\" style=\"stroke:#ff0000; fill: #ff6666;\" />\n";

   double scale = get_scale();
   lig_build::pos_t centre = get_ligand_centre();
   
   // cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
   // cairo_set_line_width(cr, 0.07 * scale * median_bond_length_);

   for (unsigned int ib=0; ib<bonds.size(); ib++) {
      int idx_1 = bonds[ib].get_atom_1_index();
      int idx_2 = bonds[ib].get_atom_2_index();
      if ((idx_1 != UNASSIGNED_INDEX) && (idx_2 != UNASSIGNED_INDEX)) {
	 lig_build::bond_t::bond_type_t bt = bonds[ib].get_bond_type();

	 std::pair<bool, bool> shorten = shorten_flags(ib);

	 bool at_1_in_ring_flag = in_ring_p(idx_1);
	 bool at_2_in_ring_flag = in_ring_p(idx_2);

	 std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > other_connections_to_first_atom =
	    make_other_connections_to_first_atom_info(ib);
	 std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > other_connections_to_second_atom =
	    make_other_connections_to_second_atom_info(ib);

         std::string bond_string = bonds[ib].draw_bond(atoms[idx_1], atoms[idx_2],
                                                       at_1_in_ring_flag,
                                                       at_2_in_ring_flag,
                                                       bt,
                                                       shorten.first, shorten.second,
                                                       other_connections_to_first_atom,
                                                       other_connections_to_second_atom,
                                                       centre, scale);
         s += make_bond_comment(ib, bonds[ib]);
         s += bond_string;
      }
   }

   s += std::string("<!-- Atom Labels -->\n");
   for (unsigned int iat=0; iat<atoms.size(); iat++) {
      std::string ele = atoms[iat].element;
      std::vector<unsigned int> local_bonds = bonds_having_atom_with_atom_index(iat);
      bool gl_flag = false;
      if (ele != "C") {
	 lig_build::atom_id_info_t atom_id_info = make_atom_id_by_using_bonds(iat, ele, local_bonds, gl_flag);
	 atoms[iat].set_atom_id(atom_id_info.atom_id); // quick hack
	 if (false)
	    std::cout << "in render(): atom_index " << iat << " with charge "
		      << atoms[iat].charge << " made atom_id_info "
		      << atom_id_info << std::endl;
	 s += atoms[iat].make_text_item(atom_id_info, centre, scale, median_bond_length_);
      } else {
#if 0
	 // a super-atom carbon
	 if (local_bonds.size() == 1) {
	    atoms[iat].set_colour(cr);
	    lig_build::atom_id_info_t atom_id_info =
	       make_atom_id_by_using_bonds(iat, ele, local_bonds, gl_flag);
	    atoms[iat].make_text_item(cr, atom_id_info, centre, scale, median_bond_length_);
	 }
#endif
      }
   }
   

   s += svg_footer;
   return s;
}


double
svg_molecule_t::get_scale() const {

   std::pair<lig_build::pos_t, lig_build::pos_t> ext = ligand_extents();
   double delta_x = ext.second.x - ext.first.x;
   double delta_y = ext.second.y - ext.first.y;
   double delta = (delta_x > delta_y) ? delta_x : delta_y;
   // The scale transforms from molecule coordinates to the 0->1 square
   // (centering is also applied of course)
   // What makes this tricky is that atom labels (e.g. "NH2") can go beyond
   // the limits of the atoms. So we need to make space for them by making
   // the molecules a bit smaller that exact fit to the box. How much smaller?
   // Not clear. A bit less than 0.75? BHH falls off the edge at 0.75 (with
   // scale 0.09).
   // So let's make the scale limit a bit smaller (was 0.09, now 0.089).
   // Also, we don't want massive "zoomed in" representation of tiny molecules
   // (e.g. ALM).
   double scale_lim = 0.089; // heuristic
   double scale = scale_lim;
   // we could be more clever here by looking at the atoms on the edge of the molecle
   // if they are non-carbon, reduce the scale a bit. 0.75 is a safe value that
   // allows non-C edge atom representation.
   if (delta > 1)
      scale = 0.74/delta; // was 0.75 (0.75 cuts atom labels on occassion e.g. 1386433)
   if (scale > scale_lim)
      scale = scale_lim;

   std::cout << "get_scale() returns " << scale << std::endl;
   return scale;
}


#endif // MAKE_ENHANCED_LIGAND_TOOLS
