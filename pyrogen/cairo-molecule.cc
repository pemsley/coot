
#include "cairo-molecule.hh"

// this needs to be here (also).  It is in wmolecule.cc and hence the library also.
// But if this is not here I get unresovled symbol for this destructor when compiling
// this exectuable on the mac (clang).
// template<class cairo_atom_t, class cairo_bond_t> lig_build::molecule_t<cairo_atom_t, cairo_bond_t>::~molecule_t() {}

// coot::cairo_molecule_t::~cairo_molecule_t() { }


#ifdef MAKE_ENHANCED_LIGAND_TOOLS
void
coot::cairo_molecule_t::import_rdkit_mol(RDKit::ROMol *rdkm, int iconf) {

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

      if (n_mol_atoms > 0) {
	 double centre_x = sum_x/double(n_mol_atoms);
	 double centre_y = sum_y/double(n_mol_atoms);
	 // centre_correction = lig_build::pos_t(centre_x, centre_y);
	 // mol_in_min_y = min_y;
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
	 cairo_atom_t mol_at(pos, element, charge);

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
	    const cairo_atom_t &cat1 = atoms[idx_1];
	    const cairo_atom_t &cat2 = atoms[idx_2];
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

	    // needed?
	    std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > empty;

	    // do we need to pass idx_1 and idx_2, neighbour info?
// 	    cairo_bond_t bond(idx_1, idx_2, cat1, cat2, shorten_first, shorten_second,
// 			      bt, empty, empty);
	    cairo_bond_t bond(idx_1, idx_2, bt);
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
   // std::cout << "median_bond_length: " << median_bond_length_ << std::endl;
}
#endif // MAKE_ENHANCED_LIGAND_TOOLS

void
coot::cairo_atom_t::set_colour(cairo_t *cr) const {

   if (element == "C")  { cairo_set_source_rgb(cr, 0.1, 0.1, 0.1); } else {
   if (element == "O")  { cairo_set_source_rgb(cr, 0.8, 0.0, 0.0); } else {
   if (element == "N")  { cairo_set_source_rgb(cr, 0.2, 0.2, 0.8); } else {
   if (element == "S")  { cairo_set_source_rgb(cr, 0.6, 0.4, 0.2); } else {
   if (element == "F")  { cairo_set_source_rgb(cr, 0.0, 0.5, 0.0); } else {
   if (element == "Cl") { cairo_set_source_rgb(cr, 0.0, 0.5, 0.0); } else {
   if (element == "Br") { cairo_set_source_rgb(cr, 0.5, 0.2, 0.0); } else {
   if (element == "I")  { cairo_set_source_rgb(cr, 0.3, 0.0, 0.3); } else {
   if (element == "P")  { cairo_set_source_rgb(cr, 0.8, 0.5, 0.0); } else { // orange
   if (element == "Fe") { cairo_set_source_rgb(cr, 0.6, 0.3, 0.0); } else { // dark orange
   if (element == "H")  { cairo_set_source_rgb(cr, 0.5, 0.5, 0.5); } else { // not white

      cairo_set_source_rgb(cr, 0.7, 0.3, 0.9);
   }}}}}}}}}}}
}

void
coot::cairo_atom_t::make_text_item(cairo_t *cr,
				   const lig_build::atom_id_info_t &atom_id_info,
				   const lig_build::pos_t &centre, double scale,
				   double median_bond_length) const {

   for (unsigned int i=0; i<atom_id_info.n_offsets(); i++) {

      cairo_set_font_size(cr, 0.44 * scale * median_bond_length);
      lig_build::pos_t p = cairo_molecule_t::mol_coords_to_cairo_coords(atom_position, centre, scale);
      p += atom_id_info.offsets[i].tweak * scale * 0.030 * median_bond_length;

      // should these positions depend on the median_bond_length_?

      if (atom_id_info[i].text_pos_offset == lig_build::offset_text_t::UP)
	 p.y -= 0.36 * scale * median_bond_length;
      if (atom_id_info[i].text_pos_offset == lig_build::offset_text_t::DOWN)
	 p.y += 0.36 * scale * median_bond_length;

      if (atom_id_info.size_hint == -1)
	 cairo_set_font_size(cr, 0.44 * scale * 0.7 * median_bond_length);

      if (atom_id_info.offsets[i].subscript) {
	 p.y += 0.2 * scale * median_bond_length;
	 cairo_set_font_size(cr, 0.66 * scale * 0.533 * median_bond_length);
      }
      if (atom_id_info.offsets[i].superscript) {
	 cairo_set_font_size(cr, 0.66 * scale * 0.533 * median_bond_length);
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
	    cairo_text_extents_t te;
	    // we need to "centre" the first letter of the text, ie, the N, not the NH.
	    std::string t0(txt.substr(0,1));
	    cairo_text_extents(cr, t0.c_str(), &te);
	    cairo_move_to(cr, p.x, p.y);
	    cairo_rel_move_to(cr,
			      - te.x_bearing - te.width / 2,
			      - te.y_bearing - te.height / 2);
	    if (false)
	       std::cout << "show text \"" << txt << "\" at " << p << " for t0 "
			 << t0 << " with move-rel "
			 << te.x_bearing << " - " << te.width  << "/2 "
			 << te.y_bearing << " - " << te.height << "/2 "
			 << std::endl;
	    cairo_show_text(cr, txt.c_str());
	    cairo_stroke(cr);
	 } else {
	    std::cout << "oops empty text!" << std::endl;
	 }
      }
   }
}

// shorten_first and shorten_second are set depending on the element of the atoms and number of bonds.
void
coot::cairo_bond_t::draw_bond(cairo_t *cr,
			      const coot::cairo_atom_t &at_1,
			      const coot::cairo_atom_t &at_2,
			      bool at_1_in_ring_flag, bool at_2_in_ring_flag,
			      lig_build::bond_t::bond_type_t bt,
			      bool shorten_first, bool shorten_second,
			      const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_first_atom,
			      const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
			      const lig_build::pos_t &centre, double scale) {

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

   // for case
   lig_build::pos_t p1;
   lig_build::pos_t p2;

   // std::cout << "draw_bond for bt " << bt << " between " << at_1 << " and " << at_2 << std::endl;

   switch (bt) {
   case lig_build::bond_t::SINGLE_BOND:
   case lig_build::bond_t::SINGLE_OR_DOUBLE:
   case lig_build::bond_t::SINGLE_OR_AROMATIC:
   case lig_build::bond_t::DELOC_ONE_AND_A_HALF:
   case lig_build::bond_t::BOND_ANY:
      {
	 // push down the decision making for the shortening
	 bool at_1_is_singleton = false;
	 bool at_2_is_singleton = false;
	 if (other_connections_to_first_atom.size()  == 0) at_1_is_singleton = true;
	 if (other_connections_to_second_atom.size() == 0) at_2_is_singleton = true;
	 std::pair<lig_build::pos_t, lig_build::pos_t> bp =
	    coords_for_single_bond(at_1, at_2, at_1_is_singleton, at_2_is_singleton);
	 pos_1 = bp.first;
	 pos_2 = bp.second;
	 p1 = cairo_molecule_t::mol_coords_to_cairo_coords(pos_1, centre, scale);
	 p2 = cairo_molecule_t::mol_coords_to_cairo_coords(pos_2, centre, scale);
	 cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
	 cairo_move_to(cr, p1.x, p1.y);
	 cairo_line_to(cr, p2.x, p2.y);
	 cairo_stroke(cr);
      }
      break;

   case lig_build::bond_t::DOUBLE_BOND:

   case lig_build::bond_t::DOUBLE_OR_AROMATIC:
      {
	 if (have_centre_pos()) {
	    draw_double_in_ring_bond(cr, pos_1_in, pos_2_in, shorten_first, shorten_second,
				     centre, scale);
	 } else {
	    draw_double_bond(cr, at_1, at_2,
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
	 bool dashed_inner = true;
	 draw_double_in_ring_bond(cr, pos_1_in, pos_2_in, shorten_first, shorten_second,
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

	 lig_build::pos_t sc_p1 = cairo_molecule_t::mol_coords_to_cairo_coords(p1, centre, scale);
	 lig_build::pos_t sc_p2 = cairo_molecule_t::mol_coords_to_cairo_coords(p2, centre, scale);
	 lig_build::pos_t sc_p3 = cairo_molecule_t::mol_coords_to_cairo_coords(p3, centre, scale);
	 lig_build::pos_t sc_p4 = cairo_molecule_t::mol_coords_to_cairo_coords(p4, centre, scale);
	 lig_build::pos_t sc_p5 = cairo_molecule_t::mol_coords_to_cairo_coords(p5, centre, scale);
	 lig_build::pos_t sc_p6 = cairo_molecule_t::mol_coords_to_cairo_coords(p6, centre, scale);

	 cairo_move_to(cr, sc_p1.x, sc_p1.y);
	 cairo_line_to(cr, sc_p2.x, sc_p2.y);
	 cairo_stroke(cr);
	 cairo_move_to(cr, sc_p3.x, sc_p3.y);
	 cairo_line_to(cr, sc_p4.x, sc_p4.y);
	 cairo_stroke(cr);
	 cairo_move_to(cr, sc_p5.x, sc_p5.y);
	 cairo_line_to(cr, sc_p6.x, sc_p6.y);
	 cairo_stroke(cr);
      }
      break;

   case IN_BOND:
      {
	 // set of lines
	 std::vector<std::pair<lig_build::pos_t, lig_build::pos_t> > vp =
	    lig_build::pos_t::make_wedge_in_bond(pos_1, pos_2);
	 if (vp.size()) {
	    cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
	    for (unsigned int i=0; i<vp.size(); i++) { 

	       lig_build::pos_t p1 = cairo_molecule_t::mol_coords_to_cairo_coords(vp[i].first,  centre, scale);
	       lig_build::pos_t p2 = cairo_molecule_t::mol_coords_to_cairo_coords(vp[i].second, centre, scale);

	       cairo_move_to(cr, p1.x, p1.y);
	       cairo_line_to(cr, p2.x, p2.y);
	    }
	    cairo_stroke(cr);
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
	       draw_sheared_or_darted_wedge_bond(cr, pos_1, pos_2, other_connections_to_second_atom,
						 centre, scale);
	       done_darted = true;
	    }
	 }

	 if (! done_darted) {
	    // filled shape (normal wedge)
	    std::vector<lig_build::pos_t> v =
	       lig_build::pos_t::make_wedge_out_bond(pos_1, pos_2);

	    lig_build::pos_t p = cairo_molecule_t::mol_coords_to_cairo_coords(v[0], centre, scale);
	    cairo_move_to(cr, p.x, p.y);
	    cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
	    for (unsigned int i=1; i<v.size(); i++) {
	       lig_build::pos_t p_i = cairo_molecule_t::mol_coords_to_cairo_coords(v[i], centre, scale);
	       cairo_line_to(cr, p_i.x, p_i.y);
	    }
	    cairo_close_path(cr);
	    cairo_fill(cr);
	    cairo_stroke(cr);
	 }
      }
      break;
   case BOND_UNDEFINED:
      break;
   }
}


void
coot::cairo_bond_t::draw_sheared_or_darted_wedge_bond(cairo_t *cr,
						      const lig_build::pos_t &pos_1,
						      const lig_build::pos_t &pos_2,
						      const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
						      const lig_build::pos_t &centre,
						      double scale) const {
   std::vector<lig_build::pos_t> v =
      coords_for_sheared_or_darted_wedge_bond(pos_1, pos_2,
					      other_connections_to_second_atom);

   // if there is a triple bond connected to the second atom, we want to draw an ordinary
   // bond between the atom postions as well as the wedge (because the wedge is shortened
   // for aesthetic reasons).
   //
   if (other_connections_to_second_atom.size() == 1) {
      const lig_build::pos_t  &third_atom_pos = other_connections_to_second_atom[0].first.atom_position;
      const lig_build::bond_t &third_bond     = other_connections_to_second_atom[0].second;
      if (third_bond.get_bond_type() == lig_build::bond_t::TRIPLE_BOND) {
	 lig_build::pos_t p1 = cairo_molecule_t::mol_coords_to_cairo_coords(pos_1, centre, scale);
	 lig_build::pos_t p2 = cairo_molecule_t::mol_coords_to_cairo_coords(pos_2, centre, scale);
	 cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
	 cairo_move_to(cr, p1.x, p1.y);
	 cairo_line_to(cr, p2.x, p2.y);
	 cairo_stroke(cr);
      }
   }

   // darts have 5 points, sheared bonds have 4 points. We don't want to stroke the
   // darts (I think).
   //
   if (v.size() == 4) {
      // stroked shape
      //
      lig_build::pos_t p = cairo_molecule_t::mol_coords_to_cairo_coords(v[0], centre, scale);
      cairo_move_to(cr, p.x, p.y);
      for (unsigned int i=1; i<v.size(); i++) {
	 lig_build::pos_t p_i = cairo_molecule_t::mol_coords_to_cairo_coords(v[i], centre, scale);
	 cairo_line_to(cr, p_i.x, p_i.y);
      }
      cairo_close_path(cr);
      cairo_stroke(cr);
   }

   // filled shape
   //
   lig_build::pos_t p = cairo_molecule_t::mol_coords_to_cairo_coords(v[0], centre, scale);
   cairo_move_to(cr, p.x, p.y);
   for (unsigned int i=1; i<v.size(); i++) {
      lig_build::pos_t p_i = cairo_molecule_t::mol_coords_to_cairo_coords(v[i], centre, scale);
      cairo_line_to(cr, p_i.x, p_i.y);
   }
   cairo_close_path(cr);
   cairo_fill(cr);
   cairo_stroke(cr);
   
}


// static
lig_build::pos_t
coot::cairo_molecule_t::cairo_molecule_t::mol_coords_to_cairo_coords(const lig_build::pos_t &pos_1,
						   const lig_build::pos_t &centre,
						   double scale) {

   lig_build::pos_t p1 = (pos_1 - centre) * scale;
   p1.y = -p1.y; // canvas is upside down c.f. normal/real-world/molecule coordinates
   p1 += lig_build::pos_t(0.5,0.5);
   return p1;

}


// This is for a double bond in a ring
//
// pos_1 and pos_2 are the coordinates of the atoms, not the pre-shorted
// coordinates
//
void
coot::cairo_bond_t::draw_double_in_ring_bond(cairo_t *cr,
					     const lig_build::pos_t &pos_1_in,
					     const lig_build::pos_t &pos_2_in,
					     bool shorten_first,
					     bool shorten_second,
					     const lig_build::pos_t &centre,
					     double scale, bool dashed_inner) {
   lig_build::pos_t pos_1 = pos_1_in;
   lig_build::pos_t pos_2 = pos_2_in;

   double shorten_fraction = 0.74;
   if (shorten_first)
      pos_1 = lig_build::pos_t::fraction_point(pos_2_in, pos_1_in, shorten_fraction);
   if (shorten_second)
      pos_2 = lig_build::pos_t::fraction_point(pos_1_in, pos_2_in, shorten_fraction);

   // outside "normal" bond
   lig_build::pos_t p1 = cairo_molecule_t::mol_coords_to_cairo_coords(pos_1, centre, scale);
   lig_build::pos_t p2 = cairo_molecule_t::mol_coords_to_cairo_coords(pos_2, centre, scale);

   // for inner bond
   std::pair<lig_build::pos_t, lig_build::pos_t> p = 
      make_double_aromatic_short_stick(pos_1_in, pos_2_in, shorten_first, shorten_second);

   cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
   cairo_move_to(cr, p1.x, p1.y);
   cairo_line_to(cr, p2.x, p2.y);
   cairo_stroke(cr);

   if (dashed_inner) {
      double dashlength = 0.015; // 0.01 is also fine
      cairo_set_dash(cr, &dashlength, 1, 0);
   }
   p1 = cairo_molecule_t::mol_coords_to_cairo_coords(p.first,  centre, scale);
   p2 = cairo_molecule_t::mol_coords_to_cairo_coords(p.second, centre, scale);
   cairo_move_to(cr, p1.x, p1.y);
   cairo_line_to(cr, p2.x, p2.y);
   cairo_stroke(cr);
   if (dashed_inner)
      cairo_set_dash(cr, NULL, 0, 0); // restore
}

// not in-ring
//
// shorten_first and shorten_second are set depending on the element of the atoms and number of bonds.
void
coot::cairo_bond_t::draw_double_bond(cairo_t *cr,
				     const lig_build::atom_t &at_1,
				     const lig_build::atom_t &at_2,
				     bool shorten_first, bool shorten_second,
				     const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_first_atom,
				     const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
				     const lig_build::pos_t &centre, double scale) {

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

      // boring case :-)

      std::pair<std::pair<lig_build::pos_t, lig_build::pos_t>, std::pair<lig_build::pos_t, lig_build::pos_t> > p =
	 make_double_bond(at_1.atom_position, at_2.atom_position, shorten_first, shorten_second);

      lig_build::pos_t p1 = cairo_molecule_t::mol_coords_to_cairo_coords(p.first.first,  centre, scale);
      lig_build::pos_t p2 = cairo_molecule_t::mol_coords_to_cairo_coords(p.first.second, centre, scale);

      cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
      cairo_move_to(cr, p1.x, p1.y);
      cairo_line_to(cr, p2.x, p2.y);
      cairo_stroke(cr);

      p1 = cairo_molecule_t::mol_coords_to_cairo_coords(p.second.first,  centre, scale);
      p2 = cairo_molecule_t::mol_coords_to_cairo_coords(p.second.second, centre, scale);

      cairo_move_to(cr, p1.x, p1.y);
      cairo_line_to(cr, p2.x, p2.y);
      cairo_stroke(cr);

   } else {
      // elegant case, offset, shorten a bond
      std::pair<std::pair<lig_build::pos_t, lig_build::pos_t>, std::pair<lig_build::pos_t, lig_build::pos_t> > bonds =
	 make_double_bond(at_1.atom_position, at_2.atom_position, shorten_first, shorten_second,
			  other_connections_to_first_atom, other_connections_to_second_atom);

      lig_build::pos_t p1 = cairo_molecule_t::mol_coords_to_cairo_coords(bonds.first.first,  centre, scale);
      lig_build::pos_t p2 = cairo_molecule_t::mol_coords_to_cairo_coords(bonds.first.second, centre, scale);
      cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
      cairo_move_to(cr, p1.x, p1.y);
      cairo_line_to(cr, p2.x, p2.y);
      p1 = cairo_molecule_t::mol_coords_to_cairo_coords(bonds.second.first,  centre, scale);
      p2 = cairo_molecule_t::mol_coords_to_cairo_coords(bonds.second.second, centre, scale);
      cairo_move_to(cr, p1.x, p1.y);
      cairo_line_to(cr, p2.x, p2.y);
      cairo_stroke(cr);
   }
}


double
coot::cairo_molecule_t::get_scale() const {
   lig_build::pos_t centre = get_ligand_centre();
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

   return scale;
}


// not const because it changes atom ids.
void
coot::cairo_molecule_t::render(cairo_t *cr) {

   bool debug = true;
   double scale = get_scale();
   lig_build::pos_t centre = get_ligand_centre();

   // ------------ the font size and the line width depend on the size of the
   //              extents (and here, scale), the smaller the scale the smaller
   //              the font
   //
   // cairo_set_source_rgb (cr, 0, 0, 0);
   // cairo_rectangle (cr, 0.25, 0.25, 0.5, 0.5);
   // cairo_set_source_rgb(cr, 0, 0.5, 0);
   // cairo_rectangle (cr, 0.01, 0.01, 0.98, 0.98);
   // cairo_stroke(cr);
   cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);

   // std::cout << "in render: scale: " << scale << std::endl;

   cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
   cairo_set_line_width(cr, 0.07 * scale * median_bond_length_);
   cairo_font_extents_t fe;
   cairo_text_extents_t te;
   cairo_select_font_face (cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
   cairo_font_extents (cr, &fe);

   for (unsigned int ib=0; ib<bonds.size(); ib++) {
      int idx_1 = bonds[ib].get_atom_1_index();
      int idx_2 = bonds[ib].get_atom_2_index();
      if ((idx_1 != UNASSIGNED_INDEX) && (idx_2 != UNASSIGNED_INDEX)) {
	 lig_build::bond_t::bond_type_t bt = bonds[ib].get_bond_type();

	 std::pair<bool, bool> shorten = shorten_flags(ib);

	 lig_build::pos_t pos_1 =  atoms[idx_1].atom_position;
	 lig_build::pos_t pos_2 =  atoms[idx_2].atom_position;
	 // c.f. canvas_item_for_bond

	 bool at_1_in_ring_flag = in_ring_p(idx_1);
	 bool at_2_in_ring_flag = in_ring_p(idx_2);

	 std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > other_connections_to_first_atom =
	    make_other_connections_to_first_atom_info(ib);
	 std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > other_connections_to_second_atom =
	    make_other_connections_to_second_atom_info(ib);

	 bonds[ib].draw_bond(cr, atoms[idx_1], atoms[idx_2],
			     at_1_in_ring_flag,
			     at_2_in_ring_flag,
			     bt,
			     shorten.first, shorten.second,
			     other_connections_to_first_atom,
			     other_connections_to_second_atom,
			     centre, scale);
      }
   }

   for (unsigned int iat=0; iat<atoms.size(); iat++) {
      std::string ele = atoms[iat].element;
      std::vector<unsigned int> local_bonds = bonds_having_atom_with_atom_index(iat);
      bool gl_flag = false;
      if (ele != "C") {
	 lig_build::atom_id_info_t atom_id_info =
	    make_atom_id_by_using_bonds(iat, ele, local_bonds, gl_flag);
	 atoms[iat].set_atom_id(atom_id_info.atom_id); // quick hack
	 if (false)
	    std::cout << "in render(): atom_index " << iat << " with charge "
		      << atoms[iat].charge << " made atom_id_info "
		      << atom_id_info << std::endl;
	 atoms[iat].set_colour(cr);
	 atoms[iat].make_text_item(cr, atom_id_info, centre, scale, median_bond_length_);
      } else {
	 // a super-atom carbon
	 if (local_bonds.size() == 1) {
	    atoms[iat].set_colour(cr);
	    lig_build::atom_id_info_t atom_id_info =
	       make_atom_id_by_using_bonds(iat, ele, local_bonds, gl_flag);
	    atoms[iat].make_text_item(cr, atom_id_info, centre, scale, median_bond_length_);
	 }
      }
   }
   if (false)
      debug_box(cr);
}

void
coot::cairo_molecule_t::debug_box(cairo_t *cr) {

   cairo_set_line_width(cr, 0.02);

   cairo_set_source_rgb(cr, 1, 0, 0); // red
   cairo_move_to(cr, 0.01, 0.01);
   cairo_line_to(cr, 0.01, 0.99);
   cairo_stroke(cr);

   cairo_set_source_rgb(cr, 0.6, 0.4, 0); // darker red
   cairo_move_to(cr, 0.01, 0.99);
   cairo_line_to(cr, 0.99, 0.99);
   cairo_stroke(cr);

   cairo_set_source_rgb(cr, 0, 1, 1); // cyan
   cairo_move_to(cr, 0.99, 0.99);
   cairo_line_to(cr, 0.99, 0.01);
   cairo_stroke(cr);

   cairo_set_source_rgb(cr, 0.3, 0.4, 0.5); // dark cyan
   cairo_move_to(cr, 0.99, 0.01);
   cairo_line_to(cr, 0.01, 0.01);
   cairo_stroke(cr);

   cairo_rectangle (cr, 0.25, 0.25, 0.5, 0.5);
   // cairo_set_source_rgb(cr, 0, 0.5, 0);
   // ncairo_rectangle (cr, 0.02, 0.02, 0.98, 0.98);
   cairo_stroke(cr);
}

void
coot::cairo_molecule_t::render_to_file(const std::string &png_file_name, unsigned int npx,
				       const std::pair<bool, colour_holder> &bg_col) {

   cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, npx, npx);
   cairo_t *cr = cairo_create(surface);

   // Drawing space is 0->1 with 0,0 top left
   // if cairo_image_surface_create called with 240,240, then cairo_scale called with
   // 120,120 puts the image, half-size in top left quadrant
   cairo_scale(cr, npx, npx);
   // background colour:
   if (bg_col.first) {
      cairo_set_source_rgb(cr, bg_col.second.red, bg_col.second.green, bg_col.second.blue);
      cairo_paint (cr);
   }
   render(cr);

   /* Write output and clean up */
   cairo_surface_write_to_png(surface, png_file_name.c_str());
   cairo_destroy(cr);
   cairo_surface_destroy(surface);

}

// static
cairo_status_t
coot::cairo_molecule_t::png_stream_writer(void *closure_in,
					  const unsigned char *data,
					  unsigned int length) {

   std::string *s_ptr = static_cast<std::string *>(closure_in);
   *s_ptr += std::string(reinterpret_cast<const char *>(data), length); // it's safe!
   return CAIRO_STATUS_SUCCESS;
}

std::string
coot::cairo_molecule_t::render_to_png_string(const std::vector<unsigned int> &atom_highlight_list,
					     const std::vector<unsigned int> &bond_highlight_list,
					     bool use_highlight_bond_indices_flag,
					     unsigned int npx) {

   bool set_background = false; // this (or a non-None background colour) should be passed
   std::string s;
   s.reserve(12000);

   cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, npx, npx);
   cairo_t *cr = cairo_create(surface);
   cairo_scale(cr, npx, npx);
   double scale = get_scale();
   lig_build::pos_t centre = get_ligand_centre();

   // background colour:
   if (set_background) {
      cairo_set_source_rgb (cr, 1, 1, 0.9);
      cairo_paint (cr);
   }

   draw_atom_highlights(cr, centre, scale, atom_highlight_list,
			bond_highlight_list, use_highlight_bond_indices_flag);
   render(cr);
   cairo_surface_write_to_png_stream(surface, png_stream_writer, reinterpret_cast<void *> (&s));
   cairo_destroy(cr);
   cairo_surface_destroy(surface);
   return s;
}

#if CAIRO_HAS_SVG_SURFACE
#include <cairo/cairo-svg.h>
#endif // CAIRO_HAS_SVG_SURFACE

std::string
coot::cairo_molecule_t::render_to_svg_string(const std::vector<unsigned int> &atom_highlight_list,
					     const std::vector<unsigned int> &bond_highlight_list,
					     bool use_highlight_bond_indices_flag,
					     unsigned int npx) {

   // consider consolidating this and the png string version

   std::string s;

#ifdef CAIRO_HAS_SVG_SURFACE
   s.reserve(12000);

   cairo_surface_t *surface = cairo_svg_surface_create_for_stream(png_stream_writer, reinterpret_cast<void *>(&s),
								  npx, npx);
   cairo_t *cr = cairo_create(surface);
   cairo_scale(cr, npx, npx);
   double scale = get_scale();
   lig_build::pos_t centre = get_ligand_centre();

   draw_atom_highlights(cr, centre, scale, atom_highlight_list,
			bond_highlight_list, use_highlight_bond_indices_flag);
   render(cr);

   cairo_destroy(cr);
   cairo_surface_destroy(surface);
#endif // CAIRO_HAS_SVG_SURFACE
   return s;
}


void
coot::cairo_molecule_t::set_highlight_colour(cairo_t *cr, unsigned int idx) {

   // cairo_set_source_rgba(cr, 0.8, 0.2, 0.2, 0.5); // not sure that it needs to be transparent
                                                     // if it's drawn first

   cairo_set_source_rgb(cr, 0.9, 0.7, 0.7);

}

std::vector<unsigned int>
coot::cairo_molecule_t::find_bonds_for_atoms(const std::vector<unsigned int> &highlight_atom_indices) const {

   std::vector<unsigned int> v;
   for (std::size_t ib=0; ib<bonds.size(); ib++) {
      unsigned int idx_1 = bonds[ib].get_atom_1_index();
      unsigned int idx_2 = bonds[ib].get_atom_2_index();
      if (std::find(highlight_atom_indices.begin(),
		    highlight_atom_indices.end(), idx_1) != highlight_atom_indices.end()) {
	 if (std::find(highlight_atom_indices.begin(),
		       highlight_atom_indices.end(), idx_2) != highlight_atom_indices.end()) {
	    v.push_back(ib);
	 }
      }
   }
   return v;
}

void
coot::cairo_molecule_t::draw_atom_highlights(cairo_t *cr,
					     const lig_build::pos_t &centre,
					     double scale,
					     const std::vector<unsigned int> &highlight_atom_indices,
					     const std::vector<unsigned int> &highlight_bond_indices,
					     bool use_highlight_bond_indices_flag) {

   // if use_highlight_bond_indices_flag then use the bonds in the vector
   // possibly none
   // else
   // work out the bonds yourself.

   // what is the line width now?
   cairo_set_line_width(cr, 0.07 * scale * median_bond_length_);
   unsigned int n_atoms = get_number_of_atoms_including_hydrogens();
   for (std::size_t i=0; i<highlight_atom_indices.size(); i++) {
      unsigned int idx = highlight_atom_indices[i];
      if (idx < n_atoms) {
	 set_highlight_colour(cr, i);

	 lig_build::pos_t pos_at = atoms[idx].atom_position;
	 lig_build::pos_t p = cairo_molecule_t::mol_coords_to_cairo_coords(pos_at, centre, scale);

	 cairo_new_sub_path(cr);
	 cairo_arc(cr, p.x, p.y, 0.03, 0, 2 * M_PI);
	 cairo_close_path(cr);
	 cairo_fill(cr);
	 cairo_stroke(cr);
	 // cairo_stroke_preserve(cr);
      }
   }
   std::vector<unsigned int> bond_indices;
   if (! use_highlight_bond_indices_flag) {
      bond_indices = find_bonds_for_atoms(highlight_atom_indices);
   } else {
      bond_indices = highlight_bond_indices;
   }
   for (std::size_t ib=0; ib<bond_indices.size(); ib++) {
      lig_build::pos_t pos_1 = atoms[bonds[bond_indices[ib]].get_atom_1_index()].atom_position;
      lig_build::pos_t pos_2 = atoms[bonds[bond_indices[ib]].get_atom_2_index()].atom_position;
      // set the colour
      lig_build::pos_t p1 = cairo_molecule_t::mol_coords_to_cairo_coords(pos_1, centre, scale);
      lig_build::pos_t p2 = cairo_molecule_t::mol_coords_to_cairo_coords(pos_2, centre, scale);
      std::cout << "line: " << p1 << " " << p2 << std::endl;
      cairo_set_line_width(cr, 0.3 * scale * median_bond_length_);
      cairo_move_to(cr, p1.x, p1.y);
      cairo_line_to(cr, p2.x, p2.y);
      cairo_stroke(cr);
   }
}


// npx is a option arg, default 300
// background_colour is an option arg, default Null
void
coot::cairo_png_depict_from_mmcif(const std::string &mmcif_file_name,
				  const std::string &comp_id,
				  const std::string &png_file_name,
				  unsigned int npx,
				  PyObject *background_colour_py) {

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

   int rn = 42;
   coot::protein_geometry geom;
   geom.set_verbose(false);
   geom.init_refmac_mon_lib(mmcif_file_name, rn);

   bool dictionary_only = false;

   if (dictionary_only) {
      if (geom.have_dictionary_for_residue_type_no_dynamic_add(comp_id)) {
	 std::pair<bool, coot::dictionary_residue_restraints_t> dp =
	    geom.get_monomer_restraints(comp_id, coot::protein_geometry::IMOL_ENC_ANY);
	 if (dp.first) {
	    const coot::dictionary_residue_restraints_t &d = dp.second;
	    try {
	       RDKit::RWMol rdkm = coot::rdkit_mol(d);
	       RDKit::MolOps::removeHs(rdkm);
	       // coot::remove_non_polar_Hs(&rdkm); either is good.
	       RDKit::MolOps::Kekulize(rdkm);
	       int iconf = RDDepict::compute2DCoords(rdkm, NULL, true);
	       RDKit::Conformer &conf = rdkm.getConformer(iconf);
	       RDKit::WedgeMolBonds(rdkm, &conf);
	       coot::cairo_molecule_t mol(&rdkm, iconf);
	       mol.render_to_file(png_file_name, npx);
	    }
	    catch (...) {
	    }
	 }
      }

   } else {

      try {
	 if (geom.have_dictionary_for_residue_type_no_dynamic_add(comp_id)) {
	    std::pair<bool, coot::dictionary_residue_restraints_t> dp =
	       geom.get_monomer_restraints(comp_id, coot::protein_geometry::IMOL_ENC_ANY);
	    if (dp.first) {

	       const coot::dictionary_residue_restraints_t &rest = dp.second;
	       bool idealized = false;
	       bool try_autoload_if_needed = false;
	       mmdb::Residue *r = geom.get_residue(comp_id, idealized, try_autoload_if_needed);
	       if (r) {
		  bool undelocalize_flag = true;
		  RDKit::RWMol mol_rw = coot::rdkit_mol(r, rest, "", undelocalize_flag);

		  // coot::remove_non_polar_Hs(&rdkm); either is good.
		  RDKit::MolOps::removeHs(mol_rw);
		  RDKit::MolOps::Kekulize(mol_rw);
		  int iconf = RDDepict::compute2DCoords(mol_rw, NULL, true);
		  RDKit::Conformer &conf = mol_rw.getConformer(iconf);
		  RDKit::WedgeMolBonds(mol_rw, &conf);
		  coot::cairo_molecule_t mol(&mol_rw, iconf);
		  std::pair<bool, colour_holder> bg_col;
		  bg_col.first = false;
		  if (background_colour_py) {
		     if  (PyString_Check(background_colour_py)) {
			std::string s = PyString_AsString(background_colour_py);
			bg_col.second = colour_holder(s);
			bg_col.first = true;
		     } else {
			// std::cout << "cairo_png_depict(): No background colour " << std::endl;
		     }
		  }
		  mol.render_to_file(png_file_name, npx, bg_col);
	       }
	    }
	 }
      }
      catch (...) {
      }
      
   }
#endif

}


#ifdef MAKE_ENHANCED_LIGAND_TOOLS
// highlightAtoms default arg null pointer.
// draw_atom_highlights is optional arg, default true
//
// If highlight_atom_indices comes from a GetSubstructMatch, then
// the result will be a tuple.  We should allow that too.
std::string
coot::cairo_png_string_from_mol(RDKit::ROMol *m, int iconf,
				PyObject *highlight_atom_list,
				PyObject *highlight_bond_list,
				PyObject *highlight_atom_colours_dict,
				PyObject *highlight_bond_colours_dict,
				unsigned int npx) {

   bool png_vs_svg_mode = true;
   return cairo_image_string_from_mol(m, iconf, highlight_atom_list, highlight_bond_list,
				      highlight_atom_colours_dict, highlight_bond_colours_dict,
				      png_vs_svg_mode, npx);
}
#endif // MAKE_ENHANCED_LIGAND_TOOLS

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
std::string
coot::cairo_svg_string_from_mol(RDKit::ROMol *m, int iconf,
				PyObject *highlight_atom_list,
				PyObject *highlight_bond_list,
				PyObject *highlight_atom_colours_dict,
				PyObject *highlight_bond_colours_dict,
				unsigned int npx) {

   bool png_vs_svg_mode = false; // make an svg
   return cairo_image_string_from_mol(m, iconf, highlight_atom_list, highlight_bond_list,
				      highlight_atom_colours_dict, highlight_bond_colours_dict,
				      png_vs_svg_mode, npx);
}
#endif // MAKE_ENHANCED_LIGAND_TOOLS

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
// highlightAtoms default arg null pointer.
// draw_atom_highlights is optional arg, default true
//
// If highlight_atom_indices comes from a GetSubstructMatch, then
// the result will be a tuple.  We should allow that too.
std::string
coot::cairo_image_string_from_mol(RDKit::ROMol *m, int iconf,
				  PyObject *highlight_atom_list,
				  PyObject *highlight_bond_list,
				  PyObject *highlight_atom_colours_dict,
				  PyObject *highlight_bond_colours_dict,
				  bool png_vs_svg_mode,
				  unsigned int npx) {

   // we need to distinguish between highlight_bond_list being
   // an empty list (which was set by the user)
   // and None (default value)

   std::string s;
   int n_confs = m->getNumConformers();

   if (n_confs == 0) {
      std::cout << "WARNING:: molecule has no conformers" << std::endl;
   } else {
      if (iconf == -1) {
	 // use the most recent one
	 iconf = n_confs -1;
      }

      if (highlight_atom_colours_dict) {
	 if (PyDict_Check(highlight_atom_colours_dict)) {
	    unsigned int l = PyDict_Size(highlight_atom_colours_dict);
	    if (l > 0) {
	       for (std::size_t i=0; i<l; i++) {
		  // PyObject *obj_py = PyDict_GetItem(highlight_atom_colours_dict, i);
	       }
	    }
	 }
      }

      RDKit::Conformer &conf = m->getConformer(iconf);
      RDKit::WedgeMolBonds(*m, &conf);
      coot::cairo_molecule_t mol(m, iconf);
      std::vector<unsigned int> highlight_atom_indices;
      std::vector<unsigned int> highlight_bond_indices;
      bool use_highlight_bond_indices_flag = false;
      if (highlight_atom_list) {
	 if (PyList_Check(highlight_atom_list)) {
	    unsigned int l = PyList_Size(highlight_atom_list);
	    if (l > 0) {
	       for (std::size_t i=0; i<l; i++) {
		  PyObject *obj_py = PyList_GetItem(highlight_atom_list, i);
		  if PyInt_Check(obj_py) {
		     long item = PyInt_AsLong(obj_py);
		     if (item >= 0) {
			highlight_atom_indices.push_back(static_cast<unsigned int>(item));
		     }
		  }
	       }
	    }
	 } else {
	    if (PyTuple_Check(highlight_atom_list)) {
	       unsigned int l = PyTuple_Size(highlight_atom_list);
	       if (l > 0) {
		  for (std::size_t i=0; i<l; i++) {
		     PyObject *obj_py = PyTuple_GetItem(highlight_atom_list, i);
		     if PyInt_Check(obj_py) {
			long item = PyInt_AsLong(obj_py);
			if (item >= 0) {
			   highlight_atom_indices.push_back(static_cast<unsigned int>(item));
			}
		     }
		  }
	       }
	    }
	 }
      }
      if (highlight_bond_list) {
	 if (PyList_Check(highlight_bond_list)) {
	    // user has set
	    use_highlight_bond_indices_flag = true;
	    unsigned int l = PyList_Size(highlight_bond_list);
	    if (l > 0) {
	       for (std::size_t i=0; i<l; i++) {
		  PyObject *obj_py = PyList_GetItem(highlight_bond_list, i);
		  if PyInt_Check(obj_py) {
		     long item = PyInt_AsLong(obj_py);
		     if (item >= 0) {
			highlight_bond_indices.push_back(static_cast<unsigned int>(item));
		     }
		  }
	       }
	    }
	 }
      }
      if (png_vs_svg_mode)
	 s = mol.render_to_png_string(highlight_atom_indices, highlight_bond_indices,
				      use_highlight_bond_indices_flag, npx);
      else
	 s = mol.render_to_svg_string(highlight_atom_indices, highlight_bond_indices,
				      use_highlight_bond_indices_flag, npx);
   }
   return s;
}
#endif
