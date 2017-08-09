
#include "cairo-molecule.hh"

// this needs to be here (also).  It is in wmolecule.cc and hence the library also.
// But if this is not here I get unresovled symbol for this destructor when compiling
// this exectuable on the mac (clang).
template<class cairo_atom_t, class cairo_bond_t> lig_build::molecule_t<cairo_atom_t, cairo_bond_t>::~molecule_t() {}

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
	 RDKit::ATOM_SPTR at_p = (*rdkm)[iat];
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
	 RDKit::ATOM_SPTR at_p = (*rdkm)[iat];
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
}
#endif // MAKE_ENHANCED_LIGAND_TOOLS

void
coot::cairo_atom_t::set_colour(cairo_t *cr) const {

   if (element == "O")  cairo_set_source_rgb(cr, 0.8, 0.0, 0);
   if (element == "N")  cairo_set_source_rgb(cr, 0.2, 0.2, 0.8);
   if (element == "S")  cairo_set_source_rgb(cr, 0.5, 0.5, 0);
   if (element == "F") cairo_set_source_rgb (cr, 0,   0.5, 0);
   if (element == "Cl") cairo_set_source_rgb(cr, 0,   0.5, 0);
   if (element == "Br") cairo_set_source_rgb(cr, 0.5, 0.2, 0);
   if (element == "I")  cairo_set_source_rgb(cr, 0.3, 0.0, 0.3);
   if (element == "C")  cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
}

void
coot::cairo_atom_t::make_text_item(cairo_t *cr,
				   const lig_build::atom_id_info_t &atom_id_info,
				   const lig_build::pos_t &centre, double scale) const {

   for (unsigned int i=0; i<atom_id_info.n_offsets(); i++) {

      cairo_set_font_size(cr, 0.66 * scale);
      lig_build::pos_t p = cairo_molecule_t::mol_coords_to_cairo_coords(atom_position, centre, scale);
      p += atom_id_info.offsets[i].tweak * scale * 0.045;

      if (atom_id_info[i].text_pos_offset == lig_build::offset_text_t::UP)
	 p.y -= 0.55 * scale;
      if (atom_id_info[i].text_pos_offset == lig_build::offset_text_t::DOWN)
	 p.y += 0.55 * scale;

      if (atom_id_info.size_hint == -1)
	 cairo_set_font_size(cr, 0.66 * scale * 0.7);

      if (atom_id_info.offsets[i].subscript) {
	 p.y += 0.3 * scale;
	 cairo_set_font_size(cr, 0.66 * scale * 0.8);
      }
      if (atom_id_info.offsets[i].superscript) {
	 cairo_set_font_size(cr, 0.66 * scale * 0.8);
	 p.y -= 0.3 * scale;
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

void
coot::cairo_bond_t::draw_bond(cairo_t *cr,
			      const lig_build::pos_t &pos_1_in,
			      const lig_build::pos_t &pos_2_in,
			      bool shorten_first, bool shorten_second,
			      bool at_1_in_ring_flag, bool at_2_in_ring_flag,
			      lig_build::bond_t::bond_type_t bt,
			      const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_first_atom,
			      const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
			      const lig_build::pos_t &centre, double scale) {

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

   // std::cout << "draw_bond for bt " << bt << std::endl;

   switch (bt) {
   case lig_build::bond_t::SINGLE_BOND:
   case lig_build::bond_t::SINGLE_OR_DOUBLE:
   case lig_build::bond_t::SINGLE_OR_AROMATIC:
   case lig_build::bond_t::DELOC_ONE_AND_A_HALF:
   case lig_build::bond_t::BOND_ANY:
      {
	 p1 = cairo_molecule_t::mol_coords_to_cairo_coords(pos_1, centre, scale);
	 p2 = cairo_molecule_t::mol_coords_to_cairo_coords(pos_2, centre, scale);
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
	    draw_double_bond(cr, pos_1, pos_2,
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
	    //
	    if (! at_2_in_ring_flag) {
	       draw_sheared_or_darted_wedge_bond(cr, pos_1, pos_2, other_connections_to_second_atom, centre, scale);
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

   // filled shape
   //
   p = cairo_molecule_t::mol_coords_to_cairo_coords(v[0], centre, scale);
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

   cairo_move_to(cr, p1.x, p1.y);
   cairo_line_to(cr, p2.x, p2.y);
   cairo_stroke(cr);

   if (dashed_inner) {
      double dashlength = 0.008;
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
void
coot::cairo_bond_t::draw_double_bond(cairo_t *cr,
				     const lig_build::pos_t &pos_1,
				     const lig_build::pos_t &pos_2,
				     bool shorten_first, bool shorten_second,
				     const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_first_atom,
				     const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
				     const lig_build::pos_t &centre, double scale) {

   std::pair<std::pair<lig_build::pos_t, lig_build::pos_t>, std::pair<lig_build::pos_t, lig_build::pos_t> > p = make_double_bond(pos_1, pos_2);

   lig_build::pos_t p1 = cairo_molecule_t::mol_coords_to_cairo_coords(p.first.first,  centre, scale);
   lig_build::pos_t p2 = cairo_molecule_t::mol_coords_to_cairo_coords(p.first.second, centre, scale);

   if ((other_connections_to_second_atom.size() == 0) ||
       (other_connections_to_first_atom.size()  == 0)) {

      // boring case :-)
      //
      cairo_move_to(cr, p1.x, p1.y);
      cairo_line_to(cr, p2.x, p2.y);
      cairo_stroke(cr);

      p1 = cairo_molecule_t::mol_coords_to_cairo_coords(p.second.first,  centre, scale);
      p2 = cairo_molecule_t::mol_coords_to_cairo_coords(p.second.second, centre, scale);

      cairo_move_to(cr, p1.x, p1.y);
      cairo_line_to(cr, p2.x, p2.y);
      cairo_stroke(cr);

   } else {
      // elegant case, shorten a bond
      std::pair<std::pair<lig_build::pos_t, lig_build::pos_t>, std::pair<lig_build::pos_t, lig_build::pos_t> > bonds =
	 make_double_bond(pos_1, pos_2, shorten_first, shorten_second, other_connections_to_first_atom, other_connections_to_second_atom);
      p1 = cairo_molecule_t::mol_coords_to_cairo_coords(bonds.first.first,  centre, scale);
      p2 = cairo_molecule_t::mol_coords_to_cairo_coords(bonds.first.second, centre, scale);
      cairo_move_to(cr, p1.x, p1.y);
      cairo_line_to(cr, p2.x, p2.y);
      p1 = cairo_molecule_t::mol_coords_to_cairo_coords(bonds.second.first,  centre, scale);
      p2 = cairo_molecule_t::mol_coords_to_cairo_coords(bonds.second.second, centre, scale);
      cairo_move_to(cr, p1.x, p1.y);
      cairo_line_to(cr, p2.x, p2.y);
      cairo_stroke(cr);
   }
}

// not const because it changes atom ids.
void
coot::cairo_molecule_t::render(cairo_t *cr) {

   // ------------ the font size and the line width depend on the size of the
   //              extents (and here, scale), the smaller the scale the smaller
   //              the font
   //
   // cairo_set_source_rgb (cr, 0, 0, 0);
   // cairo_rectangle (cr, 0.25, 0.25, 0.5, 0.5);
   // cairo_set_source_rgb(cr, 0, 0.5, 0);
   // cairo_rectangle (cr, 0.01, 0.01, 0.98, 0.98);
   // cairo_stroke(cr);
   cairo_set_source_rgb(cr, 0, 0, 0);

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
      scale = 0.75/delta;
   if (scale > scale_lim)
      scale = scale_lim;

   // std::cout << "scale: " << scale << std::endl;

   cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
   cairo_set_line_width(cr, 0.11 * scale);
   cairo_set_font_size(cr,  0.66 * scale);
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

	 std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > other_connections_to_first_atom =
	    make_other_connections_to_first_atom_info(ib);
	 std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > other_connections_to_second_atom =
	    make_other_connections_to_second_atom_info(ib);

	 bool at_1_in_ring_flag = in_ring_p(idx_1);
	 bool at_2_in_ring_flag = in_ring_p(idx_2);

	 bonds[ib].draw_bond(cr, pos_1, pos_2, shorten.first, shorten.second,
			     at_1_in_ring_flag,
			     at_2_in_ring_flag,
			     bt,
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
	 atoms[iat].make_text_item(cr, atom_id_info, centre, scale);
      } else {
	 // a super-atom carbon
	 if (local_bonds.size() == 1) {
	    atoms[iat].set_colour(cr);
	    lig_build::atom_id_info_t atom_id_info =
	       make_atom_id_by_using_bonds(iat, ele, local_bonds, gl_flag);
	    atoms[iat].make_text_item(cr, atom_id_info, centre, scale);
	 }
      }
   }
}

void
coot::cairo_molecule_t::render_to_file(const std::string &png_file_name, unsigned int npx) {

   cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, npx, npx);
   cairo_t *cr = cairo_create(surface);

   // Drawing space is 0->1 with 0,0 top left
   // if cairo_image_surface_create called with 240,240, then cairo_scale called with
   // 120,120 puts the image, half-size in top left quadrant
   cairo_scale(cr, npx, npx);
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
coot::cairo_molecule_t::render_to_string(unsigned int npx) {

   std::string s;
   s.reserve(12000);

   cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, npx, npx);
   cairo_t *cr = cairo_create(surface);
   cairo_scale(cr, npx, npx);
   render(cr);
   cairo_surface_write_to_png_stream(surface, png_stream_writer, (void *) &s);
   cairo_destroy(cr);
   cairo_surface_destroy(surface);
   return s;
}

void
coot::cairo_png_depict(const std::string &mmcif_file_name,
		       const std::string &comp_id,
		       const std::string png_file_name,
		       unsigned int npx) {

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
		  mol.render_to_file(png_file_name, npx);
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
std::string
coot::cairo_png_string_from_mol(RDKit::ROMol *m, int iconf, unsigned int npx) {

   std::string s;
   int n_confs = m->getNumConformers();

   if (n_confs > 0) {
      if (iconf == -1) {
	 // use the most recent one
	 iconf = n_confs -1;
      }
      RDKit::Conformer &conf = m->getConformer(iconf);
      RDKit::WedgeMolBonds(*m, &conf);
      coot::cairo_molecule_t mol(m, iconf);
      s = mol.render_to_string(npx);
   }
   return s;
#endif
}
