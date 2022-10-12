
#include "utils/coot-utils.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot_molecule.hh"
#include "ideal/pepflip.hh"
#include "rama-plot-phi-psi.hh"

bool
coot::molecule_t::is_valid_model_molecule() const {

   bool status = false;
   if (atom_sel.mol)
      status = true;
   return status;

}

bool
coot::molecule_t::is_valid_map_molecule() const {

   bool status = false;
   if (! xmap.is_null()) {
      status = true;
   }
   return status;
}

std::pair<bool, coot::residue_spec_t>
coot::molecule_t::cid_to_residue_spec(const std::string &cid) {

   bool status = false;
   coot::residue_spec_t rs;
   if (atom_sel.mol) {
      int selHnd = atom_sel.mol->NewSelection(); // d
      mmdb::Residue **SelResidues;
      int nSelResidues = 0;
      atom_sel.mol->Select(selHnd, mmdb::STYPE_CHAIN, cid.c_str(), mmdb::SKEY_NEW);
      atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);
      atom_sel.mol->DeleteSelection(selHnd);
      if (nSelResidues > 0) {
         mmdb::Residue *residue_p = SelResidues[0];
         coot::residue_spec_t rs_inner(residue_p);
         rs = rs_inner;
         status = true;
      }
   }
   return std::make_pair(status, rs);
}

int coot::molecule_t::flipPeptide(const coot::residue_spec_t &rs, const std::string &alt_conf) {

   int result = coot::pepflip(atom_sel.mol, rs.chain_id, rs.res_no, rs.ins_code, alt_conf);
   return result;

}

// public - because currently making bonds is not done on molecule construction
void
coot::molecule_t::make_bonds(coot::protein_geometry *geom, coot::rotamer_probability_tables *rot_prob_tables_p) {

   bonds_box_type = coot::COLOUR_BY_CHAIN_BONDS;
   make_bonds_type_checked(geom, rot_prob_tables_p, __FUNCTION__);

   std::cout << "debug:: in molecule_t::make_bonds() " << bonds_box.n_bonds() << " bonds " << bonds_box.n_atoms() << " atoms "
             << std::endl;
}



// private
void
coot::molecule_t::makebonds(coot::protein_geometry *geom, coot::rotamer_probability_tables *rotamer_tables_p, std::set<int> &no_bonds_to_these_atoms) {

   bool force_rebond = true;
   bool do_rotamer_markup = true; // pass this
   make_colour_by_chain_bonds(geom, no_bonds_to_these_atoms, true, false, do_rotamer_markup, rotamer_tables_p, force_rebond);

}


void
coot::molecule_t::make_colour_by_chain_bonds(coot::protein_geometry *geom,
                                             const std::set<int> &no_bonds_to_these_atoms,
                                             bool change_c_only_flag,
                                             bool goodsell_mode,
                                             bool do_rota_markup,
                                             coot::rotamer_probability_tables *tables_p,
                                             bool force_rebonding) {

   // We don't want to rebond if we don't have to (i.e the mode requested is the current mode)
   // so check the previous value of bonds_box_type so that we can know if it can be skipped.

   bool draw_hydrogens_flag = true; // pass this
   bool draw_missing_loops_flag = true; // pass this

   Bond_lines_container bonds(geom, no_bonds_to_these_atoms, draw_hydrogens_flag);

   std::cout << "debug:: in make_colour_by_chain_bonds() adding tables_p " << tables_p << std::endl;

   bonds.add_rotamer_tables(tables_p);
   bonds.do_colour_by_chain_bonds(atom_sel, false, imol_no, draw_hydrogens_flag,
                                  draw_missing_loops_flag,
                                  change_c_only_flag, goodsell_mode, do_rota_markup);

   // std::cout << "------------------- calling make_graphical_bonds_no_thinning() " << std::endl;

   bonds_box = bonds.make_graphical_bonds_no_thinning(); // make_graphical_bonds() is pretty
                                                         // stupid when it comes to thining.

   // bonds_box = bonds.make_graphical_bonds(); // make_graphical_bonds() is pretty
                                                // stupid when it comes to thining.

   // testing previous values of bonds_box_type
   if (bonds_box_type != coot::COLOUR_BY_CHAIN_BONDS)
      force_rebonding = true;

   if (goodsell_mode)
      if (bonds_box_type != coot::COLOUR_BY_CHAIN_GOODSELL)
         force_rebonding = true;

   bonds_box_type = coot::COLOUR_BY_CHAIN_BONDS;

   if (goodsell_mode)
      bonds_box_type = coot::COLOUR_BY_CHAIN_GOODSELL;

   // 20221011-PE Hmm... is this needed in this API? I don't think so
   //
   // if (force_rebonding)
   //    make_glsl_bonds_type_checked(__FUNCTION__);

}

void
coot::molecule_t::make_ca_bonds() {

}


void
coot::molecule_t::make_bonds_type_checked(coot::protein_geometry *geom_p, coot::rotamer_probability_tables *rotamer_probability_tables_p, const char *caller) {

   bool draw_missing_loops_flag = false; // pass this
   bool rotate_colour_map_on_read_pdb_c_only_flag = true; // pass this or make class data item

   bool debug = false;

   // Note caller can be 0 (e.g. with clang) - so be aware of that when debugging.

   std::string caller_s("NULL");
   if (caller) caller_s = std::string(caller);

   bool is_intermediate_atoms_molecule = false; // 20221005-PE IMPORT-HACK
   if (debug)
      std::cout << "debug:: plain make_bonds_type_checked() --------start--------- called by "
                << caller_s << "() with is_intermediate_atoms_molecule: " << is_intermediate_atoms_molecule
                << std::endl;
   if (debug)
      std::cout << "--------- make_bonds_type_checked() called with bonds_box_type "
                << bonds_box_type << " vs "
                << "NORMAL_BONDS " << coot::NORMAL_BONDS << " "
                << "BONDS_NO_HYDROGENS " << coot::BONDS_NO_HYDROGENS << " "
                << "COLOUR_BY_CHAIN_BONDS " << coot::COLOUR_BY_CHAIN_BONDS << " "
                << "COLOUR_BY_MOLECULE_BONDS " << coot::COLOUR_BY_MOLECULE_BONDS << " "
                << "CA_BONDS " << coot::CA_BONDS << " "
                << "CA_BONDS_PLUS_LIGANDS " << coot::CA_BONDS_PLUS_LIGANDS << " "
                << "COLOUR_BY_USER_DEFINED_COLOURS___BONDS " << coot::COLOUR_BY_USER_DEFINED_COLOURS____BONDS << " "
                << std::endl;

   // Delete this in due course
   // graphics_info_t g; // urgh!  (But the best solution?)

   bool force_rebonding = true; // if we get here, this must be true (?)

   // coot::protein_geometry *geom_p = g.Geom_p();

   std::set<int> dummy;

   if (bonds_box_type == coot::NORMAL_BONDS)
      makebonds(geom_p, nullptr, dummy);

   if (bonds_box_type == coot::BONDS_NO_HYDROGENS)
      makebonds(geom_p, nullptr, dummy);
   if (bonds_box_type == coot::CA_BONDS)
      make_ca_bonds();
   if (bonds_box_type == coot::COLOUR_BY_CHAIN_BONDS || bonds_box_type == coot::COLOUR_BY_CHAIN_GOODSELL) {
      // Baah, we have to use the static in graphics_info_t here as it
      // is not a per-molecule property.
      std::set<int> s;
      bool goodsell_mode = false;
      if (bonds_box_type == coot::COLOUR_BY_CHAIN_GOODSELL)
         goodsell_mode = true;
      bool do_rota_markup = true;

      make_colour_by_chain_bonds(geom_p, s, rotate_colour_map_on_read_pdb_c_only_flag, goodsell_mode, do_rota_markup, rotamer_probability_tables_p, force_rebonding);
   }

#if 0 // not implemenented yet
   if (bonds_box_type == coot::COLOUR_BY_MOLECULE_BONDS)
      make_colour_by_molecule_bonds(force_rebonding);
   if (bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS)
      make_ca_plus_ligands_bonds(g.Geom_p());
   if (bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS_AND_SIDECHAINS)
      make_ca_plus_ligands_and_sidechains_bonds(g.Geom_p());
   if (bonds_box_type == coot::BONDS_NO_WATERS)
      bonds_no_waters_representation();
   if (bonds_box_type == coot::BONDS_SEC_STRUCT_COLOUR)
      bonds_sec_struct_representation();
   if (bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR)
      ca_plus_ligands_sec_struct_representation(g.Geom_p());
   if (bonds_box_type == coot::COLOUR_BY_RAINBOW_BONDS)
      ca_plus_ligands_rainbow_representation(g.Geom_p());
   if (bonds_box_type == coot::COLOUR_BY_OCCUPANCY_BONDS)
      occupancy_representation();
   if (bonds_box_type == coot::COLOUR_BY_B_FACTOR_BONDS)
      b_factor_representation();
   if (bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS_B_FACTOR_COLOUR)
      b_factor_representation_as_cas();
   if (bonds_box_type == coot::COLOUR_BY_USER_DEFINED_COLOURS____BONDS)
      user_defined_colours_representation(g.Geom_p(), true, g.draw_missing_loops_flag); // hack,
                                                             // because we need to remeber somehow
                                                             // if this was called with all-atom or CA-only.
                                                             // See c-interface.cc
                                                             // graphics_to_user_defined_atom_colours_representation()
                                                             // Perhaps we need two functions
                                                             // user_defined_colours_representation_all()
                                                             // user_defined_colours_representation_Calpha() [+ ligands]
   if (bonds_box_type == coot::COLOUR_BY_USER_DEFINED_COLOURS_CA_BONDS)
      user_defined_colours_representation(geom_p, false, draw_missing_loops_flag); // hack,

#endif


#if 0 // 20221005-PE not sure what these are
   // all these will need to be changed or removed
   update_additional_representations(glci, g.Geom_p());
   update_fixed_atom_positions();
   update_ghosts();
   update_extra_restraints_representation();
#endif

   if (debug) {
      std::cout << "debug:: -------------- make_bonds_type_checked() done " << std::endl;
   }
}

std::vector<std::pair<coot::Cartesian, coot::util::phi_psi_t> >
coot::molecule_t::ramachandran_validation() const {

   auto get_HA_unit_vector = [] (mmdb::Residue *r) {
      bool status = false;
      coot::Cartesian dir;
      mmdb::Atom *CA = r->GetAtom(" CA ");
      mmdb::Atom *C  = r->GetAtom(" C  ");
      mmdb::Atom *N  = r->GetAtom(" N  ");
      mmdb::Atom *CB = r->GetAtom(" CB ");

      if (CA && C && N && CB) {
         coot::Cartesian ca_pos(CA->x, CA->y, CA->z);
         coot::Cartesian  c_pos( C->x,  C->y,  C->z);
         coot::Cartesian  n_pos( N->x,  N->y,  N->z);
         coot::Cartesian cb_pos(CB->x, CB->y, CB->z);
         coot::Cartesian dir_1 = ca_pos - c_pos;
         coot::Cartesian dir_2 = ca_pos - n_pos;
         coot::Cartesian dir_3 = ca_pos - cb_pos;
         coot::Cartesian r = dir_1 + dir_2 + dir_3;
         dir = r.unit();
         status = true;
      } else {
         if (CA && C && N) {
            coot::Cartesian ca_pos(CA->x, CA->y, CA->z);
            coot::Cartesian  c_pos( C->x,  C->y,  C->z);
            coot::Cartesian  n_pos( N->x,  N->y,  N->z);
            coot::Cartesian dir_1 = ca_pos - c_pos;
            coot::Cartesian dir_2 = ca_pos - n_pos;
            coot::Cartesian r = dir_1 + dir_2;
            dir = r.unit();
            status = true;
         }
      }
      return std::make_pair(status, dir);
   };

   std::vector<std::pair<coot::Cartesian, coot::util::phi_psi_t> > v;

   float rama_ball_pos_offset_scale = 0.6;

   rama_plot::phi_psis_for_model_t ppm(atom_sel.mol);
   // This: std::map<coot::residue_spec_t, phi_psi_t> ppm.phi_psi  is now filled

   std::map<residue_spec_t, rama_plot::phi_psi_t>::const_iterator it;
   for (it=ppm.phi_psi.begin(); it!=ppm.phi_psi.end(); ++it) {
      const auto &phi_psi(it->second);
      mmdb::Residue *rp = phi_psi.residue_prev;
      mmdb::Residue *rt = phi_psi.residue_this;
      mmdb::Residue *rn = phi_psi.residue_next;
      if (rp && rt && rn) {
         mmdb::Atom *at = rt->GetAtom(" CA "); // 20221006-PE alt-confs another day
         if (at) {
            coot::Cartesian pos(at->x, at->y, at->z);
            std::pair<bool, coot::Cartesian> hav = get_HA_unit_vector(rt);
            coot::Cartesian offset(0,0,rama_ball_pos_offset_scale);
            if (hav.first) offset = hav.second * rama_ball_pos_offset_scale;
            coot::util::phi_psi_t cupp(rp, rt, rn);
            auto p = std::make_pair(pos + offset, cupp);
            v.push_back(p);
         }
      }
   }

   return v;
}


// returns either the specified atom or null if not found
mmdb::Atom *
coot::molecule_t::get_atom(const coot::atom_spec_t &atom_spec) const {

   mmdb::Atom *at = coot::util::get_atom(atom_spec, atom_sel.mol);
   return at;
}


#include "utils/dodec.hh"

// returns either the specified residue or null if not found
mmdb::Residue *
coot::molecule_t::get_residue(const coot::residue_spec_t &residue_spec) const {

   mmdb::Residue *r = coot::util::get_residue(residue_spec, atom_sel.mol);
   return r;

}

coot::simple_mesh_t
coot::molecule_t::get_rotamer_dodecs(coot::protein_geometry *geom_p,
                                     coot::rotamer_probability_tables *rpt) {

   // THis function is an API version of:
   //
   // void
   // Mesh::make_graphical_bonds_rotamer_dodecs(const graphical_bonds_container &gbc,
   // const glm::vec3 &screen_up_dir)

   simple_mesh_t m;

   // use bonds_box

   auto cartesian_to_glm = [] (const coot::Cartesian &c) {
      return glm::vec3(c.x(), c.y(), c.z()); };

   auto clipper_to_glm = [] (const clipper::Coord_orth &c) {
      return glm::vec3(c.x(), c.y(), c.z()); };

   auto clipper_to_cartesian = [] (const clipper::Coord_orth &c) {
      return Cartesian(c.x(), c.y(), c.z()); };

   auto colour_holder_to_glm = [] (const coot::colour_holder &ch) {
                                  return glm::vec4(ch.red, ch.green, ch.blue, 1.0f);
                               };

   auto get_HA_unit_vector = [] (mmdb::Residue *r) {
      bool status = false;
      coot::Cartesian dir;
      mmdb::Atom *CA = r->GetAtom(" CA ");
      mmdb::Atom *C  = r->GetAtom(" C  ");
      mmdb::Atom *N  = r->GetAtom(" N  ");
      mmdb::Atom *CB = r->GetAtom(" CB ");

      if (CA && C && N && CB) {
         coot::Cartesian ca_pos(CA->x, CA->y, CA->z);
         coot::Cartesian  c_pos( C->x,  C->y,  C->z);
         coot::Cartesian  n_pos( N->x,  N->y,  N->z);
         coot::Cartesian cb_pos(CB->x, CB->y, CB->z);
         coot::Cartesian dir_1 = ca_pos - c_pos;
         coot::Cartesian dir_2 = ca_pos - n_pos;
         coot::Cartesian dir_3 = ca_pos - cb_pos;
         coot::Cartesian r = dir_1 + dir_2 + dir_3;
         dir = r.unit();
         status = true;
      } else {
         if (CA && C && N) {
            coot::Cartesian ca_pos(CA->x, CA->y, CA->z);
            coot::Cartesian  c_pos( C->x,  C->y,  C->z);
            coot::Cartesian  n_pos( N->x,  N->y,  N->z);
            coot::Cartesian dir_1 = ca_pos - c_pos;
            coot::Cartesian dir_2 = ca_pos - n_pos;
            coot::Cartesian r = dir_1 + dir_2;
            dir = r.unit();
            status = true;
         }
      }
      return std::make_pair(status, dir);
   };

   std::set<int> dummy;
   bool do_rota_markup = true;
   make_colour_by_chain_bonds(geom_p, dummy, true, false, do_rota_markup, rpt, true);

   if (bonds_box.n_rotamer_markups > 0) {

      auto &vertices = m.vertices;
      auto &triangles = m.triangles;

      glm::vec4 col(0.6, 0.2, 0.8, 1.0); // starting colour
      dodec d;
      std::vector<clipper::Coord_orth> coords = d.coords();
      std::vector<glm::vec3> dodec_postions(coords.size());
      for (unsigned int i=0; i<coords.size(); i++)
         dodec_postions[i] = clipper_to_glm(coords[i]);

      std::vector<coot::api::vnc_vertex> dodec_vertices;
      std::vector<g_triangle> dodec_triangles;
      dodec_triangles.reserve(36);

      for (unsigned int iface=0; iface<12; iface++) {

         std::vector<coot::api::vnc_vertex> face_verts;
         std::vector<g_triangle> face_triangles;
         face_triangles.reserve(3);

         std::vector<unsigned int> indices_for_face = d.face(iface);
         glm::vec3 ns(0,0,0);
         for (unsigned int j=0; j<5; j++)
            ns += dodec_postions[indices_for_face[j]];
         glm::vec3 normal = glm::normalize(ns);

         for (unsigned int j=0; j<5; j++) {
            glm::vec3 &pos = dodec_postions[indices_for_face[j]];
            coot::api::vnc_vertex v(0.5f * pos, normal, col);
            face_verts.push_back(v);
         }

         face_triangles.push_back(g_triangle(0,1,2));
         face_triangles.push_back(g_triangle(0,2,3));
         face_triangles.push_back(g_triangle(0,3,4));

         unsigned int idx_base = dodec_vertices.size();
         unsigned int idx_tri_base = dodec_triangles.size();
         dodec_vertices.insert(dodec_vertices.end(), face_verts.begin(), face_verts.end());
         dodec_triangles.insert(dodec_triangles.end(), face_triangles.begin(), face_triangles.end());
         for (unsigned int jj=idx_tri_base; jj<dodec_triangles.size(); jj++)
            dodec_triangles[jj].rebase(idx_base);
      }

      // now there is a dodec at the origin, dodec_vertices and dodec_triangle

      // let's make copies of that and move them around to the residues

      double rama_ball_pos_offset_scale = 1.2; // may need tweaking
      for (int i=0; i<bonds_box.n_rotamer_markups; i++) {
         const rotamer_markup_container_t &rm = bonds_box.rotamer_markups[i];
         const residue_spec_t &residue_spec = rm.spec;
         mmdb::Residue *residue_p = get_residue(residue_spec);
         Cartesian offset(0,0,rama_ball_pos_offset_scale);
         if (residue_p) {
            std::pair<bool, coot::Cartesian> hav = get_HA_unit_vector(residue_p);
            if (hav.first) offset = hav.second; // * rama_ball_pos_offset_scale;
         }

         glm::vec3 atom_pos = cartesian_to_glm(rm.pos) + cartesian_to_glm(offset);

         std::vector<coot::api::vnc_vertex> this_dodec_vertices = dodec_vertices; // at the origin to start

         // now move it.
         for (unsigned int j=0; j<dodec_vertices.size(); j++) {
            this_dodec_vertices[j].pos  += atom_pos;
            this_dodec_vertices[j].color = colour_holder_to_glm(rm.col);
         }
         unsigned int idx_base = vertices.size();
         unsigned int idx_tri_base = triangles.size();
         vertices.insert(vertices.end(), this_dodec_vertices.begin(), this_dodec_vertices.end());
         triangles.insert(triangles.end(), dodec_triangles.begin(), dodec_triangles.end());
         for (unsigned int jj=idx_tri_base; jj<triangles.size(); jj++)
            triangles[jj].rebase(idx_base);
      }
   }
   return m;
}

int
coot::molecule_t::auto_fit_rotamer(const std::string &chain_id, int res_no, const std::string &ins_code,
                                   const clipper::Xmap<float> &xmap) {

   int status = 0;

   return status;
}
