
#include <iomanip>

#include "molecules_container.hh"
#include "ideal/pepflip.hh"
#include "coot-utils/coot-map-utils.hh"

#include "coords/Bond_lines.h"
#include "oct.hh"

std::atomic<bool> molecules_container_t::on_going_updating_map_lock(false);


bool
molecules_container_t::is_valid_model_molecule(int imol) const {
   bool status = false;
   if (imol >= 0) {
      int ms = molecules.size();
      if (imol < ms) {
         status = molecules[imol].is_valid_model_molecule();
      }
   }
   return status;
}

bool
molecules_container_t::is_valid_map_molecule(int imol) const {
   bool status = false;
   if (imol >= 0) {
      int ms = molecules.size();
      if (imol < ms) {
         status = molecules[imol].is_valid_map_molecule();
      }
   }
   return status;
}

void
molecules_container_t::display_molecule_names_table() const {

   for (unsigned int imol=0; imol<molecules.size(); imol++) {
      std::cout << imol << " " << std::setw(40) << molecules[imol].get_name() << std::endl;
   }

}

coot::atom_spec_t
molecules_container_t::atom_cid_to_atom_spec(int imol, const std::string &cid) const {

   coot::atom_spec_t spec;
   if (is_valid_model_molecule(imol)) {
      auto p = molecules[imol].cid_to_atom_spec(cid);
      if (p.first) {
         spec = p.second;
      } else {
         std::cout << "WARNING:: molecule_class_info_t::atom_cid_to_atom_spec() no matching atom " << cid << std::endl;
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return spec;
}


coot::residue_spec_t
molecules_container_t::residue_cid_to_residue_spec(int imol, const std::string &cid) const   {

   coot::residue_spec_t spec;

   if (is_valid_model_molecule(imol)) {
      auto p = molecules[imol].cid_to_residue_spec(cid);
      if (p.first) {
         spec = p.second;
      } else {
         std::cout << "WARNING:: molecule_class_info_t::residue_cid_to_residue_spec() no matching residue " << cid << std::endl;
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return spec;
}


void
molecules_container_t::geometry_init_standard() {
   geom.init_standard();
}


int
molecules_container_t::undo(int imol) {
   int status = 0;
   if (is_valid_model_molecule(imol)) {
      status = molecules[imol].undo();
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}

int
molecules_container_t::redo(int imol) {
   int status = 0;
   if (is_valid_model_molecule(imol)) {
      status = molecules[imol].redo();
         } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}




int
molecules_container_t::flip_peptide(int imol, const coot::residue_spec_t &rs, const std::string &alt_conf) {

   int result = 0;
   if (is_valid_model_molecule(imol)) {
      result = molecules[imol].flip_peptide(rs, alt_conf);
   }
   return result;
}

int
molecules_container_t::flip_peptide_using_cid(int imol, const std::string &cid, const std::string &alt_conf) {

   int result = 0;
   if (is_valid_model_molecule(imol)) {
      auto &m = molecules[imol];
      std::pair<bool, coot::residue_spec_t> rs = m.cid_to_residue_spec(cid);
      if (rs.first)
         result = molecules[imol].flip_peptide(rs.second, alt_conf);
   }
   return result;
}


int
molecules_container_t::read_pdb(const std::string &file_name) {

   int status = -1;
   atom_selection_container_t asc = get_atom_selection(file_name);
   if (asc.read_success) {
      // 20221011-PE this constructor doesn't (at the moment) call make_bonds(). I
      // don't know if that is a good idea.
      int imol = molecules.size();
      coot::molecule_t m = coot::molecule_t(asc, imol, file_name);
      // m.make_bonds(&geom, &rot_prob_tables); // where does this go? Here or as a container function?
      molecules.push_back(m);
      status = imol;
   }
   return status;
}

int
molecules_container_t::write_coordinates(int imol, const std::string &file_name) const {
   int status = 0;
   if (is_valid_model_molecule(imol)) {
      status = molecules[imol].write_coordinates(file_name);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}


int
molecules_container_t::read_mtz(const std::string &file_name,
                                const std::string &f, const std::string &phi, const std::string &weight,
                                bool use_weight, bool is_a_difference_map) {

   int imol = -1; // currently unset
   int imol_in_hope = molecules.size();

   std::string name_in = file_name + std::string(" ") + std::string(f) + std::string(" ") + std::string(phi);
   coot::molecule_t m(name_in, imol_in_hope);
   bool status = coot::util::map_fill_from_mtz(&m.xmap, file_name, f, phi, weight, use_weight, is_a_difference_map);
   if (is_a_difference_map)
      m.set_map_is_difference_map(true);
   if (status) {
      molecules.push_back(m);
      std::cout << "DEBUG:: in read_mtz() imol map: " << imol << " diff-map-status: " << is_a_difference_map << std::endl;
      imol = imol_in_hope;
   }
   return imol;
}


coot::validation_information_t
molecules_container_t::density_fit_analysis(int imol_model, int imol_map) {

   coot::validation_information_t r;
   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map)) {
         // fill these
         mmdb::PResidue *SelResidues = 0;
         int nSelResidues = 0;

         auto atom_sel = molecules[imol_model].atom_sel;
         int selHnd = atom_sel.mol->NewSelection(); // yes, it's deleted.
         int imod = 1; // multiple models don't work on validation graphs

         atom_sel.mol->Select(selHnd, mmdb::STYPE_RESIDUE, imod,
                              "*", // chain_id
                              mmdb::ANY_RES, "*",
                              mmdb::ANY_RES, "*",
                              "*",  // residue name
                              "*",  // Residue must contain this atom name?
                              "*",  // Residue must contain this Element?
                              "*",  // altLocs
                              mmdb::SKEY_NEW // selection key
                              );
         atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

         for (int ir=0; ir<nSelResidues; ir++) {
            mmdb::Residue *residue_p = SelResidues[ir];
            coot::residue_spec_t res_spec(residue_p);
            mmdb::PAtom *residue_atoms=0;
            int n_residue_atoms;
            residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
            double residue_density_score =
               coot::util::map_score(residue_atoms, n_residue_atoms, molecules[imol_map].xmap, 1);
            std::string l = res_spec.label();
            std::string atom_name = coot::util::intelligent_this_residue_mmdb_atom(residue_p)->GetAtomName();
            const std::string &chain_id = res_spec.chain_id;
            int this_resno = res_spec.res_no;
            coot::atom_spec_t atom_spec(chain_id, this_resno, res_spec.ins_code, atom_name, "");
            coot::residue_validation_information_t rvi(res_spec, atom_spec, residue_density_score, l);
            r.add_residue_valiation_informtion(rvi, chain_id);
         }
         atom_sel.mol->DeleteSelection(selHnd);
      }
   }
   return r;
}

#include "vertex.hh" // neeeded?

coot::simple_mesh_t
molecules_container_t::test_origin_cube() const {

   coot::simple_mesh_t mesh;

   std::vector<coot::api::vnc_vertex> vertices;
   std::vector<g_triangle> triangles;

   glm::vec4 c(0.5, 0.2, 0.5, 1.0); // colour

   // bottom
   coot::api::vnc_vertex v0(glm::vec3(0, 0, 0), glm::vec3(0,0,-1), c); vertices.push_back(v0);
   coot::api::vnc_vertex v1(glm::vec3(1, 0, 0), glm::vec3(0,0,-1), c); vertices.push_back(v1);
   coot::api::vnc_vertex v2(glm::vec3(0, 1, 0), glm::vec3(0,0,-1), c); vertices.push_back(v2);
   coot::api::vnc_vertex v3(glm::vec3(1, 1, 0), glm::vec3(0,0,-1), c); vertices.push_back(v3);

   // top
   coot::api::vnc_vertex v4(glm::vec3(0, 0, 1), glm::vec3(0,0,1), c); vertices.push_back(v4);
   coot::api::vnc_vertex v5(glm::vec3(1, 0, 1), glm::vec3(0,0,1), c); vertices.push_back(v5);
   coot::api::vnc_vertex v6(glm::vec3(0, 1, 1), glm::vec3(0,0,1), c); vertices.push_back(v6);
   coot::api::vnc_vertex v7(glm::vec3(1, 1, 1), glm::vec3(0,0,1), c); vertices.push_back(v7);

   // left
   coot::api::vnc_vertex v8 (glm::vec3(0, 0, 0), glm::vec3(-1,0,0), c); vertices.push_back(v8);
   coot::api::vnc_vertex v9 (glm::vec3(0, 1, 0), glm::vec3(-1,0,0), c); vertices.push_back(v9);
   coot::api::vnc_vertex v10(glm::vec3(0, 0, 1), glm::vec3(-1,0,0), c); vertices.push_back(v10);
   coot::api::vnc_vertex v11(glm::vec3(0, 1, 1), glm::vec3(-1,0,0), c); vertices.push_back(v11);

   // right
   coot::api::vnc_vertex v12(glm::vec3(1, 0, 0), glm::vec3(1,0,0), c); vertices.push_back(v12);
   coot::api::vnc_vertex v13(glm::vec3(1, 1, 0), glm::vec3(1,0,0), c); vertices.push_back(v13);
   coot::api::vnc_vertex v14(glm::vec3(1, 0, 1), glm::vec3(1,0,0), c); vertices.push_back(v14);
   coot::api::vnc_vertex v15(glm::vec3(1, 1, 1), glm::vec3(1,0,0), c); vertices.push_back(v15);

   // front
   coot::api::vnc_vertex v16(glm::vec3(0, 0, 0), glm::vec3(0,-1,0), c); vertices.push_back(v16);
   coot::api::vnc_vertex v17(glm::vec3(1, 0, 0), glm::vec3(0,-1,0), c); vertices.push_back(v17);
   coot::api::vnc_vertex v18(glm::vec3(0, 0, 1), glm::vec3(0,-1,0), c); vertices.push_back(v18);
   coot::api::vnc_vertex v19(glm::vec3(1, 0, 1), glm::vec3(0,-1,0), c); vertices.push_back(v19);

   // back
   coot::api::vnc_vertex v20(glm::vec3(0, 1, 0), glm::vec3(0,1,0), c); vertices.push_back(v20);
   coot::api::vnc_vertex v21(glm::vec3(1, 1, 0), glm::vec3(0,1,0), c); vertices.push_back(v21);
   coot::api::vnc_vertex v22(glm::vec3(0, 1, 1), glm::vec3(0,1,0), c); vertices.push_back(v22);
   coot::api::vnc_vertex v23(glm::vec3(1, 1, 1), glm::vec3(0,1,0), c); vertices.push_back(v23);

   triangles.push_back(g_triangle( 0, 1, 2));
   triangles.push_back(g_triangle( 1, 3, 2));
   triangles.push_back(g_triangle( 4, 5, 6));
   triangles.push_back(g_triangle( 5, 7, 6));
   triangles.push_back(g_triangle( 8, 9,10));
   triangles.push_back(g_triangle( 9,11,10));
   triangles.push_back(g_triangle(12,13,14));
   triangles.push_back(g_triangle(13,15,14));
   triangles.push_back(g_triangle(16,17,18));
   triangles.push_back(g_triangle(17,19,18));
   triangles.push_back(g_triangle(20,21,22));
   triangles.push_back(g_triangle(21,23,22));

   for (auto &vertex : vertices) {
      vertex.pos *= 10.0;
      vertex.pos += glm::vec3(-5.0, -5.0, -5.0);
   }

   coot::simple_mesh_t m(vertices, triangles);
   // m.translate(glm::vec3(-0.5, -0.5, -0.5));
   return m;
}

std::vector<std::pair<coot::Cartesian, coot::util::phi_psi_t> >
molecules_container_t::ramachandran_validation(int imol) const {

   std::vector<std::pair<coot::Cartesian, coot::util::phi_psi_t> > v;
   if (is_valid_model_molecule(imol))
      v = molecules[imol].ramachandran_validation();

   return v;
}


coot::simple_mesh_t
molecules_container_t::ramachandran_validation_markup_mesh(int imol) const {

   // this function should be pushed into the coot::molecule_t class
   // (which means that the mesh will be copied)

   unsigned int num_subdivisions = 2;  // pass this
   float rama_ball_radius = 0.5;

   auto prob_raw_to_colour_rotation = [] (float prob) {
                                         if (prob > 0.5) prob = 0.5; // 0.4 and 2.5 f(for q) might be better (not tested)
                                         // good probabilities have q = 0
                                         // bad probabilities have q 0.66
                                         double q = (1.0 - 2.0 * prob);
                                         q = pow(q, 20.0);
                                         return q;
   };

   auto cartesian_to_glm = [] (const coot::Cartesian &c) {
                              return glm::vec3(c.x(), c.y(), c.z());
   };

   auto phi_psi_probability = [] (const coot::util::phi_psi_t &phi_psi, const ramachandrans_container_t &rc) {

      const clipper::Ramachandran *rama = &rc.rama;

      if (phi_psi.residue_name() == "PRO") rama = &rc.rama_pro;
      if (phi_psi.residue_name() == "GLY") rama = &rc.rama_gly;

      // if (phi_psi.residue_name() == "ILE" || phi_psi.residue_name() == "VAL" ) rama = &rc.rama_ileval;
      // if (phi_psi.is_pre_pro())
      // if (phi_psi.residue_name() != "GLY")
      // rama = &rc.rama_pre_pro;

      double rama_prob = rama->probability(clipper::Util::d2rad(phi_psi.phi()),
                                           clipper::Util::d2rad(phi_psi.psi()));
      return rama_prob;
   };

   auto test_ramachandran_probabilities = [] (const ramachandrans_container_t &rc) {

      std::vector<const clipper::Ramachandran *> ramas = { &rc.rama, &rc.rama_gly, &rc.rama_pro, &rc.rama_non_gly_pro };

      for (unsigned int ir=0; ir<ramas.size(); ir++) {
         for (unsigned int i=0; i<10; i++) {
            for (unsigned int j=0; j<10; j++) {
               double phi = static_cast<double>(i * 36.0) - 180.0;
               double psi = static_cast<double>(j * 36.0) - 180.0;
               double p = rc.rama.probability(phi, psi);
               std::cout << ir << "   "
                         << std::setw(10) << phi << " " << std::setw(10) << psi << " "
                         << std::setw(10) << p << std::endl;
            }
         }
      }
   };

   // test_ramachandran_probabilities(ramachandrans_container); // don't use rama_pre_pro without CLIPPER_HAS_TOP8000

   coot::simple_mesh_t mesh;
   if (is_valid_model_molecule(imol)) {

      std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > octaball = tessellate_octasphere(num_subdivisions);

      std::vector<std::pair<coot::Cartesian, coot::util::phi_psi_t> > ramachandran_goodness_spots =
         ramachandran_validation(imol);
      // now convert positions into meshes of balls
      int n_ramachandran_goodness_spots = ramachandran_goodness_spots.size();
      for (int i=0; i<n_ramachandran_goodness_spots; i++) {
         const coot::Cartesian &position = ramachandran_goodness_spots[i].first;
         const coot::util::phi_psi_t &phi_psi = ramachandran_goodness_spots[i].second;
         const float &prob_raw = phi_psi_probability(phi_psi, ramachandrans_container);
         double q = prob_raw_to_colour_rotation(prob_raw);
         coot::colour_holder col = coot::colour_holder(q, 0.0, 1.0, false, std::string(""));
         glm::vec3 ball_position = cartesian_to_glm(position);
         unsigned int idx_base = mesh.vertices.size();
         unsigned int idx_tri_base = mesh.triangles.size();
         for (unsigned int ibv=0; ibv<octaball.first.size(); ibv++) {
            glm::vec4 col_v4(col.red, col.green, col.blue, 1.0f);
            const glm::vec3 &vertex_position = octaball.first[ibv];
            coot::api::vnc_vertex vertex(ball_position + rama_ball_radius * vertex_position, vertex_position, col_v4);
            mesh.vertices.push_back(vertex);
         }
         std::vector<g_triangle> octaball_triangles = octaball.second;
         mesh.triangles.insert(mesh.triangles.end(), octaball_triangles.begin(), octaball_triangles.end());

         for (unsigned int k=idx_tri_base; k<mesh.triangles.size(); k++)
            mesh.triangles[k].rebase(idx_base);
      }
   }
   return mesh;
}


mmdb::Atom *
molecules_container_t::get_atom(int imol, const coot::atom_spec_t &atom_spec) const {

   mmdb::Atom *r = nullptr;
   if (is_valid_model_molecule(imol)) {
      return molecules[imol].get_atom(atom_spec);
   }
   return r;
}

mmdb::Residue *
molecules_container_t::get_residue(int imol, const coot::residue_spec_t &residue_spec) const {

   mmdb::Residue *r = nullptr;
   if (is_valid_model_molecule(imol)) {
      return molecules[imol].get_residue(residue_spec);
   }
   return r;
}


int
molecules_container_t::writeMap(int imol, const std::string &file_name) const {

   int status= 0;
   if (is_valid_map_molecule(imol)) {
      status = molecules[imol].writeMap(file_name);
   }
   return status;

}


coot::simple_mesh_t
molecules_container_t::get_map_contours_mesh(int imol, double position_x, double position_y, double position_z,
                                             float radius, float contour_level) {

   coot::simple_mesh_t mesh;
   try {
      std::cout << ".................. here A in get_map_contours_mesh() " << imol << std::endl;
      if (is_valid_map_molecule(imol)) {
         std::cout << ".................. here B in get_map_contours_mesh() " << imol << std::endl;
         clipper::Coord_orth position(position_x, position_y, position_z);
         mesh = molecules[imol].get_map_contours_mesh(position, radius, contour_level);
      }
   }
   catch (...) {
      std::cout << "An error occured in " << __FUNCTION__<< "() - this should not happen " << std::endl;
   }
   return mesh;
}


// get the rotamer dodecs for the model
coot::simple_mesh_t
molecules_container_t::get_rotamer_dodecs(int imol) {
   coot::simple_mesh_t m;
   if (is_valid_model_molecule(imol))
      return molecules[imol].get_rotamer_dodecs(&geom, &rot_prob_tables);
   return m;
}


int
molecules_container_t::auto_fit_rotamer(int imol,
                                        const std::string &chain_id, int res_no, const std::string &ins_code,
                                        const std::string &alt_conf,
                                        int imol_map) {
   int status = 0;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_map)) {
         const clipper::Xmap<float> &xmap = molecules[imol_map].xmap;
         std::cout << "debug:: mc::auto_fit_rotamer() calling the coot_molecule version with "
                   << chain_id << " " << res_no << std::endl;
         status = molecules[imol].auto_fit_rotamer(chain_id, res_no, ins_code, alt_conf, xmap, geom);
      } else {
         std::cout << "debug:: mc::auto_fit_rotamer() not a valid map index " << imol_map << std::endl;
      }
   } else {
      std::cout << "debug:: mc::auto_fit_rotamer() not a valid model molecule " << imol << std::endl;
   }
   return status;
}


int
molecules_container_t::delete_atom(int imol,
                                   const std::string &chain_id, int res_no, const std::string &ins_code,
                                   const std::string &atom_name, const std::string &alt_conf) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec(chain_id, res_no, ins_code, atom_name, alt_conf);
      status = molecules[imol].delete_atom(atom_spec);
   }
   return status;
}

int
molecules_container_t::delete_atom_using_cid(int imol, const std::string &cid) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec = atom_cid_to_atom_spec(imol, cid);
      status = molecules[imol].delete_atom(atom_spec);
   }
   return status;
}



int
molecules_container_t::delete_residue(int imol,
                                      const std::string &chain_id, int res_no, const std::string &ins_code) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t residue_spec(chain_id, res_no, ins_code);
      status = molecules[imol].delete_residue(residue_spec);
   }
   return status;
}


int
molecules_container_t::delete_residue_using_cid(int imol, const std::string &cid) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t residue_spec = residue_cid_to_residue_spec(imol, cid);
      status = molecules[imol].delete_residue(residue_spec);
   }
   return status;
}



int
molecules_container_t::load_dictionary_file(const std::string &monomer_cif_file_name) {

   int status = 0;

   int read_number = 44;
   geom.init_refmac_mon_lib(monomer_cif_file_name, read_number);
   return status;
}

std::vector<std::string>
molecules_container_t::non_standard_residue_types_in_model(int imol) const {
   std::vector<std::string> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].non_standard_residue_types_in_model();
   }
   return v;
}

float
molecules_container_t::get_map_rmsd_approx(int imol) const {
   float rmsd = -1.1;
   if (is_valid_map_molecule(imol)) {
      rmsd = molecules[imol].get_map_rmsd_approx();
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return rmsd;
}


std::vector<coot::molecule_t::difference_map_peaks_info_t>
molecules_container_t::difference_map_peaks(int imol_map, int imol_protein, float n_rmsd) const {

   std::vector<coot::molecule_t::difference_map_peaks_info_t> v;
   if (is_valid_model_molecule(imol_protein)) {
      if (is_valid_map_molecule(imol_map)) {
         mmdb::Manager *m = mol(imol_protein);
         molecules[imol_map].difference_map_peaks(m, n_rmsd);
      } else {
         std::cout << "debug:: " << __FUNCTION__ << "(): not a valid map molecule " << imol_map << std::endl;
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol_protein << std::endl;
   }
   return v;
}



// return a useful message if the addition did not work
std::pair<int, std::string>
molecules_container_t::add_terminal_residue(int imol, const std::string &chain_id, int res_no, const std::string &ins_code) {

   int status = 0;
   std::string message;
   if (is_valid_model_molecule(imol)) {
      // status = molecules[imol].add_terminal_residue(chain_id, res_no, ins_code);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }

   return std::make_pair(status, message);

}

void
molecules_container_t::associate_data_mtz_file_with_map(int imol_map, const std::string &data_mtz_file_name,
                                                        const std::string &f_col, const std::string &sigf_col,
                                                        const std::string &free_r_col) {
   if (is_valid_map_molecule(imol_map)) {
      // 20221018-PE if free_r_col is not valid then Coot will (currently) crash on the structure factor calculation
      molecules[imol_map].associate_data_mtz_file_with_map(data_mtz_file_name, f_col, sigf_col, free_r_col);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid map molecule " << imol_map << std::endl;
   }
}

/*! \brief Calculate structure factors from the model and update the given difference
           map accordingly */

// copied from:
// void
// graphics_info_t::sfcalc_genmap(int imol_model,
//                                int imol_map_with_data_attached,
//                                int imol_updating_difference_map) {
void
molecules_container_t::sfcalc_genmap(int imol_model,
                                     int imol_map_with_data_attached,
                                     int imol_updating_difference_map) {

   // I am keen for this function to be fast - so that it can be used with cryo-EM structures
   //
   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map_with_data_attached)) {
         if (true) {
            if (is_valid_map_molecule(imol_updating_difference_map)) {
               if (molecules[imol_updating_difference_map].is_difference_map_p()) {
                  clipper::Xmap<float> *xmap_p = &molecules[imol_updating_difference_map].xmap;
                  try {
                     if (! on_going_updating_map_lock) {
                        on_going_updating_map_lock = true;
                        molecules[imol_map_with_data_attached].fill_fobs_sigfobs();
                        const clipper::HKL_data<clipper::data32::F_sigF> *fobs_data =
                           molecules[imol_map_with_data_attached].get_original_fobs_sigfobs();
                        const clipper::HKL_data<clipper::data32::Flag> *free_flag =
                           molecules[imol_map_with_data_attached].get_original_rfree_flags();
                        if (fobs_data && free_flag) {
                           molecules[imol_model].sfcalc_genmap(*fobs_data, *free_flag, xmap_p);
                        } else {
                           std::cout << "sfcalc_genmap() either fobs_data or free_flag were not set " << std::endl;
                        }
                        on_going_updating_map_lock = false;
                     } else {
                        std::cout << "DEBUG:: on_going_updating_map_lock was set! - aborting map update." << std::endl;
                     }
                  }
                  catch (const std::runtime_error &rte) {
                     std::cout << rte.what() << std::endl;
                  }
               } else {
                  std::cout << "sfcalc_genmap() not a valid difference map " << imol_updating_difference_map << std::endl;
               }
            } else {
               std::cout << "sfcalc_genmap() not a valid map (diff) " << imol_updating_difference_map << std::endl;
            }
         }
      } else {
         std::cout << "sfcalc_genmap() not a valid map " << imol_map_with_data_attached << std::endl;
      }
   } else {
      std::cout << "sfcalc_genmap() not a valid model " << imol_model << std::endl;
   }
}


coot::util::sfcalc_genmap_stats_t
molecules_container_t::sfcalc_genmaps_using_bulk_solvent(int imol_model,
                                                         int imol_map_2fofc,  // this map should have the data attached.
                                                         int imol_map_fofc) {
   coot::util::sfcalc_genmap_stats_t stats;
   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map_2fofc)) {
         if (is_valid_map_molecule(imol_map_fofc)) {
            if (molecules[imol_map_fofc].is_difference_map_p()) {
               try {
                  if (! on_going_updating_map_lock) {
                     on_going_updating_map_lock = true;
                     molecules[imol_map_2fofc].fill_fobs_sigfobs();

                     // 20210815-PE used to be const reference (get_original_fobs_sigfobs() function changed too)
                     // const clipper::HKL_data<clipper::data32::F_sigF> &fobs_data = molecules[imol_map_with_data_attached].get_original_fobs_sigfobs();
                     // const clipper::HKL_data<clipper::data32::Flag> &free_flag = molecules[imol_map_with_data_attached].get_original_rfree_flags();
                     // now the full object (40us for RNAse test).
                     // 20210815-PE OK, the const reference was not the problem. But we will leave it as it is now, for now.
                     //
                     clipper::HKL_data<clipper::data32::F_sigF> *fobs_data_p = molecules[imol_map_2fofc].get_original_fobs_sigfobs();
                     clipper::HKL_data<clipper::data32::Flag>   *free_flag_p = molecules[imol_map_2fofc].get_original_rfree_flags();

                     if (fobs_data_p && free_flag_p) {

                        if (true) { // sanity check data

                           const clipper::HKL_info &hkls_check = fobs_data_p->base_hkl_info();
                           const clipper::Spacegroup &spgr_check = hkls_check.spacegroup();
                           const clipper::Cell &cell_check = fobs_data_p->base_cell();
                           const clipper::HKL_sampling &sampling_check = fobs_data_p->hkl_sampling();

                           std::cout << "DEBUG:: in sfcalc_genmaps_using_bulk_solvent() imol_map_with_data_attached "
                                     << imol_map_2fofc << std::endl;

                           std::cout << "DEBUG:: Sanity check in graphics_info_t:sfcalc_genmaps_using_bulk_solvent(): HKL_info: "
                                     << "base_cell: " << cell_check.format() << " "
                                     << "spacegroup: " << spgr_check.symbol_xhm() << " "
                                     << "sampling is null: " << sampling_check.is_null() << " "
                                     << "resolution: " << hkls_check.resolution().limit() << " "
                                     << "invsqreslim: " << hkls_check.resolution().invresolsq_limit() << " "
                                     << "num_reflections: " << hkls_check.num_reflections()
                                     << std::endl;
                        }

                        clipper::Xmap<float> &xmap_2fofc = molecules[imol_map_2fofc].xmap;
                        clipper::Xmap<float> &xmap_fofc  = molecules[imol_map_fofc].xmap;
                        stats = molecules[imol_model].sfcalc_genmaps_using_bulk_solvent(*fobs_data_p, *free_flag_p, &xmap_2fofc, &xmap_fofc);

                     } else {
                        std::cout << "ERROR:: null data pointer in graphics_info_t::sfcalc_genmaps_using_bulk_solvent() " << std::endl;
                     }
                     on_going_updating_map_lock = false;
                  }
               }
               catch (const std::runtime_error &rte) {
                  std::cout << rte.what() << std::endl;
               }
            }
         }
      }
   }
   return stats;
}

int
molecules_container_t::gru_points_total() const { // the sum of all the gru ponts accumulated
   return gru_points_t::total(gru_point_history);
}

int
molecules_container_t::calculate_new_gru_points(int imol_diff_map) {

   float rmsd = get_map_rmsd_approx(imol_diff_map);
   if (! gru_point_history.empty()) {
      const gru_points_t &prev = gru_point_history.back();
      gru_points_t new_points(rmsd, prev);
      gru_point_history.push_back(new_points);
      return new_points.map_gru_points_delta;
   } else {
      gru_points_t prev = gru_points_t(rmsd);
      gru_points_t new_points(rmsd, prev);
      gru_point_history.push_back(new_points);
      return new_points.map_gru_points_delta;
   }
}
