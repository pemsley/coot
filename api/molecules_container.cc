
#include <iomanip>

#include "molecules_container.hh"
#include "ideal/pepflip.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-map-utils.hh"

#include "coords/Bond_lines.h"
#include "oct.hh"

// statics
std::atomic<bool> molecules_container_t::restraints_lock(false);
std::atomic<bool> molecules_container_t::on_going_updating_map_lock(false);
ctpl::thread_pool molecules_container_t::static_thread_pool(4); // or so
std::string molecules_container_t::restraints_locking_function_name; // I don't know why this needs to be static
std::vector<atom_pull_info_t> molecules_container_t::atom_pulls;
// 20221018-PE not sure that this needs to be static.
clipper::Xmap<float> *molecules_container_t::dummy_xmap = new clipper::Xmap<float>;

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

int
molecules_container_t::close_molecule(int imol) {

   int status = 0;
   int ms = molecules.size(); // type change
   if (imol < ms) {
      if (imol >= 0) {
         molecules[imol].close_yourself();
         status = 1;
      }
   }
   return status;
}

std::string
molecules_container_t::get_molecule_name(int imol) const {

   int ms = molecules.size();
   if (imol < ms)
      if (imol >= 0)
         return molecules[imol].get_name();
   return std::string("");

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
molecules_container_t::flip_peptide(int imol, const coot::atom_spec_t &as, const std::string &alt_conf) {

   int result = 0;
   if (is_valid_model_molecule(imol)) {
      result = molecules[imol].flip_peptide(as, alt_conf);
   }
   return result;
}

int
molecules_container_t::flip_peptide_using_cid(int imol, const std::string &atom_cid, const std::string &alt_conf) {

   int result = 0;
   if (is_valid_model_molecule(imol)) {
      auto &m = molecules[imol];
      std::pair<bool, coot::atom_spec_t> as = m.cid_to_atom_spec(atom_cid);
      if (as.first) {
         const auto &atom_spec = as.second;
         result = molecules[imol].flip_peptide(atom_spec, alt_conf); // N check in here
      }
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
molecules_container_t::get_monomer_from_dictionary(const std::string &comp_id,
                                                   bool idealised_flag) {

   int istat = -1; // unfound molecule

   int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
   mmdb::Manager *mol = geom.mol_from_dictionary(comp_id, imol_enc, idealised_flag);
   if (mol) {
      int imol = molecules.size();
      std::string name = comp_id;
      name += "_from_dict";
      // graphics_info_t::molecules[imol].install_model(imol, asc, g.Geom_p(), name, 1);
      // move_molecule_to_screen_centre_internal(imol);
      atom_selection_container_t asc = make_asc(mol);
      coot::molecule_t m = coot::molecule_t(asc, imol, name);
      molecules.push_back(m);
      istat = imol;
   } else {
      std::cout << "WARNING:: Null mol from mol_from_dictionary() with comp_id " << comp_id << " "
		<< idealised_flag << std::endl;
   }
   return istat;
}


int
molecules_container_t::get_monomer(const std::string &comp_id) {

   int imol = get_monomer_from_dictionary(comp_id, 1); // idealized
   return imol;
}

// 20221030-PE nice to have one day
// int
// molecules_container_t::get_monomer_molecule_by_network(const std::string &text) {

//    int imol = -1;
//    return imol;
// }


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
   bool status = coot::util::map_fill_from_mtz(&m.xmap, file_name, f, phi, weight, use_weight, map_sampling_rate);
   if (is_a_difference_map)
      m.set_map_is_difference_map(true);
   if (status) {
      molecules.push_back(m);
      imol = imol_in_hope;
      std::cout << "DEBUG:: in read_mtz() " << file_name << " imol map: " << imol
                << " diff-map-status: " << is_a_difference_map << std::endl;
   }
   return imol;
}

#include "clipper-ccp4-map-file-wrapper.hh"
#include "coot-utils/slurp-map.hh"

int
molecules_container_t::read_ccp4_map(const std::string &file_name, bool is_a_difference_map) {

   int imol = -1; // currently unset
   int imol_in_hope = molecules.size();
   bool done = false;

   if (coot::util::is_basic_em_map_file(file_name)) {

      // fill xmap
      bool check_only = false;
      coot::molecule_t m(file_name, imol_in_hope);
      clipper::Xmap<float> &xmap = m.xmap;
      done = coot::util::slurp_fill_xmap_from_map_file(file_name, &xmap, check_only);
      try {
         // what's the point of this now?
         clipper_map_file_wrapper file;
         file.open_read(file_name);
         // set_is_em_map(file); // sets is_em_map_cached_flag
         // em = is_em_map_cached_flag;
      }
      catch (const clipper::Message_base &exc) {
      std::cout << "WARNING:: failed to open " << file_name << std::endl;
      // bad_read = true;
      }

      if (done)
         molecules.push_back(m);
   }

   if (! done) {
         std::cout << "INFO:: attempting to read CCP4 map: " << file_name << std::endl;
         // clipper::CCP4MAPfile file;
         clipper_map_file_wrapper w_file;
         try {
            w_file.open_read(file_name);

            // em = set_is_em_map(file);

            bool use_xmap = true; // not an nxmap
            if (true) {
               clipper::Grid_sampling fgs = w_file.grid_sampling();
               clipper::Cell fcell = w_file.cell();
               double vol = fcell.volume();
               if (vol < 1.0) {
                  std::cout << "WARNING:: non-sane unit cell volume " << vol << " - skip read"
                            << std::endl;
                  // bad_read = true;
               } else {
                  try {
                     clipper::CCP4MAPfile file;
                     file.open_read(file_name);
                     clipper::Xmap<float> xmap;
                     file.import_xmap(xmap);
                  }
                  catch (const clipper::Message_generic &exc) {
                     std::cout << "WARNING:: failed to read " << file_name
                               << " Bad ASU (inconsistant gridding?)." << std::endl;
                     // bad_read = true;
                  }
               }
            }
         } catch (const clipper::Message_base &exc) {
            std::cout << "WARNING:: failed to open " << file_name << std::endl;
            // bad_read = true;
         }
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

// returns either the specified atom or null if not found
mmdb::Atom *
molecules_container_t::get_atom_using_cid(int imol, const std::string &cid) const {

   mmdb::Atom *at = nullptr;
   if (is_valid_model_molecule(imol)) {
      std::pair<bool, coot::atom_spec_t> p = molecules[imol].cid_to_atom_spec(cid);
      if (p.first)
         at = molecules[imol].get_atom(p.second);
   }
   return at;
}

// returns either the specified residue or null if not found
mmdb::Residue *
molecules_container_t::get_residue_using_cid(int imol, const std::string &cid) const {
   mmdb::Residue *residue_p = nullptr;
   if (is_valid_model_molecule(imol)) {
      std::pair<bool, coot::residue_spec_t> p = molecules[imol].cid_to_residue_spec(cid);
      if (p.first)
         residue_p = molecules[imol].get_residue(p.second);
   }
   return residue_p;
}


int
molecules_container_t::move_molecule_to_new_centre(int imol, float x, float y, float z) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::Cartesian new_centre(x,y,z);
      status = molecules[imol].move_molecule_to_new_centre(new_centre);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}

coot::Cartesian
molecules_container_t::get_molecule_centre(int imol) const {

   coot::Cartesian c;
   if (is_valid_model_molecule(imol)) {
      c = molecules[imol].get_molecule_centre();
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return c;
}


int
molecules_container_t::writeMap(int imol, const std::string &file_name) const {

   int status= 0;
   if (is_valid_map_molecule(imol)) {
      status = molecules[imol].writeMap(file_name);
   }
   return status;

}

// Mode is "COLOUR-BY-CHAIN-AND-DICTIONARY" or "CA+LIGANDS"
coot::simple_mesh_t
molecules_container_t::get_bonds_mesh(int imol, const std::string &mode) {

   auto tp_0 = std::chrono::high_resolution_clock::now();

   coot::simple_mesh_t sm;
   if (is_valid_model_molecule(imol)) {
      sm = molecules[imol].get_bonds_mesh(mode, &geom);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   auto tp_1 = std::chrono::high_resolution_clock::now();
   if (show_timings) {
      auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
      std::cout << "---------- timings: for get_bonds_mesh(): : " << d10 << " milliseconds " << std::endl;
   }

   return sm;
}


coot::simple_mesh_t
molecules_container_t::get_map_contours_mesh(int imol, double position_x, double position_y, double position_z,
                                             float radius, float contour_level) {

   auto tp_0 = std::chrono::high_resolution_clock::now();
   coot::simple_mesh_t mesh;
   try {
      if (is_valid_map_molecule(imol)) {
         clipper::Coord_orth position(position_x, position_y, position_z);
         mesh = molecules[imol].get_map_contours_mesh(position, radius, contour_level);
      }
   }
   catch (...) {
      std::cout << "An error occured in " << __FUNCTION__<< "() - this should not happen " << std::endl;
   }
   auto tp_1 = std::chrono::high_resolution_clock::now();
   if (show_timings) {
      auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
      std::cout << "---------- timings: for get_map_contours_mesh(): : " << d10 << " milliseconds " << std::endl;
   }
   return mesh;
}


// get the rotamer dodecs for the model
coot::simple_mesh_t
molecules_container_t::get_rotamer_dodecs(int imol) {
   coot::simple_mesh_t m;
   if (is_valid_model_molecule(imol)) {
      return molecules[imol].get_rotamer_dodecs(&geom, &rot_prob_tables);
   } else {
      std::cout << "WARNING:: in " << __FUNCTION__ << "() imol " << imol << " was not a valid model molecule " << std::endl;
   }
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
molecules_container_t::delete_residue_using_cid(int imol, const std::string &residue_cid) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t residue_spec = residue_cid_to_residue_spec(imol, residue_cid);
      status = molecules[imol].delete_residue(residue_spec);
   }
   return status;
}

int
molecules_container_t::delete_residue_atoms_using_cid(int imol, const std::string &atom_cid) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec = atom_cid_to_atom_spec(imol, atom_cid);
      coot::residue_spec_t residue_spec(atom_spec);
      status = molecules[imol].delete_residue(residue_spec);
   }
   return status;
}

int
molecules_container_t::delete_residue_atoms_with_alt_conf(int imol, const std::string &chain_id,
                                                          int res_no, const std::string &ins_code,
                                                          const std::string &alt_conf) {
  int status = 0;
  return status;
}



int
molecules_container_t::delete_chain_using_cid(int imol, const std::string &cid) {

   int status = 0;
   if (is_valid_model_molecule(imol))
      status = molecules[imol].delete_chain_using_atom_cid(cid);
   return status;
}



//where scope in ["ATOM","WATER", "RESIDUE","CHAIN","MOLECULE"]
int
molecules_container_t::delete_using_cid(int imol, const std::string &cid, const std::string &scope) {

   int status = 0;
   if (scope == "ATOM")
      status = delete_atom_using_cid(imol, cid);
   if (scope == "RESIDUE")
      status = delete_residue_atoms_using_cid(imol, cid);
   if (scope == "CHAIN")
      status = delete_chain_using_cid(imol, cid);
   if (scope == "MOLECULE")
      status = close_molecule(imol);
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
         mmdb::Manager *m = get_mol(imol_protein);
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
molecules_container_t::add_terminal_residue_directly(int imol, const std::string &chain_id, int res_no, const std::string &ins_code) {

   std::string new_res_type = "ALA";
   int status = 0;
   std::string message;

   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_refinement_map)) {
         clipper::Xmap<float> &xmap = molecules[imol_refinement_map].xmap;
         coot::residue_spec_t residue_spec(chain_id, res_no, ins_code);
         std::pair<int, std::string> m = molecules[imol].add_terminal_residue_directly(residue_spec, new_res_type,
                                                                                       geom, xmap);
         status  = m.first;
         message = m.second;
      } else {
         std::cout << "debug:: " << __FUNCTION__ << "(): not a valid map molecule " << imol_refinement_map << std::endl;
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return std::make_pair(status, message);
}

// 20221023-PE return an int for now so that I can write the binding
int
molecules_container_t::add_terminal_residue_directly_using_cid(int imol, const std::string &cid) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec = atom_cid_to_atom_spec(imol, cid);
      if (! atom_spec.empty()) {
         auto p = add_terminal_residue_directly(imol, atom_spec.chain_id, atom_spec.res_no, atom_spec.ins_code);
         status = p.first;
      }
   }
   return status;
}


// reset the gru_points (calls reset_the_gru_points()), updates the maps (using internal/clipper SFC)
// so, update your contour lines meshes after calling this function.
int
molecules_container_t::connect_updating_maps(int imol_model, int imol_map_2fofc, int imol_map_fofc) {
   int status = 0;

   return status;
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


// static
void
molecules_container_t::thread_for_refinement_loop_threaded() {

   // I think that there is a race condition here
   // check_and_warn_inverted_chirals_and_cis_peptides()
   // get called several times when the refine loop ends
   // (with success?).

   bool use_graphics_interface_flag = false;
   bool refinement_immediate_replacement_flag = true;

#if 0 // 20221018-PE this might not be the right thing

   if (restraints_lock) {
      if (false)
         std::cout << "debug:: thread_for_refinement_loop_threaded() restraints locked by "
                   << restraints_locking_function_name << std::endl;
      return;
   } else {

      if (use_graphics_interface_flag) {

         if (!refinement_immediate_replacement_flag) {

            // if there's not a refinement redraw function already running start up a new one.
            if (threaded_refinement_redraw_timeout_fn_id == -1) {
               GSourceFunc cb = GSourceFunc(regenerate_intermediate_atoms_bonds_timeout_function_and_draw);
	       // int id = gtk_timeout_add(15, cb, NULL);

               int timeout_ms = 15;
               timeout_ms = 30; // 20220503-PE try this value
	       int id = g_timeout_add(timeout_ms, cb, NULL);
               threaded_refinement_redraw_timeout_fn_id = id;
            }
         }
      }

      continue_threaded_refinement_loop = true;
      std::thread r(refinement_loop_threaded);
      r.detach();
   }
#endif

}


int
molecules_container_t::refine_direct(int imol, std::vector<mmdb::Residue *> rv, const std::string &alt_loc,
                                     mmdb::Manager *mol) {

   bool make_trans_peptide_restraints = true;
   bool do_rama_plot_restraints = false;

   int status =  0;
   std::vector<coot::atom_spec_t> fixed_atom_specs;
   std::vector<std::pair<bool,mmdb::Residue *> > local_residues;
   for (const auto &r : rv)
      local_residues.push_back(std::make_pair(false, r));

   if (is_valid_map_molecule(imol_refinement_map)) {
      clipper::Xmap<float> &xmap = molecules[imol_refinement_map].xmap;
      std::vector<mmdb::Link> links;
      coot::restraints_container_t restraints(local_residues,
                                              links,
                                              geom,
                                              mol,
                                              fixed_atom_specs, &xmap);
      
      if (refinement_is_quiet)
         restraints.set_quiet_reporting();

      std::cout << "DEBUG:: using restraints with map_weight " << map_weight << std::endl;
      restraints.add_map(map_weight);
      coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
      flags = coot::TYPICAL_RESTRAINTS;
      coot::pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;

      int n_threads = 4; // coot::get_max_number_of_threads();
      ctpl::thread_pool thread_pool(n_threads);
      restraints.thread_pool(&thread_pool, n_threads);
      
      int imol = 0; // dummy
      restraints.make_restraints(imol, geom, flags, 1, make_trans_peptide_restraints,
                                 1.0, do_rama_plot_restraints, true, true, false, pseudos);
      int nsteps_max = 4000;
      short int print_chi_sq_flag = 1;
      restraints.minimize(flags, nsteps_max, print_chi_sq_flag);
      coot::geometry_distortion_info_container_t gd = restraints.geometric_distortions();
      if (! refinement_is_quiet)
         gd.print();
      
   } else {
      std::cout << "WARNING:: refinement map " << imol_refinement_map << " is not a valid map" << std::endl;
   }
   return status;
}

int
molecules_container_t::refine_residues_using_atom_cid(int imol, const std::string &cid, const std::string &mode) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_refinement_map)) {
         coot::atom_spec_t spec = atom_cid_to_atom_spec(imol, cid);
         status = refine_residues(imol, spec.chain_id, spec.res_no, spec.ins_code, spec.alt_conf, mode);
      } else {
         std::cout << "Not a valid map molecule " << imol_refinement_map << std::endl;
      }
   } else {
      std::cout << "Not a valid model molecule " << imol << std::endl;
   }
   return status;
}



int
molecules_container_t::refine_residues(int imol, const std::string &chain_id, int res_no, const std::string &ins_code,
                                       const std::string &alt_conf, const std::string &mode) {
   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t residue_spec(chain_id, res_no, ins_code);
      mmdb::Manager *mol = get_mol(imol);
      std::vector<mmdb::Residue *> rv = molecules[imol].select_residues(residue_spec, mode);
      if (! rv.empty()) {
         status = refine_direct(imol, rv, alt_conf, mol);
      } else {
         std::cout << "WARNING:: in refine_residues() - empty residues." << std::endl;
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}

int
molecules_container_t::refine_residue_range(int imol, const std::string &chain_id, int res_no_start, int res_no_end) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      std::vector<mmdb::Residue *> rv = molecules[imol].select_residues(chain_id, res_no_start, res_no_end);
      if (! rv.empty()) {
         std::string alt_conf = "";
         mmdb::Manager *mol = get_mol(imol);
         status = refine_direct(imol, rv, alt_conf, mol);
      } else {
         std::cout << "WARNING:: in refine_residues() - empty residues." << std::endl;
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}



coot::refinement_results_t
molecules_container_t::refine_residues_vec(int imol,
                                           const std::vector<mmdb::Residue *> &residues,
                                           const std::string &alt_conf,
                                           mmdb::Manager *mol) {
   bool use_map_flag = true;
   if (false)
      std::cout << "INFO:: refine_residues_vec() with altconf \""
		<< alt_conf << "\"" << std::endl;

   coot::refinement_results_t rr = generate_molecule_and_refine(imol, residues, alt_conf, mol, use_map_flag);
   return rr;
}

// return -1 on failure to find a residue for insertion index
//
int
molecules_container_t::find_serial_number_for_insert(int seqnum_new,
                                                     const std::string &ins_code_for_new,
                                                     mmdb::Chain *chain_p) const {

   int iserial_no = -1;
   if (chain_p) {
      int current_diff = 999999;
      int nres = chain_p->GetNumberOfResidues();
      for (int ires=0; ires<nres; ires++) { // ires is a serial number
	 mmdb::Residue *residue = chain_p->GetResidue(ires);

	 // we are looking for the smallest negative diff:
	 //
	 int diff = residue->GetSeqNum() - seqnum_new;
	 if ( (diff > 0) && (diff < current_diff) ) {
	    iserial_no = ires;
	    current_diff = diff;
	 } else {
	    if (diff == 0) {
	       std::string ins_code_this = residue->GetInsCode();
	       if (ins_code_this > ins_code_for_new) {
		  iserial_no = ires;
		  break;
	       }
	    }
	 }
      }
   }
   return iserial_no;
}


#include "coords/mmdb-extras.h"

std::pair<mmdb::Manager *, std::vector<mmdb::Residue *> >
molecules_container_t::create_mmdbmanager_from_res_vector(const std::vector<mmdb::Residue *> &residues,
                                                          int imol,
                                                          mmdb::Manager *mol_in,
                                                          std::string alt_conf) {

   // returned entities
   mmdb::Manager *new_mol = 0;
   std::vector<mmdb::Residue *> rv; // gets checked

   float dist_crit = 5.0;
   bool debug = false;

   if (debug) {
      std::cout << "############ starting create_mmdbmanager_from_res_vector() with these "
		<< " residues " << std::endl;
      for (std::size_t ii=0; ii<residues.size(); ii++)
	 std::cout << "   " << coot::residue_spec_t(residues[ii])  << std::endl;
      int udd_atom_index_handle = mol_in->GetUDDHandle(mmdb::UDR_ATOM, "atom index");
      std::cout << "############ udd for atom index from seeding molecule " << udd_atom_index_handle
		<< std::endl;
      for (std::size_t ii=0; ii<residues.size(); ii++) {
	 mmdb::Residue *residue_p = residues[ii];
	 mmdb::Atom **residue_atoms = 0;
	 int n_residue_atoms;
	 residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
	 for (int iat=0; iat<n_residue_atoms; iat++) {
	    mmdb::Atom *at = residue_atoms[iat];
	    int idx = -1;
	    at->GetUDData(udd_atom_index_handle, idx);
	    std::cout << "#### input residue atom " << coot::atom_spec_t(at) << " had udd index "
		      << idx << std::endl;
	 }
      }
   }

   int n_flanker = 0; // a info/debugging counter

   if (residues.size() > 0) {

      // Also add the index of the reference residue (the one in molecules[imol].atom_selection.mol)
      // to the molecule that we are construction here. So that we can properly link
      // the residues in restraints_container (there we rather need to know the references indices,
      // not the indices from the fragment molecule)
      //

      std::pair<bool,std::string> use_alt_conf(false, "");
      if (! alt_conf.empty())
	 use_alt_conf = std::pair<bool, std::string> (true, alt_conf);

      std::cout << "----------------- in create_mmdbmanager_from_res_vector() alt_conf is "
                << "\"" << alt_conf << "\"" << std::endl;
      std::cout << "----------------- in create_mmdbmanager_from_res_vector() use_alt_conf is "
                << use_alt_conf.first << "\"" << use_alt_conf.second << "\"" << std::endl;

      std::pair<bool, mmdb::Manager *> n_mol_1 =
	 coot::util::create_mmdbmanager_from_residue_vector(residues, mol_in, use_alt_conf);

      // check that first is sane, so indent all this lot (when it works)

      if (n_mol_1.first) {

	 int index_from_reference_residue_handle =
	    n_mol_1.second->GetUDDHandle(mmdb::UDR_RESIDUE, "index from reference residue");

	 if (false) { // debug
	    int imod = 1;
	    mmdb::Model *model_p = n_mol_1.second->GetModel(imod);
	    if (model_p) {
	       int n_chains = model_p->GetNumberOfChains();
	       for (int ichain=0; ichain<n_chains; ichain++) {
		  mmdb::Chain *chain_p = model_p->GetChain(ichain);
		  int nres = chain_p->GetNumberOfResidues();
		  for (int ires=0; ires<nres; ires++) {
		     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
		     int n_atoms = residue_p->GetNumberOfAtoms();
		     for (int iat=0; iat<n_atoms; iat++) {
			mmdb::Atom *at = residue_p->GetAtom(iat);
			int idx = -1;
			at->GetUDData(index_from_reference_residue_handle, idx);
			std::cout << "   create_mmdbmanager_from_residue_vector() returns this mol atom "
				  << iat << " " << coot::atom_spec_t(at) << " with idx " << idx << std::endl;
		     }
		  }
	       }
	    }
	 }

	 new_mol = n_mol_1.second;
	 mmdb::Model *model_p = new_mol->GetModel(1);

	 // how many (free) residues were added to that model? (add them to rv)
	 //
	 int n_chains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<n_chains; ichain++) {
	    mmdb::Chain *chain_p = model_p->GetChain(ichain);
	    int nres = chain_p->GetNumberOfResidues();
	    for (int ires=0; ires<nres; ires++) {
	       mmdb::Residue *residue_p = chain_p->GetResidue(ires);
	       rv.push_back(residue_p);
	    }
	 }

	 if (false) {
	    for (std::size_t ir=0; ir<rv.size(); ir++) {
	       mmdb::Residue *r = rv[ir];
	       std::cout << "Moving Residue " << coot::residue_spec_t(r) << std::endl;
	       mmdb::Atom **residue_atoms = 0;
	       int n_residue_atoms;
	       r->GetAtomTable(residue_atoms, n_residue_atoms);
	       for (int iat=0; iat<n_residue_atoms; iat++) {
		  mmdb::Atom *at = residue_atoms[iat];
		  std::cout << "    " << coot::atom_spec_t(at) << std::endl;
	       }
	    }
	 }

	 short int whole_res_flag = 0;
	 int atom_index_udd_handle = molecules[imol].atom_sel.UDDAtomIndexHandle;

	 // Now the flanking residues:
	 //
	 std::vector<mmdb::Residue *> flankers_in_reference_mol;
	 flankers_in_reference_mol.reserve(32); // say

	 // find the residues that are close to the residues of
	 // residues that are not part of residues
	 //
	 // We don't have quite the function that we need in coot-utils,
	 // so we need to munge residues in to local_residues:
	 std::vector<std::pair<bool, mmdb::Residue *> > local_residues;
	 local_residues.resize(residues.size());
	 for (std::size_t ires=0; ires<residues.size(); ires++)
	    local_residues[ires] = std::pair<bool, mmdb::Residue *>(false, residues[ires]);
	 std::map<mmdb::Residue *, std::set<mmdb::Residue *> > rnr =
	    coot::residues_near_residues(local_residues, mol_in, dist_crit);
	 // now fill @var{flankers_in_reference_mol} from rnr, avoiding residues
	 // already in @var{residues}.
	 std::map<mmdb::Residue *, std::set<mmdb::Residue *> >::const_iterator it;
	 for (it=rnr.begin(); it!=rnr.end(); ++it) {
	    const std::set<mmdb::Residue *> &s = it->second;
	    std::set<mmdb::Residue *>::const_iterator its;
	    for (its=s.begin(); its!=s.end(); ++its) {
	       mmdb::Residue *tres = *its;
	       if (std::find(residues.begin(), residues.end(), tres) == residues.end())
		  if (std::find(flankers_in_reference_mol.begin(), flankers_in_reference_mol.end(), tres) == flankers_in_reference_mol.end())
		     flankers_in_reference_mol.push_back(tres);
	    }
	 }

	 // So we have a vector of residues that were flankers in the
	 // reference molecule, we need to add copies of those to
	 // new_mol (making sure that they go into the correct chain).
	 //
	 if (false) { // debug
	    std::cout << "debug:: ############ Found " << flankers_in_reference_mol.size()
		      << " flanking residues" << std::endl;

	    for (unsigned int ires=0; ires<flankers_in_reference_mol.size(); ires++)
	       std::cout << "     #### flankers_in_reference_mol: " << ires << " "
			 << coot::residue_spec_t(flankers_in_reference_mol[ires]) << std::endl;
	 }


	 for (unsigned int ires=0; ires<flankers_in_reference_mol.size(); ires++) {
	    mmdb::Residue *r;

	    std::string ref_res_chain_id = flankers_in_reference_mol[ires]->GetChainID();

	    mmdb::Chain *chain_p = NULL;
	    int n_new_mol_chains = model_p->GetNumberOfChains();
	    for (int ich=0; ich<n_new_mol_chains; ich++) {
	       if (ref_res_chain_id == model_p->GetChain(ich)->GetChainID()) {
		  chain_p = model_p->GetChain(ich);
		  break;
	       }
	    }

	    if (! chain_p) {
	       // Add a new one then.
	       chain_p = new mmdb::Chain;
	       chain_p->SetChainID(ref_res_chain_id.c_str());
	       model_p->AddChain(chain_p);
	    }

	    if (false)
	       std::cout << "debug:: flankers_in_reference_mol " << ires << " "
			 << coot::residue_spec_t(flankers_in_reference_mol[ires]) << " "
			 << "had index " << flankers_in_reference_mol[ires]->index
			 << std::endl;

            // get rid of this function at some stage
            bool embed_in_chain = false;
	    r = coot::deep_copy_this_residue_old_style(flankers_in_reference_mol[ires],
					               alt_conf, whole_res_flag,
					               atom_index_udd_handle, embed_in_chain);

	    if (r) {

	       r->PutUDData(index_from_reference_residue_handle, flankers_in_reference_mol[ires]->index);

	       // copy over the atom indices. UDDAtomIndexHandle in mol_n becomes UDDOldAtomIndexHandle
	       // indices in the returned molecule

	       int sni = find_serial_number_for_insert(r->GetSeqNum(), r->GetInsCode(), chain_p);

	       if (false) { // debug
		  mmdb::Atom **residue_atoms = 0;
		  int n_residue_atoms;
		  std::cout << "Flanker Residue " << coot::residue_spec_t(r) << std::endl;
		  r->GetAtomTable(residue_atoms, n_residue_atoms);
		  for (int iat=0; iat<n_residue_atoms; iat++) {
		     mmdb::Atom *at = residue_atoms[iat];
		     std::cout << "    " << coot::atom_spec_t(at) << std::endl;
		  }
	       }

	       if (sni == -1)
		  chain_p->AddResidue(r); // at the end
	       else
		  chain_p->InsResidue(r, sni);
	       r->seqNum = flankers_in_reference_mol[ires]->GetSeqNum();
	       r->SetResName(flankers_in_reference_mol[ires]->GetResName());
	       n_flanker++;

	       if (false)
		  std::cout << "debug:: create_mmdbmanager_from_residue_vector() inserted/added flanker "
			    << coot::residue_spec_t(r) << std::endl;

	    }
	 }

	 // super-critical for correct peptide bonding in refinement!
	 //
	 coot::util::pdbcleanup_serial_residue_numbers(new_mol);

	 if (debug) {
	    int imod = 1;
	    mmdb::Model *model_p = new_mol->GetModel(imod);
	    if (model_p) {
	       int n_chains = model_p->GetNumberOfChains();
	       for (int ichain=0; ichain<n_chains; ichain++) {
		  mmdb::Chain *chain_p = model_p->GetChain(ichain);
		  int nres = chain_p->GetNumberOfResidues();
		  for (int ires=0; ires<nres; ires++) {
		     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
		     std::cout << "create_mmdb..  ^^^ " << coot::residue_spec_t(residue_p) << " "
			       << residue_p << " index " << residue_p->index
			       << std::endl;
		  }
	       }
	    }
	 }

	 if (debug)
	    std::cout << "DEBUG:: in create_mmdbmanager_from_res_vector: " << rv.size()
		      << " free residues and " << n_flanker << " flankers" << std::endl;
      }
   }

   return std::pair <mmdb::Manager *, std::vector<mmdb::Residue *> > (new_mol, rv);
}



std::string
molecules_container_t::adjust_refinement_residue_name(const std::string &resname) const {

   std::string r = resname;
   if (resname == "UNK") r = "ALA"; // hack for KC/buccaneer.
   if (resname.length() > 2)
      if (resname[2] == ' ')
	 r = resname.substr(0,2);
   return r;
}


// Return 0 (first) if any of the residues don't have a dictionary
// entry and a list of the residue type that don't have restraints.
//
std::pair<int, std::vector<std::string> >
molecules_container_t::check_dictionary_for_residue_restraints(int imol, mmdb::PResidue *SelResidues, int nSelResidues) {

   int status;
   bool status_OK = 1; // pass, by default
   std::vector<std::string> res_name_vec;

   for (int ires=0; ires<nSelResidues; ires++) {
      std::string resn(SelResidues[ires]->GetResName());
      std::string resname = adjust_refinement_residue_name(resn);
      status = geom.have_dictionary_for_residue_type(resname, imol, cif_dictionary_read_number);
      cif_dictionary_read_number++;
      if (! status) {
	 status_OK = 0;
	 res_name_vec.push_back(resname);
      }

      if (0)
	 std::cout << "DEBUG:: have_dictionary_for_residues() on residue "
		   << ires << " of " << nSelResidues << ", "
		   << resname << " returns "
		   << status << std::endl;
      cif_dictionary_read_number++;
   }
   return std::pair<int, std::vector<std::string> > (status_OK, res_name_vec);
}

std::pair<int, std::vector<std::string> >
molecules_container_t::check_dictionary_for_residue_restraints(int imol, const std::vector<mmdb::Residue *> &residues) {

   std::vector<std::string> res_name_vec;
   std::pair<int, std::vector<std::string> > r(0, res_name_vec);
   for (unsigned int i=0; i<residues.size(); i++) {
      std::string resname = adjust_refinement_residue_name(residues[i]->GetResName());
      int status = geom.have_dictionary_for_residue_type(resname, imol, cif_dictionary_read_number);
      if (! status) {
	 r.first = 0;
	 r.second.push_back(resname);
      }
      cif_dictionary_read_number++; // not sure why this is needed.
   }
   return r;
}


std::vector<std::pair<mmdb::Residue *, std::vector<coot::dict_torsion_restraint_t> > >
molecules_container_t::make_rotamer_torsions(const std::vector<std::pair<bool, mmdb::Residue *> > &local_residues) const {

   std::vector<std::pair<mmdb::Residue *, std::vector<coot::dict_torsion_restraint_t> > > v;
   for (unsigned int i=0; i<local_residues.size(); i++) {
      if (! local_residues[i].first) {
         mmdb::Residue *residue_p = local_residues[i].second;
         std::string rn(residue_p->GetResName());
         if (coot::util::is_standard_amino_acid_name(rn)) {
            std::string alt_conf; // run through them all, ideally.
            coot::rotamer rot(residue_p, alt_conf, 1);
            coot::closest_rotamer_info_t cri = rot.get_closest_rotamer(rn);
            if (cri.residue_chi_angles.size() > 0) {
               std::vector<coot::dict_torsion_restraint_t> dictionary_vec;
               std::vector<std::vector<std::string> > rotamer_atom_names = rot.rotamer_atoms(rn);

               if (cri.residue_chi_angles.size() != rotamer_atom_names.size()) {

                  std::cout << "-------------- mismatch for " << coot::residue_spec_t(residue_p) << " "
                            << cri.residue_chi_angles.size() << " "  << rotamer_atom_names.size()
                            << " ---------------" << std::endl;
               } else {

                  for (unsigned int ichi=0; ichi<cri.residue_chi_angles.size(); ichi++) {
                     // we have to convert chi angles to atom names
                     double esd = 3.0; // 20210315-PE was 10.0. I want them tighter than that.
                     int per = 1;
                     std::string id = "chi " + coot::util::int_to_string(cri.residue_chi_angles[ichi].first);
                     const std::string &atom_name_1 = rotamer_atom_names[ichi][0];
                     const std::string &atom_name_2 = rotamer_atom_names[ichi][1];
                     const std::string &atom_name_3 = rotamer_atom_names[ichi][2];
                     const std::string &atom_name_4 = rotamer_atom_names[ichi][3];
                     double torsion = cri.residue_chi_angles[ichi].second;
                     coot::dict_torsion_restraint_t dr(id, atom_name_1, atom_name_2, atom_name_3, atom_name_4,
                                                       torsion, esd, per);
                     dictionary_vec.push_back(dr);
                  }

                  if (dictionary_vec.size() > 0) {
                     std::pair<mmdb::Residue *, std::vector<coot::dict_torsion_restraint_t> > p(residue_p, dictionary_vec);
                     v.push_back(p);
                  }
               }
            }
         }
      }
   }
   return v;
}



atom_selection_container_t
molecules_container_t::make_moving_atoms_asc(mmdb::Manager *residues_mol,
                                             const std::vector<mmdb::Residue *> &residues) const {

   // This also rebonds the imol_moving_atoms molecule

   atom_selection_container_t local_moving_atoms_asc;
   local_moving_atoms_asc.UDDAtomIndexHandle = -1;
   local_moving_atoms_asc.UDDOldAtomIndexHandle = residues_mol->GetUDDHandle(mmdb::UDR_ATOM, "old atom index");

   int SelHnd = residues_mol->NewSelection();

   for (unsigned int ir=0; ir<residues.size(); ir++) {
      const char *chain_id = residues[ir]->GetChainID();
      const char *inscode = residues[ir]->GetInsCode();
      int resno = residues[ir]->GetSeqNum();
      residues_mol->Select(SelHnd, mmdb::STYPE_ATOM,
			   0, chain_id,
			   resno, // starting resno, an int
			   inscode, // any insertion code
			   resno, // ending resno
			   inscode, // ending insertion code
			   "*", // any residue name
			   "*", // atom name
			   "*", // elements
			   "*",  // alt loc.
			   mmdb::SKEY_OR);
   }

   local_moving_atoms_asc.mol = residues_mol;
   local_moving_atoms_asc.SelectionHandle = SelHnd;
   residues_mol->GetSelIndex(local_moving_atoms_asc.SelectionHandle,
			     local_moving_atoms_asc.atom_selection,
			     local_moving_atoms_asc.n_selected_atoms);


   if (true) {
      std::cout << "returning an atom selection for all moving atoms "
		<< local_moving_atoms_asc.n_selected_atoms << " atoms "
		<< std::endl;
   }

   // This new block added so that we don't draw atoms in the "static" molecule when we have the
   // corresponding atoms in the moving atoms.
   //
#if 0 // 20221018-PE there is no drawing at the momment
   const atom_selection_container_t &imol_asc = molecules[imol_moving_atoms].atom_sel;
   std::set<int> atom_set = coot::atom_indices_in_other_molecule(imol_asc, local_moving_atoms_asc);

   if (false) { // debug atoms in other molecule
      std::set<int>::const_iterator it;
      for(it=atom_set.begin(); it!=atom_set.end(); it++) {
	 int idx = *it;
	 mmdb::Atom *at = imol_asc.atom_selection[idx];
	 coot::atom_spec_t as(at);
	 std::cout << " this is a moving atom: " << idx << " " << as << std::endl;
      }
   }

   if (false) { // debug old atom index
      for (int i=0; i<local_moving_atoms_asc.n_selected_atoms; i++) {
	 mmdb::Atom *at = local_moving_atoms_asc.atom_selection[i];
	 coot::atom_spec_t as(at);
	 int idx = -1;
	 at->GetUDData(local_moving_atoms_asc.UDDOldAtomIndexHandle, idx);
	 std::cout << "DEBUG:: in make_moving_atoms_asc " << as << " idx " << idx << std::endl;
      }
   }
   // now rebond molecule imol without bonds to atoms in atom_set
   if (atom_set.size() > 0) {
      if (regenerate_bonds_needs_make_bonds_type_checked_flag) {
         molecules[imol_moving_atoms].make_bonds_type_checked(atom_set, __FUNCTION__);
      }
   }
#endif

   return local_moving_atoms_asc;
}

// static
void
molecules_container_t::all_atom_pulls_off() {
   for (std::size_t i=0; i<atom_pulls.size(); i++)
      atom_pulls[i].off();
   atom_pulls.clear();
}


// return the state of having found restraints.
bool
molecules_container_t::make_last_restraints(const std::vector<std::pair<bool,mmdb::Residue *> > &local_residues,
				      const std::vector<mmdb::Link> &links,
				      const coot::protein_geometry &geom,
				      mmdb::Manager *mol_for_residue_selection,
				      const std::vector<coot::atom_spec_t> &fixed_atom_specs,
				      coot::restraint_usage_Flags flags,
				      bool use_map_flag,
				      const clipper::Xmap<float> *xmap_p) {

   bool do_torsion_restraints = true; // make this a data member
   double torsion_restraints_weight = 10.0;
   bool convert_dictionary_planes_to_improper_dihedrals_flag = false;
   double geometry_vs_map_weight = 25.5;
   bool do_trans_peptide_restraints = true;
   double rama_plot_restraints_weight = 20.0;
   bool do_rama_restraints = false;
   bool make_auto_h_bond_restraints_flag = false;
   coot::pseudo_restraint_bond_type pseudo_bonds_type = coot::NO_PSEUDO_BONDS;
   bool use_harmonic_approximation_for_NBCs = false;
   double pull_restraint_neighbour_displacement_max_radius = 1.0;
   double lennard_jones_epsilon = 1.0;
   int restraints_rama_type = 1;
   bool do_rotamer_restraints = false;
   double geman_mcclure_alpha = 0.1;
   bool do_numerical_gradients =  false;
   bool draw_gl_ramachandran_plot_flag = false;
   bool use_graphics_interface_flag = false;


   if (last_restraints) {
      std::cout << "----------------------------------------------" << std::endl;
      std::cout << "----------------------------------------------" << std::endl;
      std::cout << "    ERROR:: A: last_restraints not cleared up " << std::endl;
      std::cout << "----------------------------------------------" << std::endl;
      std::cout << "----------------------------------------------" << std::endl;
   }

   if (false) { // these are the passed residues, nothing more.
      std::cout << "debug:: on construction of restraints_container_t local_residues: "
		<< std::endl;
      for (std::size_t jj=0; jj<local_residues.size(); jj++) {
	 std::cout << "   " << coot::residue_spec_t(local_residues[jj].second)
		   << " is fixed: " << local_residues[jj].first << std::endl;
      }
   }

   // moving_atoms_extra_restraints_representation.clear();
   continue_threaded_refinement_loop = true; // no longer set in refinement_loop_threaded()

   // the refinment of torsion seems a bit confused? If it's in flags, why does it need an flag
   // of its own? I suspect that it doesn't. For now I will keep it (as it was).
   //
   bool do_residue_internal_torsions = false;
   if (do_torsion_restraints) {
      do_residue_internal_torsions = 1;
   }

   last_restraints = new
      coot::restraints_container_t(local_residues,
				   links,
				   geom,
				   mol_for_residue_selection,
				   fixed_atom_specs, xmap_p);

   std::cout << "debug:: on creation last_restraints is " << last_restraints << std::endl;

   last_restraints->set_torsion_restraints_weight(torsion_restraints_weight);

   if (convert_dictionary_planes_to_improper_dihedrals_flag) {
      last_restraints->set_convert_plane_restraints_to_improper_dihedral_restraints(true);
   }

   // This seems not to work yet.
   // last_restraints->set_dist_crit_for_bonded_pairs(9.0);

   if (use_map_flag)
      last_restraints->add_map(geometry_vs_map_weight);

   unsigned int n_threads = coot::get_max_number_of_threads();
   if (n_threads > 0)
      last_restraints->thread_pool(&static_thread_pool, n_threads);

   all_atom_pulls_off();
   particles_have_been_shown_already_for_this_round_flag = false;

   // elsewhere do this:
   // gtk_widget_remove_tick_callback(glareas[0], wait_for_hooray_refinement_tick_id);

   // moving_atoms_visited_residues.clear(); // this is used for HUD label colour

   int n_restraints = last_restraints->make_restraints(imol_moving_atoms,
						       geom, flags,
						       do_residue_internal_torsions,
						       do_trans_peptide_restraints,
						       rama_plot_restraints_weight,
						       do_rama_restraints,
						       true, true, make_auto_h_bond_restraints_flag,
						       pseudo_bonds_type);
                                                       // link and flank args default true

   if (use_harmonic_approximation_for_NBCs) {
      std::cout << "INFO:: using soft harmonic restraints for NBC" << std::endl;
      last_restraints->set_use_harmonic_approximations_for_nbcs(true);
   }

   if (pull_restraint_neighbour_displacement_max_radius > 1.99) {
      last_restraints->set_use_proportional_editing(true);
      last_restraints->pull_restraint_neighbour_displacement_max_radius =
         pull_restraint_neighbour_displacement_max_radius;
   }

   last_restraints->set_geman_mcclure_alpha(geman_mcclure_alpha);
   last_restraints->set_lennard_jones_epsilon(lennard_jones_epsilon);
   last_restraints->set_rama_type(restraints_rama_type);
   last_restraints->set_rama_plot_weight(rama_plot_restraints_weight); // >2? danger of non-convergence
                                                                       // if planar peptide restraints are used
   // Oh, I see... it's not just the non-Bonded contacts of the hydrogens.
   // It's the planes, chiral and angles too. Possibly bonds too.
   // How about marking non-H atoms in restraints that contain H atoms as
   // "invisible"? i.e. non-H atoms are not influenced by the positions of the
   // Hydrogen atoms (but Hydrogen atoms *are* influenced by the positions of the
   // non-Hydrogen atoms). This seems like a lot of work. Might be easier to turn
   // off angle restraints for H-X-X (but not H-X-H) in the first instance, that
   // should go most of the way to what "invisible" atoms would do, I imagine.
   // is_H_non_bonded_contact should be renamed to is_H_turn_offable_restraint
   // or something.
   //
   // last_restraints->set_apply_H_non_bonded_contacts(false);

   if (do_rotamer_restraints) {
      std::vector<std::pair<mmdb::Residue *, std::vector<coot::dict_torsion_restraint_t> > > rotamer_torsions = make_rotamer_torsions(local_residues);
      std::cout << "debug:: calling add_or_replace_torsion_restraints_with_closest_rotamer_restraints() from make_last_restraints() " << std::endl;
      last_restraints->add_or_replace_torsion_restraints_with_closest_rotamer_restraints(rotamer_torsions);
   }

   if (molecules[imol_moving_atoms].extra_restraints.has_restraints()) {
      std::cout << "debug:: calling add_extra_restraints() from make_last_restraints() " << std::endl;
      last_restraints->add_extra_restraints(imol_moving_atoms, "user-defined from make_last_restraints()",
                                            molecules[imol_moving_atoms].extra_restraints, geom);
   }

   if (do_numerical_gradients)
      last_restraints->set_do_numerical_gradients();

   bool found_restraints_flag = false;

   if (last_restraints->size() > 0) {

      last_restraints->analyze_for_bad_restraints();
      thread_for_refinement_loop_threaded();
      found_restraints_flag = true;
      // rr.found_restraints_flag = true;
      draw_gl_ramachandran_plot_flag = true;

      // are you looking for conditionally_wait_for_refinement_to_finish() ?

      if (refinement_immediate_replacement_flag) {
         // wait until refinement finishes
         while (restraints_lock) {
            std::this_thread::sleep_for(std::chrono::milliseconds(7));
            std::cout << "INFO:: make_last_restraints() [immediate] restraints locked by "
                      << restraints_locking_function_name << std::endl;
         }
      }

   } else {
      continue_threaded_refinement_loop = false;
      if (use_graphics_interface_flag) {
         // GtkWidget *widget = create_no_restraints_info_dialog();
         // GtkWidget *widget = widget_from_builder("no_restraints_info_dialog");
         // gtk_widget_show(widget);
      }
   }

   return found_restraints_flag;
}


// simple mmdb::Residue * interface to refinement.  20081216
coot::refinement_results_t
molecules_container_t::generate_molecule_and_refine(int imol,  // needed for UDD Atom handle transfer
                                                    const std::vector<mmdb::Residue *> &residues_in,
                                                    const std::string &alt_conf,
                                                    mmdb::Manager *mol,
                                                    bool use_map_flag) {

   // 20221018-PE make a function in the class
   auto set_refinement_flags = [] () {
      return coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_AND_CHIRALS;
   };
   int cif_dictionary_read_number = 44; // make this a class member

   bool do_torsion_restraints = true;
   bool do_rama_restraints = false; // or true?
   bool moving_atoms_have_hydrogens_displayed = false;


   coot::refinement_results_t rr(0, GSL_CONTINUE, "");

   if (is_valid_map_molecule(imol_refinement_map) || (! use_map_flag)) {
      // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
      coot::restraint_usage_Flags flags = set_refinement_flags();
      bool do_residue_internal_torsions = false;
      if (do_torsion_restraints) {
	 do_residue_internal_torsions = 1;
	 flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_AND_CHIRALS;
      }

      if (do_rama_restraints)
	 // flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_RAMA;
	 flags = coot::ALL_RESTRAINTS;

      std::vector<coot::atom_spec_t> fixed_atom_specs = molecules[imol].get_fixed_atoms();

      // refinement goes a bit wonky if there are multiple occurrances of the same residue
      // in input residue vector, so let's filter out duplicates here
      //
      std::vector<mmdb::Residue *> residues;
      std::set<mmdb::Residue *> residues_set;
      std::set<mmdb::Residue *>::const_iterator it;
      for (std::size_t i=0; i<residues_in.size(); i++)
	 residues_set.insert(residues_in[i]);
      residues.reserve(residues_set.size());
      for(it=residues_set.begin(); it!=residues_set.end(); ++it)
	 residues.push_back(*it);

      // OK, so the passed residues are the residues in the graphics_info_t::molecules[imol]
      // molecule.  We need to do 2 things:
      //
      // convert the mmdb::Residue *s of the passed residues to the mmdb::Residue *s of residues mol
      //
      // and
      //
      // in create_mmdbmanager_from_res_vector() make sure that that contains the flanking atoms.
      // The create_mmdbmanager_from_res_vector() from this class is used, not coot::util
      //
      // The flanking atoms are fixed the passed residues are not fixed.
      // Keep a clear head.

      std::vector<std::string> residue_types = coot::util::residue_types_in_residue_vec(residues);
      // use try_dynamic_add()
      bool have_restraints = geom.have_dictionary_for_residue_types(residue_types, imol, cif_dictionary_read_number);
      cif_dictionary_read_number += residue_types.size();

      if (have_restraints) {

	 std::string residues_alt_conf = alt_conf;
	 imol_moving_atoms = imol;
	 std::pair<mmdb::Manager *, std::vector<mmdb::Residue *> > residues_mol_and_res_vec =
	    create_mmdbmanager_from_res_vector(residues, imol, mol, residues_alt_conf);

	 if (true) { // debug
	    mmdb::Manager *residues_mol = residues_mol_and_res_vec.first;
	    int imod = 1;
	    mmdb::Model *model_p = residues_mol->GetModel(imod);
	    if (model_p) {
	       int n_chains = model_p->GetNumberOfChains();
	       for (int ichain=0; ichain<n_chains; ichain++) {
		  mmdb::Chain *chain_p = model_p->GetChain(ichain);
		  std::cout << "DEBUG:: in generate_molecule_and_refine() residues_mol_and_res_vec mol: chain: "
                            << chain_p->GetChainID() << std::endl;
		  int nres = chain_p->GetNumberOfResidues();
		  for (int ires=0; ires<nres; ires++) {
		     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
		     std::cout << "DEBUG:: in generate_molecule_and_refine() residues_mol_and_res_vec mol:   residue "
			       << coot::residue_spec_t(residue_p) << " residue "
			       << residue_p << " chain " << residue_p->chain << " index "
			       << residue_p->index << std::endl;
		  }
	       }
	    }
	 }

	 // We only want to act on these new residues and molecule, if
	 // there is something there.
	 //
	 if (residues_mol_and_res_vec.first) {

	    // Now we want to do an atom name check.  This stops exploding residues.
	    //
	    bool check_hydrogens_too_flag = false;
	    std::pair<bool, std::vector<std::pair<mmdb::Residue *, std::vector<std::string> > > >
	       icheck_atoms = geom.atoms_match_dictionary(imol, residues, check_hydrogens_too_flag, false);

	    if (! icheck_atoms.first) {

               std::cout << "WARNING:: non-matching atoms! " << std::endl;

	    } else {

	       moving_atoms_have_hydrogens_displayed = true;
	       if (! molecules[imol].hydrogen_atom_should_be_drawn())
		  moving_atoms_have_hydrogens_displayed = false;

	       atom_selection_container_t local_moving_atoms_asc =
		  make_moving_atoms_asc(residues_mol_and_res_vec.first, residues);

	       // 20221018-PE make_moving_atoms_graphics_object(imol, local_moving_atoms_asc); not today!

               int n_cis = coot::util::count_cis_peptides(local_moving_atoms_asc.mol);
               // moving_atoms_n_cis_peptides = n_cis; // 20221018-PE not today

	       std::vector<std::pair<bool,mmdb::Residue *> > local_residues;  // not fixed.
	       for (unsigned int i=0; i<residues_mol_and_res_vec.second.size(); i++)
		  local_residues.push_back(std::pair<bool, mmdb::Residue *>(0, residues_mol_and_res_vec.second[i]));

	       moving_atoms_asc_type = NEW_COORDS_REPLACE;

	       int imol_for_map = imol_refinement_map;
	       clipper::Xmap<float> *xmap_p = dummy_xmap;

	       if (is_valid_map_molecule(imol_for_map))
		  xmap_p = &molecules[imol_for_map].xmap;

	       bool found_restraints_flag = make_last_restraints(local_residues,
								 local_moving_atoms_asc.links,
								 geom,
								 residues_mol_and_res_vec.first,
								 fixed_atom_specs,
								 flags, use_map_flag, xmap_p);

               if (last_restraints) {
                  // 20220423-PE I can't do this here because setup_minimize() has not been called yet
                  // rr = last_restraints->get_refinement_results();
               }
	       rr.found_restraints_flag = found_restraints_flag;

	    }
	 }
      } else {
	 // we didn't have restraints for everything.
	 //
	 std::pair<int, std::vector<std::string> > icheck =
	    check_dictionary_for_residue_restraints(imol, residues);
	 if (icheck.first == 0) {
	    std::cout << "WARNING:: <some info here about missing residue types> " << std::endl;
	 }
      }
   }
   return rr;
}


int
molecules_container_t::mutate(int imol, const std::string &cid, const std::string &new_residue_type) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec = atom_cid_to_atom_spec(imol, cid);
      coot::residue_spec_t residue_spec(atom_spec);
      status = molecules[imol].mutate(residue_spec, new_residue_type);
      // qstd::cout << "mutate status " << status << std::endl;
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}

#include "coot-utils/blob-line.hh"

std::pair<bool, clipper::Coord_orth>
molecules_container_t::go_to_blob(float x1, float y1, float z1, float x2, float y2, float z2, float contour_level) {

   std::pair<bool, clipper::Coord_orth> p;

   clipper::Coord_orth p1(x1,y1,z1);
   clipper::Coord_orth p2(x2,y2,z2);

   // iterate through all the maps (another day)

   if (is_valid_map_molecule(imol_refinement_map)) {
      const clipper::Xmap<float> &xmap = molecules[imol_refinement_map].xmap;
      std::pair<bool, clipper::Coord_orth> pp = coot::find_peak_along_line_favour_front(p1, p2, contour_level, xmap);
      p = pp;
   }
   return p;
}


int
molecules_container_t::side_chain_180(int imol, const std::string &atom_cid) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec = atom_cid_to_atom_spec(imol, atom_cid);
      coot::residue_spec_t residue_spec(atom_spec);
      status = molecules[imol].side_chain_180(residue_spec, atom_spec.alt_conf, &geom);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;

}

std::string
molecules_container_t::jed_flip(int imol, const std::string &atom_cid, bool invert_selection) {

   std::string message;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec = atom_cid_to_atom_spec(imol, atom_cid);
      coot::residue_spec_t res_spec(atom_spec);
      std::string atom_name = atom_spec.atom_name;
      std::string alt_conf = "";
      message = molecules[imol].jed_flip(res_spec, atom_name, alt_conf, invert_selection, &geom);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return message;
}


void
molecules_container_t::coot_all_atom_contact_dots_instanced(mmdb::Manager *mol, int imol) {

   // 20221025-PE fill me later.

}

int
molecules_container_t::add_waters(int imol_model, int imol_map) {

   int status = 0;

   // 20221025-PE Fill me later

   return status;

}


int
molecules_container_t::delete_side_chain(int imol, const std::string &chain_id, int res_no, const std::string &ins_code) {

   int status = 0;

   // 20221025-PE Fill me later
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t res_spec(chain_id, res_no, ins_code);
      // molecules[imol].delete_side_chain(res_spec);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}

int
molecules_container_t::fill_side_chain(int imol, const std::string &chain_id, int res_no, const std::string &ins_code) {

   int status = 0;

   // 20221025-PE Fill me later
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t res_spec(chain_id, res_no, ins_code);
      // molecules[imol].fill_side_chain(res_spec);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}


std::vector<std::string>
molecules_container_t::chains_in_model(int imol) const {

   std::vector<std::string> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].chains_in_model();
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}

std::vector<std::pair<coot::residue_spec_t, std::string> >
molecules_container_t::get_single_letter_codes_for_chain(int imol, const std::string &chain_id) const {

   std::vector<std::pair<coot::residue_spec_t, std::string> > v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].get_single_letter_codes_for_chain(chain_id);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}
