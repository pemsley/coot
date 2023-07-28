
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>  // to_string()
#include <glm/gtx/rotate_vector.hpp>
#include <glm/ext.hpp>

#include "coot_molecule.hh"
#include "coot-utils/oct.hh"
#include "coot-utils/cylinder.hh"

void
coot::molecule_t::setup_cylinder_clashes(instanced_mesh_t &im, const atom_overlaps_dots_container_t &c,
                                         float ball_size, unsigned int num_subdivisions,
                                         const std::string &molecule_name_stub) const {

   auto coord_orth_to_glm = [] (const clipper::Coord_orth &co) {
      return glm::vec3(co.x(), co.y(), co.z());
   };

   if (c.clashes.size() > 0) {
      std::string clashes_name = molecule_name_stub + std::string(" clashes");

      instanced_geometry_t ig_empty;
      im.add(ig_empty);
      instanced_geometry_t &ig = im.geom.back();

      // -------------------vertices and triangles  --------------

      glm::vec3 start(0,0,0);
      glm::vec3 end(0,0,1);
      auto start_end = std::make_pair(start, end);
      float h = 1.0;
      unsigned int n_slices = 32;
      cylinder cyl(start_end, 1.0, 1.0, h, n_slices, 2);
      float z_scale = 0.37;
      // I use 1.1 here so that the clash markup is a bit fatter
      // than a typical ball.
      float unstubby_cap_factor = 0.3;
      bool cigar_mode = true;
      if (cigar_mode) {
         z_scale = 0.25;
         unstubby_cap_factor = 0.92;
      }
      cyl.set_unstubby_rounded_cap_factor(unstubby_cap_factor);
      cyl.add_octahemisphere_start_cap();
      cyl.add_octahemisphere_end_cap();
      if (cigar_mode)
         cyl.z_translate(-0.2);


      ig.vertices.resize(cyl.vertices.size());
      for (unsigned int i=0; i<cyl.vertices.size(); i++)
         ig.vertices[i] = api::vn_vertex(cyl.vertices[i].pos, cyl.vertices[i].pos);
      ig.triangles = cyl.triangles;

      // -------------------instancing --------------

      glm::vec4 colour = colour_holder_to_glm(colour_holder("#ff59c9"));
      for (unsigned int i=0; i<c.clashes.size(); i++) {
         std::pair<glm::vec3, glm::vec3> pos_pair_clash(glm::vec3(coord_orth_to_glm(c.clashes[i].first)),
                                                        glm::vec3(coord_orth_to_glm(c.clashes[i].second)));
         const glm::vec3 &start_pos  = pos_pair_clash.first;
         const glm::vec3 &finish_pos = pos_pair_clash.second;
         // glm::vec3 b = finish_pos - start_pos;
         glm::vec3 b = start_pos - finish_pos;
         glm::vec3 normalized_bond_orientation(glm::normalize(b));
         glm::mat4 ori = glm::orientation(normalized_bond_orientation, glm::vec3(0.0, 0.0, 1.0));
         glm::vec3 sc(1.1 * ball_size, 1.1 * ball_size, z_scale);
         glm::mat4 unit(1.0);
         glm::mat4 mt_1 = glm::translate(unit, start);
         glm::mat4 mt_2 = mt_1 * ori;
         glm::mat4 mt_3 = glm::scale(mt_2, sc);
         ig.instancing_data_B.push_back(instancing_data_type_B_t(start_pos, colour, sc, ori));
      }
   }
}

void
coot::molecule_t::setup_dots(instanced_mesh_t &im,
                             const atom_overlaps_dots_container_t &c,
                             float ball_size, unsigned int num_subdivisions,
                             const std::string &molecule_name_stub) const {

   instanced_geometry_t ig_empty;
   im.add(ig_empty);
   instanced_geometry_t &ig = im.geom.back();

   // -------------------vertices and triangles  --------------

   std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > octaphere_geom =
      tessellate_octasphere(num_subdivisions);
   ig.vertices.resize(octaphere_geom.first.size());
   for (unsigned int i=0; i<octaphere_geom.first.size(); i++)
      ig.vertices[i] = api::vn_vertex(octaphere_geom.first[i], octaphere_geom.first[i]);
   ig.triangles = octaphere_geom.second;

   // -------------------instancing --------------

   // setup the colour map
   std::map<std::string, coot::colour_holder> colour_map;
   colour_map["blue"      ] = colour_holder_from_colour_name("blue");
   colour_map["sky"       ] = colour_holder_from_colour_name("sky");
   colour_map["sea"       ] = colour_holder_from_colour_name("sea");
   colour_map["greentint" ] = colour_holder_from_colour_name("greentint");
   colour_map["darkpurple"] = colour_holder_from_colour_name("darkpurple");
   colour_map["green"     ] = colour_holder_from_colour_name("green");
   colour_map["orange"    ] = colour_holder_from_colour_name("orange");
   colour_map["orangered" ] = colour_holder_from_colour_name("orangered");
   colour_map["yellow"    ] = colour_holder_from_colour_name("yellow");
   colour_map["yellowtint"] = colour_holder_from_colour_name("yellowtint");
   colour_map["red"       ] = colour_holder_from_colour_name("red");
   colour_map["#55dd55"   ] = colour_holder_from_colour_name("#55dd55");
   colour_map["hotpink"   ] = colour_holder_from_colour_name("hotpink");
   colour_map["grey"      ] = colour_holder_from_colour_name("grey");
   colour_map["magenta"   ] = colour_holder_from_colour_name("magenta");
   colour_map["royalblue" ] = colour_holder_from_colour_name("royalblue");

   std::unordered_map<std::string, std::vector<coot::atom_overlaps_dots_container_t::dot_t> >::const_iterator it;
   for (it=c.dots.begin(); it!=c.dots.end(); ++it) {
      // float specular_strength = 0.5; //  default
      const std::string &type = it->first;
      const std::vector<coot::atom_overlaps_dots_container_t::dot_t> &v = it->second;
      float point_size = ball_size;
      if (type == "vdw-surface") point_size = 0.06; // 0.03 seems too small
      // if (type == "vdw-surface") specular_strength= 0.1; // dull, reduces zoomed out speckles
      std::string mesh_name = molecule_name_stub + type; // instanced_geometry_t doesn't have a name holder
      ig.name = mesh_name + std::string(" ") + type;

      glm::vec3 size(point_size, point_size, point_size);
      for (unsigned int i=0; i<v.size(); i++) {
         const atom_overlaps_dots_container_t::dot_t &dot(v[i]);
         std::map<std::string, colour_holder>::const_iterator it_inner = colour_map.find(dot.col);
         if (it_inner != colour_map.end()) {
            glm::vec4 colour = colour_holder_to_glm(it_inner->second);
            glm::vec3 position(dot.pos.x(), dot.pos.y(), dot.pos.z());
            ig.instancing_data_A.push_back(instancing_data_type_A_t(position, colour, size));
         }
      }
   }
}



//! @return the instanced mesh for the specified ligand
coot::instanced_mesh_t
coot::molecule_t::contact_dots_for_ligand(const std::string &cid, const coot::protein_geometry &geom,
                                          unsigned int num_subdivisions) const {

   // Note: std::unordered_map<std::string, std::vector<dot_t> > dots;
   // class dot_t {
   // public:
   //    double overlap;
   //    clipper::Coord_orth pos;
   //    std::string col;
   //    dot_t(double o, const std::string &col_in, const clipper::Coord_orth &pos_in) : pos(pos_in), col(col_in) {
   //       overlap = o;
   //    }
   // };

   // return all_molecule_contact_dots(geom);


   coot::instanced_mesh_t im;
   float contact_dots_density = 0.7; // 20220308-PE was 1.0
   float cdd = contact_dots_density;

   mmdb::Residue *residue_p = cid_to_residue(cid); // Make get_residue_using_cid()?

   if (residue_p) {

      mmdb::Manager *mol = atom_sel.mol;
      std::vector<mmdb::Residue *> neighbs = coot::residues_near_residue(residue_p, mol, 5);
      coot::atom_overlaps_container_t overlaps(residue_p, neighbs, mol, &geom, 0.5, 0.25);

      coot::atom_overlaps_dots_container_t c = overlaps.contact_dots_for_ligand(cdd);
      float ball_size = 0.07;
      std::string name_stub = "Molecule " + std::to_string(imol_no);

      setup_cylinder_clashes(im, c, ball_size, num_subdivisions, name_stub); //modify reference

      setup_dots(im, c, ball_size, num_subdivisions, name_stub); // modify reference

   }
   return im;
}

//! @return the instanced mesh for the specified molecule
coot::instanced_mesh_t
coot::molecule_t::all_molecule_contact_dots(const coot::protein_geometry &geom,
                                            unsigned int num_subdivisions) const {

   coot::instanced_mesh_t im;

   float contact_dots_density = 0.7; // 20220308-PE was 1.0
   float cdd = contact_dots_density;

   mmdb::Manager *mol = atom_sel.mol;
   bool ignore_waters_flag = false;
   coot::atom_overlaps_container_t overlaps(mol, &geom, ignore_waters_flag, 0.5, 0.25);

   bool make_vdw_surface_flag = false;
   coot::atom_overlaps_dots_container_t c = overlaps.all_atom_contact_dots(cdd, make_vdw_surface_flag);
   float ball_size = 0.07;
   std::string name_stub = "Molecule " + std::to_string(imol_no);

   setup_cylinder_clashes(im, c, ball_size, num_subdivisions, name_stub); //modify reference

   setup_dots(im, c, ball_size, num_subdivisions, name_stub); // modify reference

   return im;
}


std::vector<coot::geometry_distortion_info_container_t>
coot::molecule_t::geometric_distortions_from_mol(const atom_selection_container_t &asc, bool with_nbcs,
                                                 coot::protein_geometry &geom,
                                                 ctpl::thread_pool &static_thread_pool) {

   std::vector<coot::geometry_distortion_info_container_t> dcv;

   if (! asc.mol)
      return dcv;

   int n_models = asc.mol->GetNumberOfModels();

   if (n_models > 0) {

      // 20100629, crash!  Ooops, we can't run over many models
      // because geometry_graphs (and its data member `blocks' are
      // sized to the number of chains).  If we run over all models,
      // then there are too many chains for the indexing of `blocks`
      // -> crash.  So we just use the first model.

      // for (int imod=1; imod<=n_models; imod++) {
      int imod=1;
      {

         mmdb::Model *model_p = asc.mol->GetModel(imod);
         mmdb::Chain *chain_p;
         int nchains = model_p->GetNumberOfChains();
         const char *chain_id;


         for (int ichain=0; ichain<nchains; ichain++) {

            chain_p = model_p->GetChain(ichain);

            if (! chain_p->isSolventChain()) {
               chain_id = chain_p->GetChainID();

               // First make an atom selection of the residues selected to regularize.
               //
               int selHnd = asc.mol->NewSelection(); // yes, it's deleted.
               int nSelResidues;
               mmdb::PResidue *SelResidues = NULL;

               // Consider as the altconf the altconf of one of the residues (we
               // must test that the altlocs of the selected atoms to be either
               // the same as each other (A = A) or one of them is "".  We need to
               // know the mmdb syntax for "either".  Well, now I know that's ",A"
               // (for either blank or "A").
               //
               //
               //
               asc.mol->Select(selHnd, mmdb::STYPE_RESIDUE, imod,
                               chain_id,
                               mmdb::ANY_RES, "*",
                               mmdb::ANY_RES, "*",
                               "*",  // residue name
                               "*",  // Residue must contain this atom name?
                               "*",  // Residue must contain this Element?
                               "*",  // altLocs
                               mmdb::SKEY_NEW // selection key
                               );
               asc.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

               if (nSelResidues <= 0) {

                  std::cout << "ERROR:: No Residues!! This should never happen:" << std::endl;
                  std::cout << "  in geometric_distortions_from_mol()" << std::endl;

               } else { // normal

                  std::vector<coot::atom_spec_t> fixed_atom_specs;

                  // Notice that we have to make 2 atom selections, one, which includes
                  // flanking (and disulphide eventually) residues that is used for the
                  // restraints (restraints_container_t constructor) and one that is the
                  // moving atoms (which does not have flanking atoms).
                  //
                  // The restraints_container_t moves the atom of the mol that is passes to
                  // it.  This must be the same mol as the moving atoms mol so that the
                  // changed atom positions can be seen.  However (as I said) the moving
                  // atom mol should not have flanking residues shown.  So we make an asc
                  // that has the same mol as that passed to the restraints but a different
                  // atom selection (it is the atom selection that is used in the bond
                  // generation).
                  //

//                  20100210 try vector
//                   coot::restraints_container_t restraints(SelResidues, nSelResidues,
//                                                           std::string(chain_id),
//                                                           asc.mol);
                  std::vector<std::pair<bool,mmdb::Residue *> > residue_vec;
                  for (int ires=0; ires<nSelResidues; ires++)
                     residue_vec.push_back(std::pair<bool, mmdb::Residue *> (0, SelResidues[ires]));

                  std::vector<mmdb::Link> links;
                  clipper::Xmap<float> dummy_xmap;

                  coot::restraints_container_t restraints(residue_vec,
                                                          links,
                                                          geom,
                                                          asc.mol,
                                                          fixed_atom_specs,
                                                          &dummy_xmap);

                  // coot::restraint_usage_Flags flags = coot::BONDS;
                  // coot::restraint_usage_Flags flags = coot::BONDS_AND_ANGLES;
                  // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_AND_PLANES;
                  // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_TORSIONS_AND_PLANES;
                  coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_AND_NON_BONDED;
                  flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
                  flags = coot::BONDS_ANGLES_AND_PLANES;
                  flags = coot::BONDS_ANGLES_PLANES_AND_CHIRALS;

                  if (with_nbcs)
                     flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;

                  unsigned int n_threads = coot::get_max_number_of_threads();
                  if (n_threads > 0)
                     restraints.thread_pool(&static_thread_pool, n_threads);
                  short int do_residue_internal_torsions = 0;

                  //                if (do_torsion_restraints) {
                  //                   do_residue_internal_torsions = 1;
                  //                   flags = coot::BONDS_ANGLES_TORSIONS_PLANES_AND_NON_BONDED;
                  //                }

                  //                if (do_peptide_torsion_restraints)
                  //                   do_link_torsions = 1;

                  coot::pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
                  bool do_trans_peptide_restraints = false;
                  int nrestraints =
                     restraints.make_restraints(imol_no, geom,
                                                flags,
                                                do_residue_internal_torsions,
                                                do_trans_peptide_restraints,
                                                0.0, 0, false, false, false,
                                                pseudos);

                  if (nrestraints > 0) {

//                      std::cout << "DEBUG:: model " << imod << " pushing back " << nrestraints
//                                << " restraints" << std::endl;

                     dcv.push_back(restraints.geometric_distortions());

                  } else {

                     // don't give this annoying dialog if restraints
                     // have been read for the residues in this chain.
                     // e.g. a single CLs residues in a chain.
                     int cif_dictionary_read_number = 50;
                     std::vector<std::string> res_types = coot::util::residue_types_in_chain(chain_p);
                     bool hd = geom.have_dictionary_for_residue_types(res_types, imol_no,
                                                                      cif_dictionary_read_number);
                     // we need to signal to the caller that there were no restraints.
                     cif_dictionary_read_number += res_types.size();
                  }
               }
               asc.mol->DeleteSelection(selHnd);
            }
         }
      }
   }
   return dcv;
}

#include "coot-utils/shapes.hh"

//! get the mesh for ligand validation vs dictionary, coloured by badness.
//! greater then 3 standard deviations is fully red.
//! Less than 0.5 standard deviations is fully green.
coot::simple_mesh_t
coot::molecule_t::get_mesh_for_ligand_validation_vs_dictionary(const std::string &ligand_cid,
                                                               coot::protein_geometry &geom,
                                                               ctpl::thread_pool &static_thread_pool) {

   auto coord_orth_to_glm = [] (const clipper::Coord_orth &co) {
      return glm::vec3(co.x(), co.y(), co.z());
   };

   auto add_chiral_lines = [coord_orth_to_glm] (coot::simple_mesh_t &obj,
                                                clipper::Coord_orth &bl_1,
                                                clipper::Coord_orth &bl_2,
                                                clipper::Coord_orth &bl_3,
                                                clipper::Coord_orth &bl_4,
                                                float line_radius,
                                                const glm::vec4 &col,
                                                unsigned int n_slices) {
      std::pair<glm::vec3, glm::vec3> pos_pair_14(glm::vec3(coord_orth_to_glm(bl_1)),
                                                  glm::vec3(coord_orth_to_glm(bl_2)));
      std::pair<glm::vec3, glm::vec3> pos_pair_24(glm::vec3(coord_orth_to_glm(bl_2)),
                                                  glm::vec3(coord_orth_to_glm(bl_4)));
      std::pair<glm::vec3, glm::vec3> pos_pair_34(glm::vec3(coord_orth_to_glm(bl_3)),
                                                  glm::vec3(coord_orth_to_glm(bl_4)));
      float h_1c = glm::distance(pos_pair_14.first, pos_pair_14.second);
      float h_2c = glm::distance(pos_pair_24.first, pos_pair_24.second);
      float h_3c = glm::distance(pos_pair_34.first, pos_pair_34.second);
      cylinder cyl_1c(pos_pair_14, line_radius, line_radius, h_1c, col, n_slices);
      cylinder cyl_2c(pos_pair_24, line_radius, line_radius, h_2c, col, n_slices);
      cylinder cyl_3c(pos_pair_34, line_radius, line_radius, h_3c, col, n_slices);
      obj.add_submesh(simple_mesh_t(cyl_1c));
      obj.add_submesh(simple_mesh_t(cyl_2c));
      obj.add_submesh(simple_mesh_t(cyl_3c));

   };

   coot::simple_mesh_t obj;
   mmdb::Residue *residue_p = cid_to_residue(ligand_cid);
   if (residue_p) {
      mmdb::Manager *mol = coot::util::create_mmdbmanager_from_residue(residue_p);
      if (mol) {
         bool with_nbcs = true; // pass this?
         atom_selection_container_t asc = make_asc(mol);
         std::vector<coot::geometry_distortion_info_container_t> v =
            geometric_distortions_from_mol(asc, with_nbcs, geom, static_thread_pool);
         if (v.size() == 1) {
            if (v[0].geometry_distortion.size() > 1) {
               coot::geometry_distortion_info_container_t gdc = v[0];

               if (gdc.geometry_distortion.size()) {

                  std::string name = std::string("Ligand Distortion of ");
                  name += residue_p->GetChainID();
                  name += " ";
                  name += coot::util::int_to_string(residue_p->GetSeqNum());
                  name += " ";
                  name += residue_p->GetResName();

                  std::cout << "debug:: in get_mesh_for_ligand_validation_vs_dictionary with "
                            << gdc.geometry_distortion.size() << " distortions" << std::endl;

                  for (unsigned int i=0; i<gdc.geometry_distortion.size(); i++) {
                     coot::simple_restraint &rest = gdc.geometry_distortion[i].restraint;

                     if (rest.restraint_type == coot::BOND_RESTRAINT) {
                        mmdb::Atom *at_1 = residue_p->GetAtom(rest.atom_index_1);
                        mmdb::Atom *at_2 = residue_p->GetAtom(rest.atom_index_2);
                        if (at_1 && at_2) {
                           clipper::Coord_orth p1(at_1->x, at_1->y, at_1->z);
                           clipper::Coord_orth p2(at_2->x, at_2->y, at_2->z);
                           double d = sqrt((p2-p1).lengthsq());
                           double distortion = d - rest.target_value;
                           double pen_score = fabs(distortion/rest.sigma);
                           coot::colour_holder ch(pen_score, 0.1, 5, true, "");
                           ch.scale_intensity(0.5);
                           glm::vec4 col = colour_holder_to_glm(ch);
                           const unsigned int n_slices = 16;
                           float line_radius = 0.104f;  // more than 0.103 (atoms)
                           std::pair<glm::vec3, glm::vec3> pos_pair(glm::vec3(coord_orth_to_glm(p1)),
                                                                    glm::vec3(coord_orth_to_glm(p2)));

                           // obj.add_cylinder(pos_pair, ch, line_radius, n_slices, true, true,
                           // meshed_generic_display_object::ROUNDED_CAP,
                           // meshed_generic_display_object::ROUNDED_CAP);

                           cylinder cyl(pos_pair, line_radius, line_radius, d, col, n_slices);
                           obj.add_submesh(simple_mesh_t(cyl));
                        }
                     }

                     if (rest.restraint_type == coot::ANGLE_RESTRAINT) {

                        mmdb::Atom *at_1 = residue_p->GetAtom(rest.atom_index_1);
                        mmdb::Atom *at_2 = residue_p->GetAtom(rest.atom_index_2);
                        mmdb::Atom *at_3 = residue_p->GetAtom(rest.atom_index_3);
                        if (at_1 && at_2 && at_3) {
                           clipper::Coord_orth p1(at_1->x, at_1->y, at_1->z);
                           clipper::Coord_orth p2(at_2->x, at_2->y, at_2->z);
                           clipper::Coord_orth p3(at_3->x, at_3->y, at_3->z);
                           double angle_rad = clipper::Coord_orth::angle(p1, p2, p3);
                           double angle = clipper::Util::rad2d(angle_rad);
                           double distortion = fabs(angle - rest.target_value);
                           coot::colour_holder ch(distortion, 0.1, 5, true, "");
                           ch.scale_intensity(0.6);

                           try {
                              float radius = 0.66;
                              float radius_inner = 0.1; // should match bond line_radius
                              coot::arc_info_type arc_angle_info(at_1, at_2, at_3);

                              shapes::arc_t arc(arc_angle_info, radius, radius_inner, ch);
                              simple_mesh_t arc_mesh = coot::arc_mesh(arc);
                              obj.add_submesh(arc_mesh);
                           }
                           catch (const std::runtime_error &rte) {
                              std::cout << "WARNING:: " << rte.what() << std::endl;
                           }
                        }
                     }

                     if (rest.restraint_type == coot::CHIRAL_VOLUME_RESTRAINT) {
                        mmdb::Atom *at_c = residue_p->GetAtom(rest.atom_index_centre);
                        mmdb::Atom *at_1 = residue_p->GetAtom(rest.atom_index_1);
                        mmdb::Atom *at_2 = residue_p->GetAtom(rest.atom_index_2);
                        mmdb::Atom *at_3 = residue_p->GetAtom(rest.atom_index_3);
                        if (at_c && at_1 && at_2 && at_3) {
                           clipper::Coord_orth pc(at_c->x, at_c->y, at_c->z);
                           clipper::Coord_orth p1(at_1->x, at_1->y, at_1->z);
                           clipper::Coord_orth p2(at_2->x, at_2->y, at_2->z);
                           clipper::Coord_orth p3(at_3->x, at_3->y, at_3->z);
                           clipper::Coord_orth bl_1 = 0.6 * pc + 0.4 * p1;
                           clipper::Coord_orth bl_2 = 0.6 * pc + 0.4 * p2;
                           clipper::Coord_orth bl_3 = 0.6 * pc + 0.4 * p3;
                           double distortion = sqrt(fabs(gdc.geometry_distortion[i].distortion_score));
                           coot::colour_holder ch(distortion, 0.1, 5, true, "");
                           glm::vec4 col = colour_holder_to_glm(ch);

                           std::pair<glm::vec3, glm::vec3> pos_pair_12(glm::vec3(coord_orth_to_glm(bl_1)),
                                                                       glm::vec3(coord_orth_to_glm(bl_2)));
                           std::pair<glm::vec3, glm::vec3> pos_pair_13(glm::vec3(coord_orth_to_glm(bl_1)),
                                                                       glm::vec3(coord_orth_to_glm(bl_3)));
                           std::pair<glm::vec3, glm::vec3> pos_pair_23(glm::vec3(coord_orth_to_glm(bl_1)),
                                                                       glm::vec3(coord_orth_to_glm(bl_3)));
                           const unsigned int n_slices = 16;
                           float line_radius = 0.08f;
                           cylinder cyl_1(pos_pair_12, line_radius, line_radius, n_slices, true, true);
                           cylinder cyl_2(pos_pair_13, line_radius, line_radius, n_slices, true, true);
                           cylinder cyl_3(pos_pair_23, line_radius, line_radius, n_slices, true, true);
                           obj.add_submesh(simple_mesh_t(cyl_1));
                           obj.add_submesh(simple_mesh_t(cyl_2));
                           obj.add_submesh(simple_mesh_t(cyl_3));

                           // return (if possible) the atom attached to
                           // at_c that is not at_1, at_2 or at_3.
                           mmdb::Atom *at_4th = coot::chiral_4th_atom(residue_p, at_c, at_1, at_2, at_3);
                           if (at_4th) {
                              std::cout << "    " << coot::atom_spec_t(at_4th) << std::endl;
                              clipper::Coord_orth p4(at_4th->x, at_4th->y, at_4th->z);
                              clipper::Coord_orth bl_4 = 0.6 * pc + 0.4 * p4;
                              add_chiral_lines(obj, bl_1, bl_2, bl_3, bl_4, line_radius, col, n_slices); // add to obj
                           } else {
                              // make 4th tetrahedron point from the others
                              clipper::Coord_orth neighb_sum = p1 + p2 + p3;
                              clipper::Coord_orth neighb_average = 0.33333333 * neighb_sum;
                              clipper::Coord_orth dir_unit(clipper::Coord_orth(pc - neighb_average).unit());
                              clipper::Coord_orth p4(pc + 1.2 * dir_unit);
                              clipper::Coord_orth bl_4 = 0.6 * pc + 0.4 * p4;
                              add_chiral_lines(obj, bl_1, bl_2, bl_3, bl_4, line_radius, col, n_slices); // add to obj
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   return obj;
}

