
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>  // to_string()
#include <glm/gtx/rotate_vector.hpp>
#include <glm/ext.hpp>

#include "coot-utils/atom-overlaps.hh"
#include "coot_molecule.hh"
#include "coot-utils/oct.hh"
#include "coot-utils/cylinder.hh"

//! @return the instanced mesh for the specified ligand
coot::instanced_mesh_t
coot::molecule_t::contact_dots_for_ligand(const std::string &cid, const coot::protein_geometry &geom) const {

   // this function is in src/mesh-generic-display-object.cc
   auto colour_holder_to_glm = [] (const colour_holder &ch) {
      return glm::vec4(ch.red, ch.green, ch.blue, 1.0f);
   };

   auto coord_orth_to_glm = [] (const clipper::Coord_orth &co) {
      return glm::vec3(co.x(), co.y(), co.z());
   };

   auto setup_cylinder_clashes = [coord_orth_to_glm, colour_holder_to_glm]
      (instanced_mesh_t &im, const atom_overlaps_dots_container_t &c,
       float ball_size, unsigned int num_subdivisions,
       const std::string &molecule_name_stub) {

      if (c.clashes.size() > 0) {
         std::string clashes_name = molecule_name_stub + std::string(" clashes");

         instanced_geometry_t ig_empty;
         im.add(ig_empty);
         instanced_geometry_t &ig = im.geom.back();

         // -------------------vertices and triangles  --------------

         glm::vec3 start(0,0,0);
         glm::vec3 end(0,0,0);
         auto start_end = std::make_pair(start, end);
         float h = 1.0;
         unsigned int n_slices = 16;
         cylinder cyl(start_end, ball_size, ball_size, h, n_slices, 2);
         float z_scale = 0.37;
         float unstubby_cap_factor = 1.1/z_scale;
         cyl.set_unstubby_rounded_cap_factor(unstubby_cap_factor);
         cyl.add_octahemisphere_start_cap();
         cyl.add_octahemisphere_end_cap();

         ig.vertices.resize(cyl.vertices.size());
         for (unsigned int i=0; i<cyl.vertices.size(); i++)
            ig.vertices[i] = api::vn_vertex(cyl.vertices[i].pos, cyl.vertices[i].pos);
         ig.triangles = cyl.triangles;

         // -------------------instancing --------------

         glm::vec4 colour = colour_holder_to_glm(colour_holder("#ff59c9"));
         glm::vec3 size(ball_size, ball_size, ball_size);
         for (unsigned int i=0; i<c.clashes.size(); i++) {
            std::pair<glm::vec3, glm::vec3> pos_pair_clash(glm::vec3(coord_orth_to_glm(c.clashes[i].first)),
                                                           glm::vec3(coord_orth_to_glm(c.clashes[i].second)));
            const glm::vec3 &start_pos  = pos_pair_clash.first;
            const glm::vec3 &finish_pos = pos_pair_clash.second;
            glm::vec3 b = finish_pos - start_pos;
            glm::vec3 normalized_bond_orientation(glm::normalize(b));
            glm::mat4 ori = glm::orientation(normalized_bond_orientation, glm::vec3(0.0, 0.0, 1.0));
            glm::vec3 sc(1.1, 1.1, z_scale);
            glm::mat4 unit(1.0);
            glm::mat4 mt_1 = glm::translate(unit, start);
            glm::mat4 mt_2 = mt_1 * ori;
            glm::mat4 mt_3 = glm::scale(mt_2, sc);
            // mats.push_back(mt_3);
            ig.instancing_data_B.push_back(instancing_data_type_B_t(start_pos, colour, size, mt_3));// might be finish_pos
         }
      }
   };

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

   auto setup_dots = [colour_holder_to_glm] (instanced_mesh_t &im,
                                             const atom_overlaps_dots_container_t &c,
                                             float ball_size, unsigned int num_subdivisions,
                                             const std::string &molecule_name_stub) {

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
         if (type == "vdw-surface") point_size = 0.03;
         // if (type == "vdw-surface") specular_strength= 0.1; // dull, reduces zoomed out speckles
         std::string mesh_name = molecule_name_stub + type; // instanced_geometry_t doesn't have a name holder
         ig.name = mesh_name + std::string(" ") + type;

         glm::vec3 size(point_size, point_size, point_size);
         for (unsigned int i=0; i<v.size(); i++) {
            const atom_overlaps_dots_container_t::dot_t &dot(v[i]);
            std::map<std::string, coot::colour_holder>::const_iterator it = colour_map.find(dot.col);
            if (it != colour_map.end()) {
               glm::vec4 colour = colour_holder_to_glm(it->second);
               glm::vec3 position(dot.pos.x(), dot.pos.y(), dot.pos.z());
               ig.instancing_data_A.push_back(instancing_data_type_A_t(position, colour, size));
            }
         }
      }
   };

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
      float tube_radius = ball_size;
      std::string name_stub = "Molecule " + std::to_string(imol_no);
      unsigned int num_subdivisions = 3;

      setup_cylinder_clashes(im, c, ball_size, num_subdivisions, name_stub); //modify reference

      setup_dots(im, c, ball_size, num_subdivisions, name_stub); // modify reference

   }
   return im;
}

//! @return the instanced mesh for the specified molecule
coot::instanced_mesh_t
coot::molecule_t::all_molecule_contact_dots(const coot::protein_geometry &geom) const {

   coot::instanced_mesh_t im;

   return im;
}
