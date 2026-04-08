
#include <functional>

#include "translation-gizmo.hh"
#include "coot-utils/cylinder.hh"
#include "generic-vertex.hh"

// see also setup_mesh_for_translation_gizmo()

void
translation_gizmo_t::reset_scale() {
   init();
}

void
translation_gizmo_t::init() {

   gizmo_origin = coot::Cartesian(0,0,0);
   attached_to_generic_display_object_number = UNATTACHED;
   attached_to_molecule_number = UNATTACHED;

   glm::vec4 red  (0.9, 0.2, 0.2, 1.0);
   glm::vec4 green(0.2, 0.9, 0.2, 1.0);
   glm::vec4 blue (0.2, 0.2, 0.9, 1.0);

   glm::vec3 origin(0,0,0);
   glm::vec3 up_x(1,0,0);
   glm::vec3 up_y(0,1,0);
   glm::vec3 up_z(0,0,1);
   cylinder c_x(std::make_pair(origin, up_x), 0.2, 0.2, 1.4, 16);
   cylinder c_y(std::make_pair(origin, up_y), 0.2, 0.2, 1.4, 16);
   cylinder c_z(std::make_pair(origin, up_z), 0.2, 0.2, 1.4, 16);
   c_x.add_flat_start_cap();
   c_y.add_flat_start_cap();
   c_z.add_flat_start_cap();
   c_x.add_flat_end_cap();
   c_y.add_flat_end_cap();
   c_z.add_flat_end_cap();

   coot::simple_mesh_t x_axis_mesh_c0(c_x.vertices, c_x.triangles);
   coot::simple_mesh_t x_axis_mesh_c1(c_x.vertices, c_x.triangles);
   coot::simple_mesh_t x_axis_mesh_c2(c_x.vertices, c_x.triangles);
   coot::simple_mesh_t x_axis_mesh_c3(c_x.vertices, c_x.triangles);
   coot::simple_mesh_t y_axis_mesh_c0(c_y.vertices, c_y.triangles);
   coot::simple_mesh_t y_axis_mesh_c1(c_y.vertices, c_y.triangles);
   coot::simple_mesh_t y_axis_mesh_c2(c_y.vertices, c_y.triangles);
   coot::simple_mesh_t y_axis_mesh_c3(c_y.vertices, c_y.triangles);
   coot::simple_mesh_t z_axis_mesh_c0(c_z.vertices, c_z.triangles);
   coot::simple_mesh_t z_axis_mesh_c1(c_z.vertices, c_z.triangles);
   coot::simple_mesh_t z_axis_mesh_c2(c_z.vertices, c_z.triangles);
   coot::simple_mesh_t z_axis_mesh_c3(c_z.vertices, c_z.triangles);
   std::pair<glm::vec3, glm::vec3> pp_x(glm::vec3(1.0, 0.0, 0.0), glm::vec3(0.0, 0.0, 0.0));
   std::pair<glm::vec3, glm::vec3> pp_y(glm::vec3(0.0, 1.0, 0.0), glm::vec3(0.0, 0.0, 0.0));
   std::pair<glm::vec3, glm::vec3> pp_z(glm::vec3(0.0, 0.0, 1.0), glm::vec3(0.0, 0.0, 0.0));
   cylinder x_axis_cone(pp_x, 0.6, 0.0, 2.0, 16);
   cylinder y_axis_cone(pp_y, 0.6, 0.0, 2.0, 16);
   cylinder z_axis_cone(pp_z, 0.6, 0.0, 2.0, 16);
   x_axis_cone.add_flat_end_cap();
   y_axis_cone.add_flat_end_cap();
   z_axis_cone.add_flat_end_cap();
   coot::simple_mesh_t x_axis_mesh_cone(x_axis_cone.vertices, x_axis_cone.triangles);
   coot::simple_mesh_t y_axis_mesh_cone(y_axis_cone.vertices, y_axis_cone.triangles);
   coot::simple_mesh_t z_axis_mesh_cone(z_axis_cone.vertices, z_axis_cone.triangles);

   x_axis_mesh_c1.translate(  glm::vec3(2.0, 0.0, 0.0));
   x_axis_mesh_c2.translate(  glm::vec3(4.0, 0.0, 0.0));
   x_axis_mesh_c3.translate(  glm::vec3(6.0, 0.0, 0.0));
   x_axis_mesh_cone.translate(glm::vec3(8.0, 0.0, 0.0));

   y_axis_mesh_c1.translate(  glm::vec3(0.0, 2.0, 0.0));
   y_axis_mesh_c2.translate(  glm::vec3(0.0, 4.0, 0.0));
   y_axis_mesh_c3.translate(  glm::vec3(0.0, 6.0, 0.0));
   y_axis_mesh_cone.translate(glm::vec3(0.0, 8.0, 0.0));

   z_axis_mesh_c0.translate(  glm::vec3(0.0, 0.0, 1.2));
   z_axis_mesh_c1.translate(  glm::vec3(0.0, 0.0, 3.2));
   z_axis_mesh_c2.translate(  glm::vec3(0.0, 0.0, 5.2));
   z_axis_mesh_c3.translate(  glm::vec3(0.0, 0.0, 7.2));
   z_axis_mesh_cone.translate(glm::vec3(0.0, 0.0, 8.0));

   mesh.clear();

   std::vector<std::reference_wrapper<coot::simple_mesh_t> > x_sub_meshes = {
     std::ref(x_axis_mesh_c0), std::ref(x_axis_mesh_c1), std::ref(x_axis_mesh_c2),
     std::ref(x_axis_mesh_c3), std::ref(x_axis_mesh_cone) };
   std::vector<std::reference_wrapper<coot::simple_mesh_t> > y_sub_meshes = {
     std::ref(y_axis_mesh_c0), std::ref(y_axis_mesh_c1), std::ref(y_axis_mesh_c2),
     std::ref(y_axis_mesh_c3), std::ref(y_axis_mesh_cone) };
   std::vector<std::reference_wrapper<coot::simple_mesh_t> > z_sub_meshes = {
     std::ref(z_axis_mesh_c0), std::ref(z_axis_mesh_c1), std::ref(z_axis_mesh_c2),
     std::ref(z_axis_mesh_c3), std::ref(z_axis_mesh_cone) };
   for (auto &m : x_sub_meshes) { m.get().change_colour(red); }
   for (auto &m : y_sub_meshes) { m.get().change_colour(green); }
   for (auto &m : z_sub_meshes) { m.get().change_colour(blue); }
   for (auto &m : x_sub_meshes) { mesh.add_submesh(m); }
   for (auto &m : y_sub_meshes) { mesh.add_submesh(m); }
   for (auto &m : z_sub_meshes) { mesh.add_submesh(m); }

   // Now bring the whole mesh to Unit size
   mesh.scale(1.0/7.8);
   scale_factor = 1.0;

}

translation_gizmo_t::pick_info_t
translation_gizmo_t::pick(const coot::Cartesian &front, const coot::Cartesian &back) {

   pick_info_t pi(NONE);
   float sf = scale_factor;
   float dist_crit = 2.05;
   coot::Cartesian x_axis_tip(sf, 0, 0);
   coot::Cartesian y_axis_tip(0, sf, 0);
   coot::Cartesian z_axis_tip(0, 0, sf);
   x_axis_tip += gizmo_origin;
   y_axis_tip += gizmo_origin;
   z_axis_tip += gizmo_origin;
   if (x_axis_tip.within_box(front, back)) {
      float dist = x_axis_tip.distance_to_line(front, back);
      // std::cout << "dist x " << dist << std::endl;
      if (dist < dist_crit)
         pi = X_AXIS;
   }
   if (y_axis_tip.within_box(front, back)) {
      float dist = y_axis_tip.distance_to_line(front, back);
      // std::cout << "dist y " << dist << std::endl;
      if (dist < dist_crit)
         pi = Y_AXIS;
   }
   if (z_axis_tip.within_box(front, back)) {
      float dist = z_axis_tip.distance_to_line(front, back);
      // std::cout << "dist z " << dist << std::endl;
      if (dist < dist_crit)
         pi = Z_AXIS;
   }
   return pi;
}

void
translation_gizmo_t::translate(const coot::Cartesian &t) {

   gizmo_origin += t;
   glm::vec3 tt(t.x(), t.y(), t.z());
   mesh.translate(tt);
   // the caller we will need to update the Mesh
}

void
translation_gizmo_t::set_position(const coot::Cartesian &p) {

   glm::vec3 prev_origin(gizmo_origin.x(), gizmo_origin.y(), gizmo_origin.z());
   gizmo_origin = p;
   glm::vec3 pp(p.x(), p.y(), p.z());
   mesh.translate(-prev_origin);
   mesh.translate(pp);

}

void
translation_gizmo_t::scale(float sf) {

   if (sf > 0.0) {
      scale_factor *= sf;
      mesh.scale(sf);
      // we will need to update the Mesh
   }
}

void
translation_gizmo_t::set_scale_absolute(float sf_in) {

   if (sf_in > 0.0) {
      float sf = sf_in/scale_factor;
      scale_factor = sf;
      mesh.scale(sf);
      // we will need to update the Mesh
   }
}
