/*
 * api/molecules-container-raytrace.cc
 *
 * Copyright 2026 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <cmath>
#include <png.h>

#include "molecules-container.hh"
#include "coot-utils/json.hpp"

#ifdef HAVE_OSPRAY
#include <ospray/ospray_cpp.h>
#include <ospray/ospray_cpp/ext/rkcommon.h>
#endif

using json = nlohmann::json;

// ---- local helpers ----

namespace {

bool write_rgba_png(const std::string &file_name, int width, int height, const uint32_t *pixels) {

   FILE *fp = fopen(file_name.c_str(), "wb");
   if (!fp) {
      std::cout << "ERROR:: write_rgba_png(): could not open " << file_name << " for writing" << std::endl;
      return false;
   }

   png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
   if (!png_ptr) {
      fclose(fp);
      return false;
   }

   png_infop info_ptr = png_create_info_struct(png_ptr);
   if (!info_ptr) {
      png_destroy_write_struct(&png_ptr, nullptr);
      fclose(fp);
      return false;
   }

   if (setjmp(png_jmpbuf(png_ptr))) {
      png_destroy_write_struct(&png_ptr, &info_ptr);
      fclose(fp);
      return false;
   }

   png_init_io(png_ptr, fp);
   png_set_IHDR(png_ptr, info_ptr, width, height, 8,
                PNG_COLOR_TYPE_RGBA, PNG_INTERLACE_NONE,
                PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
   png_write_info(png_ptr, info_ptr);

   // OSPRay framebuffer is bottom-to-top, PNG is top-to-bottom
   std::vector<png_bytep> row_pointers(height);
   for (int y=0; y<height; y++) {
      row_pointers[y] = reinterpret_cast<png_bytep>(const_cast<uint32_t *>(&pixels[(height - 1 - y) * width]));
   }

   png_write_image(png_ptr, row_pointers.data());
   png_write_end(png_ptr, nullptr);
   png_destroy_write_struct(&png_ptr, &info_ptr);
   fclose(fp);
   return true;
}

#ifdef HAVE_OSPRAY

// Extract unique edges from a triangle mesh and build OSPRay curve geometry data.
// Each edge becomes a pair of consecutive vec4f entries in curve_vertices (x, y, z, radius)
// and a uint32 index in curve_indices pointing to the first of the pair.
void extract_wireframe_from_mesh(const coot::simple_mesh_t &mesh,
                                 float line_radius,
                                 const rkcommon::math::vec4f &colour,
                                 std::vector<rkcommon::math::vec4f> *curve_vertices_p,
                                 std::vector<uint32_t> *curve_indices_p,
                                 std::vector<rkcommon::math::vec4f> *curve_colors_p) {

   // Use a set of sorted pairs to deduplicate edges
   std::set<std::pair<unsigned int, unsigned int>> edge_set;
   for (const auto &tri : mesh.triangles) {
      unsigned int a = tri.point_id[0];
      unsigned int b = tri.point_id[1];
      unsigned int c = tri.point_id[2];
      edge_set.insert(std::minmax(a, b));
      edge_set.insert(std::minmax(b, c));
      edge_set.insert(std::minmax(a, c));
   }

   auto &cv = *curve_vertices_p;
   auto &ci = *curve_indices_p;
   auto &cc = *curve_colors_p;

   for (const auto &edge : edge_set) {
      const auto &p0 = mesh.vertices[edge.first].pos;
      const auto &p1 = mesh.vertices[edge.second].pos;
      uint32_t base_idx = cv.size();
      cv.push_back(rkcommon::math::vec4f(p0.x, p0.y, p0.z, line_radius));
      cv.push_back(rkcommon::math::vec4f(p1.x, p1.y, p1.z, line_radius));
      ci.push_back(base_idx);
      cc.push_back(colour);
      cc.push_back(colour);
   }
}

#endif // HAVE_OSPRAY

} // anonymous namespace

void
molecules_container_t::ray_trace_init() {

#ifdef HAVE_OSPRAY
   if (!ospray_is_initialized) {
      OSPError err = ospInit(nullptr, nullptr);
      if (err != OSP_NO_ERROR) {
         std::cout << "ERROR:: ray_trace_init(): ospInit() failed with error " << err << std::endl;
      } else {
         ospray_is_initialized = true;
         std::cout << "INFO:: OSPRay initialized" << std::endl;
      }
   }
#else
   std::cout << "WARNING:: ray_trace_init(): Coot was built without OSPRay support" << std::endl;
#endif
}

void
molecules_container_t::ray_trace_shutdown() {

#ifdef HAVE_OSPRAY
   if (ospray_is_initialized) {
      ospShutdown();
      ospray_is_initialized = false;
      std::cout << "INFO:: OSPRay shut down" << std::endl;
   }
#else
   std::cout << "WARNING:: ray_trace_shutdown(): Coot was built without OSPRay support" << std::endl;
#endif
}

void
molecules_container_t::ray_trace_image(const std::string &json_str) {

#ifdef HAVE_OSPRAY

   if (!ospray_is_initialized) {
      std::cout << "ERROR:: ray_trace_image(): OSPRay not initialized. Call ray_trace_init() first." << std::endl;
      return;
   }

   // ---- 0. Parse JSON parameters ----

   json j;
   try {
      j = json::parse(json_str);
   } catch (const json::parse_error &e) {
      std::cout << "ERROR:: ray_trace_image(): JSON parse error: " << e.what() << std::endl;
      return;
   }

   int img_width  = j.value("image_width",  1024);
   int img_height = j.value("image_height", 768);
   std::string output_file_stub = j.value("output_file_stub", "coot-ray-trace");
   int n_accumulation_frames = j.value("n_accumulation_frames", 16);

   rkcommon::math::vec4f bg_colour(1.0f, 1.0f, 1.0f, 1.0f);
   if (j.contains("background_colour") && j["background_colour"].is_array()) {
      auto &bc = j["background_colour"];
      if (bc.size() >= 4)
         bg_colour = rkcommon::math::vec4f(bc[0], bc[1], bc[2], bc[3]);
   }

   if (!j.contains("molecules") || !j["molecules"].is_object()) {
      std::cout << "WARNING:: ray_trace_image(): JSON must contain a \"molecules\" object" << std::endl;
      return;
   }

   // ---- 1. Generate geometry ----

   // Model meshes are triangle geometry; map meshes are wireframe (curve geometry).
   // We collect them separately and add both to the OSPRay group.

   // Combined model mesh (triangles)
   coot::simple_mesh_t combined_model_mesh;

   // Combined map wireframe (curves)
   std::vector<rkcommon::math::vec4f> curve_vertices;
   std::vector<uint32_t> curve_indices;
   std::vector<rkcommon::math::vec4f> curve_colors;

   // Bounding box across all geometry (for auto camera)
   glm::vec3 bb_min( 1e30f,  1e30f,  1e30f);
   glm::vec3 bb_max(-1e30f, -1e30f, -1e30f);
   bool have_model = false;
   glm::vec3 model_bb_min( 1e30f,  1e30f,  1e30f);
   glm::vec3 model_bb_max(-1e30f, -1e30f, -1e30f);

   for (auto &[key, mol_params] : j["molecules"].items()) {
      int imol = std::stoi(key);

      if (is_valid_model_molecule(imol)) {

         std::string colour_mode = mol_params.value("colour_mode", "COLOUR-BY-CHAIN-AND-DICTIONARY");
         float bonds_width = mol_params.value("bonds_width", 0.12f);
         float atom_ratio = mol_params.value("atom_radius_to_bond_width_ratio", 1.5f);
         int smoothness = mol_params.value("smoothness_factor", 3);
         bool draw_hydrogens = mol_params.value("draw_hydrogens", false);
         bool draw_missing_loops = mol_params.value("draw_missing_loops", false);
         bool against_dark_bg = mol_params.value("against_a_dark_background", false);

         coot::instanced_mesh_t im = molecules[imol].get_bonds_mesh_instanced(
            colour_mode, &geom, against_dark_bg,
            bonds_width, atom_ratio,
            true, 0.5f, false, false,
            smoothness, draw_hydrogens, draw_missing_loops);
         coot::simple_mesh_t sm = coot::instanced_mesh_to_simple_mesh(im);

         if (!sm.vertices.empty()) {
            for (const auto &v : sm.vertices) {
               model_bb_min = glm::min(model_bb_min, v.pos);
               model_bb_max = glm::max(model_bb_max, v.pos);
               bb_min = glm::min(bb_min, v.pos);
               bb_max = glm::max(bb_max, v.pos);
            }
            have_model = true;
            combined_model_mesh.add_submesh(sm);
            std::cout << "INFO:: ray_trace_image(): model " << imol << ": "
                      << sm.vertices.size() << " vertices, "
                      << sm.triangles.size() << " triangles" << std::endl;
         }

      } else if (is_valid_map_molecule(imol)) {

         std::string style = mol_params.value("style", "lines");

         float contour_level = mol_params.value("map_contour_level", -999.0f);
         if (contour_level < -900.0f) {
            contour_level = get_suggested_initial_contour_level(imol);
            if (contour_level < 0.0f) contour_level = 1.5f;
         }

         float map_radius = mol_params.value("map_radius", -1.0f);
         float line_width = mol_params.value("map_line_width", 0.02f);

         rkcommon::math::vec4f map_colour(0.3f, 0.5f, 0.8f, 1.0f);
         if (mol_params.contains("map_colour") && mol_params["map_colour"].is_array()) {
            auto &mc = mol_params["map_colour"];
            if (mc.size() >= 3)
               map_colour = rkcommon::math::vec4f(mc[0], mc[1], mc[2], mc.size() >= 4 ? float(mc[3]) : 1.0f);
         }

         // Determine map centre and radius from model bounding box if available
         glm::vec3 map_centre(0.0f, 0.0f, 0.0f);
         if (have_model) {
            map_centre = 0.5f * (model_bb_min + model_bb_max);
            if (map_radius < 0.0f) {
               float extent = glm::length(model_bb_max - model_bb_min);
               map_radius = extent * 0.5f + 5.0f;
            }
         } else {
            if (map_radius < 0.0f) map_radius = 20.0f;
         }

         clipper::Coord_orth position(map_centre.x, map_centre.y, map_centre.z);
         coot::simple_mesh_t map_mesh = molecules[imol].get_map_contours_mesh(
            position, map_radius, contour_level,
            map_is_contoured_using_thread_pool_flag, &thread_pool);

         std::cout << "INFO:: ray_trace_image(): map " << imol << " contoured at " << contour_level
                   << ", radius " << map_radius << ": "
                   << map_mesh.vertices.size() << " vertices, "
                   << map_mesh.triangles.size() << " triangles" << std::endl;

         if (!map_mesh.vertices.empty()) {
            // Update bounding box from map mesh
            for (const auto &v : map_mesh.vertices) {
               bb_min = glm::min(bb_min, v.pos);
               bb_max = glm::max(bb_max, v.pos);
            }

            if (style == "lines") {
               extract_wireframe_from_mesh(map_mesh, line_width, map_colour,
                                           &curve_vertices, &curve_indices, &curve_colors);
               std::cout << "INFO:: ray_trace_image(): map " << imol << " wireframe: "
                         << curve_indices.size() << " line segments" << std::endl;
            } else {
               // Solid surface - add to model mesh
               combined_model_mesh.add_submesh(map_mesh);
            }
         }

      } else {
         std::cout << "WARNING:: ray_trace_image(): molecule " << imol << " is not valid" << std::endl;
      }
   }

   bool have_triangles = !combined_model_mesh.vertices.empty() && !combined_model_mesh.triangles.empty();
   bool have_curves = !curve_vertices.empty();

   if (!have_triangles && !have_curves) {
      std::cout << "WARNING:: ray_trace_image(): no geometry to render" << std::endl;
      return;
   }

   // ---- 2. Build OSPRay geometric models ----

   std::vector<ospray::cpp::GeometricModel> geo_models;

   // Triangle mesh (models and optionally solid maps)
   if (have_triangles) {
      unsigned int n_v = combined_model_mesh.vertices.size();
      unsigned int n_t = combined_model_mesh.triangles.size();

      std::vector<rkcommon::math::vec3f> positions(n_v);
      std::vector<rkcommon::math::vec3f> normals(n_v);
      std::vector<rkcommon::math::vec4f> colors(n_v);
      for (unsigned int i=0; i<n_v; i++) {
         const auto &v = combined_model_mesh.vertices[i];
         positions[i] = rkcommon::math::vec3f(v.pos.x, v.pos.y, v.pos.z);
         normals[i]   = rkcommon::math::vec3f(v.normal.x, v.normal.y, v.normal.z);
         colors[i]    = rkcommon::math::vec4f(v.color.r, v.color.g, v.color.b, v.color.a);
      }

      std::vector<rkcommon::math::vec3ui> indices(n_t);
      for (unsigned int i=0; i<n_t; i++) {
         const auto &t = combined_model_mesh.triangles[i];
         indices[i] = rkcommon::math::vec3ui(t.point_id[0], t.point_id[1], t.point_id[2]);
      }

      ospray::cpp::Geometry mesh("mesh");
      mesh.setParam("vertex.position", ospray::cpp::CopiedData(positions));
      mesh.setParam("vertex.normal",   ospray::cpp::CopiedData(normals));
      mesh.setParam("vertex.color",    ospray::cpp::CopiedData(colors));
      mesh.setParam("index",           ospray::cpp::CopiedData(indices));
      mesh.commit();

      ospray::cpp::Material material("obj");
      material.setParam("kd", rkcommon::math::vec3f(0.8f, 0.8f, 0.8f));
      material.setParam("ks", rkcommon::math::vec3f(0.3f, 0.3f, 0.3f));
      material.setParam("ns", 20.0f);
      material.commit();

      ospray::cpp::GeometricModel model(mesh);
      model.setParam("material", material);
      model.commit();
      geo_models.push_back(model);
   }

   // Curve geometry (wireframe map)
   if (have_curves) {
      ospray::cpp::Geometry curves("curve");
      curves.setParam("vertex.position_radius", ospray::cpp::CopiedData(curve_vertices));
      curves.setParam("vertex.color",           ospray::cpp::CopiedData(curve_colors));
      curves.setParam("index",                  ospray::cpp::CopiedData(curve_indices));
      curves.setParam("type",  OSP_DISJOINT);
      curves.setParam("basis", OSP_LINEAR);
      curves.commit();

      ospray::cpp::Material curve_material("obj");
      curve_material.setParam("ks", rkcommon::math::vec3f(0.1f, 0.1f, 0.1f));
      curve_material.setParam("ns", 5.0f);
      curve_material.commit();

      ospray::cpp::GeometricModel curve_model(curves);
      curve_model.setParam("material", curve_material);
      curve_model.commit();
      geo_models.push_back(curve_model);
   }

   // ---- 3. Assemble scene ----

   ospray::cpp::Group group;
   group.setParam("geometry", ospray::cpp::CopiedData(geo_models));
   group.commit();

   ospray::cpp::Instance instance(group);
   instance.commit();

   ospray::cpp::World world;
   world.setParam("instance", ospray::cpp::CopiedData(instance));

   // Lights
   ospray::cpp::Light ambient_light("ambient");
   ambient_light.setParam("intensity", 0.9f);
   ambient_light.setParam("color", rkcommon::math::vec3f(1.0f, 1.0f, 1.0f));
   ambient_light.commit();

   ospray::cpp::Light dir_light("distant");
   dir_light.setParam("direction", rkcommon::math::vec3f(0.5f, -1.0f, 0.8f));
   dir_light.setParam("intensity", 1.0f);
   dir_light.setParam("color", rkcommon::math::vec3f(1.0f, 1.0f, 1.0f));
   dir_light.commit();

   std::vector<ospray::cpp::Light> lights = {ambient_light, dir_light};
   world.setParam("light", ospray::cpp::CopiedData(lights));
   world.commit();

   // ---- 4. Camera (auto-framing from bounding box) ----

   rkcommon::math::vec3f centre;
   centre.x = 0.5f * (bb_min.x + bb_max.x);
   centre.y = 0.5f * (bb_min.y + bb_max.y);
   centre.z = 0.5f * (bb_min.z + bb_max.z);

   float extent_x = bb_max.x - bb_min.x;
   float extent_y = bb_max.y - bb_min.y;
   float extent_z = bb_max.z - bb_min.z;
   float max_extent = std::max({extent_x, extent_y, extent_z});

   float zoom = j.value("zoom", 1.0f);
   float ortho_height = max_extent / zoom;
   float camera_distance = max_extent * 2.0f;

   // Override camera from JSON if provided (used when not rendering orthogonal views)
   if (j.contains("camera") && j["camera"].is_object()) {
      if (j["camera"].contains("height"))
         ortho_height = j["camera"]["height"];
   }

   // Build the list of views to render.
   // If "orthogonal_views" is true, render front (-Z), side (+X) and top (-Y).
   // Otherwise render a single view (from JSON camera overrides or the default front view).

   struct view_t {
      std::string suffix;
      rkcommon::math::vec3f position;
      rkcommon::math::vec3f direction;
      rkcommon::math::vec3f up;
   };

   std::vector<view_t> views;
   bool orthogonal_views = j.value("orthogonal_views", false);

   if (orthogonal_views) {
      views.push_back({"-front",
         rkcommon::math::vec3f(centre.x, centre.y, centre.z - camera_distance),
         rkcommon::math::vec3f(0.0f, 0.0f, 1.0f),
         rkcommon::math::vec3f(0.0f, 1.0f, 0.0f)});
      views.push_back({"-side",
         rkcommon::math::vec3f(centre.x + camera_distance, centre.y, centre.z),
         rkcommon::math::vec3f(-1.0f, 0.0f, 0.0f),
         rkcommon::math::vec3f(0.0f, 1.0f, 0.0f)});
      views.push_back({"-top",
         rkcommon::math::vec3f(centre.x, centre.y + camera_distance, centre.z),
         rkcommon::math::vec3f(0.0f, -1.0f, 0.0f),
         rkcommon::math::vec3f(0.0f, 0.0f, 1.0f)});
   } else {
      rkcommon::math::vec3f cam_pos(centre.x, centre.y, centre.z - camera_distance);
      rkcommon::math::vec3f cam_dir(0.0f, 0.0f, 1.0f);
      rkcommon::math::vec3f cam_up(0.0f, 1.0f, 0.0f);
      if (j.contains("camera") && j["camera"].is_object()) {
         auto &cam = j["camera"];
         if (cam.contains("position") && cam["position"].is_array() && cam["position"].size() >= 3)
            cam_pos = rkcommon::math::vec3f(cam["position"][0], cam["position"][1], cam["position"][2]);
         if (cam.contains("direction") && cam["direction"].is_array() && cam["direction"].size() >= 3)
            cam_dir = rkcommon::math::vec3f(cam["direction"][0], cam["direction"][1], cam["direction"][2]);
         if (cam.contains("up") && cam["up"].is_array() && cam["up"].size() >= 3)
            cam_up = rkcommon::math::vec3f(cam["up"][0], cam["up"][1], cam["up"][2]);
      }
      views.push_back({"", cam_pos, cam_dir, cam_up});
   }

   // Renderer
   ospray::cpp::Renderer renderer("scivis");
   renderer.setParam("aoSamples", 16);
   renderer.setParam("shadows", true);
   renderer.setParam("backgroundColor", bg_colour);
   renderer.commit();

   float aspect = static_cast<float>(img_width) / static_cast<float>(img_height);

   // ---- 5. Render each view ----

   for (const auto &view : views) {
      ospray::cpp::Camera camera("orthographic");
      camera.setParam("aspect", aspect);
      camera.setParam("position", view.position);
      camera.setParam("direction", view.direction);
      camera.setParam("up", view.up);
      camera.setParam("height", ortho_height);
      camera.commit();

      ospray::cpp::FrameBuffer framebuffer(img_width, img_height, OSP_FB_SRGBA, OSP_FB_COLOR | OSP_FB_ACCUM);
      framebuffer.clear();

      for (int i=0; i<n_accumulation_frames; i++)
         framebuffer.renderFrame(renderer, camera, world);

      std::string view_file = output_file_stub + view.suffix + ".png";
      uint32_t *fb = static_cast<uint32_t *>(framebuffer.map(OSP_FB_COLOR));
      bool png_ok = write_rgba_png(view_file, img_width, img_height, fb);
      framebuffer.unmap(fb);

      if (png_ok)
         std::cout << "INFO:: ray_trace_image(): wrote " << view_file
                   << " (" << img_width << "x" << img_height << ")" << std::endl;
      else
         std::cout << "ERROR:: ray_trace_image(): failed to write " << view_file << std::endl;
   }

#else
   std::cout << "WARNING:: ray_trace_image(): Coot was built without OSPRay support" << std::endl;
#endif
}
