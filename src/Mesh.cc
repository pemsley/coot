
#include <vector>
#include <iostream>
#include <map>
#include <set>

// Having this up here...
#include <mmdb2/mmdb_manager.h>
// Fixes this:
// In file included from /usr/include/X11/Xlib.h:44,
//                 from /usr/include/epoxy/glx.h:36,
//                 from ../../coot/src/Mesh.cc:9:
///home/paule/autobuild/Linux-penelope-pre-release-gtk3-python/include/mmdb2/mmdb_io_file.h:131:22: error: expected unqualified-id before numeric constant
//  131 |         inline bool  Success   () { return IOSuccess; }
//      |                      ^~~~~~~


#include <epoxy/gl.h>
#include "Mesh.hh"
#include "Shader.hh"
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>  // to_string()
#include <glm/gtx/rotate_vector.hpp>

#include "cylinder.hh"

#ifdef THIS_IS_HMT
// Get rid of this at some stage - needs thinking about meshed objects
#include "generic-display-object.hh"
#else
#include "old-generic-display-object.hh"
#endif

void
Mesh::init() {

   clear();
   first_time = true;
   is_instanced = false;
   is_instanced_colours = false;
   is_instanced_with_rts_matrix = false;
   use_blending = false;
   draw_this_mesh = true;
   hydrogen_bond_cylinders_angle = 0.0;
   normals_are_setup = false;
   this_mesh_is_closed = false;
   n_instances = 0;
   n_instances_allocated = 0;
   particle_draw_count = 0;
   vao = VAO_NOT_SET; // use UNSET_VAO
}

Mesh::Mesh(const std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > &indexed_vertices) {

   init();
   vertices  = indexed_vertices.first;
   triangles = indexed_vertices.second;
}

// a molecular_triangles_mesh_t is a poor man's Mesh. Why does it exist?
Mesh::Mesh(const molecular_triangles_mesh_t &mtm) {

   init();
   vertices  = mtm.vertices;
   triangles = mtm.triangles;
   name = mtm.name;
}

void
Mesh::close() {

   clear();
   draw_this_mesh = false;
   this_mesh_is_closed = true;
}

void
Mesh::import(const std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > &indexed_vertices) {

   is_instanced = false;
   is_instanced_colours = false;
   is_instanced_with_rts_matrix = false;
   use_blending = false;

   unsigned int idx_base = vertices.size();
   unsigned int idx_tri_base = triangles.size();
   vertices.insert(vertices.end(), indexed_vertices.first.begin(), indexed_vertices.first.end());
   triangles.insert(triangles.end(),
                    indexed_vertices.second.begin(),
                    indexed_vertices.second.end());
   for (unsigned int i=idx_tri_base; i<triangles.size(); i++)
      triangles[i].rebase(idx_base);
}

void
Mesh::import(const std::vector<s_generic_vertex> &gv, const std::vector<g_triangle> &indexed_vertices) {

   is_instanced = false;
   is_instanced_colours = false;
   is_instanced_with_rts_matrix = false;

   unsigned int idx_base = vertices.size();
   unsigned int idx_tri_base = triangles.size();
   vertices.insert(vertices.end(), gv.begin(), gv.end());
   triangles.insert(triangles.end(),
                    indexed_vertices.begin(),
                    indexed_vertices.end());
   for (unsigned int i=idx_tri_base; i<triangles.size(); i++)
      triangles[i].rebase(idx_base);
}

void
Mesh::import(const std::vector<position_normal_vertex> &verts, const std::vector<g_triangle> &indexed_vertices,
             const glm::vec4 &colour) {

   is_instanced = false;
   is_instanced_with_rts_matrix = false;

   unsigned int idx_base = vertices.size();
   unsigned int idx_tri_base = triangles.size();

   std::vector<s_generic_vertex> gv(verts.size());
   for (unsigned int ii=0; ii<verts.size(); ii++) {
      gv[ii].pos    = verts[ii].pos;
      gv[ii].normal = verts[ii].normal;
      gv[ii].color  = colour;
   }

   vertices.insert(vertices.end(), gv.begin(), gv.end());
   triangles.insert(triangles.end(),
                    indexed_vertices.begin(),
                    indexed_vertices.end());
   if (idx_base != 0)
      for (unsigned int i=idx_tri_base; i<triangles.size(); i++)
         triangles[i].rebase(idx_base);

   // now the caller function should call setup(shader, material) so that setup_buffers() is run

}

void
Mesh::translate_by(const glm::vec3 &t) {

   for (unsigned int ii=0; ii<vertices.size(); ii++)
      vertices[ii].pos += t;

   setup_buffers();
}



void
Mesh::debug() const {

   std::cout << "Mesh::debug() " << name << " n-vertices " << vertices.size()
             << " n-triangles " << triangles.size() << std::endl;
}


void
Mesh::setup(Shader *shader_p, const Material &material_in) {

   // Generic objects need a lot of reworking
   // because I want to "setup" only after all objects have been added - not
   // after adding every object (say 1000 spheres/dots).
   // But I don't want to draw if the object is not setup. Hmm.
   // Also adding thousands of balls is slow!
   // HOLE dots is a good example of where I should be using instancing

   // if (setup_has_been_done)
   // return;

   material = material_in;
   shader_p->Use();
   setup_buffers();
}


// maybe should not be const
void
Mesh::setup_simple_triangles(Shader *shader_p, const Material &material_in) {

   material = material_in;
   shader_p->Use();
   fill_with_simple_triangles_vertices();
   fill_with_direction_triangles();
   setup_buffers();

}

void
Mesh::setup_rama_balls(Shader *shader_p, const Material &material_in) {

   material = material_in;
   shader_p->Use();
   fill_rama_balls();
   setup_buffers();
}


void
Mesh::fill_one_dodec() {

   // These are dodec balls for the moment, not flat shaded real dodecs.
   // In a bit I will transform them to pentakis dodecs.

   dodec d;
   std::vector<clipper::Coord_orth> v = d.coords();
   vertices.resize(v.size()); // should be 20

   for (unsigned int i=0; i<v.size(); i++) {
      glm::vec3 a = glm::vec3(v[i].x(), v[i].y(), v[i].z());
      float f = static_cast<float>(i)/static_cast<float>(v.size());
      vertices[i].pos    = 0.2f * a;
      vertices[i].normal = a;
      vertices[i].color  = glm::vec4(0.3 + 0.4 * f, 0.8 - 0.8 * f, 0.1 + 0.9 * f, 1.0);
   }

   for (unsigned int i=0; i<12; i++) {
      const std::vector<unsigned int> &face = d.face(i);
      g_triangle gt_0(face[0], face[1], face[2]);
      g_triangle gt_1(face[0], face[2], face[3]);
      g_triangle gt_2(face[0], face[3], face[4]);
      triangles.push_back(gt_0);
      triangles.push_back(gt_1);
      triangles.push_back(gt_2);
   }

}

#ifdef THIS_IS_HMT
#else
#include "old-generic-display-object.hh"
#endif


void
Mesh::add_one_ball(float scale, const glm::vec3 &centre) { // i.e. a smooth-shaded pentakis dodec

   pentakis_dodec pkdd(1.05);
   // coordninate system finagalling - baah - convert dodecs to glm.
   clipper::Coord_orth centre_c(centre.x, centre.y, centre.z);
   coot::old_generic_display_object_t::pentakis_dodec_t penta_dodec(pkdd, 0.01, centre_c);
   std::vector<clipper::Coord_orth> v = penta_dodec.pkdd.d.coords();

   unsigned int vertex_index_start_base   = vertices.size();

   const std::vector<clipper::Coord_orth> &pv = penta_dodec.pkdd.pyrimid_vertices;
   for (unsigned int i=0; i<12; i++) {

      const std::vector<unsigned int> &face = penta_dodec.pkdd.d.face(i);

      // first the base point (tip of the triangles/pyrimid)
      //
      clipper::Coord_orth pvu(pv[i].unit());
      s_generic_vertex gv;
      gv.pos    = scale * (glm::vec3(pv[i].x(), pv[i].y(), pv[i].z())) + centre;
      gv.normal = glm::vec3(pvu.x(), pvu.y(), pvu.z());
      vertices.push_back(gv);

      for (unsigned int j=0; j<5; j++) {
         const clipper::Coord_orth &pt = v[face[j]];
         clipper::Coord_orth ptu(pt.unit());
         gv.pos    = scale * (glm::vec3(pt.x(),  pt.y(),  pt.z())) + centre;
         gv.normal = glm::vec3(ptu.x(), ptu.y(), ptu.z());
         vertices.push_back(gv);
      }
      unsigned int idx_base = vertex_index_start_base + i*6;
      g_triangle gt_0(idx_base, idx_base + 1, idx_base + 2);
      g_triangle gt_1(idx_base, idx_base + 2, idx_base + 3);
      g_triangle gt_2(idx_base, idx_base + 3, idx_base + 4);
      g_triangle gt_3(idx_base, idx_base + 4, idx_base + 5);
      g_triangle gt_4(idx_base, idx_base + 5, idx_base + 1);
      triangles.push_back(gt_0);
      triangles.push_back(gt_1);
      triangles.push_back(gt_2);
      triangles.push_back(gt_3);
      triangles.push_back(gt_4);
   }
}

void
Mesh::add_one_origin_ball() { // i.e. a smooth-shaded pentakis dodec

   pentakis_dodec pkdd(1.05);
   // coordninate system finagalling - baah - convert dodecs to glm.
   clipper::Coord_orth centre_c(0,0,0);
   coot::old_generic_display_object_t::pentakis_dodec_t penta_dodec(pkdd, 0.01, centre_c);
   std::vector<clipper::Coord_orth> v = penta_dodec.pkdd.d.coords();

   unsigned int vertex_index_start_base   = vertices.size();
   unsigned int triangle_index_start_base = triangles.size(); // needed?

   const std::vector<clipper::Coord_orth> &pv = penta_dodec.pkdd.pyrimid_vertices;
   for (unsigned int i=0; i<12; i++) {

      const std::vector<unsigned int> &face = penta_dodec.pkdd.d.face(i);
      float scale = 0.5;

      // first the base point (tip of the triangles/pyrimid)
      //
      clipper::Coord_orth pvu(pv[i].unit());
      s_generic_vertex gv;
      gv.pos    = scale * 0.1f * (glm::vec3(pv[i].x(), pv[i].y(), pv[i].z()));
      gv.normal = glm::vec3(pvu.x(), pvu.y(), pvu.z());
      vertices.push_back(gv);

      for (unsigned int j=0; j<5; j++) {
         const clipper::Coord_orth &pt = v[face[j]];
         clipper::Coord_orth ptu(pt.unit());
         gv.pos    = scale * 0.1f * (glm::vec3(pt.x(),  pt.y(),  pt.z()));
         gv.normal = glm::vec3(ptu.x(), ptu.y(), ptu.z());
         vertices.push_back(gv);
      }
      unsigned int idx_base = vertex_index_start_base + i*6;
      g_triangle gt_0(idx_base, idx_base + 1, idx_base + 2);
      g_triangle gt_1(idx_base, idx_base + 2, idx_base + 3);
      g_triangle gt_2(idx_base, idx_base + 3, idx_base + 4);
      g_triangle gt_3(idx_base, idx_base + 4, idx_base + 5);
      g_triangle gt_4(idx_base, idx_base + 5, idx_base + 1);
      triangles.push_back(gt_0);
      triangles.push_back(gt_1);
      triangles.push_back(gt_2);
      triangles.push_back(gt_3);
      triangles.push_back(gt_4);
   }
}


void
Mesh::add_one_origin_dodec() { // i.e. a smooth-shaded pentakis dodec

   if (false) {
      glm::vec3 p0(0,0,0);
      glm::vec3 p1(1,0,0);
      glm::vec3 p2(1,1,0);
      glm::vec3 n(0,0,0);
      glm::vec4 c(1,0,0,1);
      s_generic_vertex gv0(p0,n,c);
      s_generic_vertex gv1(p1,n,c);
      s_generic_vertex gv2(p2,n,c);
      vertices.push_back(gv0);
      vertices.push_back(gv1);
      vertices.push_back(gv2);
      triangles.push_back(g_triangle(0,1,2));
   }

   if (true) {
      float scale = 0.05;
      dodec d;
      float radius = 0.2;
      clipper::Coord_orth position(0,0,0);
      coot::old_generic_display_object_t::dodec_t dod(d, radius, position);
      // we can put a colour into dod here
   
      std::vector<clipper::Coord_orth> v = dod.d.coords();

      unsigned int vertex_index_start_base = vertices.size();
      for (unsigned int i=0; i<12; i++) {
         const std::vector<unsigned int> &face = dod.d.face(i);
         clipper::Coord_orth sum_vertex(0,0,0);
         for (unsigned int j=0; j<5; j++)
            sum_vertex += v[face[j]];
         clipper::Coord_orth face_normal(sum_vertex.unit());
         for (unsigned int j=0; j<5; j++) {
            glm::vec3 p(v[face[j]].x(), v[face[j]].y(), v[face[j]].z());
            glm::vec3 n(face_normal.x(), face_normal.y(), face_normal.z());
            glm::vec4 c(1,0,0,1);
            s_generic_vertex gv(scale * p,n,c);
            vertices.push_back(gv);
         }
         unsigned int idx_base = vertex_index_start_base + i*5;
         g_triangle gt_0(idx_base, idx_base + 1, idx_base + 2);
         g_triangle gt_1(idx_base, idx_base + 2, idx_base + 3);
         g_triangle gt_2(idx_base, idx_base + 3, idx_base + 4);
         triangles.push_back(gt_0);
         triangles.push_back(gt_1);
         triangles.push_back(gt_2);
      }
   }
}

void
Mesh::add_one_origin_cylinder(unsigned int n_slices, unsigned int n_stacks) {

   // short fat, radius 1, height 1.

   cylinder c(std::pair<glm::vec3, glm::vec3>(glm::vec3(0,0,0), glm::vec3(0,0,1)),
              1.0, 1.0, 1.0, n_slices, n_stacks);

   // vertices = std::move(c.vertices);
   // triangle_vertex_indices = std::move(c.triangle_indices_vec);

   unsigned int idx_base = vertices.size();
   unsigned int idx_tri_base = triangles.size();
   vertices.insert(vertices.end(), c.vertices.begin(), c.vertices.end());
   triangles.insert(triangles.end(),
                    c.triangle_indices_vec.begin(),
                    c.triangle_indices_vec.end());
   for (unsigned int i=idx_tri_base; i<triangles.size(); i++)
      triangles[i].rebase(idx_base);

}

#include "oct.hh"

void
Mesh::add_one_origin_octahemisphere(unsigned int num_subdivisions) {

   std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > oct =
      tessellate_hemisphere_patch(num_subdivisions);

   float turn_fraction = 0.125; // default num_subdivisions = 1
   float angle = 2.0f * M_PI * turn_fraction;
   glm::vec3 z_axis(0,0,1);
   glm::vec4 atom_colour(0.2f, 0.6f, 0.3f, 1.0f);

   vertices.resize(oct.first.size());
   for (unsigned int i=0; i<oct.first.size(); i++) {
      vertices[i].pos    = glm::rotate(oct.first[i], angle, z_axis);
      vertices[i].normal = glm::rotate(oct.first[i], angle, z_axis);
      vertices[i].color  = atom_colour;
   }
   triangles = oct.second;

}

void
Mesh::add_one_origin_octasphere(unsigned int num_subdivisions) {

   // what is this?

   std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > oct =
      tessellate_hemisphere_patch(num_subdivisions);

   float turn_fraction = 0.125; // default num_subdivisions = 1
   float angle = 2.0f * M_PI * turn_fraction;
   glm::vec3 z_axis(0,0,1);
   glm::vec4 atom_colour(0.2f, 0.6f, 0.3f, 1.0f);

   vertices.resize(oct.first.size());
   for (unsigned int i=0; i<oct.first.size(); i++) {
      vertices[i].pos    = glm::rotate(oct.first[i], angle, z_axis);
      vertices[i].normal = glm::rotate(oct.first[i], angle, z_axis);
      vertices[i].color  = atom_colour;
   }
   triangles = oct.second;

}

void
Mesh::fill_rama_balls() {

   // glm::vec3 centre(0,0,0.0);
   // add_one_ball(centre);

   for (unsigned int i=0; i<6; i++) {
      float f = static_cast<float>(i)/10.0;
      glm::vec3 c(0.44, 0.05, f - 0.2);
      float scale = 0.04;
      add_one_ball(scale, c);
   }
}

void
Mesh::setup_debugging_instancing_buffers() {

   // put these in their own file - make inst_matrices be member data

   is_instanced = true;
   is_instanced_colours = true;
   std::vector<glm::vec3> inst_trans_matrices;
   inst_trans_matrices.push_back(glm::vec3(0.25, 0.25, -0.2));
   inst_trans_matrices.push_back(glm::vec3(0.25, 0.25, -0.1));
   inst_trans_matrices.push_back(glm::vec3(0.25, 0.25,  0.0));
   inst_trans_matrices.push_back(glm::vec3(0.25, 0.25,  0.1));
   inst_trans_matrices.push_back(glm::vec3(0.25, 0.25,  0.2));
   inst_trans_matrices.push_back(glm::vec3(0.25, 0.25,  0.3));

   std::vector<glm::vec4> inst_col_matrices;
   inst_col_matrices.push_back(glm::vec4(0.8, 0.0, 0.0, 1.0));
   inst_col_matrices.push_back(glm::vec4(0.6, 0.3, 0.0, 1.0));
   inst_col_matrices.push_back(glm::vec4(0.5, 0.5, 0.0, 1.0));
   inst_col_matrices.push_back(glm::vec4(0.3, 0.6, 0.0, 1.0));
   inst_col_matrices.push_back(glm::vec4(0.0, 0.9, 0.0, 1.0));
   inst_col_matrices.push_back(glm::vec4(0.0, 0.7, 1.0, 1.0));

   if (false) { // debugging. one red ball at the origin.
      inst_trans_matrices.clear();
      inst_col_matrices.clear();
      inst_col_matrices.push_back(glm::vec4(0.8, 0.0, 0.0, 1.0));
      inst_trans_matrices.push_back(glm::vec3(0,0,0));
   }

   n_instances = inst_trans_matrices.size();
   n_instances_allocated = n_instances;

   // has the vao been genvertexarrayed before this is called?
   glBindVertexArray(vao);

   glGenBuffers(1, &inst_colour_buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, inst_colour_buffer_id);
   glBufferData(GL_ARRAY_BUFFER, n_instances * sizeof(glm::vec4), &(inst_col_matrices[0]), GL_STATIC_DRAW);
   glEnableVertexAttribArray(2);
   glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(glm::vec4), 0);
   glVertexAttribDivisor(2, 1);

   glGenBuffers(1, &inst_model_translation_buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, inst_model_translation_buffer_id);
   glBufferData(GL_ARRAY_BUFFER, n_instances * sizeof (glm::vec3), &(inst_trans_matrices[0]), GL_STATIC_DRAW);
   glEnableVertexAttribArray(3);
   glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), 0);
   glVertexAttribDivisor(3, 1);

   glBindVertexArray(0);
}

void
Mesh::setup_buffers() {

   if (triangles.empty()) return;
   if (vertices.empty()) return;

   if (first_time) {
      glGenVertexArrays (1, &vao);
      // std::cout << "############## first time: generated VAO " << vao << std::endl;
      // don't return before we set first_time = false at the end
   }

   // std::cout << "Mesh::setup_buffers() using vao " << vao << std::endl;
   glBindVertexArray(vao);
   GLenum err = glGetError();
   if (err) std::cout << "error setup_buffers() on binding vao " << vao << " error "
                      << err << std::endl;

   unsigned int n_vertices = vertices.size();

   if (first_time) {
      glGenBuffers(1, &buffer_id);
      glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
      glBufferData(GL_ARRAY_BUFFER, n_vertices * sizeof(s_generic_vertex), &(vertices[0]), GL_STATIC_DRAW);
   } else {
      glDeleteBuffers(1, &buffer_id);
      glGenBuffers(1, &buffer_id);
      glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
      glBufferData(GL_ARRAY_BUFFER, n_vertices * sizeof(s_generic_vertex), &(vertices[0]), GL_STATIC_DRAW);
   }

   // position
   glEnableVertexAttribArray(0);
   glVertexAttribPointer (0, 3, GL_FLOAT, GL_FALSE, sizeof(s_generic_vertex), 0);

   // normal
   glEnableVertexAttribArray (1);
   glVertexAttribPointer (1, 3, GL_FLOAT, GL_FALSE, sizeof(s_generic_vertex),
                          reinterpret_cast<void *>(sizeof(glm::vec3)));

   // colour
   if (! is_instanced_colours) {
      // when the colour attribute is instanced (like the model_translation), so set up the colours
      // when we do the model_translation in setup_instancing_buffers(), not here.
      glEnableVertexAttribArray(2);
      glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(s_generic_vertex),
                            reinterpret_cast<void *>(2 * sizeof(glm::vec3)));
   }

   unsigned int n_triangles = triangles.size();
   unsigned int n_bytes = n_triangles * 3 * sizeof(unsigned int);

   if (first_time) {
      glGenBuffers(1, &index_buffer_id);
      err = glGetError(); if (err) std::cout << "GL error setup_simple_triangles()\n";
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
      err = glGetError(); if (err) std::cout << "GL error setup_simple_triangles()\n";
   } else {
      glDeleteBuffers(1, &index_buffer_id);
      glGenBuffers(1, &index_buffer_id);
      err = glGetError(); if (err) std::cout << "GL error setup_simple_triangles()\n";
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
      err = glGetError(); if (err) std::cout << "GL error setup_simple_triangles()\n";
   }

   if (false)
      std::cout << "debug:: glBufferData for index buffer_id " << index_buffer_id
                << " n_triangles: " << n_triangles
                << " allocating with size: " << n_bytes << " bytes" << std::endl;

   glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_bytes, &triangles[0], GL_STATIC_DRAW);
   err = glGetError(); if (err) std::cout << "GL error setup_simple_triangles()\n";

   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   if (! is_instanced_colours) glDisableVertexAttribArray(2);
   glBindBuffer(GL_ARRAY_BUFFER, 0);
   glUseProgram(0);
   glBindVertexArray(0);

   first_time = false;

}

void
Mesh::setup_instanced_debugging_objects(Shader *shader_p, const Material &material_in) {

   material = material_in; // not currently used in the shader.
   shader_p->Use();
   unsigned int n_vert = vertices.size();
   add_one_origin_dodec(); // flat-faced for debugging
   // move the dodec
   for (unsigned int i=n_vert; i<vertices.size(); i++) {
      s_generic_vertex &v = vertices[i];
      v.pos += glm::vec3(-0.05,0.15,0);
   }
   add_one_origin_ball();
   setup_buffers(); // sets the vao
   setup_debugging_instancing_buffers();

}

void
Mesh::setup_instanced_dodecs(Shader *shader_p, const Material &material_in) {

   is_instanced = true;
   is_instanced_colours =  true; //  I think.
   material = material_in; // not currently used in the shader.
   shader_p->Use();
   add_one_origin_dodec(); // flat-faced for debugging
   setup_buffers(); // sets the vao

   setup_debugging_instancing_buffers(); // a few translated objects

}

// rts rotation, translation & scale
void
Mesh::import_and_setup_instanced_cylinders(Shader *shader_p,
                                           const Material &material_in,
                                           const std::vector<glm::mat4> &mats,
                                           const std::vector<glm::vec4> &colours) {

   GLenum err = glGetError();
   if (err) std::cout << "error import_and_setup_instanced_cylinders() -- start -- "
                      << err << std::endl;
   is_instanced = true;
   is_instanced_colours = true;
   is_instanced_with_rts_matrix = true;

   shader_p->Use();
   material = material_in; // not currently used in the shader.
   add_one_origin_cylinder(16);
   setup_buffers();

   n_instances = mats.size();
   n_instances_allocated = n_instances;
   setup_matrix_and_colour_instancing_buffers(mats, colours);
   err = glGetError(); if (err) std::cout << "error import_and_setup_instanced_cylinders() -- end -- "
                                          << err << std::endl;

}

void
Mesh::setup_rtsc_instancing(Shader *shader_p,
                            const std::vector<glm::mat4> &mats,
                            const std::vector<glm::vec4> &colours,
                            unsigned int n_instances_in,
                            const Material &material_in) {

   is_instanced = true;
   is_instanced_colours = true;
   is_instanced_with_rts_matrix = true;

   shader_p->Use();
   material = material_in;

   setup_buffers();
   n_instances = n_instances_in;
   n_instances_allocated = n_instances;

   setup_matrix_and_colour_instancing_buffers(mats, colours);
   GLenum err = glGetError(); if (err) std::cout << "   error setup_instanced_cylinders() -- end -- "
                                                 << err << std::endl;

   std::cout << "setup_rtsc_instancing(): " << vertices.size() << " vertices" << std::endl;
   std::cout << "setup_rtsc_instancing(): " << triangles.size()
             << " triangles" << std::endl;
}

void
Mesh::setup_instanced_octahemispheres(Shader *shader_p,
                                                    const Material &material_in,
                                                    const std::vector<glm::mat4> &mats,
                                                    const std::vector<glm::vec4> &colours) {

   GLenum err = glGetError(); if (err) std::cout << "   error setup_instanced_octahemispheres() "
                                                 << " -- start -- " << err << std::endl;
   is_instanced = true;
   is_instanced_colours = true;
   is_instanced_with_rts_matrix = true;

   material = material_in; // not currently used in the shader.
   shader_p->Use();
   add_one_origin_octahemisphere(2);
   setup_buffers();
   n_instances = mats.size();
   n_instances_allocated = n_instances;

   setup_matrix_and_colour_instancing_buffers(mats, colours); // maybe pass a flag MATS_AND_COLOURS
                                            // because we might have
                                            // other instanced geometry that doesn't change colour.
                                            // How about HOLE balls?

   err = glGetError(); if (err) std::cout << "   error setup_instanced_octahemispheres() -- end -- "
                                          << err << std::endl;

}


// instancing buffer for particles. Make *space* for n_particles, but set
// n_instances = 0.
void
Mesh::setup_vertex_and_instancing_buffers_for_particles(unsigned int n_particles) {

   // we want to allocate space for n_particles instances, but until the particles
   // are created, and the particle bufffer data updated, we don't want to draw
   // them, so n_instances is 0.
   //
   n_instances = 0;
   n_instances_allocated = n_particles;
   particle_draw_count = 0;

   // glm::vec3 n(0,0,1);
   // glm::vec4 c(0.8, 0.4, 0.8, 0.8);

   setup_camera_facing_polygon(5, 0.3); // calls setup_buffers() for the vertices and sets the VAO

   glBindVertexArray(vao);
   GLenum err = glGetError();
   if (err) std::cout << "GL error ####"
                      << " setup_vertex_and_instancing_buffers_for_particles() B "
                      << err << std::endl;


   // a Particle has position, velocity and colour. We need position and colour

   // I shouldn't need to make 2 buffers (ie. 2 calls to glBufferData) here!
   // Look at how the Mesh for ribbons does it.

   // instanced position

   unsigned int n_bytes = n_particles * sizeof(Particle);

   glGenBuffers(1, &inst_model_translation_buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, inst_model_translation_buffer_id);
   glBufferData(GL_ARRAY_BUFFER, n_bytes, nullptr, GL_DYNAMIC_DRAW);
   glEnableVertexAttribArray(3);
   glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(Particle), 0);
   glVertexAttribDivisor(3, 1);
   err = glGetError();
   if (err) std::cout << "GL error #####"
                      << " setup_instancing_buffers_for_particles() B "
                      << err << std::endl;

   // instanced colours - setup another buffer - extravagent.
   glGenBuffers(1, &inst_colour_buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, inst_colour_buffer_id);
   glBufferData(GL_ARRAY_BUFFER, n_particles * sizeof(Particle), nullptr, GL_DYNAMIC_DRAW);
   glEnableVertexAttribArray(4);
   glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, sizeof(Particle),
                         reinterpret_cast<void *>(2 * sizeof(glm::vec3)));
   glVertexAttribDivisor(4, 1);

   // index the quad/hex/polygon

   glGenBuffers(1, &index_buffer_id);
   err = glGetError(); if (err) std::cout << "GL error setup_instancing_buffers_for_particles()\n";
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
   err = glGetError(); if (err) std::cout << "GL error setup_instancing_buffers_for_particles()\n";
   unsigned int n_triangles = triangles.size();
   n_bytes = n_triangles * 3 * sizeof(unsigned int);
   if (false)
      std::cout << "debug:: setup_vertex_and_instancing_buffers_for_particles() "
                << "vao " << vao
                << " glBufferData for index buffer_id " << index_buffer_id
                << " n_triangles: " << n_triangles
                << " allocating with size: " << n_bytes << " bytes for triangles"
                << std::endl;
   glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_bytes, &triangles[0], GL_DYNAMIC_DRAW);
   err = glGetError(); if (err) std::cout << "GL error setup_instancing_buffers_for_particles()\n";

   err = glGetError();
   if (err) std::cout << "GL error #####"
                      << " setup_vertex_and_instancing_buffers_for_particles() --- end --- "
                      << err << std::endl;
   glBindVertexArray(0);

}


void
Mesh::setup_matrix_and_colour_instancing_buffers(const std::vector<glm::mat4> &mats,
                                                 const std::vector<glm::vec4> &colours) {

   GLenum err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() -- start -- "
                      << err << std::endl;

   n_instances = mats.size();
   n_instances_allocated = n_instances;

   glBindVertexArray(vao);
   err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() B binding-vao "
                      << err << " with vao " << vao << std::endl;


   // ------------- colours ----------------------------
   glGenBuffers(1, &inst_colour_buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, inst_colour_buffer_id);
   err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() B0 "
                      << err << std::endl;
   glBufferData(GL_ARRAY_BUFFER, n_instances * sizeof(glm::vec4), &(colours[0]), GL_DYNAMIC_DRAW);
   glEnableVertexAttribArray(3);
   err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() B1 "
                      << err << std::endl;
   glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, sizeof(glm::vec4), 0);
   err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() B2 "
                      << err << std::endl;
   glVertexAttribDivisor(3, 1);
   err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() B3 "
                      << err << std::endl;

   // ------------- orientation/position/scale matrices ----------------------------

   glGenBuffers(1, &inst_rts_buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, inst_rts_buffer_id);
   err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() C0 "
                      << err << std::endl;
   glBufferData(GL_ARRAY_BUFFER, n_instances * sizeof(glm::mat4), &(mats[0]), GL_DYNAMIC_DRAW);
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() C1 " << err << std::endl;

   glEnableVertexAttribArray(4);
   glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), 0);
   glVertexAttribDivisor(4, 1);
   err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() C3 " << err << std::endl;
   glEnableVertexAttribArray(5);
   glVertexAttribPointer(5, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), reinterpret_cast<void *>(sizeof(glm::vec4)));
   glVertexAttribDivisor(5, 1);
   err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() C4 " << err << std::endl;
   glEnableVertexAttribArray(6);
   glVertexAttribPointer(6, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), reinterpret_cast<void *>(2 * sizeof(glm::vec4)));
   glVertexAttribDivisor(6, 1);
   err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() C5 " << err << std::endl;
   glEnableVertexAttribArray(7);
   glVertexAttribPointer(7, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), reinterpret_cast<void *>(3 * sizeof(glm::vec4)));
   glVertexAttribDivisor(7, 1);
   err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() C6 " << err << std::endl;
}


void
Mesh::setup_matrix_and_colour_instancing_buffers_old(const std::vector<glm::mat4> &mats,
                                                     const std::vector<glm::vec4> &colours) {

   // this function doesn't make sense in it's current form.
   // Instead, I need to *make space* for n_mats and n_colours.
   //
   // I don't need matrices with values here - because they will be updated/replaced by
   // the tick function. So this function should be changed to pass the size of
   // matrix and colour vectors, not the actual values.

   std::cout << "----- setup_instancing_buffers(): mats size " << mats.size()
             << " colours size " << colours.size() << std::endl;

   GLenum err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() -- start -- "
                      << err << std::endl;

   n_instances = mats.size();
   n_instances_allocated = n_instances;

   std::vector<glm::mat4> inst_rts_matrices = mats;
   std::vector<glm::vec4> inst_col_matrices = colours;

   err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() A "
                      << err << std::endl;

   glBindVertexArray(vao);

   err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() B binding-vao "
                      << err << " with vao " << vao << std::endl;

   // rama balls we want to have instanced colours but hydrogen bond rotating cylinders we do not.

   // -------- colours -----------

   if (is_instanced_colours) {
      glGenBuffers(1, &inst_colour_buffer_id);
      glBindBuffer(GL_ARRAY_BUFFER, inst_colour_buffer_id);
      err = glGetError();
      if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() B0 "
                         << err << std::endl;
      std::cout << "setup_matrix_and_colour_instancing_buffers() allocating colour buffer data "
                << n_instances * sizeof(glm::vec4) << std::endl;
      glBufferData(GL_ARRAY_BUFFER, n_instances * sizeof(glm::vec4), &(inst_col_matrices[0]), GL_DYNAMIC_DRAW); // dynamic
      glEnableVertexAttribArray(2);
      err = glGetError();
      if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() B1 "
                         << err << std::endl;
      glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(glm::vec4), 0);
      err = glGetError();
      if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() B2 "
                         << err << std::endl;
      glVertexAttribDivisor(2, 1);
      err = glGetError();
      if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() B3 "
                         << err << std::endl;
   }

   // -------- rotation/translation/scale matrices -----------

   glGenBuffers(1, &inst_rts_buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, inst_rts_buffer_id);
   std::cout << "setup_matrix_and_colour_instancing_buffers() allocating matrix buffer data "
             << n_instances * 4 * sizeof(glm::mat4) << std::endl;
   glBufferData(GL_ARRAY_BUFFER, n_instances * 4 * sizeof (glm::vec4), &(inst_rts_matrices[0]), GL_DYNAMIC_DRAW); // dynamic

   err = glGetError(); if (err) std::cout << "   error setup_instancing_buffers() C1 "
                                          << err << std::endl;
   glEnableVertexAttribArray(3);
   glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), 0);
   glVertexAttribDivisor(3, 1);
   err = glGetError(); if (err) std::cout << "   error setup_instancing_buffers() C2 " << err << std::endl;
   glEnableVertexAttribArray(4);
   glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), reinterpret_cast<void *>(sizeof(glm::vec4)));
   glVertexAttribDivisor(4, 1);
   err = glGetError(); if (err) std::cout << "   error setup_instancing_buffers() C2 " << err << std::endl;
   glEnableVertexAttribArray(5);
   glVertexAttribPointer(5, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), reinterpret_cast<void *>(2 * sizeof(glm::vec4)));
   glVertexAttribDivisor(5, 1);
   err = glGetError(); if (err) std::cout << "   error setup_instancing_buffers() C3 " << err << std::endl;
   glEnableVertexAttribArray(6);
   glVertexAttribPointer(6, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), reinterpret_cast<void *>(3 * sizeof(glm::vec4)));
   glVertexAttribDivisor(6, 1);

   err = glGetError(); if (err) std::cout << "   error setup_instancing_buffers() E " << err << std::endl;
   
}


// void
// Mesh::setup_instanced_balls(Shader *shader_p, const Material &material_in) {

// }




void
Mesh::fill_with_simple_triangles_vertices() {

   float scale = 0.4;

   std::vector<s_generic_vertex> v(6);
   float z_sep = 0.4;

   v[0].pos = scale * glm::vec3(  0.0f,  0.5f,  -z_sep);
   v[1].pos = scale * glm::vec3(  0.5f, -0.36f, -z_sep);
   v[2].pos = scale * glm::vec3( -0.5f, -0.36f, -z_sep);
   v[3].pos = scale * glm::vec3(  0.0f,  0.5f,   z_sep);
   v[4].pos = scale * glm::vec3(  0.5f, -0.36f,  z_sep);
   v[5].pos = scale * glm::vec3( -0.5f, -0.36f,  z_sep);

   v[0].normal = glm::vec3( 0.2f, 0.2f,  0.9f);
   v[1].normal = glm::vec3( 0.2f, 0.9f,  0.2f);
   v[2].normal = glm::vec3( 0.9f, 0.3f,  0.1f);
   v[3].normal = glm::vec3( 0.0f, 0.9f, -0.1f);
   v[4].normal = glm::vec3( 0.9f, 0.3f, -0.2f);
   v[5].normal = glm::vec3( 0.2f, 0.6f, -0.9f);

   v[0].color = glm::vec4(0.0f, 0.0f, 0.0f, 1.f);
   v[1].color = glm::vec4(0.2f, 0.3f, 1.0f, 1.f);
   v[2].color = glm::vec4(0.5f, 0.9f, 0.2f, 1.f);
   v[3].color = glm::vec4(0.2f, 0.2f, 0.9f, 1.f);
   v[4].color = glm::vec4(0.1f, 0.9f, 0.2f, 1.f);
   v[5].color = glm::vec4(0.9f, 0.3f, 0.2f, 1.f);

   unsigned int idx_base = vertices.size();
   vertices.insert(vertices.end(), v.begin(), v.end());

   g_triangle gt_0(idx_base,   idx_base+1, idx_base+2);
   g_triangle gt_1(idx_base+3, idx_base+4, idx_base+5);

   triangles.push_back(gt_0);
   triangles.push_back(gt_1);

}

void
Mesh::fill_with_direction_triangles() {

   std::vector<s_generic_vertex> v;
   v.resize(3 * 3);
   float scale = 0.25;

   v[0].pos = scale * glm::vec3( -0.2f, 0.0f, 0.0f);
   v[1].pos = scale * glm::vec3(  0.2f, 0.0f, 0.0f);
   v[2].pos = scale * glm::vec3(  0.0f, 0.0f, 0.5f);

   v[3].pos = scale * glm::vec3( -0.2f,  0.0f, 0.0f);
   v[4].pos = scale * glm::vec3(  0.2f,  0.0f, 0.0f);
   v[5].pos = scale * glm::vec3(  0.0f,  0.5f, 0.0f);

   v[6].pos = scale * glm::vec3(  0.0f, -0.2f, 0.01f);
   v[7].pos = scale * glm::vec3(  0.0f,  0.2f, 0.01f);
   v[8].pos = scale * glm::vec3(  0.5f,  0.0f, 0.0f);

   v[0].normal = glm::vec3( 0.2f, 0.2f,  0.9f);
   v[1].normal = glm::vec3( 0.2f, 0.9f,  0.2f);
   v[2].normal = glm::vec3( 0.9f, 0.1f,  0.1f);
   v[3].normal = glm::vec3( 0.0f, 0.9f, -0.1f);
   v[4].normal = glm::vec3( 0.9f, 0.3f, -0.2f);
   v[5].normal = glm::vec3( 0.1f, 0.9f, -0.1f);
   v[6].normal = glm::vec3( 0.0f, 0.9f, -0.1f);
   v[7].normal = glm::vec3( 0.9f, 0.3f, -0.2f);
   v[8].normal = glm::vec3( 0.1f, 0.1f, -0.9f);

   v[0].color = glm::vec4(0.8f, 0.0f, 0.0f, 1.f);
   v[1].color = glm::vec4(0.8f, 0.3f, 1.0f, 1.f);
   v[2].color = glm::vec4(0.8f, 0.1f, 0.1f, 1.f);
   v[3].color = glm::vec4(0.2f, 0.8f, 0.9f, 1.f);
   v[4].color = glm::vec4(0.1f, 0.9f, 0.2f, 1.f);
   v[5].color = glm::vec4(0.1f, 0.8f, 0.1f, 1.f);
   v[6].color = glm::vec4(0.4f, 0.2f, 0.3f, 1.f);
   v[7].color = glm::vec4(0.1f, 0.4f, 0.3f, 1.f);
   v[8].color = glm::vec4(0.1f, 0.1f, 0.9f, 1.f);

   unsigned int base = vertices.size();
   vertices.insert(vertices.end(), v.begin(), v.end());
 
   g_triangle gt_0(base,base+1,base+2);
   g_triangle gt_1(base+3,base+4,base+5);
   g_triangle gt_2(base+6,base+7,base+8);
   triangles.push_back(gt_0);
   triangles.push_back(gt_1);
   triangles.push_back(gt_2);

}


void
Mesh::draw_normals(const glm::mat4 &mvp, float normal_scaling) {

   GLenum err = glGetError(); if (err) std::cout << "   error draw_normals() -- start -- "
                                                 << err << std::endl;

   if (! normals_are_setup) {
      glGenVertexArrays(1, &normals_vao);
      std::cout << "####### draw_normals() new normals_vao " << normals_vao << std::endl;
   }
   glBindVertexArray(normals_vao);

   auto vec_length = [] (const glm::vec3 &v) {
                        float s = v.x * v.x + v.y * v.y + v.z * v.z;
                        return sqrtf(s);
                     };

   if (shader_for_draw_normals.name.empty())
      shader_for_draw_normals.init("draw-normals.shader", Shader::Entity_t::GENERIC_DISPLAY_OBJECT);

   shader_for_draw_normals.Use();
   glUniformMatrix4fv(shader_for_draw_normals.mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);

   err = glGetError(); if (err) std::cout << "   error draw_normals() post mvp uniform "
                                          << err << std::endl;

   std::vector<glm::vec3> tmp_normals; // these are positions, 2 vertices per normal
   std::vector<glm::vec4> tmp_colours; // 2 colours per normal, one for each of the tmp_normals

   if (! normals_are_setup) {
      for (unsigned int i=0; i<vertices.size(); i++) {
         float f = vec_length(vertices[i].normal);
         if (false)
            std::cout << i
                      << " vertex: " << glm::to_string(vertices[i].pos)
                      << " normal: " << glm::to_string(vertices[i].normal)
                      << " length " << f << std::endl;
         tmp_normals.push_back(vertices[i].pos);
         tmp_normals.push_back(vertices[i].pos + normal_scaling * vertices[i].normal);
         tmp_colours.push_back(glm::vec4(0.4, 0.4, 0.7, 1.0));
         tmp_colours.push_back(glm::vec4(0.4, 0.4, 0.7, 1.0));
      }

      for (unsigned int i=0; i<triangles.size(); i++) {
         const glm::vec3 &p0 = vertices[triangles[i].point_id[0]].pos;
         const glm::vec3 &p1 = vertices[triangles[i].point_id[1]].pos;
         const glm::vec3 &p2 = vertices[triangles[i].point_id[2]].pos;
         glm::vec3 delta_0(0.0001, 0.0001, 0.0001);
         glm::vec3 delta_1(0.0001, 0.0001, 0.0001);
         glm::vec3 delta_2(0.0001, 0.0001, 0.0001);
         if (p0.x < 0) delta_0.x *= -1.0;
         if (p0.y < 0) delta_0.y *= -1.0;
         if (p0.z < 0) delta_0.z *= -1.0;
         if (p1.x < 0) delta_1.x *= -1.0;
         if (p1.y < 0) delta_1.y *= -1.0;
         if (p1.z < 0) delta_1.z *= -1.0;
         if (p2.x < 0) delta_2.x *= -1.0;
         if (p2.y < 0) delta_2.y *= -1.0;
         if (p2.z < 0) delta_2.z *= -1.0;
         tmp_normals.push_back(p0+delta_0);
         tmp_normals.push_back(p1+delta_1);
         tmp_normals.push_back(p0+delta_0);
         tmp_normals.push_back(p2+delta_2);
         tmp_normals.push_back(p1+delta_1);
         tmp_normals.push_back(p2+delta_2);
         tmp_colours.push_back(glm::vec4(0.5, 0.5, 0.5, 1.0));
         tmp_colours.push_back(glm::vec4(0.5, 0.5, 0.5, 1.0));
         tmp_colours.push_back(glm::vec4(0.5, 0.5, 0.5, 1.0));
         tmp_colours.push_back(glm::vec4(0.5, 0.5, 0.5, 1.0));
         tmp_colours.push_back(glm::vec4(0.5, 0.5, 0.5, 1.0));
         tmp_colours.push_back(glm::vec4(0.5, 0.5, 0.5, 1.0));
      }
   }

   unsigned int n_normals = triangles.size() * 6 + 2 * vertices.size();

   if (! normals_are_setup)
      glGenBuffers(1, &normals_buffer_id);

   glBindBuffer(GL_ARRAY_BUFFER, normals_buffer_id);
   if (! normals_are_setup)
      glBufferData(GL_ARRAY_BUFFER, n_normals * sizeof(glm::vec3), &(tmp_normals[0]), GL_STATIC_DRAW);
   glEnableVertexAttribArray(0);
   glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), 0);

   if (! normals_are_setup)
      glGenBuffers(1, &normals_colour_buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, normals_colour_buffer_id);
   if (! normals_are_setup)
      glBufferData(GL_ARRAY_BUFFER, n_normals * sizeof(glm::vec4), &(tmp_colours[0]), GL_STATIC_DRAW);
   glEnableVertexAttribArray(1);
   glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(glm::vec4), 0);
   unsigned int n_vertices = 2 * n_normals;
   glDrawArrays(GL_LINES, 0, n_vertices); // number of vertices, not number of lines
   err = glGetError(); if (err) std::cout << "   error draw_normals() post gldrawarrays "
                                          << err << std::endl;

   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);

   if (! normals_are_setup)
      normals_are_setup = true;

}


void
Mesh::draw_instanced(Shader *shader_p,
                     const glm::mat4 &mvp,
                     const glm::mat4 &view_rotation_matrix,
                     const std::map<unsigned int, lights_info_t> &lights,
                     const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                     const glm::vec4 &background_colour,
                     bool do_depth_fog) {

   // std::cout << "Mesh::draw_instanced() Mesh " << name << " -- start -- with shader " << shader_p->name << std::endl;

   if (! draw_this_mesh) return;

   unsigned int n_triangles = triangles.size();
   GLuint n_verts = 3 * n_triangles;

   if (n_triangles == 0) return;

   GLenum err = glGetError();
   if (err) std::cout << "error Mesh::draw_instanced() " << name << " " << shader_p->name
                      << " -- start -- " << err << std::endl;
   shader_p->Use();
   const std::string &shader_name = shader_p->name;
   
   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);
   err = glGetError();
   if (err) std::cout << "error:: " << shader_p->name << " draw_instanced() post mvp uniform "
                      << err << std::endl;

   glUniformMatrix4fv(shader_p->view_rotation_uniform_location, 1, GL_FALSE, &view_rotation_matrix[0][0]);
   err = glGetError();
   if (err) std::cout << "error:: Mesh::draw_instanced() " << name << " " << shader_p->name
                      << " draw_instanced() post view rotation uniform " << err << std::endl;

   std::map<unsigned int, lights_info_t>::const_iterator it;
   unsigned int light_idx = 0;
   it = lights.find(light_idx);
   if (it != lights.end())
      shader_p->setup_light(light_idx, it->second, view_rotation_matrix);
   light_idx = 1;
   it = lights.find(light_idx);
   if (it != lights.end())
      shader_p->setup_light(light_idx, it->second, view_rotation_matrix);

   shader_p->set_vec4_for_uniform("background_colour", background_colour);
   shader_p->set_bool_for_uniform("do_depth_fog", do_depth_fog);

   err = glGetError();
   if (err) std::cout << "error draw_instanced() pre-setting material " << err << std::endl;
   shader_p->set_vec4_for_uniform( "material.ambient",   material.ambient);
   shader_p->set_vec4_for_uniform( "material.diffuse",   material.diffuse);
   shader_p->set_vec4_for_uniform( "material.specular",  material.specular);
   shader_p->set_float_for_uniform("material.shininess", material.shininess);
   shader_p->set_float_for_uniform("material.specular_strength", material.specular_strength);
   err = glGetError();
   if (err) std::cout << "error draw_instanced(): " << shader_name << " pre-set eye position " << " with GL err " << err << std::endl;
   shader_p->set_vec3_for_uniform("eye_position", eye_position);
   err = glGetError();
   if (err) std::cout << "error:: Mesh::draw_instanced() " << name << " " << shader_name << " post-set eye position "
                      << " with GL err " << err << std::endl;
   err = glGetError();
   if (err) std::cout << "error:: Mesh::draw_instanced() " << shader_name << " pre-glBindVertexArray() vao " << vao
                      << " with GL err " << err << std::endl;

   if (vao == VAO_NOT_SET)
      std::cout << "ERROR:: You forgot to setup this Mesh " << name << " " << shader_p->name << std::endl;

   // std::cout << "Mesh::draw_instanced() using vao " << vao << std::endl;
   glBindVertexArray(vao);
   err = glGetError();
   if (err) std::cout << "error:: Mesh::draw_instanced() " << shader_name << " " << name
                      << " glBindVertexArray() vao " << vao << " with GL err " << err << std::endl;

   glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
   err = glGetError(); if (err) std::cout << "error draw_instanced() glBindBuffer() v " << err << std::endl;
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
   err = glGetError(); if (err) std::cout << "error draw_instanced() glBindBuffer() i " << err << std::endl;

   glEnableVertexAttribArray(0);
   glEnableVertexAttribArray(1);
   glEnableVertexAttribArray(2);
   glEnableVertexAttribArray(3);
   glEnableVertexAttribArray(4);
   glEnableVertexAttribArray(5);
   glEnableVertexAttribArray(6);
   glEnableVertexAttribArray(7);

   glBindBuffer(GL_ARRAY_BUFFER, inst_rts_buffer_id); // needed?
   err = glGetError(); if (err) std::cout << "error draw_instanced() glBindBuffer() inst_rts_buffer_id" << std::endl;


   if (false)
      std::cout << "Mesh::draw_instanced() Mesh " << name << " drawing n_verts " << n_verts << " n_instances " << n_instances
                << " with shader " << shader_p->name << std::endl;

   glDrawElementsInstanced(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr, n_instances);
   err = glGetError();
   if (err) std::cout << "error draw_instanced() glDrawElementsInstanced()"
                      << " shader: " << shader_p->name << " vao: " << vao
                      << " n_triangle_verts: " << n_verts << " n_instances: " << n_instances
                      << " with GL err " << err << std::endl;
   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   glDisableVertexAttribArray(3);
   glDisableVertexAttribArray(4);
   glDisableVertexAttribArray(5);
   glDisableVertexAttribArray(6);
   glDisableVertexAttribArray(7);
   glUseProgram(0);

}


void
Mesh::draw_particles(Shader *shader_p, const glm::mat4 &mvp, const glm::mat4 &view_rotation) {

   if (false)
      std::cout << "in draw_particles() with n_instances " << n_instances << " and n_triangles: "
                << triangles.size() << std::endl;

   // this can happen when all the particles have life 0 - and have been removed.
   if (n_instances == 0) return;
   if (triangles.empty()) return;

   particle_draw_count += 1;
   shader_p->Use();
   glBindVertexArray(vao);
   GLenum err = glGetError();
   if (err) std::cout << "error draw_particles() " << shader_p->name
                      << " glBindVertexArray() vao " << vao
                      << " with GL err " << err << std::endl;

   glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
   err = glGetError(); if (err) std::cout << "   error draw() glBindBuffer() v "
                                          << err << std::endl;
   glEnableVertexAttribArray(0); // vertex positions
   glEnableVertexAttribArray(1); // vertex normal
   glEnableVertexAttribArray(2); // vertex colours

   glBindBuffer(GL_ARRAY_BUFFER, inst_colour_buffer_id);
   glEnableVertexAttribArray(3); // instanced model/particle colours (time varying)

   glBindBuffer(GL_ARRAY_BUFFER, inst_model_translation_buffer_id);
   glEnableVertexAttribArray(4); // instanced model/particle translations (time varying)

   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);
   err = glGetError();
   if (err) std::cout << "error:: draw_particles() " << shader_p->name
                      << " draw() post mvp uniform " << err << std::endl;

   glUniformMatrix4fv(shader_p->view_rotation_uniform_location, 1, GL_FALSE, &view_rotation[0][0]);
   err = glGetError();
   if (err) std::cout << "error:: draw_particles() " << shader_p->name
                      << " draw() post view_rotation uniform " << err << std::endl;

   //
   float rotation_angle = 0.05f * static_cast<float>(particle_draw_count);

   shader_p->set_float_for_uniform("rotation_angle", rotation_angle);

   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   unsigned int n_verts = 3 * triangles.size();
   if (false)
      std::cout << "draw_particles() " << name << " with shader " << shader_p->name
                << " drawing n_instances " << n_instances << std::endl;
   glDrawElementsInstanced(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr, n_instances);
   err = glGetError();
   if (err) std::cout << "   error draw_particles() " << shader_p->name
                      << " glDrawElementsInstanced() vao " << vao
                      << " with GL err " << err << std::endl;

   glDisable(GL_BLEND);
}

void
Mesh::draw(Shader *shader_p,
           const glm::mat4 &mvp,
           const glm::mat4 &mouse_based_rotation_matrix,
           const std::map<unsigned int, lights_info_t> &lights,
           const glm::vec3 &eye_position,
           const glm::vec4 &background_colour,
           bool do_depth_fog) {

   // std::cout << "start:: Mesh::draw() " << name << " " << shader_p->name << " " << draw_this_mesh
   // << std::endl;

   if (! draw_this_mesh) return;

   unsigned int n_triangles = triangles.size();
   GLuint n_verts = 3 * n_triangles;

   // At the moment, I don't think that it's an error to come here if there are no triangles (yet)
   // Just quietly do nothing and return.
   if (n_triangles == 0) return;

   GLenum err = glGetError();
   if (err) std::cout << "error Mesh::draw() " << name << " " << shader_p->name
                      << " -- start -- " << err << std::endl;
   shader_p->Use();
   const std::string &shader_name = shader_p->name;

   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);
   err = glGetError();
   if (err) std::cout << "error:: " << shader_p->name << " draw() post mvp uniform "
                      << err << std::endl;

   glUniformMatrix4fv(shader_p->view_rotation_uniform_location, 1, GL_FALSE,
                      &mouse_based_rotation_matrix[0][0]);
   err = glGetError();
   if (err) std::cout << "error:: Mesh::draw() " << name << " " << shader_p->name
                      << " draw() post view rotation uniform " << err << std::endl;

   std::map<unsigned int, lights_info_t>::const_iterator it;
   unsigned int light_idx = 0;
   it = lights.find(light_idx);
   if (it != lights.end())
      shader_p->setup_light(light_idx, it->second, mouse_based_rotation_matrix);
   light_idx = 1;
   it = lights.find(light_idx);
   if (it != lights.end())
      shader_p->setup_light(light_idx, it->second, mouse_based_rotation_matrix);

   // add material properties, use class built-ins this time.

   // not using the built-in - hmm.
   shader_p->set_vec4_for_uniform("background_colour", background_colour);

   shader_p->set_bool_for_uniform("do_depth_fog", do_depth_fog);

   err = glGetError(); if (err) std::cout << "   error draw() pre-setting material "
                                          << err << std::endl;
   shader_p->set_vec4_for_uniform( "material.ambient",   material.ambient);
   shader_p->set_vec4_for_uniform( "material.diffuse",   material.diffuse);
   shader_p->set_vec4_for_uniform( "material.specular",  material.specular);
   shader_p->set_float_for_uniform("material.shininess", material.shininess);
   shader_p->set_float_for_uniform("material.specular_strength", material.specular_strength);

   if (false) {
      std::cout << "debug:: draw(): " << shader_p->name << " material.ambient "
                << glm::to_string(material.ambient) << std::endl;
      std::cout << "debug:: draw(): " << shader_p->name << " material.diffuse "
                << glm::to_string(material.diffuse) << std::endl;
      std::cout << "debug:: draw(): " << shader_p->name << " material.specular "
                << glm::to_string(material.specular) << std::endl;
      std::cout << "debug:: draw(): " << shader_p->name << " material.shininess "
                << material.shininess << std::endl;
      std::cout << name << " " << shader_p->name << " sent material.shininess "
                << material.shininess << std::endl;
      std::cout << name << " " << shader_p->name << " sent material.specular_strength "
                << material.specular_strength << std::endl;

   }

   err = glGetError();
   if (err) std::cout << "error draw() " << shader_name << " pre-set eye position "
                      << " with GL err " << err << std::endl;

   if (false)
      std::cout << "sending eye_position " << glm::to_string(eye_position) << std::endl;
   shader_p->set_vec3_for_uniform("eye_position", eye_position);

   err = glGetError();
   if (err) std::cout << "error:: Mesh::draw() " << name << " " << shader_name << " post-set eye position "
                      << " with GL err " << err << std::endl;

   // bind the vertices and their indices

   err = glGetError();
   if (err) std::cout << "error:: Mesh::draw() " << shader_name << " pre-glBindVertexArray() vao " << vao
                      << " with GL err " << err << std::endl;

   if (vao == VAO_NOT_SET)
      std::cout << "ERROR:: You forgot to setup this Mesh " << name << " "
                << shader_p->name << std::endl;

   // std::cout << "Mesh::draw() using vao " << vao << std::endl;
   glBindVertexArray(vao);
   err = glGetError();
   if (err) std::cout << "error:: Mesh::draw() " << shader_name << " " << name
                      << " glBindVertexArray() vao " << vao << " with GL err " << err << std::endl;

   glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
   err = glGetError(); if (err) std::cout << "   error draw() glBindBuffer() v "
                                          << err << std::endl;
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
   err = glGetError(); if (err) std::cout << "   error draw() glBindBuffer() i "
                                          << err << std::endl;

   glEnableVertexAttribArray(0);
   glEnableVertexAttribArray(1);
   glEnableVertexAttribArray(2);

   // if (is_instanced_colours) {
   // glBindBuffer(GL_ARRAY_BUFFER, inst_colour_buffer_id);
   // err = glGetError(); if (err) std::cout << "error draw() glBindBuffer() inst col "
   // << err << std::endl;
   // }

   if (is_instanced)
      glEnableVertexAttribArray(3);

   err = glGetError();
   if (err) std::cout << "error draw() glBindBuffer() Mesh::draw() post-vertex arrays "
                      << "shader " << shader_p->name << " error " << err << std::endl;

#if 0 //using a VAO means we don't need to do this (so delete it)
   if (is_instanced_with_rts_matrix) {
      glEnableVertexAttribArray(4);
      glEnableVertexAttribArray(5);
      glEnableVertexAttribArray(6);
      glBindBuffer(GL_ARRAY_BUFFER, inst_rts_buffer_id);
      err = glGetError(); if (err) std::cout << "error draw() glBindBuffer() inst rts "
                                             << err << std::endl;
   } else {
      glBindBuffer(GL_ARRAY_BUFFER, inst_model_translation_buffer_id);
      err = glGetError(); if (err) std::cout << "error draw() glBindBuffer() inst model translation buffer "
                                             << " inst_model_translation_buffer_id "
                                             << inst_model_translation_buffer_id
                                             << " shader " << shader_p->name
                                             << " error " << err << std::endl;
   }
#endif

   err = glGetError();
   if (err) std::cout << "   error draw() " << name << " pre-draw " << err << std::endl;

   if (use_blending) {
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   }

   if (is_instanced) {

      // If you are here, did you remember to use gtk_gl_area_attach_buffers(GTK_GL_AREA(di.gl_area));
      // before making a new VAO?

      if (false)
         std::cout << "debug:: Mesh::draw() instanced: " << name << " " << shader_p->name
                   << " drawing " << n_verts
                   << " triangle vertices"  << " in " << n_instances << " instances" << std::endl;
      glDrawElementsInstanced(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr, n_instances);
      err = glGetError();
      if (err) std::cout << "error draw() glDrawElementsInstanced()"
                         << " shader: " << shader_p->name
                         << " vao: " << vao
                         << " n_triangle_verts: " << n_verts
                         << " n_instances: " << n_instances
                         << " with GL err " << err << std::endl;
   } else {

      // If you are here, did you remember to use gtk_gl_area_attach_buffers(GTK_GL_AREA(di.gl_area));
      // before making a new VAO?

      if (false)
         std::cout << "debug:: Mesh::draw() " << name << " shader " << shader_p->name
                   << " vao " << vao
                   << " drawing " << n_verts << " triangle vertices"  << std::endl;

      glDrawElements(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr);
      err = glGetError();
      if (err) std::cout << "   error draw() glDrawElements()"
                         << " of Mesh \"" << name << "\""
                         << " shader: " << shader_p->name
                         << " vao " << vao
                         << " n_triangle_verts " << n_verts
                         << " with GL err " << err << std::endl;
   }

   if (use_blending) {
      glDisable(GL_BLEND);
   }
   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   if (is_instanced) glDisableVertexAttribArray(3);
   if (is_instanced) glDisableVertexAttribArray(4);
   if (is_instanced) glDisableVertexAttribArray(5);
   if (is_instanced) glDisableVertexAttribArray(6);
   glUseProgram(0);

}

// draw symmetry with lines
void
Mesh::draw_symmetry(Shader *shader_p,
                    const glm::mat4 &mvp,
                    const glm::mat4 &mouse_based_rotation_matrix,
                    const std::map<unsigned int, lights_info_t> &lights,
                    const glm::vec3 &eye_position,
                    const glm::vec4 &background_colour,
                    bool do_depth_fog) {

   if (vao == VAO_NOT_SET)
      std::cout << "ERROR:: You forgot to setup this Mesh " << name << " "
                << shader_p->name << std::endl;

   shader_p->Use();
   GLenum err = glGetError();
   if (err) std::cout << "error:: Mesh::draw_symmetry() " << shader_p->name << " " << name
                      << " use shader with GL err " << err << std::endl;

   glBindVertexArray(vao);
   err = glGetError();
   if (err) std::cout << "error:: Mesh::draw_symmetry() " << shader_p->name << " " << name
                      << " glBindVertexArray() vao " << vao << " with GL err " << err << std::endl;

   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);
   err = glGetError();
   if (err) std::cout << "error:: " << shader_p->name << " Mesh::draw_symmetry() post mvp uniform "
                      << err << std::endl;

   // is this needed?
   // glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
   err = glGetError(); if (err) std::cout << "   error draw() glBindBuffer() " << err << std::endl;
   glEnableVertexAttribArray(0);
   glEnableVertexAttribArray(1);

   shader_p->set_vec4_for_uniform("background_colour", background_colour);
   shader_p->set_bool_for_uniform("do_depth_fog", do_depth_fog);
#ifdef __APPLE__
#else
   glLineWidth(5.0);
#endif
   unsigned int n_verts = n_symmetry_atom_lines_vertices;

   unsigned int first = 0;
   if (n_verts > 1211111111) {
      n_verts = 12;
      first = 6;
   }

   glDrawArrays(GL_LINES, first, n_verts); // first and count
   err = glGetError();
   if (err) std::cout << "error:: Mesh::draw_symmetry() " << shader_p->name << " " << name
                      << " post glDrawArrays() " << vao << " with GL err " << err << std::endl;

   if (false) {
      // why does this do damage? Anyway, the "switch" is done by binding the VAO
      glDisableVertexAttribArray(0);
      glDisableVertexAttribArray(1);
   }
   glBindVertexArray(0);
}

void
Mesh::update_vertices() {

   unsigned int vs = vertices.size();
   if (vs > 0) {
      glBindVertexArray(vao); // needed?
      glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
      glBufferSubData(GL_ARRAY_BUFFER, 0, vs * sizeof(s_generic_vertex), &(vertices[0]));
   }
}

void
Mesh::update_instancing_buffer_data(const std::vector<glm::mat4> &mats,
                                    const std::vector<glm::vec4> &colours) {

   // glBufferSubData(	GLenum        target,
   //                   GLintptr      offset,
   //                   GLsizeiptr    size,
   //                   const GLvoid *data);

   unsigned int n_mats =    mats.size();
   unsigned int n_cols = colours.size();

   // No binding of the VAO?  !!!????
   // glBindVertexArray(vao); // needed?

   if (n_mats > 0) {
      glBindBuffer(GL_ARRAY_BUFFER, inst_rts_buffer_id);
      glBufferSubData(GL_ARRAY_BUFFER, 0, n_mats * 4 * sizeof(glm::vec4), &(mats[0]));
   }
   if (n_cols > 0) {
      glBindBuffer(GL_ARRAY_BUFFER, inst_colour_buffer_id);
      glBufferSubData(GL_ARRAY_BUFFER, 0, n_cols * sizeof(glm::vec4), &(colours[0]));
   }
}

void
Mesh::update_instancing_buffer_data(const std::vector<glm::mat4> &mats) {

   // glBufferSubData(	GLenum        target,
   //                   GLintptr      offset,
   //                   GLsizeiptr    size,
   //                   const GLvoid *data);

   // No binding of the VAO?  !!!????

   int n_mats =    mats.size();
   if (n_mats > n_instances_allocated) {
      std::vector<glm::vec4> dummy;
      setup_matrix_and_colour_instancing_buffers(mats, dummy);
   }

   if (n_mats > 0) {
      glBindBuffer(GL_ARRAY_BUFFER, inst_rts_buffer_id);
      glBufferSubData(GL_ARRAY_BUFFER, 0, n_mats * 4 * sizeof(glm::vec4), &(mats[0]));
   }
}

void
Mesh::update_instancing_buffer_data_for_particles(const particle_container_t &particles) {

   GLenum err = glGetError();
   if (err) std::cout << "GL error Mesh::update_instancing_buffer_data_for_particles() A0 "
                      << "binding vao " << vao << " error " << err << std::endl;

   glBindVertexArray(vao);
   err = glGetError();
   if (err) std::cout << "GL error Mesh::update_instancing_buffer_data_for_particles() A1 "
                      << "binding vao " << vao << std::endl;

   n_instances = particles.size();

   if (false)
      std::cout << "debug:: update_instancing_buffer_data() transfering " << n_instances
                << " particle/instances " << std::endl;

   if (vao == VAO_NOT_SET)
      std::cout << "You forgot to setup this Mesh " << name << std::endl;

   glBindVertexArray(vao);

   err = glGetError();
   if (err) std::cout << "GL error Mesh::update_instancing_buffer_data_for_particles() A2 "
                      << "binding vao " << vao << " error " << err << std::endl;

   // std::cout << " particle 0 position " << glm::to_string(particles.particles[0].position)
   //           << std::endl;
   glBindBuffer(GL_ARRAY_BUFFER, inst_model_translation_buffer_id);
   err = glGetError();
   if (err) std::cout << "GL error Mesh::update_instancing_buffer_data_for_particles() A3 "
                      << " vao " << vao
                      << " inst_model_translation_buffer_id " << inst_model_translation_buffer_id
                      << "\n";
   glBufferSubData(GL_ARRAY_BUFFER, 0, n_instances * sizeof(Particle), &(particles.particles[0]));
   err = glGetError();
   if (err) std::cout << "GL error Mesh::update_instancing_buffer_data_for_particles() B\n";
   glBindBuffer(GL_ARRAY_BUFFER, inst_colour_buffer_id);
   err = glGetError();
   if (err) std::cout << "GL error Mesh::update_instancing_buffer_data_for_particles() C\n";
   glBufferSubData(GL_ARRAY_BUFFER, 0, n_instances * sizeof(Particle), &(particles.particles[0]));
   err = glGetError();
   if (err) std::cout << "GL error Mesh::update_instancing_buffer_data_for_particles() D\n";

}

void
Mesh::smooth_triangles() {

   // What other vertices have the same position as this one?
   // std::map<unsigned int, std::set<unsigned int> >

   // I need to test if aa.pos == bb.pos

   // If there are lots of triangles, that will be slow
   // so instead test aa_hash == bb_hash
   // aa_hash = make_has(aa.pos)
   // which needs a hash map
   // and a glm::vec3 -> hash for position

   // but for now, do it the simple-minded way

   int n = vertices.size();
   std::map<unsigned int, std::set<unsigned int> > m;
   for (unsigned int i=0; i<triangles.size(); i++) {
      for (unsigned int j=i; j<triangles.size(); j++) {
         if (i != j) {
            // raw compare glm::vec3
            if (vertices[triangles[i].point_id[0]].pos == vertices[triangles[j].point_id[0]].pos) {
               m[triangles[i].point_id[0]].insert(triangles[j].point_id[0]);
               m[triangles[j].point_id[0]].insert(triangles[i].point_id[0]);
            }
            if (vertices[triangles[i].point_id[0]].pos == vertices[triangles[j].point_id[1]].pos) {
               m[triangles[i].point_id[0]].insert(triangles[j].point_id[1]);
               m[triangles[j].point_id[1]].insert(triangles[i].point_id[0]);
            }
            if (vertices[triangles[i].point_id[0]].pos == vertices[triangles[j].point_id[2]].pos) {
               m[triangles[i].point_id[0]].insert(triangles[j].point_id[2]);
               m[triangles[j].point_id[2]].insert(triangles[i].point_id[0]);
            }
            if (vertices[triangles[i].point_id[1]].pos == vertices[triangles[j].point_id[0]].pos) {
               m[triangles[i].point_id[1]].insert(triangles[j].point_id[0]);
               m[triangles[j].point_id[0]].insert(triangles[i].point_id[1]);
            }
            if (vertices[triangles[i].point_id[1]].pos == vertices[triangles[j].point_id[1]].pos) {
               m[triangles[i].point_id[1]].insert(triangles[j].point_id[1]);
               m[triangles[j].point_id[1]].insert(triangles[i].point_id[1]);
            }
            if (vertices[triangles[i].point_id[1]].pos == vertices[triangles[j].point_id[2]].pos) {
               m[triangles[i].point_id[1]].insert(triangles[j].point_id[2]);
               m[triangles[j].point_id[2]].insert(triangles[i].point_id[1]);
            }
            if (vertices[triangles[i].point_id[2]].pos == vertices[triangles[j].point_id[0]].pos) {
               m[triangles[i].point_id[2]].insert(triangles[j].point_id[0]);
               m[triangles[j].point_id[0]].insert(triangles[i].point_id[2]);
            }
            if (vertices[triangles[i].point_id[2]].pos == vertices[triangles[j].point_id[1]].pos) {
               m[triangles[i].point_id[2]].insert(triangles[j].point_id[1]);
               m[triangles[j].point_id[1]].insert(triangles[i].point_id[2]);
            }
            if (vertices[triangles[i].point_id[2]].pos == vertices[triangles[j].point_id[2]].pos) {
               m[triangles[i].point_id[2]].insert(triangles[j].point_id[2]);
               m[triangles[j].point_id[2]].insert(triangles[i].point_id[2]);
            }

         }
      }
   }

   float cos_min = cosf(45.0 * 2.0 * M_PI);
   std::map<unsigned int, std::set<unsigned int> >::const_iterator it;
   for (it=m.begin(); it!=m.end(); ++it) {
      glm::vec3 sum = vertices[it->first].pos;
      glm::vec3 base_normal = vertices[it->first].normal;
      unsigned int count = 1;
      std::set<unsigned int>::const_iterator it_s;
      for (it_s=it->second.begin(); it_s!=it->second.end(); ++it_s) {
         float dp = glm::dot(base_normal, vertices[*it_s].normal);
         if (dp > cos_min) {
            sum += vertices[*it_s].pos;
            count += 1;
         }
      }
      float sc = 1.0f / static_cast<float>(count);
      glm::vec3 av = sc * sum;
      vertices[it->first].pos = av;
   }

}

#include <glm/gtx/rotate_vector.hpp> // for orientation

void
Mesh::setup_camera_facing_outline() {

   if (true) {

      unsigned int idx_base = vertices.size();
      add_one_origin_cylinder();
      float half = 0.5;
      float angle = 0.5f * M_PI;
      for (unsigned int i=idx_base; i<vertices.size(); i++) {
         vertices[i].pos.x *= 0.03f;
         vertices[i].pos.y *= 0.03f;
         vertices[i].pos    = glm::rotate(vertices[i].pos,    angle, glm::vec3(1,0,0));
         vertices[i].normal = glm::rotate(vertices[i].normal, angle, glm::vec3(1,0,0));
         vertices[i].pos.x += half;
         vertices[i].pos.y += half;
         vertices[i].color = glm::vec4(0.3f, 0.4f, 0.5f, 1.0f);
      }

      idx_base = vertices.size();
      add_one_origin_cylinder();
      for (unsigned int i=idx_base; i<vertices.size(); i++) {
         vertices[i].pos.x *= 0.03f;
         vertices[i].pos.y *= 0.03f;
         vertices[i].pos    = glm::rotate(vertices[i].pos,    angle, glm::vec3(1,0,0));
         vertices[i].normal = glm::rotate(vertices[i].normal, angle, glm::vec3(1,0,0));
         vertices[i].pos.x -= half;
         vertices[i].pos.y += half;
         vertices[i].color = glm::vec4(0.3f, 0.4f, 0.5f, 1.0f);
      }

      idx_base = vertices.size();
      add_one_origin_cylinder();
      for (unsigned int i=idx_base; i<vertices.size(); i++) {
         vertices[i].pos.x *= 0.03f;
         vertices[i].pos.y *= 0.03f;
         vertices[i].pos    = glm::rotate(vertices[i].pos,    angle, glm::vec3(0,1,0));
         vertices[i].normal = glm::rotate(vertices[i].normal, angle, glm::vec3(0,1,0));
         vertices[i].pos.x += -half;
         vertices[i].pos.y += -half;
         vertices[i].color = glm::vec4(0.3f, 0.4f, 0.5f, 1.0f);
      }

      idx_base = vertices.size();
      add_one_origin_cylinder();
      for (unsigned int i=idx_base; i<vertices.size(); i++) {
         vertices[i].pos.x *= 0.03f;
         vertices[i].pos.y *= 0.03f;
         vertices[i].pos    = glm::rotate(vertices[i].pos,    angle, glm::vec3(0,1,0));
         vertices[i].normal = glm::rotate(vertices[i].normal, angle, glm::vec3(0,1,0));
         vertices[i].pos.x += -half;
         vertices[i].pos.y +=  half;
         vertices[i].color = glm::vec4(0.3f, 0.4f, 0.5f, 1.0f);
      }
   }
   setup_buffers();
}

void
Mesh::setup_camera_facing_quad() {

   float scale_x = 0.4; // pass?
   float scale_y = 0.2;

   glm::vec3 n(0,0,1);
   glm::vec4 col(1.0, 1.0, 1.0, 1.0);

   vertices.clear();
   triangles.clear();

   vertices.push_back(s_generic_vertex(glm::vec3( scale_x,  scale_y, 0.0f), n, col));
   vertices.push_back(s_generic_vertex(glm::vec3(-scale_x,  scale_y, 0.0f), n, col));
   vertices.push_back(s_generic_vertex(glm::vec3(-scale_x, -scale_y, 0.0f), n, col));
   vertices.push_back(s_generic_vertex(glm::vec3( scale_x, -scale_y, 0.0f), n, col));

   triangles.push_back(g_triangle(0,1,2));
   triangles.push_back(g_triangle(2,3,0));

   setup_buffers();

}


void
Mesh::setup_camera_facing_hex() {

   // for hexagonal coloured particles

   float s = 0.1;
   glm::vec3 n(0,0,1);
   glm::vec4 cc(0.4, 0.4, 0.4, 1.0);
   glm::vec4  c(0.4, 0.4, 0.4, 0.1);
   float ot = 0.5;
   float tt = 0.7;
   s_generic_vertex g0(s * glm::vec3(  0,   0, 0), n, cc);
   s_generic_vertex g1(s * glm::vec3( tt,  ot, 0), n, c);
   s_generic_vertex g2(s * glm::vec3( tt, -ot, 0), n, c);
   s_generic_vertex g3(s * glm::vec3(  0,  -1, 0), n, c);
   s_generic_vertex g4(s * glm::vec3(-tt, -ot, 0), n, c);
   s_generic_vertex g5(s * glm::vec3(-tt,  ot, 0), n, c);
   s_generic_vertex g6(s * glm::vec3(  0,   1, 0), n, c);
   unsigned int idx_base = vertices.size();
   unsigned int idx_tri_base = triangles.size();
   vertices.push_back(g0);
   vertices.push_back(g1);
   vertices.push_back(g2);
   vertices.push_back(g3);
   vertices.push_back(g4);
   vertices.push_back(g5);
   vertices.push_back(g6);
   triangles.push_back(g_triangle(0,1,2));
   triangles.push_back(g_triangle(0,2,3));
   triangles.push_back(g_triangle(0,3,4));
   triangles.push_back(g_triangle(0,4,5));
   triangles.push_back(g_triangle(0,5,6));
   triangles.push_back(g_triangle(0,6,1));
   if (idx_tri_base != 0)
      for (unsigned int i=idx_tri_base; i<triangles.size(); i++)
         triangles[i].rebase(idx_base);

   // setup_buffers();
 
}

void
Mesh::setup_camera_facing_polygon(unsigned int n_sides, float scale) {

   bool stellation = true; // pass this 

   float turn_per_step = 2.0f * M_PI / static_cast<float>(n_sides);
   glm::vec3 n(0,0,1);
   glm::vec4 ccol(1.0, 1.0, 1.0, 1.00);
   glm::vec4  col(0.4, 0.4, 0.4, 0.951);

   unsigned int idx_base = vertices.size();
   unsigned int idx_tri_base = triangles.size();
   vertices.push_back(s_generic_vertex(glm::vec3(0.0f, 0.0f, 0.0f), n, ccol));

   if (stellation) {
      for (unsigned int i=0; i<n_sides; i++) {
         float a1 = static_cast<float>(i) * turn_per_step;
         float a2 = (static_cast<float>(i) + 0.5f) * turn_per_step;
         float s1 = sinf(a1);
         float c1 = cosf(a1);
         float s2 = sinf(a2);
         float c2 = cosf(a2);
         glm::vec3 v1 = 1.0f * scale * glm::vec3(s1, c1, 0.0f);
         glm::vec3 v2 = 0.3f * scale * glm::vec3(s2, c2, 0.0f);
         vertices.push_back(s_generic_vertex(v1, n, col));
         vertices.push_back(s_generic_vertex(v2, n, col));
      }
      for (unsigned int idx=0; idx<n_sides; idx++) {
         unsigned int idx_this                = 2 * idx + 1;
         unsigned int idx_for_in_vertex       = 2 * idx + 2;
         unsigned int idx_for_next_out_vertex = 2 * idx + 3;
         if (idx_for_next_out_vertex == 2 * n_sides + 1) idx_for_next_out_vertex = 1;
         triangles.push_back(g_triangle(0, idx_this, idx_for_in_vertex));
         triangles.push_back(g_triangle(0, idx_for_in_vertex, idx_for_next_out_vertex));
      }
      
   } else {
      for (unsigned int i=0; i<n_sides; i++) {
         float a = static_cast<float>(i) * turn_per_step;
         float s = sinf(a);
         float c = cosf(a);
         vertices.push_back(s_generic_vertex(scale * glm::vec3(s, c, 0.0f), n, col));
      }
      for (unsigned int idx=1; idx<=n_sides; idx++) {
         unsigned int idx_next = idx + 1;
         if (idx == n_sides) idx_next = 1;
         triangles.push_back(g_triangle(0, idx, idx_next));
      }
      if (idx_tri_base != 0)
         for (unsigned int i=idx_tri_base; i<triangles.size(); i++)
            triangles[i].rebase(idx_base);
   }

   setup_buffers(); // this was not here - why not?

}


#include "cylinder-with-rotation-translation.hh"

void
Mesh::setup_hydrogen_bond_cyclinders(Shader *shader_p, const Material &material_in) {

   // call this
   // gtk_gl_area_attach_buffers(GTK_GL_AREA(graphics_info_t::glareas[0]));
   // before calling this function

   // Make some space (allocate the buffers) for the rotatation/translation instances also.
   // (the colours are not instanced)

   is_instanced_colours = false; // keep the original colours
   is_instanced = true;
   is_instanced_with_rts_matrix = true;

   material = material_in;
   shader_p->Use();

   unsigned int n_slices = 20;
   unsigned int n_stacks = 80;
   glm::vec3 start_pos(0, 0, 0);
   glm::vec3 end_pos(1, 0, 0);
   std::pair<glm::vec3, glm::vec3> pp(start_pos, end_pos);
   float height = glm::distance(start_pos, end_pos);
   float radius = 0.05;
   cylinder_with_rotation_translation c(pp, radius, radius, height, n_slices, n_stacks);
   if (true)
      c.add_spiral(); // take 2 colours

   // now convert the vertices
   std::vector<s_generic_vertex> new_vertices(c.vertices.size());
   for (unsigned int ii=0; ii<c.vertices.size(); ii++) {
      const vertex_with_rotation_translation &v = c.vertices[ii];
      s_generic_vertex gv(v.pos, v.normal, v.colour);
      // gv.color = glm::vec4(0,1,0,1);
      new_vertices[ii] = gv;
   }
   unsigned int idx_base = vertices.size();
   unsigned int idx_tri_base = triangles.size();
   vertices.insert(vertices.end(), new_vertices.begin(), new_vertices.end());
   triangles.insert(triangles.end(), c.triangle_indices_vec.begin(), c.triangle_indices_vec.end());
   for (unsigned int ii=idx_tri_base; ii<triangles.size(); ii++)
      triangles[ii].rebase(idx_base);

   setup_buffers();

   std::vector<glm::mat4> mats(1000, glm::mat4(1.0f));
   std::vector<glm::vec4> colours; //dummy
   n_instances = mats.size();
   setup_matrix_and_colour_instancing_buffers(mats, colours);
}


void
Mesh::test_cyclinders(Shader *shader_p, const Material &material_in) {

   // call this
   // gtk_gl_area_attach_buffers(GTK_GL_AREA(graphics_info_t::glareas[0]));
   // before calling this function

   is_instanced_colours = false; // keep the original colours
   is_instanced = true;
   is_instanced_with_rts_matrix = true;

   material = material_in;
   shader_p->Use();

   unsigned int n_slices = 20;
   unsigned int n_stacks = 80;
   glm::vec3 start_pos(0, 0, 0);
   glm::vec3 end_pos(1, 0, 0);
   std::pair<glm::vec3, glm::vec3> pp(start_pos, end_pos);
   float height = glm::distance(start_pos, end_pos);
   float radius = 0.05;
   cylinder_with_rotation_translation c(pp, radius, radius, height, n_slices, n_stacks);
   if (true)
      c.add_spiral(); // take 2 colours

   // now convert the vertices
   std::vector<s_generic_vertex> new_vertices(c.vertices.size());
   for (unsigned int ii=0; ii<c.vertices.size(); ii++) {
      const vertex_with_rotation_translation &v = c.vertices[ii];
      s_generic_vertex gv(v.pos, v.normal, v.colour);
      // gv.color = glm::vec4(0,1,0,1);
      new_vertices[ii] = gv;
   }
   unsigned int idx_base = vertices.size();
   unsigned int idx_tri_base = triangles.size();
   vertices.insert(vertices.end(), new_vertices.begin(), new_vertices.end());
   triangles.insert(triangles.end(), c.triangle_indices_vec.begin(), c.triangle_indices_vec.end());
   for (unsigned int ii=idx_tri_base; ii<triangles.size(); ii++)
      triangles[ii].rebase(idx_base);


   setup_buffers();

   // some test cylinders
   std::vector<glm::mat4> mats;
   std::vector<glm::vec4> colours; // empty

   for (auto i=0; i<10; i++) {
      float x = static_cast<float>(i);
      glm::mat4 m(1.0);
      glm::vec3 t(2 * x, 0, 0);
      m += glm::translate(t);
      mats.push_back(m);

      // see get_bond_matrix() in update_mats_and_colours() in the generator
   }

   glm::vec3 p1(42.08, 9.67, 14.42);
   glm::vec3 p2(40.59, 5.68, 13.24);

   glm::vec3 p3(44.88, 12.95, 8.76);
   glm::vec3 p4(46.13, 10.59, 9.97);

   float theta = 1.0;

   mats.push_back(make_hydrogen_bond_cylinder_orientation(p1, p2, theta));
   mats.push_back(make_hydrogen_bond_cylinder_orientation(p3, p4, theta));

   n_instances = mats.size();
   setup_matrix_and_colour_instancing_buffers(mats, colours);


}

// static
glm::mat4
Mesh::make_hydrogen_bond_cylinder_orientation(const glm::vec3 &p1, const glm::vec3 &p2, float theta) {

   glm::mat4 u(1.0f);
   glm::vec3 delta(p2-p1);
   float h = glm::distance(p1, p2);
   glm::mat4 sc = glm::scale(u, glm::vec3(1.0f, 1.0f, h));
   glm::mat4 rot = glm::rotate(u, theta, glm::vec3(0.0f, 0.0f, 1.0f));
   glm::vec3 normalized_bond_orientation(glm::normalize(delta));
   glm::mat4 ori = glm::orientation(normalized_bond_orientation,
                                    glm::vec3(0.0f, 0.0f, 1.0f));
   glm::mat4 t = glm::translate(u, p1);
   glm::mat4 m = t * ori * rot * sc;
   return m;
}

#include <fstream>

bool
Mesh::export_as_obj(const std::string &file_name) const {
   return export_as_obj_internal(file_name);
}


bool
Mesh::export_as_obj_internal(const std::string &file_name) const {

   bool status = true;

   std::cout << "debug:: export_as_obj_internal: n vertices:  " <<  vertices.size() << std::endl;
   std::cout << "debug:: export_as_obj_internal: n triangles: " << triangles.size() << std::endl;

   std::ofstream f(file_name.c_str());
   if (f) {
      f << "# " << name << "\n";
      f << "# " << "\n";
      f << "" << "\n";
      f << "g exported_obj\n";
      for (unsigned int i=0; i<vertices.size(); i++) {
         const s_generic_vertex &vert = vertices[i];
         f << "v " << vert.pos.x << " " << vert.pos.y << " " << vert.pos.z;
         f << " " << vert.color.r << " " << vert.color.g << " " << vert.color.b;
         f << "\n";
      }
      for (unsigned int i=0; i<vertices.size(); i++) {
         const s_generic_vertex &vert = vertices[i];
         f << "vn " << -vert.normal.x << " " << -vert.normal.y << " " << -vert.normal.z << "\n";
      }
      for (unsigned int i=0; i<triangles.size(); i++) {
         const g_triangle &tri = triangles[i];
         f << "f "
           << tri.point_id[0]+1 << "//" << tri.point_id[0]+1 << " "
           << tri.point_id[1]+1 << "//" << tri.point_id[1]+1 << " "
           << tri.point_id[2]+1 << "//" << tri.point_id[2]+1 << "\n";
      }
   } else {
      status = false;
   }
   return status;
}

bool
Mesh::export_as_obj(std::ofstream &f, unsigned int vertex_index_offset) const {

   bool status = true;
   if (f) {
      for (unsigned int i=0; i<vertices.size(); i++) {
         const s_generic_vertex &vert = vertices[i];
         f << "v " << vert.pos.x << " " << vert.pos.y << " " << vert.pos.z;
         // f << " " << vert.color.r << " " << vert.color.g << " " << vert.color.b;
         f << "\n";
      }
      for (unsigned int i=0; i<vertices.size(); i++) {
         const s_generic_vertex &vert = vertices[i];
         f << "vn " << vert.normal.x << " " << vert.normal.y << " " << vert.normal.z << "\n";
      }
      for (unsigned int i=0; i<triangles.size(); i++) {
         const g_triangle &tri = triangles[i];
         f << "f "
           << tri.point_id[0]+1+vertex_index_offset << "//" << tri.point_id[0]+1+vertex_index_offset << " "
           << tri.point_id[1]+1+vertex_index_offset << "//" << tri.point_id[1]+1+vertex_index_offset << " "
           << tri.point_id[2]+1+vertex_index_offset << "//" << tri.point_id[2]+1+vertex_index_offset << "\n";
      }
   }
   return status;
}



// We import from assimp from Model and export from Mesh - at the moment
//
#ifdef USE_ASSIMP
#include <assimp/Exporter.hpp>
#endif // USE_ASSIMP

bool
Mesh::export_as_obj_via_assimp(const std::string &file_name) const {

   unsigned int status = false;

   std::cout << "exporting to " << file_name << std::endl;

#ifdef USE_ASSIMP

   if (! vertices.empty()) {
      // this will do a copy. Is that what I want?
      aiScene scene = generate_scene();
      // now export scener to file_name;
      Assimp::Exporter ae;

      if (false) { // 20 formats - "obj" is one of them
         size_t n_formats = ae.GetExportFormatCount();
         for (size_t i=0; i<n_formats; i++) {
            const aiExportFormatDesc *fd = ae.GetExportFormatDescription(i);
            std::cout << i << " " << fd->id << " " << fd->description << std::endl;
         }
      }

      std::string format_id = "obj";

      std::cout << "----- calling ae.Export() " << std::endl;
      std::cout << "----- calling ae.Export() scene: " << &scene << std::endl;
      std::cout << "----- calling ae.Export() mMaterials " << scene.mMaterials << std::endl;
      std::cout << "----- calling ae.Export() mMaterials[0] " << scene.mMaterials[0] << std::endl;
      // Oh, scene gets copies in Export()!
      aiReturn air = ae.Export(&scene, format_id.c_str(), file_name.c_str(), 0);
      std::cout << "export status " << air << std::endl;

   }
#endif
   return status;
}

#if USE_ASSIMP
// Use make_shared and a shared or unique? pointer as the return value
aiScene
Mesh::generate_scene() const {

   aiScene scene;
   scene.mRootNode = new aiNode();

   std::cout << "debug:: scene.mNumMaterials " << scene.mNumMaterials << std::endl;

   // Materials? Maybe we have to.
   // else when Exporter::Export() does the scene copy the scene.mMaterials gets
   // set to null and GetMaterialName() fails

   scene.mNumMaterials = 1;
   scene.mMaterials = new aiMaterial *[1];
   scene.mMaterials[0] = new aiMaterial();

   aiMesh mesh;
   scene.mMeshes = new aiMesh *[1];
   scene.mNumMeshes = 1;
   scene.mMeshes[0] = new aiMesh();
   std::cout << "new aimesh mMaterialIndex " << scene.mMeshes[0]->mMaterialIndex << std::endl;
   scene.mMeshes[0]->mMaterialIndex = 0; // Hmm? Dangerous?
   scene.mRootNode->mMeshes = new unsigned int[1];
   scene.mRootNode->mMeshes[0] = 0;
   scene.mRootNode->mNumMeshes = 1;
   aiMesh *pMesh = scene.mMeshes[0];

   //  --- vertices ---

   pMesh->mVertices = new aiVector3D[vertices.size()];
   pMesh->mNumVertices = vertices.size();
   pMesh->mTextureCoords[0] = new aiVector3D[vertices.size()];       // needed?
   pMesh->mNumUVComponents[0] = vertices.size();                     // needed?
   for (unsigned int i=0; i<vertices.size(); i++) {
      const s_generic_vertex &vert = vertices[i];
      aiVector3D aiv(vert.pos.x, vert.pos.y, vert.pos.z);
      // std::cout << "vertex " << i << " " << glm::to_string(vert.pos)  << std::endl;
      pMesh->mVertices[i] = aiv;
      pMesh->mTextureCoords[0][i] = aiVector3D(0,0,0); // for now
   }

   //  --- normals ---


   //  --- triangles ---

   pMesh->mFaces = new aiFace[triangles.size()];
   pMesh->mNumFaces = triangles.size();
   for (unsigned int i=0; i<triangles.size(); i++) {
      aiFace &face = pMesh->mFaces[i];
      face.mIndices = new unsigned int[3];
      face.mIndices[0] = triangles[i].point_id[0];
      face.mIndices[1] = triangles[i].point_id[1];
      face.mIndices[2] = triangles[i].point_id[2];
   }

   return scene;
}
#endif

