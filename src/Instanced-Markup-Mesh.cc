

#include <iostream>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>  // to_string()

#include "Instanced-Markup-Mesh.hh"

#include "oct.hh"

void
Instanced_Markup_Mesh::init() {

   first_time = true;
   n_instances = 0;
   max_n_instances = 0;
   draw_this_mesh = true;
   this_mesh_is_closed = false;
   vao = VAO_NOT_SET;
}

void
Instanced_Markup_Mesh::clear() {

   vertices.clear();
   triangles.clear();

   // don't reset draw_this_mesh

   // the buffers will be reset (deleted, created) on the setup_octasphere() and setup_instancing_buffers()

}

void
Instanced_Markup_Mesh::close() {

   if (this_mesh_is_closed) {
      // nothing to do
   } else {
      clear();
      draw_this_mesh = false;
      this_mesh_is_closed = true;

      if (! first_time) {
         glDeleteBuffers(1, &buffer_id);
         glDeleteBuffers(1, &index_buffer_id);
      }
   }
}

void
Instanced_Markup_Mesh::setup_buffers() {

   if (triangles.empty()) return;
   if (vertices.empty()) return;

   if (first_time)
      glGenVertexArrays(1, &vao);

   if (vao == VAO_NOT_SET)
      std::cout << "ERROR:: in Instanced_Markup_Mesh::setup_buffers() vao not set" << std::endl;

   glBindVertexArray(vao);

   unsigned int n_vertices = vertices.size(); // 4

   if (first_time) {
      glGenBuffers(1, &buffer_id);
      glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
      glBufferData(GL_ARRAY_BUFFER, n_vertices * sizeof(Instanced_Markup_Mesh_Vertex_attrib_t),
                   &(vertices[0]), GL_DYNAMIC_DRAW);
   } else {
      glDeleteBuffers(1, &buffer_id);
      glGenBuffers(1, &buffer_id);
      glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
      glBufferData(GL_ARRAY_BUFFER, n_vertices * sizeof(Instanced_Markup_Mesh_Vertex_attrib_t),
                   &(vertices[0]), GL_DYNAMIC_DRAW);
   }

   // position (of the quad)
   glEnableVertexAttribArray(0);
   glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Instanced_Markup_Mesh_Vertex_attrib_t), 0);
   // pointyness of the quad
   glEnableVertexAttribArray(1);
   glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(Instanced_Markup_Mesh_Vertex_attrib_t),
                         reinterpret_cast<void *>(sizeof(glm::vec3)));

   unsigned int n_triangles = triangles.size();
   unsigned int n_bytes = n_triangles * 3 * sizeof(unsigned int);

   if (first_time) {
      glGenBuffers(1, &index_buffer_id);
      GLenum err = glGetError(); if (err) std::cout << "GL error HUDTextureMesh setup_buffers()\n";
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
      err = glGetError(); if (err) std::cout << "GL error HUDMesh setup_buffers()\n";
   } else {
      glDeleteBuffers(1, &index_buffer_id);
      glGenBuffers(1, &index_buffer_id);
      GLenum err = glGetError(); if (err) std::cout << "GL error HUDMesh setup_buffers()\n";
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
      err = glGetError(); if (err) std::cout << "GL error HUDMesh setup_buffers()\n";
   }

   // std::cout << "HUDMesh::setup_buffers() indices " << n_bytes << " bytes" << std::endl;
   glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_bytes, &triangles[0], GL_DYNAMIC_DRAW);
   GLenum err = glGetError(); if (err) std::cout << "GL error HUDMesh setup_simple_triangles()\n";

   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glBindBuffer(GL_ARRAY_BUFFER, 0);
   glUseProgram(0);
   glBindVertexArray(0);

   first_time = false;

   

}

void
Instanced_Markup_Mesh::setup_instancing_buffers(unsigned int max_nun_instances) {

   GLenum err = glGetError();
   if (err) std::cout << "GL ERROR:: Instanced_Markup_Mesh::setup_instancing_buffers() "
                      << name << " -- start -- " << err << std::endl;

   max_n_instances = max_nun_instances;
   n_instances = 0;

   if (vao == VAO_NOT_SET)
      std::cout << "ERROR:: in Instanced_Markup_Mesh::setup_instancing_buffers() vao not set" << std::endl;

   glBindVertexArray (vao);

   err = glGetError();
   if (err) std::cout << "GL ERROR:: Instanced_Markup_Mesh::setup_instancing_buffers() "
                      << name << " A4 " << err << std::endl;

   glEnableVertexAttribArray(2);
   glEnableVertexAttribArray(3);
   glEnableVertexAttribArray(4);
   glEnableVertexAttribArray(5);
   glEnableVertexAttribArray(6);

   // layout
   // 0 vec2 position           (vertex)
   // 1 float displacemnt       (vertex)
   // 2 vec4 colour             (instanced)
   // 3 vec3 position           (instanced)
   // 4 float size              (instanced)
   // 5 float specular_strength (instanced)
   // 6 float shininess         (instanced)

   err = glGetError();
   if (err) std::cout << "GL ERROR:: Instanced_Markup_Mesh::setup_instancing_buffers() "
                      << name << " A4 " << err << std::endl;

   glGenBuffers(1, &inst_attribs_buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, inst_attribs_buffer_id);
   glBufferData(GL_ARRAY_BUFFER, max_n_instances * sizeof(Instanced_Markup_Mesh_attrib_t),
                nullptr, GL_DYNAMIC_DRAW);
   glEnableVertexAttribArray(2);
   glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(Instanced_Markup_Mesh_attrib_t), 0);
   glVertexAttribDivisor(2, 1);

   glEnableVertexAttribArray(3);
   glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(Instanced_Markup_Mesh_attrib_t),
                         reinterpret_cast<void *>(sizeof(glm::vec4)));
   glVertexAttribDivisor(3, 1);
   
   glEnableVertexAttribArray(4);
   glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, sizeof(Instanced_Markup_Mesh_attrib_t),
                         reinterpret_cast<void *>(sizeof(glm::vec4) + sizeof(glm::vec3)));
   glVertexAttribDivisor(4, 1);
   
   glEnableVertexAttribArray(5);
   glVertexAttribPointer(5, 1, GL_FLOAT, GL_FALSE, sizeof(Instanced_Markup_Mesh_attrib_t),
                         reinterpret_cast<void *>(sizeof(glm::vec4) + sizeof(glm::vec3) + sizeof(float)));
   glVertexAttribDivisor(5, 1);
   
   glEnableVertexAttribArray(6);
   glVertexAttribPointer(6, 1, GL_FLOAT, GL_FALSE, sizeof(Instanced_Markup_Mesh_attrib_t),
                         reinterpret_cast<void *>(sizeof(glm::vec4) + sizeof(glm::vec3) + 2 * sizeof(float)));
   glVertexAttribDivisor(6, 1);

   err = glGetError();
   if (err) std::cout << "GL ERROR:: Instanced_Markup_Mesh::setup_instancing_buffers() "
                      << name << " B " << err << std::endl;

   glDisableVertexAttribArray(2);
   glDisableVertexAttribArray(3);
   glDisableVertexAttribArray(4);
   glDisableVertexAttribArray(5);
   glDisableVertexAttribArray(6);

   glBindBuffer(GL_ARRAY_BUFFER, 0);

   if (err) std::cout << "GL ERROR:: Instanced_Markup_Mesh::setup_instancing_buffers() "
                      << name << " C " << err << std::endl;
}

void
Instanced_Markup_Mesh::update_instancing_buffers(const std::vector<Instanced_Markup_Mesh_attrib_t> &balls) {

   unsigned int s = balls.size();
   n_instances = s;
   if (n_instances > max_n_instances)
      n_instances = max_n_instances;

   glBindBuffer(GL_ARRAY_BUFFER, inst_attribs_buffer_id);
   glBufferSubData(GL_ARRAY_BUFFER, 0, n_instances * sizeof(Instanced_Markup_Mesh_attrib_t), &(balls[0]));

}

#include "utils/coot-utils.hh"

void
Instanced_Markup_Mesh::setup_octasphere(unsigned int num_subdivisions) {

   // really we want a version of make_octasphere that only returns
   // vectors (glm::vec3) and triangles - radius 1.

   glm::vec3 position(0,0,0);
   glm::vec4 colour(0, 0, 0, 1);
   float radius = 1.0;

   // num_subdivisions = 3;
   //
   bool remove_redundant_vertices_flag = false; // seems not to fully work.

   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
      oct = make_octasphere(num_subdivisions, position, radius, colour, remove_redundant_vertices_flag);

   const std::vector<s_generic_vertex> &v1 = oct.first;
   const std::vector<g_triangle> &v2 = oct.second;

   vertices.resize(v1.size());

   // float irm = 1.0f / static_cast<float>(RAND_MAX);

   for (unsigned int i=0; i<v1.size(); i++) {
      vertices[i].position = v1[i].pos;
      vertices[i].displacement = 0.0; // perturbations on a sphere surface

      // sadly not - octaspheres have edges of the same position, different index.
      // (that could be fixed - by merging when we create the octasphere).
      // vertices[i].displacement = 0.6 * (1.0 - irm * static_cast<float>(coot::util::random()));
      // If we are going to do this then the normals need to be perturbed too.
      // and they will need to be part of a Instanced_Markup_Mesh_Vertex_attrib_t
      // Another time.
      // http://libnoise.sourceforge.net/tutorials/tutorial8.html
      
   }
   triangles = v2;

   std::cout << "debug:: in setup_octasphere() calling setup_buffers" << std::endl;
   setup_buffers();

}


void
Instanced_Markup_Mesh::draw(Shader *shader_p,
                            const glm::mat4 &mvp,
                            const glm::mat4 &view_rotation_matrix,
                            const std::map<unsigned int, lights_info_t> &lights,
                            const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                            const glm::vec4 &background_colour,
                            bool do_depth_fog) {

   if (! draw_this_mesh) return;

   unsigned int n_triangles = triangles.size();
   GLuint n_verts = 3 * n_triangles;

   if (n_triangles == 0) return;

   GLenum err = glGetError();
   if (err) std::cout << "error Instanced_Markup_Mesh::draw() " << name << " " << shader_p->name
                      << " -- start -- " << err << std::endl;
   shader_p->Use();
   const std::string &shader_name = shader_p->name;

   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);
   err = glGetError();
   if (err) std::cout << "error:: Instanced_Markup_Mesh::draw()" << shader_p->name
                      << " post mvp uniform " << err << std::endl;

   glUniformMatrix4fv(shader_p->view_rotation_uniform_location, 1, GL_FALSE,
                      &view_rotation_matrix[0][0]);
   err = glGetError();
   if (err) std::cout << "error:: Instanced_Markup_Mesh::draw() " << name << " " << shader_p->name
                      << " draw() post view rotation uniform " << err << std::endl;

   std::map<unsigned int, lights_info_t>::const_iterator it;
   unsigned int light_idx = 0;
   it = lights.find(light_idx);
   if (it != lights.end())
      shader_p->setup_light(light_idx, it->second, view_rotation_matrix);
   light_idx = 1;
   it = lights.find(light_idx);
   if (it != lights.end())
      shader_p->setup_light(light_idx, it->second, view_rotation_matrix);

   // add material properties, use class built-ins this time.

   // not using the build-in - hmm.
   shader_p->set_vec4_for_uniform("background_colour", background_colour);

   shader_p->set_bool_for_uniform("do_depth_fog", do_depth_fog);

   err = glGetError(); if (err) std::cout << "   error draw() pre-setting material "
                                          << err << std::endl;

   // the material is part of the attributes - not a uniform (different balls have different shininess)

   err = glGetError();
   if (err) std::cout << "error draw() " << shader_name << " pre-set eye position "
                      << " with GL err " << err << std::endl;

   if (false)
      std::cout << "sending eye_position " << glm::to_string(eye_position) << std::endl;
   shader_p->set_vec3_for_uniform("eye_position", eye_position);

   err = glGetError();
   if (err) std::cout << "error:: Mesh::draw() " << name << " " << shader_name
                      << " post-set eye position " << " with GL err " << err << std::endl;

   // bind the vertices and their indices

   err = glGetError();
   if (err) std::cout << "error:: Mesh::draw() " << shader_name
                      << " pre-glBindVertexArray() vao " << vao
                      << " with GL err " << err << std::endl;

   if (vao == 99999999)
      std::cout << "You forget to setup this mesh " << name << " "
                << shader_p->name << std::endl;

   glBindVertexArray(vao);
   err = glGetError();
   if (err) std::cout << "error:: Mesh::draw() " << shader_name << " " << name
                      << " glBindVertexArray() vao " << vao << " with GL err " << err << std::endl;

   glBindBuffer(GL_ARRAY_BUFFER, buffer_id); // needed?                       
   err = glGetError(); if (err) std::cout << "   error draw() glBindBuffer() v "
                                          << err << std::endl;
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id); // needed?                            
   err = glGetError(); if (err) std::cout << "   error draw() glBindBuffer() i "
                                          << err << std::endl;

   glEnableVertexAttribArray(0);
   glEnableVertexAttribArray(1);
   glEnableVertexAttribArray(2);
   glEnableVertexAttribArray(3);
   glEnableVertexAttribArray(4);
   glEnableVertexAttribArray(5);
   glEnableVertexAttribArray(6);

   err = glGetError();
   if (err) std::cout << "   error draw() " << name << " pre-draw " << err << std::endl;

   if (false) // are you using the correct shader? (rama-balls.shader)?
      std::cout << "debug instanced mesh draw(): " << name << " drawing " << n_verts
                << " vertices " << " n_instances " << n_instances
                << " with shader " << shader_p->name << std::endl;

   glDrawElementsInstanced(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr, n_instances);
   err = glGetError();
   if (err) std::cout << "error:: Instanced_Markup_Mesh::draw() glDrawElementsInstanced()"
                      << " shader: " << shader_p->name
                      << " vao: " << vao
                      << " n_triangle_verts: " << n_verts
                      << " n_instances: " << n_instances
                      << " with GL err " << err << std::endl;
   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   glDisableVertexAttribArray(3);
   glDisableVertexAttribArray(4);
   glDisableVertexAttribArray(5);
   glDisableVertexAttribArray(6);
   glUseProgram (0);

}
