

// THis is not compiled.

#include <epoxy/gl.h>
#include <iostream>
#include <vector>
#include <string>

#include "graphical-molecule.hh"

void
graphical_molecule::draw(Shader *shader_p,
                         const glm::mat4 &mvp,
                         const glm::mat4 &world_rotation_matrix,
                         const glm::mat4 &world_rotation_translation_matrix,
                         const std::map<unsigned int, gl_lights_info_t> &lights,
                         const glm::vec3 &eye_position) {

   // std::cout << "debug:: graphical_molecule::draw() start " << std::endl;

   GLenum err = glGetError(); if (err) std::cout << "   error draw() -- start -- "
                                                 << err << std::endl;
   shader_p->Use();
   const std::string &shader_name = shader_p->name;
   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);
   glUniformMatrix4fv(shader_p->view_rotation_uniform_location, 1, GL_FALSE, &world_rotation_matrix[0][0]);

   // std::cout << "sending view-rotation: " << glm::to_string(world_rotation_matrix) << std::endl;

   std::map<unsigned int, gl_lights_info_t>::const_iterator it;
   unsigned int light_idx = 0;
   it = lights.find(light_idx);

   // needs a shader fixup.

   // if (it != lights.end())
   // shader_p->setup_light(light_idx, it->second, world_rotation_translation_matrix);

   // add material properties, use class built-ins this time.

   err = glGetError(); if (err) std::cout << "   error draw() pre-setting material "
                                          << err << std::endl;
   shader_p->set_vec4_for_uniform( "material.ambient",   material.ambient);
   shader_p->set_vec4_for_uniform( "material.diffuse",   material.diffuse);
   shader_p->set_vec4_for_uniform( "material.specular",  material.specular);
   shader_p->set_float_for_uniform("material.shininess", material.shininess);

   if (false) {
      std::cout << "debug:: draw(): " << shader_p->name << " material.ambient " << glm::to_string(material.ambient) << std::endl;
      std::cout << "debug:: draw(): " << shader_p->name << " material.diffuse " << glm::to_string(material.diffuse) << std::endl;
      std::cout << "debug:: draw(): " << shader_p->name << " material.specular " << glm::to_string(material.specular) << std::endl;
      std::cout << "debug:: draw(): " << shader_p->name << " material.shininess " << material.shininess << std::endl;
   }

   err = glGetError();
   if (err) std::cout << "   error draw() " << shader_name << " pre-set eye position "
                      << " with GL err " << err << std::endl;
   // eye - pass this
   shader_p->set_vec3_for_uniform("eye_position", eye_position);

   err = glGetError();
   if (err) std::cout << "   error draw() " << shader_name << " post-set eye position "
                      << " with GL err " << err << std::endl;
  // bind the vertices and their indices

   err = glGetError();
   if (err) std::cout << "   error draw() " << shader_name << " pre-glBindVertexArray() vao " << vao
                      << " with GL err " << err << std::endl;
   glBindVertexArray(vao);
   err = glGetError();
   if (err) std::cout << "   error draw() " << shader_name << " glBindVertexArray() vao " << vao
                      << " with GL err " << err << std::endl;

   glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
   err = glGetError(); if (err) std::cout << "   error draw() glBindBuffer() v "
                                          << err << std::endl;
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
   err = glGetError(); if (err) std::cout << "   error draw() glBindBuffer() i "
                                          << err << std::endl;

   glEnableVertexAttribArray(0);
   glEnableVertexAttribArray(1);
   glEnableVertexAttribArray(2);
   if (is_instanced) {
      err = glGetError(); if (err) std::cout << "error draw() pre-glBindBuffer() inst colour buffer "
                                             << err << std::endl;
      glBindBuffer(GL_ARRAY_BUFFER, inst_colour_buffer_id);
      err = glGetError();
      if (err) std::cout << "error draw() glBindBuffer() inst colour " << err << std::endl;
      glEnableVertexAttribArray(3);
      glBindBuffer(GL_ARRAY_BUFFER, inst_model_translation_buffer_id);
      err = glGetError();
      if (err) std::cout << "error draw() post glBindBuffer() inst translation buffer "
                         << err << std::endl;
   }

   unsigned int n_triangles = triangle_vertex_indices.size();
   GLuint n_verts = 3 * n_triangles;

   unsigned int n_instances = 6;

   if (is_instanced) {
      glDrawElementsInstanced(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr, n_instances);
      err = glGetError();
      if (err) std::cout << "   error draw() glDrawElementsInstanced()"
                         << " shader: " << shader_p->name
                         << " vao: " << vao
                         << " n_verts: " << n_verts
                         << " with GL err " << err << std::endl;
   } else {
      // std::cout << "debug:: graphical_molecule::draw() (elements) " << n_verts << " vertices " << std::endl;
      glDrawElements(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr);
      err = glGetError();
      if (err) std::cout << "   error draw() glDrawElements()"
                         << " shader: " << shader_p->name
                         << " vao " << vao
                         << " n_verts " << n_verts
                         << " with GL err " << err << std::endl;
   }

   glDisableVertexAttribArray (0);
   glDisableVertexAttribArray (1);
   glDisableVertexAttribArray (2);
   if (is_instanced) glDisableVertexAttribArray(3);
   if (is_instanced) glDisableVertexAttribArray(4);
   glUseProgram (0);

}

void
graphical_molecule::setup_buffers() {

   // note to self: which Shader are we using when this function gets called?

   glGenVertexArrays (1, &vao);
   glBindVertexArray (vao);

   glGenBuffers(1, &buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
   std::cout << "debug:: sizeof(vertices) " << sizeof(vertices) << std::endl;
   std::cout << "debug:: sizeof(vertices[0]) " << sizeof(vertices[0]) << std::endl;
   std::cout << "debug:: sizeof(s_generic_vertex) " << sizeof(s_generic_vertex) << std::endl;
   unsigned int n_vertices = vertices.size();
   std::cout << "debug:: glBufferData for simple_triangles_buffer_id " << buffer_id
             << " allocating with size " << n_vertices * sizeof(vertices[0]) << std::endl;
   glBufferData(GL_ARRAY_BUFFER, n_vertices * sizeof(vertices[0]), &(vertices[0]), GL_STATIC_DRAW);

   // position
   glEnableVertexAttribArray(0);
   glVertexAttribPointer (0, 3, GL_FLOAT, GL_FALSE, sizeof(s_generic_vertex), 0);

   // normal
   glEnableVertexAttribArray (1);
   glVertexAttribPointer (1, 3, GL_FLOAT, GL_FALSE, sizeof(s_generic_vertex),
                          reinterpret_cast<void *>(sizeof(glm::vec3)));

   // colour
   if (! is_instanced) {
      // when the colour attribute is instanced (like the model_translation), so set up the colours
      // when we do the model_translation in setup_instancing_buffers(), not here.
      glEnableVertexAttribArray(2);
      glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(s_generic_vertex),
                            reinterpret_cast<void *>(2 * sizeof(glm::vec3)));
   }

   glGenBuffers(1, &index_buffer_id);
   GLenum err = glGetError(); if (err) std::cout << "GL error setup_simple_triangles()\n";
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
   err = glGetError(); if (err) std::cout << "GL error setup_simple_triangles()\n";
   unsigned int n_triangles = triangle_vertex_indices.size();
   unsigned int n_bytes = n_triangles * 3 * sizeof(unsigned int);
   glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_bytes, &triangle_vertex_indices[0], GL_STATIC_DRAW);
   err = glGetError(); if (err) std::cout << "GL error setup_simple_triangles()\n";

   glDisableVertexAttribArray (0);
   glDisableVertexAttribArray (1);
   if (! is_instanced) glDisableVertexAttribArray(2);
   glBindBuffer(GL_ARRAY_BUFFER, 0);
   glUseProgram(0);

   std::cout << "------------- done setup_buffer() " << std::endl;

}

void
graphical_molecule::setup_simple_triangles(Shader *shader_p, const Material &material_in) {

   material = material_in;
   shader_p->Use();
   fill_with_simple_triangles_vertices();
   // fill_with_direction_triangles();
   setup_buffers();

}

void
graphical_molecule::draw_normals(const glm::mat4 &mvp) {

   GLenum err = glGetError(); if (err) std::cout << "   error draw_normals() -- start -- "
                                                 << err << std::endl;

   if (! normals_setup)
      glGenVertexArrays(1, &vao_normals);

   glBindVertexArray(vao_normals);

   auto vec_length = [](const glm::vec3 &v) {
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

   for (unsigned int i=0; i<vertices.size(); i++) {
      float f = vec_length(vertices[i].normal);
      if (false)
         std::cout << i
                   << " vertex: " << glm::to_string(vertices[i].pos)
                   << " normal: " << glm::to_string(vertices[i].normal)
                   << " length " << f << std::endl;
      tmp_normals.push_back(vertices[i].pos);
      tmp_normals.push_back(vertices[i].pos + 0.2f * vertices[i].normal);
   }

   for (unsigned int i=0; i<triangle_vertex_indices.size(); i++) {
      const glm::vec3 &p0 = vertices[triangle_vertex_indices[i].point_id[0]].pos;
      const glm::vec3 &p1 = vertices[triangle_vertex_indices[i].point_id[1]].pos;
      const glm::vec3 &p2 = vertices[triangle_vertex_indices[i].point_id[2]].pos;
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
   }

   if (! normals_setup)
      glGenBuffers(1, &normals_buffer_id);

   glBindBuffer(GL_ARRAY_BUFFER, normals_buffer_id);
   unsigned int n_normals = tmp_normals.size();
   glBufferData(GL_ARRAY_BUFFER, n_normals * sizeof(glm::vec3), &(tmp_normals[0]), GL_STATIC_DRAW);
   glEnableVertexAttribArray(0);
   glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), 0);
   glDrawArrays(GL_LINES, 0, n_normals); // number of vertices, not number of lines
   err = glGetError(); if (err) std::cout << "   error draw_normals() post gldrawarrays "
                                          << err << std::endl;

   normals_setup = true;    
}

void
graphical_molecule::fill_with_simple_triangles_vertices() {

   vertices.resize(2 * 3);
   float scale = 0.4;

   vertices[0].pos = scale * glm::vec3(  0.f,   0.5f,  -0.2f);
   vertices[1].pos = scale * glm::vec3(  0.5f, -0.36f, -0.2f);
   vertices[2].pos = scale * glm::vec3( -0.5f, -0.36f, -0.2f);
   vertices[3].pos = scale * glm::vec3(  0.0f,  0.5f,   0.2f);
   vertices[4].pos = scale * glm::vec3(  0.5f, -0.36f,  0.2f);
   vertices[5].pos = scale * glm::vec3( -0.5f, -0.36f,  0.2f);

   vertices[0].normal = glm::vec3( 0.0f, 0.0f, 1.0f);
   vertices[1].normal = glm::vec3( 0.0f, 0.0f, 1.0f);
   vertices[2].normal = glm::vec3( 0.0f, 0.0f, 1.0f);
   vertices[3].normal = glm::vec3( 0.0f, 0.0f, 1.0f);
   vertices[4].normal = glm::vec3( 0.0f, 0.0f, 1.0f);
   vertices[5].normal = glm::vec3( 0.0f, 0.0f, 1.0f);

   vertices[0].color = glm::vec4(0.0f, 0.0f, 0.0f, 1.f);
   vertices[1].color = glm::vec4(0.2f, 0.3f, 1.0f, 1.f);
   vertices[2].color = glm::vec4(0.5f, 0.9f, 0.2f, 1.f);
   vertices[3].color = glm::vec4(0.2f, 0.2f, 0.9f, 1.f);
   vertices[4].color = glm::vec4(0.1f, 0.9f, 0.2f, 1.f);
   vertices[5].color = glm::vec4(0.9f, 0.3f, 0.2f, 1.f);

   g_triangle gt_0(0,1,2);
   g_triangle gt_1(3,4,5);
   triangle_vertex_indices.push_back(gt_0);
   triangle_vertex_indices.push_back(gt_1);

}
