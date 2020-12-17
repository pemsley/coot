
#include <iostream>
#include "LigandViewMesh.hh"
#include <glm/gtx/string_cast.hpp>  // to_string()

void
LigandViewMesh::init() {

   clear();
   first_time = true;
   draw_this_mesh = true;
   this_mesh_is_closed = false;
   vao_lines = VAO_NOT_SET;
   vao_triangles = VAO_NOT_SET;
   vao_text = VAO_NOT_SET;
}

void
LigandViewMesh::close() {

      draw_this_mesh = false;
      this_mesh_is_closed = true; // and delete the buffers if not first time,
                                  //  so don't inline this function
}

void
LigandViewMesh::clear() {

   lines_vertices.clear();
   triangles_vertices.clear();

}

void
LigandViewMesh::setup_buffers() {

   unsigned int n_lines_vertices     =     lines_vertices.size();
   unsigned int n_triangles_vertices = triangles_vertices.size();

   // -------------------------------------------- lines -------------------------------------

   if (first_time)
      glGenVertexArrays(1, &vao_lines);

   glBindVertexArray(vao_lines);
   for (unsigned int i=0; i<5; i++) {
      std::cout << i << " " << glm::to_string(lines_vertices[i]) << std::endl;
   }

   if (first_time) {
      glGenBuffers(1, &lines_buffer_id);
      glBindBuffer(GL_ARRAY_BUFFER, lines_buffer_id);
      glBufferData(GL_ARRAY_BUFFER, n_lines_vertices * sizeof(glm::vec2), &(lines_vertices[0]), GL_STATIC_DRAW);
   } else {
      glDeleteBuffers(1, &lines_buffer_id);
      glGenBuffers(1, &lines_buffer_id);
      glBindBuffer(GL_ARRAY_BUFFER, lines_buffer_id);
      glBufferData(GL_ARRAY_BUFFER, n_lines_vertices * sizeof(glm::vec2), &(lines_vertices[0]), GL_STATIC_DRAW);
   }

   // Whatever buffer is bound using glBindBuffer() affects this vertexattrib call
   glEnableVertexAttribArray(0); // position
   glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(glm::vec2), 0);

   // -------------------------------------------- triangles -------------------------------------

   if (first_time)
      glGenVertexArrays(1, &vao_triangles);

   glBindVertexArray(vao_triangles);
   if (first_time) {
      glGenBuffers(1, &triangles_buffer_id);
      glBindBuffer(GL_ARRAY_BUFFER, triangles_buffer_id);
      glBufferData(GL_ARRAY_BUFFER, n_triangles_vertices * sizeof(glm::vec2), &(triangles_vertices[0]), GL_STATIC_DRAW);
   } else {
      glDeleteBuffers(1, &triangles_buffer_id);
      glGenBuffers(1, &triangles_buffer_id);
      glBindBuffer(GL_ARRAY_BUFFER, triangles_buffer_id);
      glBufferData(GL_ARRAY_BUFFER, n_triangles_vertices * sizeof(glm::vec2), &(triangles_vertices[0]), GL_STATIC_DRAW);
   }
   // Whatever buffer is bound using glBindBuffer() affects this vertexattrib call
   glEnableVertexAttribArray(0); // position
   glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(glm::vec2), 0);

   GLenum err = glGetError();
   if (err) std::cout << "error:: LigandViewMesh::setup_buffers() " << " " << name << std::endl;

   first_time = false;
}

void
LigandViewMesh::import(const std::vector<glm::vec2> &lines_vertices_in, const std::vector<glm::vec2> &triangle_vertices_in)  {

   lines_vertices = lines_vertices_in;
   triangles_vertices = triangle_vertices_in;

   setup_buffers();
}


void
LigandViewMesh::draw(Shader *shader_p, float aspect_ratio) {

   if (vao_lines == VAO_NOT_SET) {
      // not an error necessarily
      std::cout << "LigandViewMesh::draw() vao not set yet" << std::endl;
      return;
   }

   shader_p->Use();
   shader_p->set_float_for_uniform("aspect_ratio", aspect_ratio);

   // ----------------------------- lines --------------------------------------

   glBindVertexArray(vao_lines);
   GLenum err = glGetError();
   if (err) std::cout << "error:: LigandViewMesh::draw() " << shader_p->name << " " << name
                      << " glBindVertexArray() vao_lines " << vao_lines << " with GL err " << err << std::endl;

   unsigned int n_vertices = lines_vertices.size();

   glBindBuffer(GL_ARRAY_BUFFER, lines_buffer_id);
   glEnableVertexAttribArray(0);

   // std::cout << "debug:: LigandViewMesh::draw() glDrawArrays() n_vertices " << n_vertices << std::endl;
   glDrawArrays(GL_LINES, 0, n_vertices);
   err = glGetError();
   if (err) std::cout << "error:: LigandViewMesh::draw() " << shader_p->name << " " << name
                      << " glDrawArrays" << " with GL err " << err << std::endl;

   glDisableVertexAttribArray(0);

   // ----------------------------- triangles --------------------------------------

   glBindVertexArray(vao_triangles);
   err = glGetError();
   if (err) std::cout << "error:: LigandViewMesh::draw() " << shader_p->name << " " << name
                      << " glBindVertexArray() vao_triangles " << vao_triangles
                      << " with GL err " << err << std::endl;

   n_vertices = triangles_vertices.size();

   glBindBuffer(GL_ARRAY_BUFFER, triangles_buffer_id);
   glEnableVertexAttribArray(0);

   // std::cout << "debug:: LigandViewMesh::draw() glDrawArrays() n_vertices " << n_vertices << std::endl;
   glDrawArrays(GL_TRIANGLES, 0, n_vertices);
   err = glGetError();
   if (err) std::cout << "error:: LigandViewMesh::draw() " << shader_p->name << " " << name
                      << " glDrawArrays" << " with GL err " << err << std::endl;

   glDisableVertexAttribArray(0);


   // ----------------------------- text --------------------------------------

   glBindVertexArray(vao_text);

   // ----------------------------- done --------------------------------------

   glDisableVertexAttribArray(0);
   glUseProgram(0);

}

