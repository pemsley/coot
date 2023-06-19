
#include <iostream>
#include "LigandViewMesh.hh"
#include <glm/gtx/string_cast.hpp>  // to_string()

void
LigandViewMesh::init() {

   clear();
   first_time = true;
   draw_this_mesh = true;
   this_mesh_is_closed = false;
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

   triangles_vertices.clear();

}

void
LigandViewMesh::setup_buffers() {

   unsigned int n_triangles_vertices = triangles_vertices.size();

   // -------------------------------------------- lines -------------------------------------

   if (n_triangles_vertices == 0) return;

   if (false) {
      std::cout << "debug:: LigandViewMesh::setup_buffers vao_triangles:  " << vao_triangles << std::endl;
      std::cout << "debug:: LigandViewMesh::setup_buffers first_time: " << first_time << std::endl;
   }

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
LigandViewMesh::import(const std::vector<glm::vec2> &triangle_vertices_in)  {

   triangles_vertices = triangle_vertices_in;

   setup_buffers();
}


void
LigandViewMesh::draw(Shader *shader_p, float widget_height, float widget_width) {

   // std::cout << "---------- LigandViewMesh::draw() start " << std::endl;

   shader_p->Use();
   float aspect_ratio = widget_width/widget_height;
   if (aspect_ratio < 1.0) aspect_ratio = 1.0;
   shader_p->set_float_for_uniform("aspect_ratio", aspect_ratio);

   // ----------------------------- triangles --------------------------------------

   if (vao_triangles == VAO_NOT_SET) {
      std::cout << "LigandViewMesh::draw() vao_triangles not set yet" << std::endl;
      return;
   }

   {
      glBindVertexArray(vao_triangles);
      GLenum err = glGetError();
      if (err) std::cout << "error:: LigandViewMesh::draw() " << shader_p->name << " " << name
                         << " glBindVertexArray() vao_triangles " << vao_triangles
                         << " with GL err " << err << std::endl;

      unsigned int n_vertices = triangles_vertices.size();

      glBindBuffer(GL_ARRAY_BUFFER, triangles_buffer_id); // remove this when fixed.
      glEnableVertexAttribArray(0);

      // std::cout << "debug:: LigandViewMesh::draw() triangles glDrawArrays() n_vertices " << n_vertices << std::endl;
      glDrawArrays(GL_TRIANGLES, 0, n_vertices);
      err = glGetError();
      if (err) std::cout << "error:: LigandViewMesh::draw() " << shader_p->name << " " << name
                         << " glDrawArrays" << " with GL err " << err << std::endl;

      glDisableVertexAttribArray(0);
   }


   // text not done here?
   // // ----------------------------- text --------------------------------------

   // glBindVertexArray(vao_text);

   // // ----------------------------- done --------------------------------------

   // glDisableVertexAttribArray(0);


   glUseProgram(0);

}

