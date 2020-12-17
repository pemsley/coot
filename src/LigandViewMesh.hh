
#ifndef LIGAND_VIEW_MESH_HH
#define LIGAND_VIEW_MESH_HH

#include <vector>
#include "Shader.hh"

class LigandViewMesh {

   // contains vectors for both lines and triangles. Don't use indexing to draw.   

   enum { VAO_NOT_SET = 99999999 };
   GLuint vao;
   GLuint lines_buffer_id;
   GLuint triangles_buffer_id;
   std::vector<glm::vec2> lines_vertices;
   std::vector<glm::vec2> triangles_vertices;
   // I need a container for text too.
   std::string name;
   void init();
   bool first_time;
   bool draw_this_mesh;
   bool this_mesh_is_closed;
   void setup_buffers();

public:
   LigandViewMesh() { init(); }
   explicit LigandViewMesh(const std::string &n) : name(n) { init(); }
   std::string get_name() const { return name; }
   void import(const std::vector<glm::vec2> &lines_vertices, const std::vector<glm::vec2> &triangle_vertices);
   bool get_draw_status() const { return draw_this_mesh; }
   void draw(Shader *shader_p, float aspect_ratio);
   void close();
   void clear();
};


#endif // LIGAND_VIEW_MESH_HH
