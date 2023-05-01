
#ifndef DRAW_HH
#define DRAW_HH

void setup_monkey_head();

void draw_monkey_head();

void draw_one_triangle();

void draw_molecular_triangles(GtkWidget *widget);

void setup_for_single_triangle();

void draw_single_triangle();

void gtk3_draw_molecules();
void test_gtk3_adjustment_changed(GtkAdjustment *adj, GtkWidget *window);

struct shader_program_source {
   std::string VertexSource;
   std::string FragmentSource;
};

shader_program_source parse_shader(const std::string &file_name);
unsigned int CreateShader(const std::string &vertex_shader, const std::string &fragment_shader);

// static int programID_global;
// static GLuint VertexArrayID;
// static GLuint molecules_imol_VertexArrayID[10];
// static int location_global;



#endif // DRAW_HH
