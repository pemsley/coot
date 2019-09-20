
GtkWidget *my_gtkglarea(GtkWidget *vbox);
void my_glarea_add_signals_and_events(GtkWidget *glarea);

void init_central_cube();
void draw_central_cube(GtkGLArea *glarea);

glm::vec4 new_unproject(float z);

glm::vec4 new_unproject(float mouse_x, float mouse_y, float z);

glm::mat4 get_molecule_mvp();
