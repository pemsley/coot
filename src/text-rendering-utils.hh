
#include "Shader.hh"

void  draw_hud_text(int widget_width, int widget_height, Shader &shader);
void setup_hud_text(int widget_width, int widget_height, Shader &shader);

void RenderText(Shader &shader, std::string text, GLfloat x, GLfloat y, GLfloat scale, glm::vec3 color);

void render_atom_label(Shader &shader, std::string text, glm::vec3 projected_position,
                       GLfloat scale, glm::vec3 color);

