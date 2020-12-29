   class pulse_data_t {
   public:
      int n_pulse_steps;
      int n_pulse_steps_max;
      pulse_data_t(int n1, int n2) {
         n_pulse_steps = n1;
         n_pulse_steps_max = n2;
      }
   };

      {
         auto identification_pulse_func = [] (GtkWidget *widget,
                                              GdkFrameClock *frame_clock,
                                              gpointer data) {

                                             gboolean continue_status = 1;
                                             pulse_data_t *pulse_data = reinterpret_cast<pulse_data_t *>(data);
                                             pulse_data->n_pulse_steps += 1;
                                             // std::cout << "pulse: " << pulse_data->n_pulse_steps << std::endl;
                                             if (pulse_data->n_pulse_steps > pulse_data->n_pulse_steps_max) {
                                                continue_status = 0;
                                                lines_mesh_for_identification_pulse.clear();
                                             } else {
                                                float ns = pulse_data->n_pulse_steps;
                                                lines_mesh_for_identification_pulse.update_buffers_for_pulse(ns);
                                             }
                                             graphics_draw();
                                             return gboolean(continue_status);
                                        };

         pulse_data_t *pulse_data = new pulse_data_t(0, 30);
         gpointer user_data = reinterpret_cast<void *>(pulse_data);
         glm::vec3 pos = cartesian_to_glm(current_centre);
         gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0]));
         lines_mesh_for_identification_pulse.setup_pulse(pos, &shader_for_lines_pulse);
         gtk_widget_add_tick_callback(glareas[0], identification_pulse_func, user_data, NULL);

      }
