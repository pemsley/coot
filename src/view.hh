

namespace coot {
   
   // This is for a Pymol-like feature to switch between views:
   // 
   class view_info_t {
   public:
      float zoom;
      coot::Cartesian rotation_centre;
      std::string view_name;
      std::string description;
      bool is_simple_spin_view_flag;
      bool is_action_view_flag;
      int n_spin_steps;
      float degrees_per_step;
      float quat[4];
      std::string action;
      view_info_t(float *quat_in, const coot::Cartesian &rot_centre_in,
		  float zoom_in, const std::string &view_name_in) {
	is_simple_spin_view_flag = 0;
	is_action_view_flag = 0;
	 zoom = zoom_in;
	 rotation_centre = rot_centre_in;
	 view_name = view_name_in;
	 for (int i=0; i<4; i++) 
	    quat[i] = quat_in[i];
      }
      view_info_t() {
      	is_simple_spin_view_flag = 0;
	is_action_view_flag = 0;
      }
      // a spin view 
      view_info_t(const std::string &view_name_in, int nsteps, float degrees_total) {
	is_simple_spin_view_flag = 1;
	is_action_view_flag = 0;
	view_name = view_name_in;
	n_spin_steps = nsteps;
	if (n_spin_steps > 0) 
	  degrees_per_step = degrees_total/n_spin_steps;
	else 
	  degrees_per_step = 0.5;
      }
      // an action view
      view_info_t(const std::string &view_name_in, const std::string &funct) {
	is_action_view_flag = 1;
	view_name = view_name_in;
	action = funct;
      }
      bool matches_view (const coot::view_info_t &view) const;
      float quat_length() const;
      void add_description(const std::string &descr) { 
	description = descr;
      }
      static view_info_t interpolate(const view_info_t &view1,
				     const view_info_t &view2,
				     int n_steps);
      static float dot_product(const view_info_t &view1,
			       const view_info_t &view2);

      friend std::ostream& operator<<(std::ostream &stream, 
				      view_info_t &view);
   };
   std::ostream& operator<<(std::ostream &stream, 
			    view_info_t &view);
}
