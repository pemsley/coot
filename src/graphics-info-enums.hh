
namespace coot {
   enum {NEW_COORDS_UNSET = 0,       // moving_atoms_asc_type values
	 NEW_COORDS_ADD = 1,                 // not used?
	 NEW_COORDS_REPLACE = 2,
	 NEW_COORDS_REPLACE_CHANGE_ALTCONF = 3,
	 NEW_COORDS_INSERT = 4,
         NEW_COORDS_INSERT_CHANGE_ALTCONF = 5};
   enum {STATE_SCM = 1, STATE_PYTHON = 2};
   enum undo_type { UNDO, REDO};
   enum display_mode_type {MONO_MODE=0, HARDWARE_STEREO_MODE=1,
			   SIDE_BY_SIDE_STEREO=2,  // cross-eye
			   DTI_SIDE_BY_SIDE_STEREO=3,
                           SIDE_BY_SIDE_STEREO_WALL_EYE=4,
                           ZALMAN_STEREO=5};
   enum accept_reject_text_type { CHI_SQUAREDS, CHIRAL_CENTRES};
   enum chooser_selector_type { OLD_STYLE, CHOOSER_STYLE };
   enum chooser_overwrite_type { CHOOSER_OVERWRITE, CHOOSER_OVERWRITE_PROTECT };
   enum accept_reject_dialog_type { DIALOG, DIALOG_DOCKED };
   enum accept_reject_dialog_docked_type { DIALOG_DOCKED_HIDE = 0,
                                           DIALOG_DOCKED_SHOW = 1};
   enum fixed_atom_pick_state_t { FIXED_ATOM_NO_PICK = 0,
				  FIXED_ATOM_FIX = 1,
				  FIXED_ATOM_UNFIX = 2 };
   enum ncs_matrix_type { NCS_SSM  = 0,
                          NCS_LSQ  = 1,
			  NCS_LSQ2 = 2};
   namespace model_toolbar {
     enum toolbar_position_type { RIGHT  = 0,
                                  LEFT   = 1,
                                  TOP    = 2,
                                  BOTTOM = 3};
   }
   namespace refmac {
     enum refmac_refinement_method_type { RESTRAINED     = 0 ,
					  RIGID_BODY     = 1 ,
					  RESTRAINED_TLS = 2 };
     enum refmac_phase_input_type { NO_PHASES = 0,
				    PHASE_FOM = 1 ,
				    HL        = 2 ,
                                    SAD       = 3 };
     enum refmac_use_tls_type { TLS_OFF = 0,
				TLS_ON  = 1};
     enum refmac_use_twin_type { TWIN_OFF = 0,
				 TWIN_ON  = 1};
     enum refmac_use_sad_type { SAD_OFF = 0,
				SAD_ON = 1};
     enum refmac_use_ncs_type { NCS_OFF = 0,
				NCS_ON  = 1};
     enum refmac_use_intensities_type { AMPLITUDES   = 0,
					INTENSITIES  = 1};
     enum refmac_used_mtz_file_type { MAP = 0,
				      MTZ = 1};

     class sad_atom_info_t {
     public:
	std::string atom_name;
	float fp;
	float fpp;
	float lambda;
        sad_atom_info_t(const std::string &atom_name_in, float &fp_in,
			float &fpp_in, float &lambda_in) :
           atom_name(atom_name_in), fp(fp_in), fpp(fpp_in), lambda(lambda_in) { }
     };
   }

   enum nomenclature_error_handle_type {
     AUTO_CORRECT, IGNORE, PROMPT};

   enum scripting_language_type { SCRIPT_UNSET = -1,
				  SCHEME_SCRIPT = 1,
				  PYTHON_SCRIPT = 2};

   // we can (only) use scripting_function for things that return a
   // command_arg_t (bool, int, float, string) but that should cover
   // most cases.
   coot::command_arg_t scripting_function(const std::string &function_name,
					  const std::vector<coot::command_arg_t> &args);


#ifdef DO_GEOMETRY_GRAPHS
   void set_validation_graph(int imol, coot::geometry_graph_type type, GtkWidget *dialog);
   GtkWidget *get_validation_graph(int imol, coot::geometry_graph_type type);
#endif

   class coord_orth_triple {
   public:
      clipper::Coord_orth p1;
      clipper::Coord_orth p2;
      clipper::Coord_orth p3;
   };

   class intermediate_atom_distance_t {
      coot::Cartesian static_position;
      mmdb::Atom *dynamic_atom;
      bool static_pos_filled_flag;

   public:
      intermediate_atom_distance_t() {
         dynamic_atom = 0;
         static_pos_filled_flag = 0;
      }
      explicit intermediate_atom_distance_t(const coot::Cartesian &pt) : static_position(pt) {
         dynamic_atom = 0;
         static_pos_filled_flag = 1;
      }
      explicit intermediate_atom_distance_t(mmdb::Atom *at) : dynamic_atom(at) {
         static_pos_filled_flag = 0;
      }
      void draw_dynamic_distance() const;
      bool static_position_is_filled() const { return static_pos_filled_flag; }
      bool atom_is_filled() const {
         return (dynamic_atom != 0);
      }
      void add_atom(mmdb::Atom *at) {
         dynamic_atom = at;
      }
      void add_static_point(Cartesian &pt) {
         static_position = pt;
         static_pos_filled_flag = 1;
      }
      bool filled() const {
         return (static_pos_filled_flag && dynamic_atom);
      }
   };

   class ramachandran_points_container_t {

      std::vector<std::pair<std::string, clipper::Coord_orth> > points;


   public:
      ramachandran_points_container_t() {};
      void clear() {
	 points.resize(0);
      }
      std::pair<short int, clipper::Coord_orth> get(const std::string &atom_name) const {
	 std::pair<short int, clipper::Coord_orth> v;
	 v.first = 0;
	 for (unsigned int i=0; i<points.size(); i++) {
	    if (atom_name == points[i].first) {
	       v.first = 1;
	       v.second = points[i].second;
	       break;
	    }
	 }
	 return v;
      }
      void add(const std::string &s, const clipper::Coord_orth &p) {
	 points.push_back(std::pair<std::string, clipper::Coord_orth> (s,p));
      }

   };

   class graph_rotamer_info_t {
   public:
     std::string chain_id;
     int resno;
     std::string inscode;
     float probability;
     std::string rotamer_name;
     graph_rotamer_info_t(const std::string &chain_id_in, int resno_in, const std::string &inscode_in, float prob_in, const std::string &rotamer_name_in) {
       chain_id = chain_id_in;
       resno = resno_in;
       inscode = inscode_in;
       probability = prob_in;
       rotamer_name = rotamer_name_in;
     }
   };

   // To pass rotamer info back to scripting layer, for testing.
   // Hmmmm.. confusing names, perhaps.
   //
   class rotamer_graphs_info_t {
   public:
     std::vector<graph_rotamer_info_t> info;
   };

   class diff_map_peak_helper_data {
   public:
      int ipeak;
      clipper::Coord_orth pos;
   };


   enum tube_end_t { NO_ENDS, FLAT_ENDS, ROUND_ENDS};  // use gluDisk or gluSphere.


   class console_display_commands_t {
   public:
     bool display_commands_flag;
     bool hilight_flag;
     bool hilight_colour_flag;
     int colour_prefix;
     console_display_commands_t() {
       // hilighting
       display_commands_flag = 1;
       hilight_flag = 1;
       colour_prefix = 4;
       hilight_colour_flag = 0;
     }
   };

  // for preferences
  class preference_info_t {

  public:
    int preference_type;   // e.g. PREFERENCES_bla
    int ivalue1;
    int ivalue2;
    float fvalue1;
    float fvalue2;
    float fvalue3;
  };

  class preferences_icon_info_t {

  public:
    int icon_pos;
    std::string icon_filename;
    std::string icon_text;
    std::string icon_widget;
    int show_hide_flag;
    int default_show_flag;
    preferences_icon_info_t(int icon_pos_in,
			    std::string icon_filename_in,
			    std::string icon_text_in,
			    std::string icon_widget_in,
			    int show_hide_flag_in,
			    int default_show_flag_in) {
	icon_pos = icon_pos_in;
	icon_filename = icon_filename_in;
	icon_text = icon_text_in;
	icon_widget = icon_widget_in;
	show_hide_flag = show_hide_flag_in;
	default_show_flag = default_show_flag_in;
    }
    void hide() {
	show_hide_flag = 0;
    }
    void show() {
	show_hide_flag = 1;
    }
  };

  class command_line_commands_t {
  public:
    std::vector<std::string> commands;
    bool is_python;
    command_line_commands_t() {
      is_python = 0;
    }
  };


  class ScreenVectors {
  public:
    ScreenVectors();
    Cartesian screen_x;
    Cartesian screen_y;
    Cartesian screen_z;
  };

  class saved_strand_info_t {
  public:
    coot::residue_spec_t start;
    coot::residue_spec_t end;
    int strand_idx;
    saved_strand_info_t(const coot::residue_spec_t &s, const coot::residue_spec_t &e, int idx) {
      start = s; end = e; strand_idx = idx;
    }
  };


#ifdef USE_LIBCURL
  class simple_curl_handler_t {
  public:
    CURL * c;
    std::string file_name;
    bool stop_it;
    simple_curl_handler_t(CURL *cin, std::string f) {
      file_name = f;
      c = cin;
      stop_it=0;
    }
    bool stop_is_set() {
      return stop_it;
    }
    void set_stop() {
      stop_it = 1;
    }
  };
#endif // USE_LIBCURL

} // namespace coot

enum { IN_STEREO_MONO = 0, 
       IN_STEREO_HARDWARE_STEREO=1, 
       IN_STEREO_ZALMAN_RIGHT=5, 
       IN_STEREO_ZALMAN_LEFT=6, 
       IN_STEREO_SIDE_BY_SIDE_LEFT=10,
       IN_STEREO_SIDE_BY_SIDE_RIGHT=11
};
