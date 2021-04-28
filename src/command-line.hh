
#ifndef COMMAND_LINE_HH
#define COMMAND_LINE_HH

#include <string>
#include <vector>

class command_line_data { 

public:
   std::vector<std::string> coords;
   std::vector<std::string> maps;
   std::vector<std::string> datasets;
   std::vector<std::string> auto_datasets;
   std::vector<std::string> script;
   std::vector<std::string> dictionaries;
   std::vector<std::string> command; // strings to to be evaluated
				     // from the command line
   std::vector<std::string> accession_codes;
   std::vector<std::string> emdb_codes;
   std::vector<std::string> comp_ids;
   short int hardware_stereo_flag;
   bool script_is_python_flag;
   int port;
   std::string hostname;
   std::string ccp4_project;
   std::string title;
   short int try_listener;
   short int do_graphics;
   short int small_screen_display;
   bool disable_state_script_writing;
   bool use_splash_screen;
   bool update_self;
   std::string alternate_splash_screen_file_name; 
   bool run_internal_tests_and_exit;
   bool em_mode;
   bool use_gtkbuilder;
   command_line_data() { 
     hardware_stereo_flag = 0; // default off
     port = 0;
     try_listener = 0;
     update_self = 0;
     do_graphics = 1; // use graphics by default
     disable_state_script_writing = 0; // don't disable, by default
     script_is_python_flag = 0;
     small_screen_display  = 0; // default is no small screen
     use_splash_screen = 1;
     run_internal_tests_and_exit = 0;
     em_mode = false;
     use_gtkbuilder = false;
   }
   void handle_immediate_settings();
   void roberto_pdbs(int argc, char **argv); // add any pdb files not alread added with --pdb/coords/xyzin
};


command_line_data
parse_command_line(int argc, char ** argv ); 


#ifdef __cplusplus
extern "C"
#endif
void 
handle_command_line_data(command_line_data cld); 

#endif // COMMAND_LINE_HH
