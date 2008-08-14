


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
   short int hardware_stereo_flag;
   bool script_is_python_flag;
   int port;
   std::string hostname;
   std::string ccp4_project;
   short int try_listener;
   short int do_graphics;
   bool disable_state_script_writing;
   command_line_data() { 
     hardware_stereo_flag = 0; // default off
     port = 0;
     hostname = "";
     try_listener = 0;
     do_graphics = 1; // use graphics by default
     disable_state_script_writing = 0; // don't disable, by default
     script_is_python_flag = 0;
   }
   void handle_immediate_settings(); 
};


command_line_data
parse_command_line(int argc, char ** argv ); 


void 
handle_command_line_data(command_line_data cld); 

