

namespace coot {
   class atom_attribute_setting_help_t {
   public:
     enum { UNSET, IS_FLOAT, IS_STRING, IS_INT};
      short int type;
      int i;
      float val;
      std::string s;
      atom_attribute_setting_help_t(const std::string &s_in) {
	 s = s_in;
	 type = IS_STRING;
      }
      atom_attribute_setting_help_t(float v) {
	 val = v;
	 type = IS_FLOAT;
      }
      atom_attribute_setting_help_t(int iin) {
	i = iin;
	type = IS_INT;
      }
      atom_attribute_setting_help_t() {
	 type = UNSET;
      }
   };

   class atom_attribute_setting_t {
   public: 
     atom_spec_t atom_spec;
     std::string attribute_name;
     atom_attribute_setting_help_t attribute_value;
     atom_attribute_setting_t(const std::string &chain_id_in, 
			      int resno_in, 
			      const std::string &inscode_in, 
			      const std::string &atom_name_in, 
			      const std::string &alt_conf_in, 
			      const std::string &attribute_name_in, 
			      const atom_attribute_setting_help_t &att_val) {
       atom_spec = atom_spec_t(chain_id_in, resno_in, inscode_in, atom_name_in, alt_conf_in);
       attribute_name = attribute_name_in;
       attribute_value = att_val;
     } 
   };
   
}
