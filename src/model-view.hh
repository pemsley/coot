
namespace coot { 
   class model_view_residue_button_info_t {
   public:
      model_view_residue_button_info_t(){} // for new allocator

      model_view_residue_button_info_t(const std::string &lab,
				       mmdb::Residue *res) {
	 button_label = lab;
	 residue_spec = res;
      } 
      std::string button_label;
      // mmdb::Residue *residue; No. This can go out of date.
      residue_spec_t residue_spec;
   };

   class model_view_atom_tree_item_info_t {
   public:
      // model_view_atom_tree_item_info_t() {}  // not needed?
      model_view_atom_tree_item_info_t(const std::string &label,
				       mmdb::Residue *res) {
	 button_label = label;
	 residue_spec = res;
      }
      std::string button_label;
      // mmdb::Residue *residue;
      residue_spec_t residue_spec;
   };

   class model_view_atom_tree_chain_t {
   public:
      model_view_atom_tree_chain_t() {} // for new allocator
      model_view_atom_tree_chain_t(const std::string &chain_id_in) {
	 chain_id = chain_id_in;
      }
      void add_residue(const model_view_atom_tree_item_info_t &res) {
	 tree_residue.push_back(res);
      } 
      std::vector<model_view_atom_tree_item_info_t> tree_residue;
      std::string chain_id;
   };

   // old 
   class model_view_atom_button_info_t {
   public:
      model_view_atom_button_info_t() {} // for new allocator
      model_view_atom_button_info_t(const std::string &label,
				    mmdb::Atom *atom_in) {
	 button_label = label;
	 atom = atom_in;
      }
      std::string button_label;
      mmdb::Atom *atom;
   };

}
