

#include <string>
#include <vector>
#include <mmdb2/mmdb_manager.h>

   // ------------------------------------------------------------------------
   //                  chem_mods
   // ------------------------------------------------------------------------

namespace coot {

   // a container for the data_mod_list chem_mods (used to be simply
   // chem_mod, but that name is now used as a class that contains the
   // actual chem mods (with lists of atoms, bonds, angles and so on).
   // 
   class list_chem_mod {
   public:
      std::string name;
      std::string id;
      std::string group_id;
      std::string comp_id;
      list_chem_mod(const std::string &id_in,
		    const std::string &name_in,
		    const std::string &comp_id_in,
		    const std::string &group_id_in) :
         name(name_in), id(id_in), group_id(group_id_in), comp_id(comp_id_in) {}
      friend std::ostream& operator<<(std::ostream &s, list_chem_mod mod);
   };
   std::ostream& operator<<(std::ostream &s, list_chem_mod mod);

   enum chem_mod_function_t { CHEM_MOD_FUNCTION_UNSET,
			      CHEM_MOD_FUNCTION_ADD,
			      CHEM_MOD_FUNCTION_CHANGE,
			      CHEM_MOD_FUNCTION_DELETE };
   
   class chem_mod_atom {
   public:
      chem_mod_atom(const std::string &function_in,
		    const std::string &atom_id_in,
		    const std::string &new_atom_id_in,
		    const std::string &new_type_symbol_in,
		    const std::string &new_type_energy_in,
		    mmdb::realtype new_partial_charge_in) {
	 function = CHEM_MOD_FUNCTION_UNSET;
	 if (function_in == "add")
	    function = CHEM_MOD_FUNCTION_ADD;
	 if (function_in == "delete")
	    function = CHEM_MOD_FUNCTION_DELETE;
	 if (function_in == "change")
	    function = CHEM_MOD_FUNCTION_CHANGE;
	 atom_id = atom_id_in;
	 new_atom_id = new_atom_id_in;
	 new_type_symbol = new_type_symbol_in;
	 new_type_energy = new_type_energy_in;
	 new_partial_charge = new_partial_charge_in;
      }
      chem_mod_function_t function;
      std::string atom_id;
      std::string new_atom_id;
      std::string new_type_symbol;
      std::string new_type_energy;
      mmdb::realtype new_partial_charge;
      friend std::ostream& operator<<(std::ostream &s, const chem_mod_atom &a);
   };
   std::ostream& operator<<(std::ostream &s, const chem_mod_atom &a);

   class chem_mod_tree {
      chem_mod_function_t function;
      std::string atom_id;
      std::string atom_back;
      std::string back_type;
      std::string atom_forward;
      std::string connect_type;
   public:
      chem_mod_tree (const std::string &function_in,
		     const std::string &atom_id_in,
		     const std::string &atom_back_in,
		     const std::string &back_type_in,
		     const std::string &atom_forward_in,
		     const std::string &connect_type_in) {
	 function = CHEM_MOD_FUNCTION_UNSET;
	 if (function_in == "add")
	    function = CHEM_MOD_FUNCTION_ADD;
	 if (function_in == "delete")
	    function = CHEM_MOD_FUNCTION_DELETE;
	 if (function_in == "change")
	    function = CHEM_MOD_FUNCTION_CHANGE;
	 atom_id = atom_id_in;
	 atom_back = atom_back_in;
	 back_type = back_type_in;
	 atom_forward = atom_forward_in;
	 connect_type = connect_type_in;
      }
      friend std::ostream& operator<<(std::ostream &s, const chem_mod_tree &a);
   };
   std::ostream& operator<<(std::ostream &s, const chem_mod_tree &a);

   class chem_mod_bond {
   public:
      chem_mod_bond(const std::string &function_in,
		    const std::string &atom_id_1_in,
		    const std::string &atom_id_2_in,
		    const std::string &new_type_in,
		    mmdb::realtype new_value_dist_in,
		    mmdb::realtype new_value_dist_esd_in) {
	 function = CHEM_MOD_FUNCTION_UNSET;
	 if (function_in == "add")
	    function = CHEM_MOD_FUNCTION_ADD;
	 if (function_in == "delete")
	    function = CHEM_MOD_FUNCTION_DELETE;
	 if (function_in == "change")
	    function = CHEM_MOD_FUNCTION_CHANGE;
	 atom_id_1 = atom_id_1_in;
	 atom_id_2 = atom_id_2_in;
	 new_type = new_type_in;
	 new_value_dist = new_value_dist_in;
	 new_value_dist_esd = new_value_dist_esd_in;
      }
      chem_mod_function_t function;
      std::string atom_id_1;
      std::string atom_id_2;
      std::string new_type;
      mmdb::realtype new_value_dist;
      mmdb::realtype new_value_dist_esd;
      friend std::ostream& operator<<(std::ostream &s, const chem_mod_bond &a);
   };
   std::ostream& operator<<(std::ostream &s, const chem_mod_bond &a);

   class chem_mod_angle {
   public:
      chem_mod_angle(const std::string &function_in,
		     const std::string &atom_id_1_in,
		     const std::string &atom_id_2_in,
		     const std::string &atom_id_3_in,
		     mmdb::realtype new_value_angle_in,
		     mmdb::realtype new_value_angle_esd_in) {
	 function = CHEM_MOD_FUNCTION_UNSET;
	 if (function_in == "add")
	    function = CHEM_MOD_FUNCTION_ADD;
	 if (function_in == "delete")
	    function = CHEM_MOD_FUNCTION_DELETE;
	 if (function_in == "change")
	    function = CHEM_MOD_FUNCTION_CHANGE;
	 atom_id_1 = atom_id_1_in;
	 atom_id_2 = atom_id_2_in;
	 atom_id_3 = atom_id_3_in;
	 new_value_angle = new_value_angle_in;
	 new_value_angle_esd = new_value_angle_esd_in;
      }
      chem_mod_function_t function;
      std::string atom_id_1;
      std::string atom_id_2;
      std::string atom_id_3;
      std::string new_type;
      mmdb::realtype new_value_angle;
      mmdb::realtype new_value_angle_esd;
      friend std::ostream& operator<<(std::ostream &s, const chem_mod_angle &a);
   };
   std::ostream& operator<<(std::ostream &s, const chem_mod_angle &a);

   class chem_mod_tor {
   public:
      chem_mod_tor(const std::string &function_in,
		   const std::string &atom_id_1_in,
		   const std::string &atom_id_2_in,
		   const std::string &atom_id_3_in,
		   const std::string &atom_id_4_in,
		   mmdb::realtype new_value_angle_in,
		   mmdb::realtype new_value_angle_esd_in,
		   int new_period_in) {
	 function = CHEM_MOD_FUNCTION_UNSET;
	 if (function_in == "add")
	    function = CHEM_MOD_FUNCTION_ADD;
	 if (function_in == "delete")
	    function = CHEM_MOD_FUNCTION_DELETE;
	 if (function_in == "change")
	    function = CHEM_MOD_FUNCTION_CHANGE;
	 atom_id_1 = atom_id_1_in;
	 atom_id_2 = atom_id_2_in;
	 atom_id_3 = atom_id_3_in;
	 atom_id_4 = atom_id_4_in;
	 new_value_angle = new_value_angle_in;
	 new_value_angle_esd = new_value_angle_esd_in;
	 new_period = new_period_in;
      }
      chem_mod_function_t function;
      std::string atom_id_1;
      std::string atom_id_2;
      std::string atom_id_3;
      std::string atom_id_4;
      std::string new_type;
      mmdb::realtype new_value_angle;
      mmdb::realtype new_value_angle_esd;
      int new_period;
      friend std::ostream& operator<<(std::ostream &s, const chem_mod_tor &a);
   };
   std::ostream& operator<<(std::ostream &s, const chem_mod_tor &a);

   class chem_mod_plane {
   public:
      chem_mod_plane(const std::string &plane_id_in,
		     const std::string &function_in) {
	 function = CHEM_MOD_FUNCTION_UNSET;
	 if (function_in == "add")
	    function = CHEM_MOD_FUNCTION_ADD;
	 if (function_in == "delete")
	    function = CHEM_MOD_FUNCTION_DELETE;
	 if (function_in == "change")
	    function = CHEM_MOD_FUNCTION_CHANGE;
	 plane_id = plane_id_in;
      }
      chem_mod_function_t function;
      std::string plane_id;
      std::vector<std::pair<std::string, mmdb::realtype> > atom_id_esd;
      void add_atom(const std::string &atom_id, mmdb::realtype esd) {
	 std::pair<std::string, mmdb::realtype> p(atom_id, esd);
	 atom_id_esd.push_back(p);
      }
      friend std::ostream& operator<<(std::ostream &s, const chem_mod_plane &a);
   };
   std::ostream& operator<<(std::ostream &s, const chem_mod_plane &a);

   class chem_mod_chir {
   public:
      chem_mod_chir(const std::string &function_in,
		    const std::string &atom_id_centre_in,
		    const std::string &atom_id_1_in,
		    const std::string &atom_id_2_in,
		    const std::string &atom_id_3_in,
		    int new_volume_sign_in) {
	 function = CHEM_MOD_FUNCTION_UNSET;
	 if (function_in == "add")
	    function = CHEM_MOD_FUNCTION_ADD;
	 if (function_in == "delete")
	    function = CHEM_MOD_FUNCTION_DELETE;
	 if (function_in == "change")
	    function = CHEM_MOD_FUNCTION_CHANGE;
	 atom_id_centre = atom_id_centre_in;
	 atom_id_1 = atom_id_1_in;
	 atom_id_2 = atom_id_2_in;
	 atom_id_3 = atom_id_3_in;
	 new_volume_sign = new_volume_sign_in;
      }
      chem_mod_function_t function;
      std::string atom_id_centre;
      std::string atom_id_1;
      std::string atom_id_2;
      std::string atom_id_3;
      int new_volume_sign;
      friend std::ostream& operator<<(std::ostream &s, const chem_mod_chir &a);
   };
   std::ostream& operator<<(std::ostream &s, const chem_mod_chir &a);




      class chem_mod {

      public:
	 chem_mod() {};
	 std::vector<chem_mod_atom>  atom_mods;
	 std::vector<chem_mod_tree>  tree_mods;
	 std::vector<chem_mod_bond>  bond_mods;
	 std::vector<chem_mod_angle> angle_mods;
	 std::vector<chem_mod_tor>   tor_mods;
	 std::vector<chem_mod_plane> plane_mods;
	 std::vector<chem_mod_chir>  chir_mods;
	 void add_mod_atom(const chem_mod_atom &chem_atom) {
	    atom_mods.push_back(chem_atom);
	 }
	 void add_mod_tree(const chem_mod_tree &chem_tree) {
	    tree_mods.push_back(chem_tree);
	 }
	 void add_mod_bond(const chem_mod_bond &chem_bond) {
	    bond_mods.push_back(chem_bond);
	 }
	 void add_mod_angle(const chem_mod_angle &chem_angle) {
	    angle_mods.push_back(chem_angle);
	 }
	 void add_mod_tor(const chem_mod_tor &chem_tor) {
	    tor_mods.push_back(chem_tor);
	 }
	 void add_mod_plane(const chem_mod_plane &chem_plane) {
	    plane_mods.push_back(chem_plane);
	 }
	 void add_mod_chir(const chem_mod_chir &chem_chir) {
	    chir_mods.push_back(chem_chir);
	 }
	 void add_plane_atom(const std::string &plane_id,
			     const std::string &function,
			     const std::string &atom_name,
			     double dist) {
	    bool done = false;
	    for (unsigned int iplane=0; iplane<plane_mods.size(); iplane++) {
	       if (plane_mods[iplane].plane_id == plane_id) {
		  chem_mod_function_t cmft = CHEM_MOD_FUNCTION_UNSET;
		  if (function == "add")    cmft = CHEM_MOD_FUNCTION_ADD;
		  if (function == "change") cmft = CHEM_MOD_FUNCTION_CHANGE;
		  if (function == "delete") cmft = CHEM_MOD_FUNCTION_DELETE;
		  if (plane_mods[iplane].function == cmft) {
		     plane_mods[iplane].add_atom(atom_name, dist);
		     done = true;
		     break;
		  }
	       }
	    }
	    if (! done) {
	       chem_mod_plane cmpl(plane_id, function);
	       cmpl.add_atom(atom_name, dist);
	       plane_mods.push_back(cmpl);
	    }
	 }
	 friend std::ostream& operator<<(std::ostream &s, chem_mod mod);
      };
   std::ostream& operator<<(std::ostream &s, chem_mod mod);

      class min_chem_mod {

      public:
	 min_chem_mod() {};
	 std::vector<chem_mod_atom>  atom_mods;
// 	 std::vector<chem_mod_tree>  tree_mods;
// 	 std::vector<chem_mod_bond>  bond_mods;
// 	 std::vector<chem_mod_angle> angle_mods;
// 	 std::vector<chem_mod_tor>   tor_mods;
// 	 std::vector<chem_mod_plane> plane_mods;
 	 std::vector<chem_mod_chir>  chir_mods;
	 void add_mod_atom(const chem_mod_atom &chem_atom) {
	    atom_mods.push_back(chem_atom);
	 }
#if 0
	 void add_mod_tree(const chem_mod_tree &chem_tree) {
	    // tree_mods.push_back(chem_tree);
	 }
	 void add_mod_bond(const chem_mod_bond &chem_bond) {
	    // bond_mods.push_back(chem_bond);
	 }
	 void add_mod_angle(const chem_mod_angle &chem_angle) {
	    // angle_mods.push_back(chem_angle);
	 }
	 void add_mod_tor(const chem_mod_tor &chem_tor) {
	    // tor_mods.push_back(chem_tor);
	 }
	 void add_mod_plane(const chem_mod_plane &chem_plane) {
	    // plane_mods.push_back(chem_plane);
	 }
	 void add_mod_chir(const chem_mod_chir &chem_chir) {
	    chir_mods.push_back(chem_chir);
	 }
	 void add_plane_atom(const std::string &plane_id,
			     const std::string &function,
			     const std::string &atom_name,
			     double dist) {
	 }
	 friend std::ostream& operator<<(std::ostream &s, min_chem_mod mod);
#endif
      };
   std::ostream& operator<<(std::ostream &s, min_chem_mod mod);

}
