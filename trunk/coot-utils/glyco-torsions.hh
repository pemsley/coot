
#include "utils/coot-utils.hh"
#include "coot-coord-utils.hh"

namespace coot { 

   typedef std::pair<bool, std::string> BS;

   class atom_by_torsion_base_t {

   public:
      std::string atom_name;
      std::string element;
      std::pair<bool, std::string> prior_atom_1;
      std::pair<bool, std::string> prior_atom_2;
      std::pair<bool, std::string> prior_atom_3;
   
      atom_by_torsion_base_t(const std::string &new_atom_name,
			     const std::string &ele_in,
			     const std::pair<bool, std::string> &prior_atom_1_in,
			     const std::pair<bool, std::string> &prior_atom_2_in,
			     const std::pair<bool, std::string> &prior_atom_3_in) {
	 atom_name = new_atom_name;
	 element = ele_in;
	 prior_atom_1 = prior_atom_1_in;
	 prior_atom_2 = prior_atom_2_in;
	 prior_atom_3 = prior_atom_3_in;
      }

      atom_by_torsion_base_t() {}
   }; 

   class atom_by_torsion_t : public atom_by_torsion_base_t {

      double bond_length; // to atom_1;
      double angle;       // to atom_2;
      double torsion;     // to atom_3;
      bool is_filled;
   
   public:
      atom_by_torsion_t() { is_filled = false;};
      atom_by_torsion_t(const atom_by_torsion_base_t &names,
			double bond_length_in, 
			double angle_in, 
			double torsion_in) // degrees
	 : atom_by_torsion_base_t(names) {
	 bond_length = bond_length_in;
	 angle = angle_in;
	 torsion = torsion_in;
	 is_filled = true;
      }
      atom_by_torsion_t(const atom_by_torsion_base_t &names,
			CResidue *residue_1_p,  // reference/lower
			CResidue *residue_2_p);
      
      bool filled() const { return is_filled; }
      clipper::Coord_orth pos(CResidue *base, CResidue *extending) const;
      friend std::ostream& operator<<(std::ostream &o, const atom_by_torsion_t &atb);
   };

   class link_by_torsion_base_t {
   public:
      std::vector<atom_by_torsion_base_t> atom_torsions;
      link_by_torsion_base_t() {}
      void add(const atom_by_torsion_base_t &at) {
	 atom_torsions.push_back(at);
      }
      bool filled() { return (atom_torsions.size()); }
   };

   class link_by_torsion_t : public link_by_torsion_base_t {
      void init(CResidue *ref_res_p, CResidue *ext_res_p);
   public:
      std::string new_residue_type;
      int new_res_no;
      std::vector<atom_by_torsion_t> geom_atom_torsions;
      link_by_torsion_t() {}
      link_by_torsion_t(const link_by_torsion_base_t &lbtb,
			CResidue *ref_res_p, CResidue *ext_res_p) :
	 link_by_torsion_base_t(lbtb) {
	 init(ref_res_p, ext_res_p);
      }

      link_by_torsion_t(const std::string &file_name) { read(file_name); }
   
      //    // old
      //    link_by_torsion_t(const std::string &n,
      // 		     int new_res_no_in,
      // 		     const link_by_torsion_base_t &lbtb,
      // 		     CResidue *ref_res_p, CResidue *ext_res_p)
      //       : link_by_torsion_base_t(lbtb) {
      //       new_residue_type = n;
      //       new_res_no = new_res_no_in;
      //       init(ref_res_p, ext_res_p);
      //    }
   
      // caller deletes
      CResidue *make_residue(CResidue *base_residue_p) const;
      void add(const atom_by_torsion_t &at) { geom_atom_torsions.push_back(at); }
      void print() const;
      void write(const std::string &file_name) const;
      void read(const std::string &file_name); // read what is written by write() function

      // helper function
      static std::pair<CResidue *, CResidue *> get_residue_pair(CMMDBManager *mol);

   };

   link_by_torsion_base_t get_names_for_link_type(const std::string &link_type);

   link_by_torsion_base_t pyranose_link_1_2_to_core();
   link_by_torsion_base_t pyranose_link_1_3_to_core();
   link_by_torsion_base_t pyranose_link_1_4_to_core();
   link_by_torsion_base_t pyranose_link_1_6_to_core();
   link_by_torsion_base_t pyranose_link_2_3_to_core();
   link_by_torsion_base_t asn_pyranose_link_to_core();
   
   link_by_torsion_base_t mannose_decorations(); // need GLC NAG
   link_by_torsion_base_t glucose_decorations(); 
   link_by_torsion_base_t NAG_decorations(); 
   

} // namespace coot
