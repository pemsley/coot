

#ifndef MOGUL_INTERFACE_HH
#define MOGUL_INTERFACE_HH

#include <string>
#include <vector>

#include "protein-geometry.hh"

namespace coot {

   // this can throw a std::runtime_error
   class mogul_distribution {
   public:
      mogul_distribution(const std::vector<std::string> &bits);
      mogul_distribution() {}
      float bin_start;
      float bin_end;
      float bin_width;
      int n_bins;
      std::vector<int> counts;
   };

   class mogul_item {
   public:
      enum type_t { BOND, ANGLE, TORSION, PLANE };
      // 1-based counting
      int idx_1;
      int idx_2;
      int idx_3;
      int idx_4;
      std::vector<int> plane_indices;
      int counts;
      float value, mean, median, std_dev, z;
      type_t type;
      float max_badness;
      mogul_distribution distribution;
   public:
      mogul_item() { max_badness = 5.0; }
      mogul_item(int atom_idx_1_in, int atom_idx_2_in,
		 float value_in,
		 int counts_in, float mean_in, float median_in, float std_dev_in, float z_in) {
	 type = BOND;
	 value = value_in;
	 idx_1 = atom_idx_1_in;
	 idx_2 = atom_idx_2_in;
	 counts = counts_in;
	 mean = mean_in;
	 median = median_in;
	 std_dev = std_dev_in;
	 z = z_in;
	 max_badness = 5.0; 
      }
      mogul_item(int atom_idx_1_in, int atom_idx_2_in, int atom_idx_3_in,
		 float value_in,
		 int counts_in, float mean_in, float median_in, float std_dev_in, float z_in) {
	 type = ANGLE;
	 value = value_in;
	 idx_1 = atom_idx_1_in;
	 idx_2 = atom_idx_2_in;
	 idx_3 = atom_idx_3_in;
	 counts = counts_in;
	 mean = mean_in;
	 median = median_in;
	 std_dev = std_dev_in;
	 z = z_in;
	 max_badness = 5.0; 
      }
      mogul_item(int atom_idx_1_in, int atom_idx_2_in, int atom_idx_3_in, int atom_idx_4_in,
		 float value_in,
		 int counts_in, float mean_in, float median_in, float std_dev_in, float z_in) {
	 type = TORSION;
	 value = value_in;
	 idx_1 = atom_idx_1_in;
	 idx_2 = atom_idx_2_in;
	 idx_3 = atom_idx_3_in;
	 idx_4 = atom_idx_4_in;
	 counts = counts_in;
	 mean = mean_in;
	 median = median_in;
	 std_dev = std_dev_in;
	 z = z_in;
	 max_badness = 5.0; 
      }
      void set_max_z_badness(float b) { max_badness = b; }
      std::string colour() const; // uses max_badness to return a hex colour string.
      void add_distribution(const mogul_distribution &d) { distribution = d; }
      int max_counts_in_a_bin() const {
	 int r = 0;
	 for (unsigned int i=0; i<distribution.counts.size(); i++) { 
	    if (distribution.counts[i] > r)
	       r = distribution.counts[i];
	 }
	 return r;
      }
   };

   class mogul {
      std::vector<mogul_item> items;
      // can throw a runtime_exception
      mogul_item parse_bond_line(const std::vector<std::string> &bits,
				 const std::vector<std::string> &stats_bits,
				 const std::vector<std::string> &distrubtion_bits) const;
      mogul_item parse_angle_line(const std::vector<std::string> &bits,
				  const std::vector<std::string> &stats_bits,
				  const std::vector<std::string> &distribution_bits) const;
      float max_z_badness; // passed to mogul items.
      
      std::string get_bond_type(coot::dictionary_residue_restraints_t &restraints,
				const std::string &name_1,
				const std::string &name_2) const;
   public:
      mogul() { max_z_badness = 5.0; }
      mogul(const std::string &file_name) {
	 max_z_badness = 5.0;
	 parse(file_name);
      }
      void parse(const std::string &file_name);
      unsigned int n_items() const { return items.size(); }
      const mogul_item &operator[](unsigned int item_no) const {
	 return items[item_no];
      }
      void set_max_z_badness(float b) { max_z_badness = b; }
      // the following coot:: is for swig
      coot::dictionary_residue_restraints_t make_restraints(CResidue *residue_p,
							    const std::string &comp_id,
							    const coot::protein_geometry &geom);
   };

}

#endif // MOGUL_INTERFACE_HH

