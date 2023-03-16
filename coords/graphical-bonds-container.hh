
#ifndef GRAPHICAL_BONDS_CONTAINER_HH
#define GRAPHICAL_BONDS_CONTAINER_HH

#include "graphics-line.hh"
#include "torus-description.hh"
#include "rotamer-container.hh"
#include "coot-utils/coot-rama.hh"
#include "ramachandran-container.hh"
#include "coot-utils/cis-peptide-info.hh"

// A poor man's vector. Contains the set of lines for each element (well,
// basically colour) type.
// 
template<class T> class graphical_bonds_lines_list {

 public:
   int num_lines;
   T *pair_list;
   bool thin_lines_flag;

   graphical_bonds_lines_list() {
      num_lines = 0;
      pair_list = NULL;
      thin_lines_flag = 0;
   }
};


class graphical_bonds_atom_info_t {
public:
   bool is_hydrogen_atom;
   bool is_water; // don't not display this in sticks-only mode - or the water
                  // will disappear (needs rebonding otherwise - not just a display flag)
   float radius_scale; // Waters (and perhaps metals) should have big radii, so that
                       // they are easier to see.
   coot::Cartesian position;
   mmdb::Atom *atom_p; // this should be a shared pointer I think.
                       // we don't want to be looking at this pointer
                       // if some other part of the code has deleted the atom.
   int model_number; // -1 is unset
   int atom_index;
   graphical_bonds_atom_info_t(const coot::Cartesian &pos, int atom_index_in, bool is_hydrogen_atom_in) :
      position(pos) {
      model_number = -1;
      is_hydrogen_atom = is_hydrogen_atom_in;
      is_water = false;
      atom_index = atom_index_in;
      atom_p = 0;
      radius_scale = 1.0;
   }
   graphical_bonds_atom_info_t() {
      model_number = -1;
      is_hydrogen_atom = false;
      atom_index = -1; // unset
      radius_scale = 1.0;
      atom_p = 0;
      is_water = false;
   }
   // this is a bit of a weird construction
   bool radius_for_atom_should_be_big(mmdb::Atom *atom_p) const {

      // 20190822-PE: you might like to add other tests here.
      // 20200608-PE: I did!
      mmdb::Residue *r = atom_p->GetResidue();
      if (r) {
         std::string res_name = r->GetResName();
         if (res_name == "HOH")
            return true;
         if (res_name == "CA")
            return true;
         if (res_name == "MG")
            return true;
         if (res_name == "IOD")
            return true;
         if (res_name == "CL")
            return true;
         if (res_name == "NA")
            return true;
         if (res_name == "K")
            return true;
      }
      return false;
   }

   float get_radius_scale_for_atom(mmdb::Atom *atom_p) const {

      float scale = 1.0f;
      mmdb::Residue *r = atom_p->GetResidue();
      if (r) {
         std::string ele(atom_p->element);
         if (ele == " H") return 0.5f;
         std::string res_name = r->GetResName();
         if (res_name == "HOH")
            return 2.6f;
         if (res_name == "CA")
            return 4.0f;
         if (res_name == "MG")
            return 4.0f;
         if (res_name == "IOD")
            return 4.0f;
         if (res_name == "CL")
            return 4.0f;
         if (res_name == "NA")
            return 4.0f;
         if (res_name == "K")
            return 4.0f;
      }
      return scale;
   }

   void set_radius_scale_for_atom(mmdb::Atom *at, bool make_fat_atom) {
      // 20230224-PE if there is no dictionary, then we want big fat atoms
      radius_scale = get_radius_scale_for_atom(at);
      if (make_fat_atom)
         radius_scale = 6.0;
      if (radius_scale > 6.0) radius_scale = 6.0;
   }

};

template<class T> class graphical_bonds_points_list {

public:
   unsigned int num_points;
   unsigned int current_count;

   // use a is-H-atom-flag for first
   T *points;
   
   graphical_bonds_points_list() {
      current_count = 0;
      num_points = 0;
      points = NULL;
   }

   graphical_bonds_points_list(unsigned int size) {
      current_count = 0;
      num_points = size;
      points = new T[size];
   }

   void add_point(const T &pt) {
      points[current_count] = pt;
      current_count++;
   } 
};

class graphical_bonds_cis_peptide_markup {

public:

   bool is_pre_pro_cis_peptide;
   bool is_twisted; // twisted trans
   int model_number; // -1 is unset
   coot::Cartesian pt_ca_1; 
   coot::Cartesian pt_c_1;
   coot::Cartesian pt_n_2;
   coot::Cartesian pt_ca_2;
   coot::atom_index_quad atom_index_quad;
   graphical_bonds_cis_peptide_markup(const coot::Cartesian &pt_ca_1_in,
				      const coot::Cartesian &pt_c_1_in,
				      const coot::Cartesian &pt_n_2_in,
				      const coot::Cartesian &pt_ca_2_in,
				      bool is_pre_pro_cis_peptide_in,
				      bool is_twisted_in,
				      int model_number_in) :
      pt_ca_1(pt_ca_1_in), pt_c_1(pt_c_1_in), pt_n_2(pt_n_2_in), pt_ca_2(pt_ca_2_in) {
      is_pre_pro_cis_peptide = is_pre_pro_cis_peptide_in;
      is_twisted = is_twisted_in;
      model_number = model_number_in;
   }

   void add_atom_index_quad(const coot::atom_index_quad &iq) {
      atom_index_quad = iq;
   }

   graphical_bonds_cis_peptide_markup() {
      model_number = -1;
      is_pre_pro_cis_peptide = false;
      is_twisted = false;
   } 
};

// Uses graphical_bonds_lines_list
// 
// a graphical_bonds_container is used by draw_molecule() and are
// created from a Bond_lines_container (which uses a vector).
// 
class graphical_bonds_container { 

 public:
   
   enum { NO_BOND,
	  BONDED_WITH_STANDARD_ATOM_BOND,
	  BONDED_WITH_BOND_TO_HYDROGEN,
	  BONDED_WITH_HETATM_BOND /* by dictionary */ };
   int num_colours;
   graphical_bonds_lines_list<graphics_line_t> *bonds_;

   int symmetry_has_been_created;
   graphical_bonds_lines_list<graphics_line_t> *symmetry_bonds_;

   // if the distance between CAs in a missing loop is longer than is possible given
   // the residue number difference, then we want to mark up that line with
   // dots along the line joining the residues.  This should work similarly with P-P
   // for nucleic acid chains - but I won't change the function name.
   coot::Cartesian *bad_CA_CA_dist_spots_ptr;
   coot::Cartesian *zero_occ_spots_ptr;
   coot::Cartesian *deuterium_spots_ptr;
   std::pair<coot::Cartesian, float> *ramachandran_goodness_spots_ptr;
   int n_zero_occ_spots;
   int n_bad_CA_CA_dist_spots;
   int n_deuterium_spots;
   int n_ramachandran_goodness_spots;

   // first is is-H-atom-flag
   graphical_bonds_atom_info_t *atom_centres_;
   int n_atom_centres_;
   int *atom_centres_colour_;
   std::vector<coot::torus_description_t> rings;
   int n_consolidated_atom_centres;
   graphical_bonds_points_list<graphical_bonds_atom_info_t> *consolidated_atom_centres;
   int n_cis_peptide_markups;
   graphical_bonds_cis_peptide_markup *cis_peptide_markups;
   int n_rotamer_markups;
   rotamer_markup_container_t *rotamer_markups;
   
   graphical_bonds_container() {
      num_colours = 0; 
      bonds_ = NULL;
      symmetry_has_been_created = 0; 
      symmetry_bonds_ = NULL; 
      zero_occ_spots_ptr = NULL;
      bad_CA_CA_dist_spots_ptr = NULL;
      n_bad_CA_CA_dist_spots = 0;
      n_zero_occ_spots = 0;
      deuterium_spots_ptr = NULL;
      n_deuterium_spots = 0;
      atom_centres_colour_ = NULL;
      atom_centres_ = NULL; 
      n_atom_centres_ = 0;
      n_ramachandran_goodness_spots = 0;
      ramachandran_goodness_spots_ptr = NULL;
      consolidated_atom_centres = NULL;
      n_consolidated_atom_centres = 0;
      n_cis_peptide_markups = 0;
      cis_peptide_markups = NULL;
      n_rotamer_markups = 0;
      rotamer_markups = NULL;
   }

   void clear_up() {

      if (bonds_)
	 for (int icol=0; icol<num_colours; icol++)
	    delete [] bonds_[icol].pair_list;
      if (symmetry_bonds_)
	 for (int icol=0; icol<num_colours; icol++)
	    delete [] symmetry_bonds_[icol].pair_list;

      delete [] bonds_;  // null testing part of delete
      delete [] symmetry_bonds_; 
      delete [] atom_centres_;
      delete [] atom_centres_colour_;
      bonds_ = NULL;
      symmetry_bonds_ = NULL;
      atom_centres_ = NULL;
      atom_centres_colour_ = NULL;
      if (n_zero_occ_spots) 
	 delete [] zero_occ_spots_ptr;
      if (n_deuterium_spots)
	 delete [] deuterium_spots_ptr;
      if (n_ramachandran_goodness_spots)
	 delete [] ramachandran_goodness_spots_ptr;
      n_zero_occ_spots = 0;
      n_deuterium_spots = 0;
      n_ramachandran_goodness_spots = 0;
      n_atom_centres_ = 0;
      if (consolidated_atom_centres) {
	 for (int i=0; i<n_consolidated_atom_centres; i++)
	    delete [] consolidated_atom_centres[i].points;
	 delete [] consolidated_atom_centres;
	 consolidated_atom_centres = NULL;
      }
      delete [] cis_peptide_markups;
      cis_peptide_markups = NULL;
      n_rotamer_markups = 0;
      delete [] rotamer_markups;
      rotamer_markups = NULL;
   }

   graphical_bonds_container(const std::vector<graphics_line_t> &a) { 

      std::cout << "constructing a graphical_bonds_container from a vector " 
		<< "of size " << a.size() << std::endl;

      num_colours = 1;
      
      bonds_ = new graphical_bonds_lines_list<graphics_line_t>[1]; // only 1 graphical_bonds_lines_list needed
      bonds_[0].pair_list = new graphics_line_t[(a.size())];
      bonds_[0].num_lines = a.size();

      // copy over
      for(int i=0; i<bonds_[0].num_lines; i++)
	 bonds_[0].pair_list[i] = a[i];

      symmetry_bonds_ = NULL; 
      symmetry_has_been_created = 0; 
      zero_occ_spots_ptr = NULL;
      n_zero_occ_spots = 0;
      deuterium_spots_ptr = NULL;
      n_deuterium_spots = 0;
      atom_centres_colour_ = NULL;
      atom_centres_ = NULL; 
      n_atom_centres_ = 0;
      n_ramachandran_goodness_spots = 0;
      ramachandran_goodness_spots_ptr = NULL;
      consolidated_atom_centres = NULL;
      n_consolidated_atom_centres = 0;
      n_cis_peptide_markups = 0;
      cis_peptide_markups = NULL;
   }
      
   void add_colour(const std::vector<graphics_line_t> &a);
   void add_zero_occ_spots(const std::vector<coot::Cartesian> &spots);
   void add_bad_CA_CA_dist_spots(const std::vector<coot::Cartesian> &spots);
   void add_deuterium_spots(const std::vector<coot::Cartesian> &spots);
   void add_ramachandran_goodness_spots(const std::vector<std::pair<coot::Cartesian,
					coot::util::phi_psi_t> > &spots,
					const ramachandrans_container_t &rc);
   void add_rotamer_goodness_markup(const std::vector<rotamer_markup_container_t> &ric);

   void add_atom_centres(const std::vector<graphical_bonds_atom_info_t> &centres,
			 const std::vector<int> &colours);
   bool have_rings() const { return rings.size(); }
   bool empty() const { return (bonds_ == NULL); }

   unsigned int n_bonds() const; // count them up
   unsigned int n_atoms() const;

   void add_cis_peptide_markup(const std::vector<coot::util::cis_peptide_quad_info_t> &cis_peptide_quads);
};

#endif // GRAPHICAL_BONDS_CONTAINER_HH
