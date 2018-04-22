
#include <utility>
#include <vector>
#include <clipper/core/coords.h>

#include "mini-mol/mini-mol.hh"

namespace coot {
   class ligand_score_card {
      int n_ligand_atoms; // non-H.
      int ligand_no;
      double atom_point_score; // the old ranking metric, now use get_score()
   public:
      // consider using a member function here:
      bool many_atoms_fit;
      double score_per_atom;
      std::pair<bool, double> correlation;
      std::vector<std::pair<clipper::Coord_orth, float> > scored_characteristic_low_density_points;
      
      ligand_score_card() {
	 ligand_no = -1; // unset
	 atom_point_score = 0.0;
	 many_atoms_fit = 0;
	 n_ligand_atoms = 0;
	 score_per_atom = 0.0;
	 correlation.first = false;
	 correlation.second = -1;
      }
      void set_ligand_number(int ilig) {
	 ligand_no = ilig;
      }
      void set_n_ligand_atoms(int n) {
	 n_ligand_atoms = n;
      }
      bool operator<(const ligand_score_card &other) const {
	 return (other.atom_point_score < atom_point_score);
      }
      void add(const double &d) {
	 atom_point_score += d;
      }
      double get_score() const;
      friend std::ostream& operator<<(std::ostream &s, const ligand_score_card &lsc);
   };
   std::ostream& operator<<(std::ostream &s, const ligand_score_card &lsc);

   class scored_ligand_eraser {
   public:
      float max_correl;
      scored_ligand_eraser(float max_correl_in) {
	 max_correl = max_correl_in;
      }
      // return shall-we-delete? status
      bool operator() (const std::pair<minimol::molecule, ligand_score_card> &sl) {
	 if (! sl.second.correlation.first) 
	    return true;
	 if (sl.second.correlation.second < max_correl)
	    return true;
	 return false;
      }
   };
}
