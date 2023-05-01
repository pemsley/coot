

#ifndef STACK_AND_PAIR_HH
#define STACK_AND_PAIR_HH

#include <mmdb2/mmdb_manager.h>
#include <vector>
#include <string>

#include <clipper/core/coords.h>

#include "geometry/residue-and-atom-specs.hh"
#include "geometry/protein-geometry.hh"


namespace coot {

   // it is sensible to stack chain by chain, but pairing needs to be checked across chains (also).

   class stack_and_pair {
      double angle_crit;
      std::map<mmdb::Residue *, clipper::Coord_orth> normal_map;
      std::pair<bool, clipper::Coord_orth> get_base_normal(mmdb::Residue *residue_p) const;
      int mark_donors_and_acceptors(mmdb::Manager *mol, int selection_handle, const protein_geometry &geom);
      std::pair<bool,clipper::Coord_orth> get_base_centre(mmdb::Residue *residue_this) const;
      std::set<std::string> base_atom_name_set;
      void init();
      std::vector<std::string> get_base_atom_names(mmdb::Residue *residue_p) const;
   public:
      stack_and_pair() {}
      stack_and_pair(mmdb::Manager *mol, const std::vector<std::pair<bool,mmdb::Residue *> > &residues_vec);
      stack_and_pair(mmdb::Manager *mol, int selection_handle);

      // suitable to construct a parallel plane restraint (with no alt-confs at the moment though)
      class stacked_planes_info_t {
      public:
	 stacked_planes_info_t(mmdb::Residue *r1, mmdb::Residue *r2,
			       const std::vector<std::string> &an1,
			       const std::vector<std::string> &an2) : res_1(r1), res_2(r2),
								      atom_names_1(an1),
								      atom_names_2(an2) {}
	 mmdb::Residue *res_1;
	 mmdb::Residue *res_2;
	 std::vector<std::string> atom_names_1;
	 std::vector<std::string> atom_names_2;
      };

      // base pairing
      class paired_residues_info_t {
      public:
	 paired_residues_info_t(mmdb::Residue *r1, mmdb::Residue *r2,
				const std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > &atom_pair_vec_in) :
	    res_1(r1), res_2(r2), atom_pair_vec(atom_pair_vec_in) {}
	 mmdb::Residue *res_1;
	 mmdb::Residue *res_2;
	 std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > atom_pair_vec;
      };

      // don't use the selection
      std::vector<stacked_planes_info_t> stacked_residues(mmdb::Manager *mol);

      // Watson & Crick, Wobble, Reverse Wobble?
      // If this function knows that all of the residues of mol are moving, then
      // it can do the selection more quickly
      std::vector<paired_residues_info_t>
      paired_residues(mmdb::Manager *mol,
		      const std::vector<std::pair<bool, mmdb::Residue *> > &residues_vec,
		      bool residues_are_all_moving_flag,
		      const protein_geometry &geom);

      std::vector<std::pair<residue_spec_t, residue_spec_t> >
      paired_residue_specs(mmdb::Manager *mol,
			   const std::vector<std::pair<bool, mmdb::Residue *> > &residues_vec);

      bool contains_nucleic_acid(mmdb::Atom **SelAtom, int nselatom);
   
      std::map<mmdb::Residue *, clipper::Coord_orth> calculate_residue_normals(mmdb::Atom **SelAtom,
									       int n_sel_atoms);

      std::map<mmdb::Residue *, clipper::Coord_orth> calculate_residue_normals(const std::vector<std::pair<bool,mmdb::Residue *> > &residues_vec) const;
      

      bool similar_normals(mmdb::Residue *res_1, mmdb::Residue *res_2,
			   const std::map<mmdb::Residue *, clipper::Coord_orth> &normal_map) const;

   };



}

#endif
