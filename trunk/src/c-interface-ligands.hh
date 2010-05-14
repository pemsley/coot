
#include <map>

#include "mmdb_manager.h"
#include "clipper/core/coords.h"
#include "coot-coord-utils.hh"
#include "protein-geometry.hh"


namespace coot { 

   class fle_residues_helper_t {
   public:
      bool is_set;
      clipper::Coord_orth centre;
      residue_spec_t spec;
      std::string residue_name;
      fle_residues_helper_t() { is_set = 0; }
      fle_residues_helper_t(const clipper::Coord_orth &pt,
			    const residue_spec_t &spec_in,
			    const std::string &res_name_in) {
	 centre = pt;
	 spec = spec_in;
	 residue_name = res_name_in;
	 is_set = 1;
      }
   };

   class pi_stacking_instance_t {
   public:
      enum { PI_PI_STACKING, PI_CATION_STACKING };
      CResidue *res;
      int type; // pi-pi or pi-cation
      std::vector<std::string> ligand_ring_atom_names;
      pi_stacking_instance_t(CResidue *res_in, int type_in,
			     const std::vector<std::string> &ring_atoms) {
	 res = res_in;
	 type = type_in;
	 ligand_ring_atom_names = ring_atoms;
      }
   };

   class pi_stacking_container_t {
      float get_pi_overlap(CResidue *res, clipper::Coord_orth &pt) const;
      std::pair<clipper::Coord_orth, clipper::Coord_orth>
      get_ring_pi_centre_points(const std::vector<std::string> &ring_atom_names,
				CResidue *res_ref) const;
      
      // can throw an exception if not enough points found in pts.
      std::pair<clipper::Coord_orth, clipper::Coord_orth>
      ring_centre_and_normal(const std::vector<clipper::Coord_orth> &pts) const;

      // TRP has 2 rings, so we have to return a vector
      // 
      std::vector<std::vector<std::string> >
      ring_atom_names(const std::string &residue_name) const;
      
      float overlap_of_pi_spheres(const clipper::Coord_orth &pt1,
				  const clipper::Coord_orth &pt2) const;
      
   public:
      // a vector of residues and types
      std::vector<pi_stacking_instance_t> stackings;
      pi_stacking_container_t (const dictionary_residue_restraints_t &monomer_restraints,
			       const std::vector<CResidue *> &filtered_residues,
			       CResidue *res_ref);
   };


   void write_solvent_accessibilities(const std::vector<std::pair<coot::atom_spec_t, float> > &sav,
				      CResidue *reference_residue);

   // The bonds from the protein to the ligand which contain
   // ligand-atom-name residue-spec and bond type (acceptor/donor).
   // These (ligand atom names) will have to be mapped to x y position
   // of the flat ligand.
   // 
   class fle_ligand_bond_t {
   public:
      enum { H_BOND_DONOR_MAINCHAIN,
	     H_BOND_DONOR_SIDECHAIN,
	     H_BOND_ACCEPTOR_MAINCHAIN, 
	     H_BOND_ACCEPTOR_SIDECHAIN,
	     METAL_CONTACT_BOND,
	     BOND_OTHER };  // must sync this to lbg.hh
      std::string ligand_atom_name;
      int bond_type; // acceptor/donor
      residue_spec_t res_spec;
      double bond_length;
      fle_ligand_bond_t(const std::string &name_in,
			int bond_type_in,
			double bl_in,
			residue_spec_t spec_in) {
	 ligand_atom_name = name_in;
	 bond_type = bond_type_in;
	 res_spec = spec_in;
	 bond_length = bl_in;
      }
   };
   std::map<std::string, std::string> make_flat_ligand_name_map(CResidue *flat_res);

   std::vector<fle_ligand_bond_t> get_fle_ligand_bonds(CResidue *res_ref,
						       const std::vector<CResidue *> &residues,
						       const std::map<std::string, std::string> &name_map);
   void write_fle_centres(const std::vector<fle_residues_helper_t> &v,
			  const std::vector<coot::fle_ligand_bond_t> &bonds_to_ligand,
			  const std::vector<coot::solvent_exposure_difference_helper_t> &sed,
			  const pi_stacking_container_t &stacking,
			  CResidue *res_flat);

   bool is_a_metal(CResidue *res);


}
