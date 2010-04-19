
#include <map>

#include "mmdb_manager.h"
#include "clipper/core/coords.h"
#include "coot-coord-utils.hh"


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
			  CResidue *res_flat);


}
