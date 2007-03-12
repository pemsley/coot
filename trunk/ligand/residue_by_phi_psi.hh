
#include "ligand.hh"
#include "clipper/core/ramachandran.h"


namespace coot { 

   class phi_psi_pair { 
   public:
      phi_psi_pair(float a, float b) { 
	 phi = a;
	 psi = b;
      }
      float phi;
      float psi;
   }; 

  class residue_by_phi_psi : public ligand { 

     const CMMDBManager *mol;
     int ires_terminus;
     std::string chain_id;
     std::string residue_type;
     std::string terminus_type;
     float rama_max;
     clipper::Ramachandran rama;


     // This is not const because GetAtom() fo CResidue is not a const
     // function.  Argh!
     // 
     CResidue *residue_p; // the residue of the last atom (we
                          // clicked on an atom of it).

     phi_psi_pair get_phi_psi_by_random() const;
     void init_phi_psi_plot(); 

     minimol::residue 
     construct_next_res_from_rama_angles(float phi, float psi, int seqno,
				const clipper::Coord_orth &previous_n,
				const clipper::Coord_orth &previous_ca,
				const clipper::Coord_orth &previous_c) const; 
     minimol::residue 
     construct_prev_res_from_rama_angles(float phi, float psi, int seqno,
				const clipper::Coord_orth &next_n,
				const clipper::Coord_orth &next_ca,
				const clipper::Coord_orth &next_c) const; 

     
     minimol::residue construct_joining_res(const phi_psi_pair &pp,
						  int seqno,
						  const clipper::Coord_orth &next_n,
						  const clipper::Coord_orth &next_ca,
						  const clipper::Coord_orth &next_c) const;
	
     std::vector<clipper::Coord_orth> get_connecting_residue_atoms() const; 
     minimol::fragment fit_terminal_residue_generic(int n_trials, int offset, 
							  short int do_rigid_body_refinement);

     minimol::fragment
     make_2_res_joining_frag(const std::string &chain_id,
			     const phi_psi_pair &pp1,
			     const phi_psi_pair &pp2,
			     int seqnum,
			     int offset, // + or - 1
			     const clipper::Coord_orth &next_n,
			     const clipper::Coord_orth &next_ca,
			     const clipper::Coord_orth &next_c) const;

     void debug_compare_check(const coot::minimol::residue &mres,
			      std::vector<minimol::atom *> atoms_p);

     
  public:
     residue_by_phi_psi(const CMMDBManager *mol_in,
			const std::string &terminus_type, // "N", or "C"
			const CResidue *res_p,
			const std::string &chain_id, 
			const std::string &res_type);

     minimol::molecule best_fit_phi_psi(int n_trials, short int do_rigid_body_refinement,
					int add_other_residue_flag);

  };

  minimol::residue 
  build_N_terminal_ALA(float phi, float psi, int seqno,
			const clipper::Coord_orth &previous_n,
			const clipper::Coord_orth &previous_ca,
			const clipper::Coord_orth &previous_c); 
   minimol::residue 
   build_C_terminal_ALA(float phi, float psi, int seqno,
			const clipper::Coord_orth &next_n,
			const clipper::Coord_orth &next_ca,
			const clipper::Coord_orth &next_c); 


} // namespace coot
