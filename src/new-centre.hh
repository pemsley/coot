
namespace coot { 

   typedef enum { NO_LIGANDS=0, NORMAL_CASE=1, SINGLE_LIGAND_NO_MOVEMENT } new_ligand_position_type;

   class new_centre_info_t {
   public:
      new_ligand_position_type type;
      clipper::Coord_orth position;
      std::string info_string;
      residue_spec_t residue_spec;
      new_centre_info_t(new_ligand_position_type pt,
			               const clipper::Coord_orth &pos,
			               const residue_spec_t &res_spec_in) {
         type = pt;
         position = pos;
         residue_spec = res_spec_in;
      } 
   };

}
