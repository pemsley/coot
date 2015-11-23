
namespace coot { 

   class animated_ligand_interactions_t : public fle_ligand_bond_t { 
   public:
      animated_ligand_interactions_t(const fle_ligand_bond_t &lb) :
	 fle_ligand_bond_t(lb) { }
      void draw(mmdb::Manager *mol,
		const gl_context_info_t &gl_info,
		const long &start_time) const;
   };

   // encapsulate this into molecule_class_info_t?
   //
   // trivial helper class to get specs and distance for atoms when a
   // link is made.
   // 
   class dict_link_info_t {
      bool check_for_order_switch(mmdb::Residue *residue_ref,
				  mmdb::Residue *residue_new,
				  const std::string &link_type,
				  const protein_geometry &geom) const;
   public:
      // this can throw a std::runtime_error
      dict_link_info_t (mmdb::Residue *residue_ref, mmdb::Residue *residue_new,
			const std::string &link_type, const protein_geometry &geom);
      atom_spec_t spec_ref;
      atom_spec_t spec_new;
      double dist;
   };
}
