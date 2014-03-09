
namespace coot {

   class pi_stacking_instance_t {
   public:

      // CATION_PI_STACKING sets ligand_cationic_atom_name, not the
      // ligand_ring_atom_names vector.
      //
      enum stacking_t {
	 NO_STACKING,
	 PI_PI_STACKING,
	 PI_CATION_STACKING, // for cations on the protein residues (ligand pi)
	 CATION_PI_STACKING, // for cations on the ligand (protein TRY, PRO, TRP)
      };
      CResidue *res;
      stacking_t type; // pi-pi or pi-cation
      std::vector<std::string> ligand_ring_atom_names;
      float overlap_score; 
      std::string ligand_cationic_atom_name; // for cations on the ligand
      
      pi_stacking_instance_t(CResidue *res_in, stacking_t type_in,
			     const std::vector<std::string> &ring_atoms) {
	 res = res_in;
	 type = type_in;
	 ligand_ring_atom_names = ring_atoms;
      }
      
      // and the constructor for CATION_PI_STACKING
      // 
      pi_stacking_instance_t(CResidue *residue_in,
			     const std::string &ligand_atom_name_in) {
	 type = CATION_PI_STACKING;
	 res = residue_in;
	 ligand_cationic_atom_name = ligand_atom_name_in;
	 overlap_score = 0;
      }
      friend std::ostream& operator<< (std::ostream& s, const pi_stacking_instance_t &spec);
   };
   std::ostream& operator<< (std::ostream& s, const pi_stacking_instance_t &spec);

   class pi_stacking_container_t {
   private: 
      // can throw an exception
      std::pair<float, pi_stacking_instance_t::stacking_t>
      get_pi_overlap_to_ligand_ring(CResidue *res, const clipper::Coord_orth &pt) const;

      float get_pi_overlap_to_ligand_cation(CResidue *res, const clipper::Coord_orth &pt) const;
      
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
				  const clipper::Coord_orth &pt2,
				  const double &m1_pt_1, const double &m2_pt_1,
				  const double &m1_pt_2, const double &m2_pt_2) const;

      float
      overlap_of_cation_pi(const clipper::Coord_orth &ligand_pi_point,
			   const clipper::Coord_orth &cation_atom_point) const;

      std::vector<clipper::Coord_orth> get_cation_atom_positions(CResidue *res) const;
      // by search through res_ref
      std::vector<std::pair<std::string, clipper::Coord_orth> >
	 get_ligand_cations(CResidue *res, const coot::dictionary_residue_restraints_t &monomer_restraints) const;

      std::vector<std::vector<std::string> >
      get_aromatic_ring_list(const dictionary_residue_restraints_t &monomer_restraints) const;

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
      std::vector<std::vector<std::string> > get_aromatic_ring_list(const dictionary_residue_restraints_t &monomer_restraints,
								    const RDKit::ROMol &mol) const;
      std::vector<std::vector<std::string> > get_aromatic_ring_list(const RDKit::ROMol &mol) const;
#endif // MAKE_ENHANCED_LIGAND_TOOLS

      void init(const dictionary_residue_restraints_t &monomer_restraints,
		const std::vector<CResidue *> &filtered_residues,
		CResidue *res_ref,
		const std::vector<std::vector<std::string> > &aromatic_ring_list_atom_names);

      
   public:
      // a vector of residues and types
      std::vector<pi_stacking_instance_t> stackings;
      pi_stacking_container_t (const dictionary_residue_restraints_t &monomer_restraints,
			       const std::vector<CResidue *> &filtered_residues,
			       CResidue *res_ref);

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
      pi_stacking_container_t (const dictionary_residue_restraints_t &monomer_restraints,
			       const std::vector<CResidue *> &filtered_residues,
			       CResidue *res_ref,
			       const RDKit::ROMol &mol);
#endif // MAKE_ENHANCED_LIGAND_TOOLS      
   };

}
