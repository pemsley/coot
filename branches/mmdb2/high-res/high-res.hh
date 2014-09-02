
#include "mini-mol/mini-mol.hh"

namespace coot { 

   class high_res { 
      coot::minimol::molecule globular_molecule;
      std::pair<clipper::Coord_orth, CMMDBManager *> get_middle_pos(const coot::minimol::molecule &mol) const;
      void fill_globular_protein(const coot::minimol::molecule &mol, 
				 const clipper::Coord_orth &target,
				 CMMDBManager *mmdb_mol);
      void fill_globular_protein_by_fragments(const coot::minimol::molecule &mol, 
					      const clipper::Coord_orth &target,
					      CMMDBManager *mmdb_mol);
      void make_trees();
      void mark_neighbours(int iatom, int igroup,
			   const std::string &atom_name,
			   const std::vector<std::vector<int> > &neighbours,
			   PPCAtom atom_selection, int uddhandle);
      coot::minimol::molecule
      filter_on_groups(const std::vector<std::vector<int> > &groups,
		       CMMDBManager *mol,
		       PPCAtom atom_selection,
		       int n_selected_atom) const;

      static bool fragment_sorter(const coot::minimol::fragment &a,
				  const coot::minimol::fragment &b);

   public:
      high_res(const coot::minimol::molecule &m);
      // use fill_globular_protein_by_fragments:
      high_res(const coot::minimol::molecule &m, int iflag); 
      high_res(int i);
      // This one is for Kevin, who wanted to pass the program a
      // centre, not let it calculate it itself.
      high_res(const coot::minimol::molecule &m,
	       const clipper::Coord_orth &given_centre);
      void buccafilter_neighbours(); // slim down the multiple copies in space
			  // and leave only one averaged trace
      void buccafilter(); // slim down the multiple copies in space
			   // and leave only first fragment of
			   // overlapping fragments (fragments sorted
			   // by length first).
      void add_os();
      void add_cbetas();
      void output_pdb(const std::string &s) const;
   }; 

}
