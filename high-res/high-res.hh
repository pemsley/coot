
#include "mini-mol/mini-mol.hh"

namespace coot { 

   class high_res {
      minimol::molecule globular_molecule;
      std::pair<clipper::Coord_orth, mmdb::Manager *> get_middle_pos(const minimol::molecule &mol) const;
      void fill_globular_protein(const minimol::molecule &mol,
				 const clipper::Coord_orth &target,
				 mmdb::Manager *mmdb_mol);
      void fill_globular_protein_by_fragments(const minimol::molecule &mol,
					      const clipper::Coord_orth &target,
					      mmdb::Manager *mmdb_mol);
      void make_trees();
      void mark_neighbours(int iatom, int igroup,
			   const std::string &atom_name,
			   const std::vector<std::vector<int> > &neighbours,
			   mmdb::PPAtom atom_selection, int uddhandle);
      minimol::molecule
      filter_on_groups(const std::vector<std::vector<int> > &groups,
		       mmdb::Manager *mol,
		       mmdb::PPAtom atom_selection,
		       int n_selected_atom) const;

      static bool fragment_sorter(const minimol::fragment &a,
				  const minimol::fragment &b);

   public:
      high_res(const minimol::molecule &m);
      // use fill_globular_protein_by_fragments:
      high_res(const minimol::molecule &m, int iflag);
      high_res(int i);
      // This one is for Kevin, who wanted to pass the program a
      // centre, not let it calculate it itself.
      high_res(const minimol::molecule &m,
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
