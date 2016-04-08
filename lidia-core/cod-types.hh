
#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#ifndef COD_TYPES_HH
#define COD_TYPES_HH

#include <string>
#include "use-rdkit.hh"

namespace cod {

   class third_neighbour_info_t {
   public:
      RDKit::Atom *atom_p;
      std::string ele;
      unsigned int degree;
      third_neighbour_info_t() {
	 degree = 0;
      }
      third_neighbour_info_t(RDKit::Atom *atom_p_in,
			     const std::string &e,
			     unsigned int d) {
	 atom_p = atom_p_in;
	 ele = e;
	 degree = d;
      }
      bool operator==(const third_neighbour_info_t &t) const {
	 return (t.atom_p == atom_p);
      }
      bool operator<(const third_neighbour_info_t &t) const {
	 return (t.atom_p < atom_p);
      }
   };

   // can throw a std::runtime_error
   //
   // return a vector with as many items as there are atoms in rdkit_mol.
   // (throw an error if it can't do so).
   //
   // rdkit_mol is not const because there is no const beginAtoms() operator.
   std::vector<std::string> get_cod_atom_types(RDKit::ROMol &rdkit_mol,
					       bool add_name_as_property=false);

   // nb_level = 0 mean this atom neighbour info - the whole thing
   // 
   // can throw a std::runtime_error
   //
   std::pair<std::string, std::list<third_neighbour_info_t> >
   get_cod_atom_type(RDKit::Atom *atom_parent_p, // the parent of atom_p (if any)
		     RDKit::Atom *atom_p,
		     RDKit::Atom *atom_base_p,
		     const RDKit::ROMol &rdkit_mol,
		     int nb_level=0);
   
   cod::third_neighbour_info_t
   get_cod_nb_3_type(RDKit::Atom *atom_parent_p, // the parent of atom_p
				 RDKit::Atom *atom_p,
				 RDKit::Atom *atom_base_p,
				 const RDKit::ROMol &rdkit_mol);

   bool check_for_3rd_nb_info(RDKit::Atom *atom_parent_p,
			      RDKit::Atom *atom_p,
			      RDKit::Atom *atom_base_p,
			      const RDKit::ROMol &rdkm);

   std::string make_cod_type(RDKit::Atom *base_atom_p,
			     const std::string &atom_ele,
			     const std::vector<std::string> &neighbour_types,
			     const std::list<third_neighbour_info_t> &tniv,
			     int level);
   // which calls:
   std::string make_cod_3rd_neighb_info_type(const std::list<third_neighbour_info_t> &tniv);


   std::vector<std::string> sort_neighbours(const std::vector<std::string> &neighbours_in,
					    int level);
   bool neighbour_sorter(const std::string &a, const std::string &b);

   // we need info for neighbours (count) and aromatic status
   class ring_info_t {
      std::vector<int> atom_indices; // mirrors the atomRings() data (type) in RDKit
      bool aromaticity;
   public:
      ring_info_t(const std::vector<int> &ai) { atom_indices = ai; aromaticity = false; }
      ring_info_t(const std::vector<int> &ai, bool aromaticity_in)
      { atom_indices = ai; aromaticity = aromaticity_in; }
      void add_atom(unsigned int i) { atom_indices.push_back(i); }
      void set_aromaticity(bool arom) { aromaticity = arom; }
      unsigned int size() const { return atom_indices.size(); }
      bool get_aromaticity() const { return aromaticity; }
   };

   void handle_bigger_rings_from_fused_rings(RDKit::ROMol &rdkm,
					     const std::vector<std::vector<int> > &fused_rings);
   bool is_ring_member(unsigned int iat,   const std::vector<std::vector<int> > &fused_rings);

   std::vector<std::vector<int> > trace_path(unsigned int idx,
					     const std::map<int, std::vector<int> > &bond_map,
					     unsigned int n_max_bonds);
   std::vector<std::vector<int> > 
   trace_path(unsigned int idx,
	      std::vector<int> in_path_indices,
	      unsigned int target_idx,
	      const std::map<int, std::vector<int> > &bond_map,
	      unsigned int level);

   class bond_table_record_t {
   public:
      std::string cod_type_1;
      std::string cod_type_2;
      double mean;
      double std_dev;
      unsigned int count;
      bond_table_record_t(const std::string &cod_type_1_in,
			  const std::string &cod_type_2_in,
			  const double &mean_in,
			  const double &std_dev_in,
			  unsigned int count_in) {
	 cod_type_1 = cod_type_1_in;
	 cod_type_2 = cod_type_2_in;
	 mean = mean_in;
	 std_dev = std_dev_in;
	 count = count_in;
      }
      bool operator<(const bond_table_record_t &btr) const {

	 if (cod_type_1 == btr.cod_type_1) {
	    return (cod_type_2 < btr.cod_type_2);
	 } else {
	    return (cod_type_1 < btr.cod_type_1);
	 }
      }
      friend std::ostream &operator<<(std::ostream &s, const bond_table_record_t &btr);
   };
   std::ostream &operator<<(std::ostream &s, const bond_table_record_t &brc);

   class bond_record_container_t {
   public:
      bond_record_container_t() {}
      std::vector<bond_table_record_t> bonds;
      void add(const bond_table_record_t &rec) {
	 bonds.push_back(rec);
      }
      void add_table(const bond_record_container_t &brc) {
	 for (unsigned int i=0; i<brc.bonds.size(); i++) { 
	    bonds.push_back(brc.bonds[i]);
	 }
      }
      void sort() { std::sort(bonds.begin(), bonds.end()); }
      unsigned int size() { return bonds.size(); }
      bool write(const std::string file_name) const;
   };

   
   bond_record_container_t read_acedrg_table(const std::string &file_name);

   // return the consolidated table
   // 
   bond_record_container_t
   read_acedrg_table_dir(const std::string &dir_name);
   
} 

#endif // COD_TYPES_HH

#endif // MAKE_ENHANCED_LIGAND_TOOLS
