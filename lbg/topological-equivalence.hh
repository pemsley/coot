
// -----------------------------------------------------------------
//                   chirality
// -----------------------------------------------------------------

class topological_equivalence_t {

   // internal copy of input params.
   // 
   std::vector<widgeted_atom_t> atoms;
   std::vector<widgeted_bond_t> bonds;
   
   std::vector<bool> unique; // is the atom index marked as unique? initially all 0.
   std::vector<int> isn;     // invariant-sequence-numbers
   std::map<std::string, std::vector<int> > atom_map;

   // the number of different EC values in the molecules.
   int n_extended_connectivity(const std::vector<long int> &equivalent_classes) const;

   bool continue_ec_calculations_p(const std::vector<long int> &curr_eqv, 
				   const std::vector<long int> &prev_eqv);
   
   // fiddles with unique, return true if at least one unique was assigned.
   bool assign_uniques(const std::vector<long int> &extended_connectivity);

   void assign_invariant_sequence_number(const std::vector<long int> &curr_ec);

   // return the next isn index to be used.
   int assign_invariant_sequence_number(const std::vector<long int> &curr_ec,
					const std::vector<std::pair<int, int> > &atom_index,
					int next_index);

   bool atoms_have_unassigned_isn_p() const;

   // return a flag to let us know that it was done.
   bool mark_isn(int atom_index, int i_s_n); 

   // old
   // fiddles with unique
   bool identified_unique_p(const std::vector<long int> &curr_eqv, 
			    const std::vector<long int> &prev_eqv);

   std::vector<long int> assign_initial_topo_indices();
   std::vector<long int> assign_topo_indices(const std::vector<long int> &prev_eqv,
					     int round);

   // Return a list of atom indices that are connected to 3 or 4 other
   // atoms (return the indices of those other atoms too.
   // 
   std::vector<std::pair<int, std::vector<int> > > tetrahedral_atoms() const;

      
public:
   topological_equivalence_t(const std::vector<widgeted_atom_t> &atoms,
			     const std::vector<widgeted_bond_t> &bonds);

   std::vector<std::string> chiral_centres() const;

};

