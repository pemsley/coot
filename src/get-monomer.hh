
#include <string>


/*! \brief import libcheck monomer give the 3-letter code. 

@return the new molecule number, if not -1 (error). */
int get_monomer(const char *three_letter_code);

/* Use the protein geometry dictionary to retrieve a set of
   coordinates quickly from cif data read in from the RCSB's Chemical
   Component Library.  There are no restraints from this method
   though. */
int get_monomer_from_dictionary(const char *three_letter_code, int idealised_flag);


// return a new molecule number
int get_monomer_molecule_by_network_and_dict_gen(const std::string &text);
