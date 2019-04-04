
#include <string>


/*! \file
  \brief Coot Scripting Interface - Alignment utilities
*/

//!  \brief associate an alignment with a chain in a model molecule
//!
//! The pir_alignment is a string (with newlines)
//!
void associate_pir_alignment(int imol, std::string chain_id, std::string pir_alignment);

//!  \brief associate an alignment in a file with a chain in a model molecule
//!
void associate_pir_alignment_from_file(int imol, std::string chain_id, std::string pir_alignment_file_name);

//!  \brief apply the mutations of the associated alignment
//!
void apply_pir_alignment(int imol, std::string chain_id);

