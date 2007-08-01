
#include <iostream>
#include "bfkurt.hh"

CMMDBManager *get_atom_selection(std::string pdb_name); 


int main(int argc, char **argv) { 

   if (argc > 1) { 
      std::string pdb_file_name = argv[1];
      CMMDBManager *mol = get_atom_selection(pdb_file_name);
      if (! mol) { 
	 std::cout << "Failed to read pdb file\n";
      } else {
	 coot_extras::b_factor_analysis bfa(mol, 0);
	 short int only_questionables = 1;
	 bfa.write_table("bfactan.xml", pdb_file_name, only_questionables);
      }
   } else { 
      std::cout << "Usage: " << argv[0] << " pdb-file-name\n";
   }
   return 0;
} 


CMMDBManager *
get_atom_selection(std::string pdb_name) {

   int err;
   CMMDBManager* MMDBManager;

   // Needed for the error message printing: 
   // MMDBManager->GetInputBuffer(S, lcount);
   // Used by reference and as a pointer.  Grimness indeed.
   int  error_count;
   char error_buf[500];

   //   Make routine initializations
   //
   InitMatType();

   MMDBManager = new CMMDBManager;
   
   MMDBManager->SetFlag ( MMDBF_IgnoreBlankLines |
			  MMDBF_IgnoreDuplSeqNum |
			  MMDBF_IgnoreNonCoorPDBErrors );
   
   std::cout << "Reading coordinate file: " << pdb_name.c_str() << "\n";
   err = MMDBManager->ReadCoorFile((char *)pdb_name.c_str());
   
   if (err) {
      // does_file_exist(pdb_name.c_str());
      std::cout << "There was an error reading " << pdb_name.c_str() << ". \n";
      std::cout << "ERROR " << err << " READ: "
		<< GetErrorDescription(err) << std::endl;
      //
      // This makes my stomach churn too. Sorry.
      // 
      MMDBManager->GetInputBuffer(error_buf, error_count);
      if (error_count >= 0) { 
	 std::cout << "         LINE #" << error_count << "\n     "
		   << error_buf << std::endl << std::endl;
      } else {
	 if (error_count == -1) { 
	    std::cout << "       CIF ITEM: " << error_buf
		      << std::endl << std::endl;
	 }
      }
      //
   } else {
      // we read the coordinate file OK.
      //
      switch (MMDBManager->GetFileType())  {
      case MMDB_FILE_PDB    :  std::cout << " PDB"         ;
	 break;
      case MMDB_FILE_CIF    :  std::cout << " mmCIF"       ; 
	 break;
      case MMDB_FILE_Binary :  std::cout << " MMDB binary" ;
	 break;
      default:
	 std::cout << " Unknown (report as a bug!)\n";
      }
      std::cout << " file " << pdb_name.c_str() << " has been read.\n";

      // atom_selection_container.read_error_message = NULL; // its a string
   }
   
    std::cout << "Spacegroup: " << MMDBManager->GetSpaceGroup() << "\n";
    return MMDBManager;
}
