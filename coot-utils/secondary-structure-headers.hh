
#include <mmdb2/mmdb_manager.h>

namespace coot {

   class access_model : public mmdb::Model {
   public:
      access_model(mmdb::Model model_in) : mmdb::Model(model_in) { }
      void add_sheets(mmdb::Sheets *sheets_in) {
	 sheets.nSheets = sheets_in->nSheets;
	 sheets.sheet   = sheets_in->sheet;
      }
      void add_helix(mmdb::Helix *helix) {
	 // helices is of type SSContainer
	 helices.AddData(helix);
      }
   };

   class secondary_structure_header_records {

   public:

      class helix_info_t {
      public:
	 mmdb::Residue *start_res;
	 mmdb::Residue *end_res;
	 unsigned int length;
	 helix_info_t(mmdb::Residue *r1, mmdb::Residue *r2, unsigned int l) : start_res(r1), end_res(r2), length(l) {}
      };


   private:

      // the order of all the sheets
      std::vector<std::vector<unsigned int> >
      get_sheet_order(mmdb::Manager *mol,
		      mmdb::Model *model_p,
		      const std::vector<std::vector<mmdb::Residue *> > &strands_with_residues);

      void make_sheets(mmdb::Manager *mol,
		       mmdb::Model *model_p,
		       const std::vector<std::vector<mmdb::Residue *> > &strands_with_residues);

      void make_helices(mmdb::Manager *mol,
			mmdb::Model *model_p,
			const std::vector<helix_info_t> &helices);

      void add_secondary_structure_header_records(mmdb::Model *model_p,
	       const std::vector<std::vector<mmdb::Residue *> > &helices_with_residues,
						  const std::vector<std::vector<mmdb::Residue *> > &strands_with_residues);

      std::string sheet_index_to_sheet_id(unsigned int idx) {
	 return std::string("A");
      }
      std::string helix_index_to_helix_id(unsigned int idx) {
	 // C++-11:
	 // return std::to_string(idx);
	 // old style C++
	 char buff[20];
	 std::sprintf(buff, "%d", idx+1);
	 return std::string(buff);
      }

   public:
      secondary_structure_header_records(mmdb::Manager *mol, bool needsCalcSecStructure);
   };
}
