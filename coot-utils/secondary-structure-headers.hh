
#ifndef SECONDARY_STRUCTURE_HEADERS_HH
#define SECONDARY_STRUCTURE_HEADERS_HH

#include <cstdio>
#include <string.h>

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
	 // helices protected and is of type SSContainer
	 if (helix) {
	    helices.AddData(helix);
	 }
      }
      ~access_model() {}
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

      class strand_relation_t {
      public:
	 unsigned int strand_idx;
	 enum sense_t { FIRST, PARALLEL, ANTI_PARALLEL, NO_RESULT };
	 sense_t sense;
	 strand_relation_t(unsigned int idx, sense_t s) : strand_idx(idx), sense(s) {}
	 bool operator==(const strand_relation_t &sr_in) const {
	    return (sr_in.strand_idx == strand_idx);
	 }
	 bool operator<(const strand_relation_t &sr_in) const {
	    return (sr_in.strand_idx < strand_idx);
	 }
	 static sense_t get_strand_sense(const std::vector<mmdb::Residue *> &strand_1,
					 const std::vector<mmdb::Residue *> &strand_2);
	 static int sense_to_pdb_sense(const sense_t &s) {
	    if (s == FIRST)    return 0;
	    if (s == PARALLEL) return 1;
	    if (s == ANTI_PARALLEL) return -1;
	    return -2; // Haha.
	 }
      };


   private:

      // the order of all the sheets
      std::vector<std::vector<strand_relation_t> >
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
         char c = 'A';
         c += idx;
	 std::string r(1,c);
         return r;
      }
      std::string helix_index_to_helix_id(unsigned int idx) {
	 // C++-11:
	 // return std::to_string(idx);
	 // old style C++
	 char buff[20];
	 std::snprintf(buff, 5, "%d", idx+1);
	 return std::string(buff);
      }

   public:
      secondary_structure_header_records() {};
      secondary_structure_header_records(mmdb::Manager *mol, bool needsCalcSecStructure);
   };
}
#endif // SECONDARY_STRUCTURE_HEADERS_HH
