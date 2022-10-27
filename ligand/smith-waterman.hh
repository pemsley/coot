
#ifndef SMITH_WATERMAN_HH
#define SMITH_WATERMAN_HH

#include <vector>
#include <string>
#include <map>

#include <mmdb2/mmdb_manager.h>
#include <clipper/core/xmap.h>
#include "geometry/protein-geometry.hh"

// delte this on merge
std::string coot_util_single_letter_to_3_letter_code(char code);

namespace sm_wat {

   // ---------------------------------------------------------------------
   //                 the testing/devel interface
   // ---------------------------------------------------------------------
   
   float H_i_j(int i, int j,
            const std::vector<std::vector<std::pair<bool, float> > > &H,
            const std::string &seq_1,
            const std::string &seq_2);

   void print_H(const std::vector<std::vector<std::pair<bool, float> > > &H);
   void fill_scoring_matrix(std::vector<std::vector<std::pair<bool, float> > > &H,
                            const std::string &seq_1,
                            const std::string &seq_2);
   float H_i_j(int i, int j,
               const std::vector<std::vector<std::pair<bool, float> > > &H,
               const std::string &seq_1,
               const std::string &seq_2);

   float score_with_method_1(int i, int j,
                             const std::vector<std::vector<std::pair<bool, float> > > &H,
                             const std::string &seq_1,
                             const std::string &seq_2);
   float score_with_method_2(int i, int j,
                             const std::vector<std::vector<std::pair<bool, float> > > &H,
                             const std::string &seq_1,
                             const std::string &seq_2);
   float score_with_method_3(int i, int j,
                             const std::vector<std::vector<std::pair<bool, float> > > &H,
                             const std::string &seq_1,
                             const std::string &seq_2);
   
   std::vector<std::vector<std::pair<bool, float> > > construct_H(const std::string &seq_1, 
                                                                  const std::string &seq_2);
   float score_against_type(char c, const std::string &type);
   float score_against_type(int idx, const std::string &seq, const std::string &type);

   std::pair<std::vector<int>, std::vector<int> > backtrack_seq(const std::vector<std::vector<std::pair<bool, float> > > &H);

   // ---------------------------------------------------------------------
   //                 the real interface
   // ---------------------------------------------------------------------

   class cell_t {
   public:
      int i;
      int j;
      cell_t() { i=-1; j=-1; }
      cell_t(const int &ii, const int &jj) : i(ii), j(jj) {}
      bool operator==(const cell_t &other) const {
         if (other.i == i)
            if (other.j == j)
               return true;
         return false;
      }
      bool non_zero_indexed() const {
         if (i > 0)
            if (j > 0)
               return true;
         return false;
      }
   };

   std::vector<cell_t> smith_waterman(const std::string &sequence,
                                      const std::vector<std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > > > &scored_residues);

   std::vector<std::vector<std::pair<bool, float> > > construct_H(const std::string &target_sequence);
   float s(char a, const std::map<std::string, std::pair<std::string, double> > &b);
   float W_gap_sequence(int k);
   float W_gap_residues(int k);

   float score_with_method_1(int seq_idx,
                             int types_idx,
                             const std::vector<std::vector<std::pair<bool, float> > > &H,
                             const std::string &target_sequence,
                             const std::vector<std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > > > &scored_residues);

   float score_with_method_2(int seq_idx,
                             int types_idx,
                             const std::vector<std::vector<std::pair<bool, float> > > &H,
                             const std::string &target_sequence,
                             const std::vector<std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > > > &scored_residues);

   float score_with_method_3(int seq_idx,
                             int types_idx,
                             const std::vector<std::vector<std::pair<bool, float> > > &H,
                             const std::string &target_sequence,
                             const std::vector<std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > > > &scored_residues);

   void fill_scoring_matrix(std::vector<std::vector<std::pair<bool, float> > > &H,
                            const std::string &target_sequence,
                            const std::vector<std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > > > &scored_residues);
   
   float H_i_j(int seq_idx, int types_idx, const std::vector<std::vector<std::pair<bool, float> > > &H,
               const std::string &target_sequence,
               const std::vector<std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > > > &scored_residues);

   std::pair<float, std::vector<sm_wat::cell_t> >
   backtrack(const std::vector<std::vector<std::pair<bool, float> > > &H);

   // backtrack non-top (the passed cell vector is the trace(s) from better scoring traces)
   std::pair<float, std::vector<sm_wat::cell_t> >
   backtrack_others(const std::vector<std::vector<std::pair<bool, float> > > &H,
                    const std::vector<cell_t> &cells_from_better_traces);

   std::vector<std::vector<std::pair<bool, float> > > construct_H(const std::string &target_sequence,
                                                                  const std::vector<std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > > > &scored_residues);

   void print_alignment(const std::vector<cell_t> &indexed_sequences, // currently reverse order
                        const std::string &sequence,
                        const std::vector<std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > > > &scored_residues);


   std::vector<std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > > >
   get_side_chain_density_scores_for_residues(const std::vector<mmdb::Residue *> &a_run_of_residues,
                                              const clipper::Xmap<float> &xmap);
   

   void apply_alignment_to_model(const std::vector<sm_wat::cell_t> &alignment,
                                 const std::string &target_sequence,
                                 const std::vector<std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > > > &scored_residues);


   // this is the API function that you want
   void align_and_mutate_and_backrub(mmdb::Manager *mol, const std::string &seq, const clipper::Xmap<float> &xmap,
                                     const coot::protein_geometry &pg);
   
}



#endif // SMITH_WATERMAN_HH
