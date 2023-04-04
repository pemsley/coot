
#include <fstream>
#include <iomanip>
#include <thread>
#include <chrono>

#include <boost/math/distributions/skew_normal.hpp>
#include <mmdb2/mmdb_tables.h>  // for mmdb::Get1LetterCode()

#include "compat/coot-sysdep.h"
#include "analysis/stats.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "utils/coot-utils.hh"
#include "side-chain-densities.hh"
#include "richardson-rotamer.hh"

coot::density_box_t::density_box_t(float *density_box_in,
                                   mmdb::Residue *residue_p_in,
                                   int n_steps_in) {
   init();
   density_box = density_box_in;
   residue_p = residue_p_in;
   n_steps = n_steps_in;
   mean=0;
   var = 0;
}


void
coot::side_chain_densities::get_results_addition_lock() {

   bool unlocked = false;
   while (! results_addition_lock.compare_exchange_weak(unlocked, true)) {
      std::this_thread::sleep_for(std::chrono::nanoseconds(100));
      unlocked = false;
   }
}

void
coot::side_chain_densities::release_results_addition_lock() {

   results_addition_lock = false;
}

coot::side_chain_densities::side_chain_densities() {

   int n_steps = 5;
   float grid_box_radius = 5.0; // half the width/length of the box (not diagonal)

   std::string n_steps_str = std::to_string(n_steps);
   std::string grid_box_radius_str = util::float_to_string_using_dec_pl(grid_box_radius, 1);
   std::string pdd = package_data_dir(); // xxx/share/coot
   std::string fn = "useable-grid-points-nstep=" + n_steps_str + ",box_radius=" +
      grid_box_radius_str + "-charybdis.data";
   std::string dir_1 = util::append_dir_dir(pdd,   "data");
   std::string dir_2 = util::append_dir_dir(dir_1, "assign-side-chains");
   std::string pathed_file_name = util::append_dir_file(dir_2, fn);
   std::string side_chain_data_sub_dir = "side-chain-data";

   std::string pathed_scd_dir_name =  util::append_dir_file(dir_2, side_chain_data_sub_dir);

   init(n_steps, grid_box_radius, pathed_file_name);
   set_data_dir(pathed_scd_dir_name);

}

void
coot::density_box_t::self_normalize() {

   int n = 2 * n_steps + 1;
   int nnn = n * n * n;
   double sum = 0;
   int n_grid_points = 0;
   for (int i=0; i<nnn; i++) {
      if (density_box[i] > 0.0) {
         sum += density_box[i];
         n_grid_points++;
      }
   }
   if (n_grid_points > 0) {
      double av = sum/static_cast<double>(n_grid_points);
      double sc = 1.0/av;
      for (int i=0; i<nnn; i++)
         if (density_box[i] > -1000)
            density_box[i] *= sc;
   }
}

bool
coot::side_chain_densities::like_the_others(const std::map<int, std::string> &chain,
                                            const std::vector<std::map<int, std::string> > &other_chains) const {

   bool is_like_the_others = false;
   unsigned int n_chain = chain.size();
   for (std::size_t i=0; i<other_chains.size(); i++) {
      const std::map<int, std::string> &other_chain = other_chains[i];
      unsigned int n = n_chain;
      if (other_chains.size() < n)
         n = other_chains.size();
      unsigned int n_match = 0;
      std::map<int, std::string>::const_iterator it;
      for (it=chain.begin(); it!=chain.end(); ++it) {
         const int &key = it->first;
         const std::string &val = it->second;

         std::map<int, std::string>::const_iterator it_other = other_chain.find(key);
         if (it_other != other_chain.end())
            if (val == it_other->second)
               n_match++;
      }

      float frac = static_cast<float>(n_match) / static_cast<float>(n);
      if (frac > 0.7) {
         is_like_the_others = true;
         break;
      }
   }

   return is_like_the_others;
}

std::map<int, std::string>
coot::side_chain_densities::make_sequence_for_chain(mmdb::Chain *chain_p) const {

   std::map<int, std::string> m;
   int n_residues = chain_p->GetNumberOfResidues();
   for (int i=0; i<n_residues; i++) {
      mmdb::Residue *residue_p = chain_p->GetResidue(i);
      if (residue_p) {
         int res_no = residue_p->GetSeqNum();
         std::string res_name = residue_p->GetResName();
         m[res_no] = res_name;
      }
   }
   return m;

}

void coot::side_chain_densities::proc_mol(const std::string &id,
                                          mmdb::Manager *mol,
                                          const clipper::Xmap<float> &xmap) {

   // don't sample chains that have a sequence similar to chains that we've
   // processed before
   //
   std::vector<std::map<int, std::string> > done_chains;

   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         if (chain_p) {
            std::map<int, std::string> sequence_for_chain = make_sequence_for_chain(chain_p);
            if (! like_the_others(sequence_for_chain, done_chains))
               proc_chain(id, chain_p, xmap);
            done_chains.push_back(sequence_for_chain);
         }
      }
   }

   normalize_density_boxes(id);
   write_density_boxes();
   for (std::size_t i=0; i<density_boxes.size(); i++)
      density_boxes[i].clear();
}

// this function is not used?
void
coot::side_chain_densities::fill_useable_grid_points_vector(const std::string &file_name) {

   if (! file_name.empty()) {
      std::ifstream f(file_name.c_str());
      if (f) {
         std::string line;
         while (std::getline(f, line)) {
            std::vector<std::string> words = util::split_string_no_blanks(line);
            if (words.size() == 1) {
               int idx = util::string_to_int(words[0]);
               useable_grid_points.insert(idx);
            }
         }
      } else {
         std::cout << "ERROR:: side_chain_densities::fill_useable_grid_points_vector file name not found "
                   << file_name << std::endl;
      }
   }

}

bool
coot::side_chain_densities::is_close_to_atoms(const std::vector<std::pair<double, clipper::Coord_orth> > &atom_positions,
                                              const clipper::Coord_orth &test_position) const {

   bool is_close = false;
   for (std::size_t i=0; i<atom_positions.size(); i++) {
      double max_atom_radius_sqrd = atom_positions[i].first;
      double d = (atom_positions[i].second - test_position).lengthsq();
      if (d < max_atom_radius_sqrd) {
         is_close = true;
         break;
      }
   }
   return is_close;

}

std::string coot::side_chain_densities::get_rotamer_name(mmdb::Residue *res) const {
   std::string alt_conf;
   mmdb::Manager *mol = 0;
   richardson_rotamer rr(res, alt_conf, mol, 0.0, 1);
   rotamer_probability_info_t prob = rr.probability_of_this_rotamer();
   std::string rotamer_name = util::remove_whitespace(prob.rotamer_name);
   return rotamer_name;

}

clipper::Coord_orth
coot::side_chain_densities::make_pt_in_grid(int ix, int iy, int iz, const float &step_size,
                                            const std::vector<clipper::Coord_orth> &axes) const {

   clipper::Coord_orth pt(0,0,0);
   pt += ix * step_size * axes[0];
   pt += iy * step_size * axes[1];
   pt += iz * step_size * axes[2];

   return pt;

}

double
coot::side_chain_densities::get_relabun(const std::string &res_name) {

   if (relabun.empty()) {
      // vertibrates
      relabun["ALA"] = 7.4; relabun["ARG"] = 4.2; relabun["ASN"] = 4.4; relabun["ASP"] = 5.9; relabun["CYS"] = 3.3;
      relabun["GLU"] = 5.8; relabun["GLN"] = 3.7; relabun["GLY"] = 7.4; relabun["HIS"] = 2.9; relabun["ILE"] = 3.8;
      relabun["LEU"] = 7.6; relabun["LYS"] = 7.2; relabun["MET"] = 1.8; relabun["PHE"] = 4.0; relabun["PRO"] = 5.0;
      relabun["SER"] = 8.1; relabun["THR"] = 6.2; relabun["TRP"] = 1.3; relabun["TYR"] = 3.3; relabun["VAL"] = 6.8;
      relabun["MSE"] = 1.8;
   }

   const std::map<std::string, double>::const_iterator it = relabun.find(res_name);
   if (it == relabun.end())
      return 2.0;
   else
      return it->second;
};

// return the "guessed" sequence
std::string
coot::side_chain_densities::guess_the_sequence(mmdb::Manager *mol,
                                               const std::string &chain_id,
                                               int resno_start, int resno_end,
                                               const clipper::Xmap<float> &xmap,
                                               bool verbose_output_mode) {

   std::vector<std::pair<mmdb::Residue *, std::string> > best_guess; // a bit of fun

   bool all_chain = false;
   if (resno_end < resno_start) all_chain = true;

   // What is the probability of each rotamer at each residue?
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         if (chain_p) {
            std::string this_chain_id(chain_p->GetChainID());
            if (this_chain_id == chain_id) {
               int n_residues = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<n_residues; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     int res_no = residue_p->GetSeqNum();
                     std::string res_name = residue_p->GetResName();
                     bool do_it = all_chain;
                     if (res_no >= resno_start) {
                        if (res_no <= resno_end) {
                           do_it = true;
                        }
                     }

                     if (do_it) {
                        std::map<std::string, std::pair<std::string, double> > likelihood_map =
                           likelihood_of_each_rotamer_at_this_residue(residue_p, xmap, false, verbose_output_mode);
                        std::map<std::string, std::pair<std::string, double> >::const_iterator it;
                        double best_score = -999999999999999.9;
                        std::string best_type;
                        std::string best_rotamer;
                        // first get the best_score
                        for (it=likelihood_map.begin(); it!=likelihood_map.end(); ++it) {
                           const std::string &res_name_this = it->first;
                           double score = it->second.second + 3.0 * log(get_relabun(res_name_this)); // hacketty hack!
                           if (score > best_score) {
                              best_score = score;
                              best_type = res_name_this;
                              best_rotamer =  it->second.first;
                           }
                        }
                        std::cout << "   " << coot::residue_spec_t(residue_p) << " " << residue_p->GetResName() << std::endl;

                        for (it=likelihood_map.begin(); it!=likelihood_map.end(); ++it) {
                           std::string markup_string;
                           const std::string &res_name_this = it->first;
                           double score = it->second.second + 3.0 * log(get_relabun(res_name_this)); // hacketty hack!

                           {
                              if (best_score < 0) {
                                 if (score > 1.1 * best_score) markup_string = " ooo";
                              } else {
                                 if (score > 0)
                                    if (score > 0.5 * best_score)
                                       markup_string = " ooo";
                              }
                              if (score == best_score) markup_string = " ***";
                           }

                           std::cout << "debug:: guess_the_sequence(): " << residue_spec_t(residue_p) << " " << res_name_this << " "
                                     << std::setw(7) << it->second.first << " "
                                     << std::fixed << std::right << std::setw(8) << std::setprecision(2)
                                     << it->second.second << " "
                                     << std::fixed << std::right << std::setw(5) << std::setprecision(2)
                                     << 10.0 * log(get_relabun(res_name_this)) << markup_string << std::endl;
                        }
                        std::pair<mmdb::Residue *, std::string> p(residue_p, best_type);
                        best_guess.push_back(p);
                     }
                  }
               }
            }
         }
      }
   }

   std::string best_guess_sequence;
   for (std::size_t i=0; i<best_guess.size(); i++) {
      mmdb::Residue *res = best_guess[i].first;
      std::string best_guess_for_residue = best_guess[i].second;
      if (verbose_output_mode) // this is for devel/analysis/testing, where the answer is known
         std::cout << res->GetChainID() << " " << res->GetSeqNum()
                   << " real-type: "  << res->GetResName()
                   << " best-guess: " << best_guess_for_residue << std::endl;
      char r;
      mmdb::Get1LetterCode(best_guess_for_residue.c_str(), &r);
      best_guess_sequence += r;
   }

   // std::cout << "best guess sequence:\n" << best_guess_sequence << std::endl;

   return best_guess_sequence;
}

std::map<std::string, std::pair<std::string, double> >
coot::side_chain_densities::likelihood_of_each_rotamer_at_this_residue(mmdb::Residue *residue_p,
                                                                       const clipper::Xmap<float> &xmap,
                                                                       bool limit_to_correct_rotamers_only,
                                                                       bool verbose_output_mode) {

   auto print_probability_map = [residue_p] (const std::map<std::string, std::pair<std::string, double> > &probability_map) {
      std::map<std::string, std::pair<std::string, double> >::const_iterator it;
      std::cout << "            " << coot::residue_spec_t(residue_p) << " " << residue_p->GetResName() << std::endl;
      for (it=probability_map.begin(); it!=probability_map.end(); ++it) {
         const auto &key = it->first;
         const auto &value = it->second;
         std::cout << "debug:: likelihood_of_each_rotamer_at_this_residue() key " << key << " " << value.first << " ll: " << value.second << "\n";
      }
   };

   std::map<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > >::const_iterator itc =
      best_score_for_res_type_cache.find(residue_p);
   if (itc != best_score_for_res_type_cache.end()) {
      return itc->second;
   } else {
      std::map<std::string, std::pair<std::string, double> > s =
         get_rotamer_likelihoods(residue_p, xmap, limit_to_correct_rotamers_only, verbose_output_mode);
      if (false)
         print_probability_map(s);
      best_score_for_res_type_cache[residue_p] = s;
      return s;
   }
}

std::vector<mmdb::Residue *>
coot::side_chain_densities::make_a_run_of_residues(mmdb::Manager *mol, const std::string &chain_id,
                                                   int resno_start, int resno_end) const {

   std::vector<mmdb::Residue *> a_run_of_residues;

   // What is the probability of each rotamer at each residue?
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         if (chain_p) {
            std::string this_chain_id(chain_p->GetChainID());
            if (this_chain_id == chain_id) {
               int n_residues = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<n_residues; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     int res_no = residue_p->GetSeqNum();
                     if (res_no >= resno_start) {
                        if (res_no <= resno_end) {
                           a_run_of_residues.push_back(residue_p);
                        }
                     }
                  }
               }
            }
         }
      }
   }
   return a_run_of_residues;
}


// return an error messages or a vector of residues
std::pair<std::string, std::vector<mmdb::Residue *> >
coot::side_chain_densities::setup_test_sequence(mmdb::Manager *mol,
                                                const std::string &chain_id, int resno_start, int resno_end,
                                                const clipper::Xmap<float> &xmap) {

   // What is the probability of each rotamer at each residue?

   auto calculate_CB_ideal_pos = [] (mmdb::Residue *residue_p) {
                                    mmdb::Atom *cb_atom = residue_p->GetAtom(" CB ");
                                    bool status = false;
                                    clipper::Coord_orth pos;
                                    if (! cb_atom) {
                                       auto CA = residue_p->GetAtom(" CA ");
                                       auto C  = residue_p->GetAtom(" C  ");
                                       auto N  = residue_p->GetAtom(" N  ");
                                       if (CA) {
                                          if (C) {
                                             if (N) {
                                                clipper::Coord_orth N_pos = co(N);
                                                clipper::Coord_orth C_pos = co(C);
                                                clipper::Coord_orth CA_pos = co(CA);
                                                clipper::Coord_orth C_to_N = N_pos - C_pos;
                                                clipper::Coord_orth C_to_N_mid_point(0.5 * (N_pos + C_pos));
                                                clipper::Coord_orth CA_to_CN_mid_point = C_to_N_mid_point - CA_pos;
                                                clipper::Coord_orth CA_to_CN_mid_point_uv(CA_to_CN_mid_point.unit());
                                                clipper::Coord_orth perp(clipper::Coord_orth::cross(C_to_N, CA_to_CN_mid_point));
                                                clipper::Coord_orth perp_uv(perp.unit());
                                                // guess and fiddle these - good enough
                                                clipper::Coord_orth CB_pos(CA_pos + 1.21 * perp_uv - 0.95 * CA_to_CN_mid_point_uv);
                                                pos = CB_pos;
                                                status =  true;
                                             }
                                          }
                                       }
                                    }
                                    return std::make_pair(status, pos);
                                 };
   

   std::vector<mmdb::Residue *> a_run_of_residues = make_a_run_of_residues(mol, chain_id, resno_start, resno_end);

   if (! a_run_of_residues.empty())
      fill_residue_blocks(a_run_of_residues, xmap); // return fast if already filled

   std::string error_message; // no error message
   for (auto &residue_p : a_run_of_residues) {
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      bool found_N  = false;
      bool found_CA = false;
      bool found_C  = false;
      bool found_O  = false;
      bool found_CB = false;
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            std::string atom_name(at->GetAtomName());
            if (atom_name == " N  ") found_N  = true;
            if (atom_name == " CA ") found_CA = true;
            if (atom_name == " CB ") found_CB = true;
            if (atom_name == " C  ") found_C  = true;
            if (atom_name == " O  ") found_O  = true;
         }
      }
      if (found_N && found_C && found_O && found_CA && found_CB) {
         // happy path
      } else {
         if (found_N && found_C && found_O && found_CA) {
            // invent a CB and add it into the residue
            auto cb_pos = calculate_CB_ideal_pos(residue_p);
            if (cb_pos.first) {
               mmdb::Atom *cb = new mmdb::Atom;
               cb->SetCoordinates(cb_pos.second.x(), cb_pos.second.y(), cb_pos.second.z(), 1.0f, 20.0f);
               cb->SetAtomName(" CB ");
               cb->SetElementName(" C");
               residue_p->AddAtom(cb);
               mol->FinishStructEdit();
            }
         } else {
            std::string chain_id = residue_p->GetChainID();
            int res_no = residue_p->GetSeqNum();
            error_message = "ERROR:: setup-test_sequence(): missing main-chain atom in residue " + chain_id + std::string(" ") + std::to_string(res_no);
         }
      }
   }

   return std::pair<std::string, std::vector<mmdb::Residue *> >(error_message, a_run_of_residues);
}

void
coot::side_chain_densities::setup_likelihood_of_each_rotamer_at_every_residue(const std::vector<mmdb::Residue *> &a_run_of_residues,
                                                                              const clipper::Xmap<float> &xmap) {

   if (! a_run_of_residues.empty()) {
      int n_residues = a_run_of_residues.size();
      for (int i=0; i<n_residues; i++) {
         mmdb::Residue *residue_p = a_run_of_residues[i];
         std::map<std::string, std::pair<std::string, double> > likelihood_map =
            likelihood_of_each_rotamer_at_this_residue(residue_p, xmap); // why is xmap needed here?
                                                                         // residues block should be filled
                                                                         // in setup_test_sequence().
                                                                         // Remove arg xmap from compare_block_vs_rotamer()
                                                                         // and fixup other functions to match
      }
   }
   
}

void coot::side_chain_densities::test_sequence(const std::vector<mmdb::Residue *> &a_run_of_residues,
                                               const clipper::Xmap<float> &xmap,
                                               const std::string &sequence_name,
                                               const std::string &sequence,
                                               bool print_slider_results) {

   auto make_pdb_reference_sequence = [] (const std::vector<mmdb::Residue *> &a_run_of_residues) {
                                         std::string true_sequence;
                                         int n_residues = a_run_of_residues.size();
                                         for (int i=0; i<n_residues; i++) {
                                            mmdb::Residue *residue_p = a_run_of_residues[i];
                                            std::string res_name(residue_p->GetResName());
                                            char r;
                                            mmdb::Get1LetterCode(res_name.c_str(), &r);
                                            true_sequence += r;
                                         }
                                         return true_sequence;
                                      };

   std::map<char, std::string> code_cache;
   auto fill_code_cache = [] (std::map<char, std::string> *code_cache_p) {
                             (*code_cache_p)['G'] = std::string("GLY");
                             (*code_cache_p)['A'] = std::string("ALA");
                             (*code_cache_p)['V'] = std::string("VAL");
                             (*code_cache_p)['S'] = std::string("SER");
                             (*code_cache_p)['N'] = std::string("ASN");
                             (*code_cache_p)['P'] = std::string("PRO");
                             (*code_cache_p)['D'] = std::string("ASP");
                             (*code_cache_p)['C'] = std::string("CYS");
                             (*code_cache_p)['Q'] = std::string("GLN");
                             (*code_cache_p)['E'] = std::string("GLU");
                             (*code_cache_p)['H'] = std::string("HIS");
                             (*code_cache_p)['I'] = std::string("ILE");
                             (*code_cache_p)['L'] = std::string("LEU");
                             (*code_cache_p)['K'] = std::string("LYS");
                             (*code_cache_p)['M'] = std::string("MET");
                             (*code_cache_p)['F'] = std::string("PHE");
                             (*code_cache_p)['T'] = std::string("THR");
                             (*code_cache_p)['W'] = std::string("TRP");
                             (*code_cache_p)['Y'] = std::string("TYR");
                             (*code_cache_p)['R'] = std::string("ARG");
                             (*code_cache_p)['X'] = std::string("");
                          };
   fill_code_cache(&code_cache);
   auto get_res_type_from_single_letter_code = [code_cache] (char c) {
                                                  std::map<char, std::string>::const_iterator it = code_cache.find(c);
                                                  if (it != code_cache.end()) {
                                                     return std::cref(it->second);
                                                  } else {
                                                     std::runtime_error message("miss");
                                                     throw(message);
                                                  }
                                               };

   std::vector<results_t> results;
   std::string gene_name = sequence_name;
   std::vector<std::string> parts = util::split_string_no_blanks(gene_name);
   if (parts.size() > 0) gene_name = parts[0];

   if (! a_run_of_residues.empty()) {
      int n_residues = a_run_of_residues.size();
      std::vector<std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > > > scored_residues(n_residues);
      for (int i=0; i<n_residues; i++) {
         mmdb::Residue *residue_p = a_run_of_residues[i];
         scored_residues[i] =
            std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > >(residue_p, likelihood_of_each_rotamer_at_this_residue(residue_p, xmap));
      }
      
      // insta-fail when the protein sequence for test is shorter than the model.
      if (sequence.length() < a_run_of_residues.size()) return;

      std::string sequence_from_pdb_model = make_pdb_reference_sequence(a_run_of_residues);

      if (print_slider_results)
         std::cout << "----------------- slider ------------ " << std::endl;

      int sequence_length = sequence.length();
      int offset_max = sequence.length() - n_residues; // n_residues is the size of a run of residues

      if (print_slider_results) {
         std::cout << "INFO:: testing sequence " << sequence << std::endl;
         std::cout << "INFO::   model sequence " << sequence_from_pdb_model << std::endl;
         std::cout << "INFO:: offset_max " << offset_max << std::endl;
      }

      for (int offset=0; offset<=offset_max; offset++) {
         // std::cout << "trying loop with offset " << offset << std::endl;
         int n_scored_residues = 0;
         double sum_score = 0;
         std::string running_sequence;
         for (int ires=0; ires<n_residues; ires++) {

            if (false) // Xs in sequences result in sequence count mismatches
               std::cout << "   compare offset" << offset << " ires " << ires << " i+o = "<< (ires+offset)
                         << " and sequence_length " << sequence_length << std::endl;
            if ((ires+offset) < sequence_length) {
               mmdb::Residue *residue_p = scored_residues[ires].first;
               const std::map<std::string, std::pair<std::string, double> > &scored_map = scored_residues[ires].second;

               if (false) { // debug
                  std::map<std::string, std::pair<std::string, double> >::const_iterator it_debug;
                  for (it_debug=scored_map.begin();
                       it_debug!=scored_map.end();
                       ++it_debug) {
                     std::cout << "   " << it_debug->first << " : " << it_debug->second.first << " " << it_debug->second.second << ", ";
                  }
                  if (! scored_map.empty()) std::cout << "\n" << std::endl;
               }

               char letter = sequence[ires+offset];
               try {
                  const std::string &res_type = get_res_type_from_single_letter_code(letter);
                  if (false)
                     std::cout << "   ----------- debug:: offset " << offset << " ires " << ires
                               << " res_type from " << letter << " is \"" << res_type << "\"" << std::endl;
                  if (! res_type.empty()) {
                     std::map<std::string, std::pair<std::string, double> >::const_iterator it = scored_map.find(res_type);
                     if (it != scored_map.end()) {
                        const double &score = it->second.second;
                        running_sequence += letter;
                        sum_score += score;
                        n_scored_residues++;
                        if (false)
                           std::cout << "debug:: adding "
                                     << std::right << std::fixed << std::setw(9) << std::setprecision(4)
                                     << score << " for " << letter << " " << "ires " << ires << " "
                                     << residue_spec_t(residue_p) << std::endl;
                     } else {
                        // res_type was not in the scored_map map. Because it couldn't be scored. The input model didn't
                        // have a CB perhaps?. The scoring function should complain!
                        if (false) { // debug
                           std::cout << "   DEBUG:: test_sequence(): Failed to find " << res_type
                                     << " in this map: which has size " << scored_map.size()
                                     << " for ires " << ires << " offset " << offset << std::endl;
                           std::map<std::string, std::pair<std::string, double> >::const_iterator it_debug;
                           std::cout << "   This was the map: POINT B, it has size " << scored_map.size() << std::endl;
                           for (it_debug=scored_map.begin(); it_debug!=scored_map.end(); ++it_debug) {
                              std::cout << "   " << it_debug->first << " : " << it_debug->second.first << " " << it_debug->second.second << ", ";
                           }
                           if (! scored_map.empty()) std::cout << std::endl;
                        }
                        running_sequence += '.';
                     }
                     // std::cout << "   running sequence is now " << running_sequence << std::endl;
                  } else {
                     if (false) // this happens a lot
                        std::cout << "WARNING:: empty residue type for sequence letter " << letter << std::endl;
                     running_sequence += '-';
                     // no point in continuing because if we don't add to n_scored_residues, then
                     // the test below (n_scored_residues == n_residues) will never be true.
                     break;
                  }
               }
               catch (const std::runtime_error &rte) {
                  // there was an X in the sequence or some other unknown single residue type
                  break;
               }
            }
         }

         if (n_scored_residues == n_residues) {
            results.push_back(results_t(offset, sum_score, n_scored_residues, running_sequence, sequence_name, sequence_from_pdb_model));
            if (print_slider_results)
               std::cout << "INFO:: offset " << offset << " sum_score " << std::setw(10) << sum_score
                         << " n_scored_residues " << n_scored_residues << " " << running_sequence
                         << " gene-name " << gene_name
                         << " model-sequence " << sequence_from_pdb_model << std::endl;
         } else {
            if (print_slider_results) // happens a lot due to Xs in sequence
               std::cout << "INFO:: Failed to push back a result because " << n_scored_residues << " != " << n_residues << std::endl;
         }
      }
      // auto tp_1 = std::chrono::high_resolution_clock::now();
      // auto d10 = std::chrono::duration_cast<std::chrono::microseconds>(tp_1 - tp_0).count();
      // std::cout << "TIMINGS:: test_sequence() " << d10 << " microseconds" << std::endl;
   }
   get_results_addition_lock();
   // std::cout << "storing results of size " << results.size() << " for sequence with name " << sequence_name << std::endl;
   results_container[sequence_name] = results;
   release_results_addition_lock();
}

#include "analysis/stats.hh"

coot::side_chain_densities::results_t
coot::side_chain_densities::get_result(bool only_probable, bool print_sequencing_solutions_flag) const {

   std::cout << "------ Here in get_result() " << only_probable << " " << print_sequencing_solutions_flag << std::endl;

   auto print_sequencing_solutions = [] (const std::map<std::string, std::vector<results_t> > &results_container) {
                                        unsigned int n_top = 20;

                                        if (results_container.empty()) return;

                                        std::cout << "INFO:: print_sequencing_solutions() given " << results_container.begin()->second.size()
                                                  << " results for the first sequence" << std::endl;

                                        std::vector<results_t> all_results;
                                        std::map<std::string, std::vector<results_t> >::const_iterator it;
                                        for (it=results_container.begin(); it!=results_container.end(); ++it) {
                                           const auto &results = it->second;
                                           for (unsigned int i=0; i<results.size(); i++) {
                                              auto &result(results[i]);
                                              all_results.push_back(result);
                                           }
                                        }
                                        auto results_sorter = [] (const results_t &r1, const results_t &r2) {
                                                                 return (r2.sum_score < r1.sum_score);
                                                              };
                                        std::sort(all_results.begin(), all_results.end(), results_sorter);
                                        if (all_results.size() < n_top) n_top = all_results.size();
                                        for (unsigned int i=0; i<n_top; i++) {
                                           std::cout << std::setw(2) << i << " " << std::setw(10) << all_results[i].sum_score << " "
                                                     << all_results[i].sequence << std::endl;
                                        }
                                     };

   double best_score = -9e10;
   results_t best_results;

   std::map<std::string, std::vector<results_t> >::const_iterator it;
   for (it=results_container.begin(); it!=results_container.end(); ++it) {
      const std::vector<results_t> &v = it->second;
      for (unsigned int i=0; i<v.size(); i++) {
         const results_t &result = v[i];
         if (result.sum_score > best_score) {
            best_score = result.sum_score;
            best_results = result;
         }
      }
   }

   if (only_probable) {
      // let's do some statistics - is the top hit much better than the rest?
      std::vector<double> scores;
      for (it=results_container.begin(); it!=results_container.end(); ++it) {
         const std::vector<results_t> &v = it->second;
         for (unsigned int i=0; i<v.size(); i++) {
            scores.push_back(v[i].sum_score);
         }
      }
      if (scores.size() > 2) {
         // do we have a clear solution? If not, clear the returned sequence
         std::sort(scores.begin(), scores.end());
         std::reverse(scores.begin(), scores.end());
         double top_score = scores[0];
         double next_best_score = scores[1];
         double top_solution_delta = top_score - next_best_score;
         unsigned int n_data = 21;
         if (scores.size() < n_data) n_data = scores.size();
         stats::single data;
         for (unsigned int i=1; i<n_data; i++)
            data.add(scores[i]);
         double var = data.variance();
         double sd = std::sqrt(var);
         std::cout << "INFO:: get_result(): top_solution_delta: " << top_solution_delta << " vs s.d. others: " << sd << std::endl;
         std::cout << "INFO:: get_result(): ratio of top score delta to std dev others: " << top_solution_delta/sd;
         if (top_solution_delta > 3.0 * sd) {
            std::cout << " ###" << std::endl; // give us a hash marker to let us know it's a good one
            // happy path
         } else {
            std::cout << std::endl;
            // clear the solution
            std::cout << "No clear solution found" << std::endl;
            best_results.sequence.clear();
         }

         if (print_sequencing_solutions_flag) {
            print_sequencing_solutions(results_container);
         }
      } else {
         std::cout << "WARNING:: in get_result() only found " << scores.size() << " scored results" << std::endl;
      }
   } else {
      if (print_sequencing_solutions_flag)
         print_sequencing_solutions(results_container);
   }

   return best_results;
}

bool
coot::side_chain_densities::test_grid_point_to_coords_interconversion() const {

   bool success = true;
   int n_per_side = 2 * n_steps + 1;
   for (int ix= -n_steps; ix<=n_steps; ix++) {
      for (int iy= -n_steps; iy<=n_steps; iy++) {
         for (int iz= -n_steps; iz<=n_steps; iz++) {
            int idx =
               (ix + n_steps) * n_per_side * n_per_side +
               (iy + n_steps) * n_per_side +
               (iz + n_steps);

            std::tuple<int, int, int> t = grid_index_to_grid(idx);
            if (std::get<0>(t) == ix)
               if (std::get<1>(t) == iy)
                  if (std::get<2>(t) == iz)
                     if (false) // noise - if the function works :-)
                        std::cout << "PASS " << std::setw(2) << ix << " " << std::setw(2) << iy << " " << std::setw(2) << iz << " idx " << idx
                                  << " decode " << std::get<0>(t) << " " << std::get<1>(t) << " " << std::get<2>(t)
                                  << std::endl;

            if (std::get<0>(t) != ix) {
               std::cout << "FAIL " << std::setw(2) << ix << " " << std::setw(2) << iy << " " << std::setw(2) << iz << " idx " << idx
                         << " decode " << std::get<0>(t) << " " << std::get<1>(t) << " " << std::get<2>(t)
                         << std::endl;
               success = false;
            }
            if (std::get<1>(t) != iy) {
               std::cout << "FAIL " << std::setw(2) << ix << " " << std::setw(2) << iy << " " << std::setw(2) << iz << " idx " << idx
                         << " decode " << std::get<0>(t) << " " << std::get<1>(t) << " " << std::get<2>(t)
                         << std::endl;
               success = false;
            }
            if (std::get<2>(t) != iz) {
               std::cout << "FAIL " << std::setw(2) << ix << " " << std::setw(2) << iy << " " << std::setw(2) << iz << " idx " << idx
                         << " decode " << std::get<0>(t) << " " << std::get<1>(t) << " " << std::get<2>(t)
                         << std::endl;
               success = false;
            }
         }
      }
   }

   return success;

}


// gen_pts_file_name is optional argument for generating initial usable grid points.
// Calling function must clear the returned density_box_t
//
coot::density_box_t
coot::side_chain_densities::sample_map(mmdb::Residue *residue_this_p,
                                       mmdb::Residue *residue_next_p,
                                       mode_t mode,
                                       const clipper::Coord_orth &cb_pt,
                                       const std::vector<clipper::Coord_orth> &axes,
                                       const clipper::Xmap<float> &xmap,
                                       std::string gen_pts_file_name) const {

   bool gen_usable_points_flag = false;
   if (mode == GEN_USABLE_POINTS) gen_usable_points_flag = true;

   // std::cout << "In sample_map() Here 1 " << residue_this_p << " " << residue_next_p
   //           << " gen_usable_points_flag " << gen_usable_points_flag << std::endl;

   // 3 modes:
   //
   // (1) Generating the grid point list
   // (2) Sampling the maps to provide data for the distributions
   // (3) Sampling the user's map to generate a grid of data to test against the distributions.

   // This function is a bit ugly, because it used to generate grid points
   // and to sample the map - it is not split into 2 or 3 functions - which would
   // be cleaner, because I want only one version of the code that converts
   // from grid indices to 3D points.

   float step_size = grid_box_radius/static_cast<float>(n_steps);
   int n_per_side = 2 * n_steps + 1;
   int n_box_vol = n_per_side * n_per_side * n_per_side;

   if (! residue_this_p) return density_box_t(0,0,0);
   if (axes.empty())     return density_box_t(0,0,0);

   std::string res_name = residue_this_p->GetResName();
   std::string rot_name = get_rotamer_name(residue_this_p); // doesn't matter for SAMPLE_FOR_RESIDUE mode

   if (mode == SAMPLE_FOR_DB) {
      if (rot_name.empty()) {
         if (res_name == "GLY") {
            rot_name = "pseudo";
         } else {
            return density_box_t(0,0,0);
         }
      }
   }

   clipper::Coord_orth ca_pt(-1,-1,-1);

   // std::cout << "In sample_map() Here 2 " << residue_this_p << " " << residue_next_p << std::endl;

   int n_atoms = residue_this_p->GetNumberOfAtoms();
   std::vector<clipper::Coord_orth> residue_atom_positions;
   std::vector<std::pair<double, clipper::Coord_orth> > main_chain_atom_positions;
   residue_atom_positions.reserve(n_atoms);
   for (int i=0; i<n_atoms; i++) {
      mmdb::Atom *at = residue_this_p->GetAtom(i);
      if (! at->isTer()) {
         clipper::Coord_orth pos = co(at);
         residue_atom_positions.push_back(pos);
         std::string atom_name(at->name);
         if (atom_name == " N  " || atom_name == " C  " || atom_name == " O  " ||
             atom_name == " H  " || atom_name == " CA ") {
            double r = 2.8;
            if (atom_name == " CA ")
               r = 1.6;
            std::pair<double, clipper::Coord_orth> p(r, pos);
            main_chain_atom_positions.push_back(p);
         }
         if (atom_name == " CA ") ca_pt = co(at);
      }
   }
   if (gen_usable_points_flag) {
      // need to reject points around N of next residue too
      n_atoms = residue_next_p->GetNumberOfAtoms();
      for (int i=0; i<n_atoms; i++) {
         mmdb::Atom *at = residue_next_p->GetAtom(i);
         if (! at->isTer()) {
            std::string atom_name(at->name);
            if (atom_name == " N  ") {
               clipper::Coord_orth pos = co(at);
               std::pair<double, clipper::Coord_orth> p(3.0, pos);
               main_chain_atom_positions.push_back(p);
            }
         }
      }
   }

   stats::single block_stats;

   // std::cout << "In sample_map() Here 3 " << residue_this_p << " " << residue_next_p << std::endl;

   std::ofstream f; // for generating the (index) side chain points file
   if (gen_usable_points_flag)
      f.open(gen_pts_file_name.c_str());

   // use shared pointer
   float *density_box = new float[n_box_vol];

   float unset_value = -1001.1;
   for (int i=0; i<n_box_vol; i++) density_box[i] = unset_value;

   // make main_chain_clashing_points - these should be the same for every residue
   // i.e. use something like this to make a reference list once
   // the grid size and n_steps is optimized.
   //
   // Note:  exclude grid points that are within 3(?)A of C, O, N and next N
   //
   std::set<int> main_chain_clashing_points;
   if (gen_usable_points_flag) {
      for (int ix= -n_steps; ix<=n_steps; ix++) {
         for (int iy= -n_steps; iy<=n_steps; iy++) {
            for (int iz= -n_steps; iz<=n_steps; iz++) {
               int idx =
                  (ix + n_steps) * n_per_side * n_per_side +
                  (iy + n_steps) * n_per_side +
                  (iz + n_steps);
               clipper::Coord_orth pt_in_grid = make_pt_in_grid(ix, iy, iz, step_size, axes);
               clipper::Coord_orth pt_grid_point = cb_pt + pt_in_grid;
               clipper::Coord_orth ca_to_cb = cb_pt - ca_pt;
               clipper::Coord_orth ca_to_pt = pt_grid_point - ca_pt;
               double dp = clipper::Coord_orth::dot(ca_to_cb, ca_to_pt);
               if (dp < 0.0)
                  main_chain_clashing_points.insert(idx);
            }
         }
      }
   }

   // std::cout << "In sample_map() Here 4 " << residue_this_p << " " << residue_next_p << std::endl;
   
   for (int ix= -n_steps; ix<=n_steps; ix++) {
      for (int iy= -n_steps; iy<=n_steps; iy++) {
         for (int iz= -n_steps; iz<=n_steps; iz++) {
            int idx =
               (ix + n_steps) * n_per_side * n_per_side +
               (iy + n_steps) * n_per_side +
               (iz + n_steps);
            // note: when generating useable points, this set is empty
            if (gen_usable_points_flag || useable_grid_points.find(idx) != useable_grid_points.end()) {
               clipper::Coord_orth pt_in_grid = make_pt_in_grid(ix, iy, iz, step_size, axes);
               clipper::Coord_orth pt_grid_point = cb_pt + pt_in_grid;
               if (gen_usable_points_flag) {
                  if (! is_close_to_atoms(main_chain_atom_positions, pt_grid_point)) {
                     if (in_sphere(pt_grid_point, cb_pt, grid_box_radius)) {
                        if (main_chain_clashing_points.find(idx) == main_chain_clashing_points.end()) {
                           f << "setting grid point " << idx << " at "
                             << pt_grid_point.x() <<  " "
                             << pt_grid_point.y() <<  " "
                             << pt_grid_point.z() << std::endl;
                        }
                     }
                  }
               } else {
                  // more usual, the decision about acceptable grid points have already
                  // been made (some time ago).
                  float dv = util::density_at_point_by_linear_interpolation(xmap, pt_grid_point);
                  density_box[idx] = dv;
                  block_stats.add(dv);
                  if (false)
                     std::cout << "debug:: in sample_map(): for block_stats " << pt_grid_point.format()
                               <<  " " << idx << " " << dv << "\n";
               }
            }
         }
      }
   }

   if (false)
      std::cout << "debug: In sample_map(): block stats " << coot::residue_spec_t(residue_this_p)
                << " size " << block_stats.size()
                << " mean " << block_stats.mean() << " sd " << sqrt(block_stats.variance()) << std::endl;


   density_box_t db(density_box, residue_this_p, n_steps); // density_box is a float *.
   db.mean = block_stats.mean();
   db.var  = block_stats.variance();

   // return negative values on failure
   std::tuple<double, double, double> ca_stats = get_stats_around_ca(residue_this_p, axes,
                                                                     0.5 * step_size, xmap);
   if (std::get<1>(ca_stats) > 0)
      db.set_around_ca_stats(std::get<0>(ca_stats), std::get<1>(ca_stats), std::get<2>(ca_stats));

   if (false)
      std::cout << "in sample_map() returning db: " << db.mean << " " << db.var << std::endl;

   return db;
}

// return negative values on failure
std::tuple<double, double, double>
coot::side_chain_densities::get_stats_around_ca(mmdb::Residue *residue_this_p,
                                                const std::vector<clipper::Coord_orth> &axes,
                                                float step_size,
                                                const clipper::Xmap<float> &xmap) const {

   double mean = 0;
   double var = 0;
   double mean_of_positives = 0;

   int n_atoms = residue_this_p->GetNumberOfAtoms();
   mmdb::Atom *ca_at = 0;
   for (int i=0; i<n_atoms; i++) {
      mmdb::Atom *at = residue_this_p->GetAtom(i);
      std::string atom_name(at->GetAtomName());
      if (atom_name == " CA ") {
         ca_at = at;
         break;
      }
   }
   if (ca_at) {

      stats::single s;
      stats::single s_above_zero;
      clipper::Coord_orth ca_pos = co(ca_at);
      int ilim = n_steps * n_steps;
      for (int ix= -n_steps; ix<=n_steps; ix++) {
         for (int iy= -n_steps; iy<=n_steps; iy++) {
            for (int iz= -n_steps; iz<=n_steps; iz++) {
               if ((ix * ix + iy * iy + iz * iz) <= ilim) {
                  clipper::Coord_orth pt_in_grid = make_pt_in_grid(ix, iy, iz, step_size, axes);
                  clipper::Coord_orth pt_grid_point = ca_pos + pt_in_grid;
                  float dv = util::density_at_point_by_linear_interpolation(xmap, pt_grid_point);
                  s.add(dv);
                  if (dv > 0.0)
                     s_above_zero.add(dv);
               }
            }
         }
      }

      mean = s.mean();
      var = s.variance();
      mean_of_positives = s_above_zero.mean();

   } else {
      // failure status
      mean = -1;
      var = -1;
      mean_of_positives = -1;
   }

   return std::tuple<double, double, double> (mean, var, mean_of_positives);
}


void
coot::side_chain_densities::normalize_density_boxes(const std::string &id) {

   // We are normalizing the boxes from our map/model (not reference data)

   // normalize_density_boxes_v2(id);

   normalize_density_boxes_v3(id); // Use CA density for scaling
}

void
coot::side_chain_densities::normalize_density_boxes_v1(const std::string &id) {

   // hacketty-hack scaling

   float sum = 0;
   float sum_sq = 0;
   int n_grid_pts = 0;

   for (std::size_t i=0; i<density_boxes.size(); i++) {
      int nnn = density_boxes[i].nnn();
      const density_box_t &db = density_boxes[i];
      for (int j=0; j<nnn; j++) {
         if (db[j] > 0.0) {
            sum += db[j];
            sum_sq += db[j] * db[j];
            n_grid_pts++;
         }
      }
   }

   if (n_grid_pts > 0) {
      float mean = sum/static_cast<float>(n_grid_pts);
      float var = sum_sq/static_cast<float>(n_grid_pts) - mean * mean;
      float scale_factor = 1.0/mean;
      std::cout << "Dataset from " << id << " mean " << mean << " scale_factor "
                << scale_factor << std::endl;
      for (std::size_t i=0; i<density_boxes.size(); i++) {
         density_box_t &db = density_boxes[i];
         db.scale_by(scale_factor); // don't scale "below zero" points
      }
   }
}

void
coot::side_chain_densities::normalize_density_boxes_v2(const std::string &id) {

   // Make the RMSd be 1.0 then the the mean to zero

   for (std::size_t i=0; i<density_boxes.size(); i++) {
      int n_grid_pts = 0;
      float sum = 0;
      float sum_sq = 0;
      int nnn = density_boxes[i].nnn();
      density_box_t &db = density_boxes[i];
      for (int j=0; j<nnn; j++) {
         if (db[j] > -1000.0) {
            sum += db[j];
            sum_sq += db[j] * db[j];
            n_grid_pts++;
         }
      }

      if (n_grid_pts > 0) {
         float mean = sum/static_cast<float>(n_grid_pts);
         float var = sum_sq/static_cast<float>(n_grid_pts) - mean * mean;
         if (var > 0.0) {
            float sd = sqrt(var);
            float scale_factor = 1.0/sd;
            for (int j=0; j<nnn; j++) {
               if (db[j] > -1000.0) {
                  db.density_box[j] *= scale_factor;
               }
            }
            sum = 0;
            for (int j=0; j<nnn; j++) {
               if (db[j] > -1000.0) {
                  sum += db[j];
               }
            }
            mean = sum/static_cast<float>(n_grid_pts);
            for (int j=0; j<nnn; j++) {
               if (db[j] > -1000.0) {
                  db.density_box[j] -= mean;
               }
            }
         }
      }
   }
}

void
coot::side_chain_densities::normalize_density_boxes_v3(const std::string &id) {

   for (std::size_t i=0; i<density_boxes.size(); i++) {
      density_box_t &db = density_boxes[i];
      db.normalize_using_ca_stats();
   }
}

void
coot::density_box_t::normalize_using_ca_stats() {

   if (! density_box) return; // clang scan-build fixup. So that
   if (var_around_ca > 0) {
      int n = nnn();
      if (mean_of_positives_around_ca > 0.0) {
         if (false)
            std::cout << "Debug: normalizing with mean of CA positives "
                      << mean_of_positives_around_ca << " for "
                      << residue_spec_t(residue_p) << std::endl;
         float sf = 0.995 * 1.0/mean_of_positives_around_ca;
         unsigned int n_scaled = 0;
         for (int j=0; j<n; j++) {
            if (density_box[j] > -1000.0) {
               density_box[j] *= sf;
               n_scaled++;
            }
         }
         // std::cout << "DEBUG:: n_scaled " << n_scaled << std::endl;
      } else {
         // mark as a baddie - perhaps this
         // should have its own flag?
         var_around_ca = -1.0;
         is_weird = true;
      }
   } else {
      std::string rn;
      is_weird = true;
      if (residue_p)
         rn = residue_p->GetResName();
      std::cout << "ERROR:: Failed variance in normalize_using_ca_stats() "
                << residue_spec_t(residue_p) << " " << rn << std::endl;
   }

}



void
coot::side_chain_densities::write_density_boxes() const {

   for (std::size_t i=0; i<density_boxes.size(); i++) {
      const density_box_t &db = density_boxes[i];
      // don't write boxes with weird/non-filled stats
      if (db.var_around_ca > 0.0)
         if (! db.is_weird)
            write_density_box(db, id);
   }
}

std::pair<float, float>
coot::density_box_t::mean_and_variance() const {

   // of non-masked points above the zero

   stats::single s;
   int n = 2 * n_steps + 1;
   int count = 0;
   int nnn = n * n * n;
   for (int i=0; i<nnn; i++) {
      const float &d = density_box[i];
      if (d > 0.0) {
         s.add(d);
         count++;
      }
   }

   float mean = 99999999999.9;
   float var = 0;
   if (count > 0) {
      mean = s.mean();
      var = s.variance();
   }
   return std::pair<float, float> (mean, var);
}

void
coot::side_chain_densities::write_density_box(const coot::density_box_t &db, const std::string &id) const {

   float *density_box = db.density_box;
   int n_steps = db.n_steps;
   mmdb::Residue *residue_p = db.residue_p;

   if (! residue_p) return;
   std::string res_name = residue_p->GetResName();
   std::string rot_name = get_rotamer_name(residue_p);
   std::string dir = "side-chain-data";

   if (!rot_name.empty()) {
      std::string rot_dir = dir + "/" + res_name + "/" + rot_name;
      std::string file_name = rot_dir + "/" + id + "-" + residue_p->GetChainID() + util::int_to_string(residue_p->GetSeqNum()) + ".tab";

      if (! file_exists(rot_dir))
         util::create_directory(rot_dir);
      // std::cout << "write_density_box() to filename " << file_name << std::endl;

      std::ofstream f(file_name.c_str());
      if (f) {
         int n_per_side = 2 * n_steps + 1;
         int n_box_vol = n_per_side * n_per_side * n_per_side;
         for (int i=0; i<n_box_vol; i++) {
            float v = density_box[i];
            if (clipper::Util::is_nan(v)) {
               std::cout << "ERROR:: " << file_name << " found a nan " << i << std::endl;
               exit(1);
            } else {
               f << v << " ";
               if (i%n_per_side == 0)
                  f << "\n";
            }
         }
         f << "\n";
      } else {
         std::cout << "WARNING:: cannot open file " << file_name << std::endl;
      }
   } else {
      // these need investigation - Coot bug?
      bool verbose = false;
      if (verbose)
         std::cout << "WARNING:: no rotamer name for " << id << " " << residue_spec_t(residue_p)
                   << " " << res_name << std::endl;
   }
}

// stores density boxes - change the name of the function
void
coot::side_chain_densities::proc_chain(const std::string &id, mmdb::Chain *chain_p,
                                       const clipper::Xmap<float> &xmap) {

   int n_residues = chain_p->GetNumberOfResidues();
   int last_res = n_residues - 2;
   for (int ires=1; ires<=last_res; ires++) {
      mmdb::Residue *t = chain_p->GetResidue(ires);
      if (t) {

        std::string rn = t->GetResName();
        if (rn == "UNK") continue;

        if (! util::is_standard_amino_acid_name(rn)) continue;

         // don't forget that ALA are useful to search, as is GLY, but
         // that will need a special function to find an imaginary CB position
         //
         std::pair<clipper::Coord_orth, std::vector<clipper::Coord_orth> > cb_pos_and_axes = get_residue_axes(t);
         const std::vector<clipper::Coord_orth> &axes = cb_pos_and_axes.second;
         if (! axes.empty()) {
            // sample_masked density around 8A around CB
            clipper::Coord_orth cb_pt = cb_pos_and_axes.first;
            // bool gen_flag = false;
            mode_t mode = SAMPLE_FOR_DB;
            density_box_t db = sample_map(t, 0, mode, cb_pt, axes, xmap);
            if (! db.empty()) {
               // std::cout << "Storing density box for residue " << residue_spec_t(t) << " "
               //           << residue_spec_t(db.residue_p) << std::endl;
               store_density_box(db); // push back to density_boxes vector
            }
         }
      }
   }
}

#include "coot-utils/c-beta-deviations.hh"

std::pair<clipper::Coord_orth, std::vector<clipper::Coord_orth> >
coot::side_chain_densities::get_residue_axes_type_GLY(mmdb::Residue *this_residue) const {

   mmdb::Atom *CA_at = 0;
   mmdb::Atom *N_at  = 0;
   mmdb::Atom *C_at  = 0;
   int n_atoms = this_residue->GetNumberOfAtoms();
   for (int i=0; i<n_atoms; i++) {
      mmdb::Atom *at = this_residue->GetAtom(i);
      std::string atom_name = at->name;
      std::string alt_loc = at->altLoc;
      if (! at->isTer()) {
         if (alt_loc.empty()) {
            if (atom_name == " CA ") CA_at = at;
            if (atom_name == " N  ")  N_at = at;
            if (atom_name == " C  ")  C_at = at;
         }
      }
   }
   if (N_at && CA_at && C_at) {
      atom_quad q(N_at, CA_at, C_at, 0); // only 3 atoms used
      clipper::Coord_orth cb_pos = make_CB_ideal_pos(q, "ALA");
      clipper::Coord_orth ca_pos = co(CA_at);
      clipper::Coord_orth c_pos = co(C_at);
      clipper::Coord_orth n_pos = co(N_at);
      std::vector<clipper::Coord_orth> axes = make_axes(ca_pos, cb_pos, c_pos, n_pos);
      return std::pair<clipper::Coord_orth, std::vector<clipper::Coord_orth> > (cb_pos, axes);
   } else {
      std::cout << "ERROR:: BAD GLY " << residue_spec_t(this_residue) << std::endl;
      std::vector<clipper::Coord_orth> axes;
      clipper::Coord_orth cb_pos(0,0,0);
      return std::pair<clipper::Coord_orth, std::vector<clipper::Coord_orth> > (cb_pos, axes);
   }
}

std::vector<clipper::Coord_orth>
coot::side_chain_densities::make_axes(const clipper::Coord_orth &pt_ca_this,
                                      const clipper::Coord_orth &pt_cb_this,
                                      const clipper::Coord_orth &pt_c_this,
                                      const clipper::Coord_orth &pt_n_this) const {

   std::vector<clipper::Coord_orth> v;
   // Here add a bit of noise to pt_cb_this to "sample"
   // density from non-ideal orientation

   clipper::Coord_orth axis_1((pt_cb_this - pt_ca_this).unit());
   clipper::Coord_orth nc(pt_c_this - pt_n_this);
   clipper::Coord_orth nc_uv(nc.unit());
   clipper::Coord_orth cp_1(clipper::Coord_orth::cross(axis_1, nc_uv));
   clipper::Coord_orth cp_2(clipper::Coord_orth::cross(cp_1, axis_1));

   clipper::Coord_orth axis_2(cp_1.unit());
   clipper::Coord_orth axis_3(cp_2.unit());

   // double dp = clipper::Coord_orth::dot(nc_uv, axis_3);
   // double theta = acos(dp);
   // std::cout << residue_spec_t(this_residue) << " "
   // << clipper::Util::rad2d(theta) << std::endl;
   v.push_back(axis_1);
   v.push_back(axis_2);
   v.push_back(axis_3);

   return v;
}


std::pair<clipper::Coord_orth, std::vector<clipper::Coord_orth> >
coot::side_chain_densities::get_residue_axes(mmdb::Residue *residue_p) const {

   // only look for atoms with alt conf ""

   std::string res_name(residue_p->GetResName());
   if (res_name == "GLY") return get_residue_axes_type_GLY(residue_p);

   // OK, we have a CB (presumably)

   std::vector<clipper::Coord_orth> v;
   clipper::Coord_orth cb_pos(0,0,0);

   mmdb::Atom *CA_at = 0;
   mmdb::Atom *CB_at = 0;
   mmdb::Atom *N_at  = 0;
   mmdb::Atom *C_at  = 0;
   mmdb::Atom *O_at  = 0;

   int n_atoms = residue_p->GetNumberOfAtoms();
   for (int i=0; i<n_atoms; i++) {
      mmdb::Atom *at = residue_p->GetAtom(i);
      std::string atom_name = at->name;
      std::string alt_loc = at->altLoc;
      if (! at->isTer()) {
         if (alt_loc.empty()) {
            if (atom_name == " CA ") CA_at = at;
            if (atom_name == " CB ") CB_at = at;
            if (atom_name == " N  ")  N_at = at;
            if (atom_name == " C  ")  C_at = at;
            if (atom_name == " O  ")  O_at = at;
         }
      }
   }

   if (CA_at && CB_at && N_at && C_at && O_at) {

         clipper::Coord_orth pt_ca_this = co(CA_at);
         clipper::Coord_orth pt_cb_this = co(CB_at);
         clipper::Coord_orth pt_c_this = co(C_at);
         clipper::Coord_orth pt_n_this = co(N_at);

         v = make_axes(pt_ca_this, pt_cb_this, pt_c_this, pt_n_this);
         cb_pos = co(CB_at);
   } else {

      if (CA_at && N_at && C_at && O_at) {
         return get_residue_axes_type_GLY(residue_p);
      } else {

         std::string missing_atoms;
         if (! CA_at) missing_atoms += " CA ";
         if (! CB_at) missing_atoms += " CB ";
         if (! C_at) missing_atoms += " C ";
         if (! N_at) missing_atoms += " N ";
         std::cout << "WARNING:: in get_residue_axes(): missing atom(s) " << missing_atoms << " for "
                   << coot::residue_spec_t(residue_p) << " " << std::endl;
      }
   }
   std::pair<clipper::Coord_orth, std::vector<clipper::Coord_orth> > p(cb_pos, v);
   return p;
}


// This function is not used - it could be deleted.
void
coot::side_chain_densities::gen_useable_grid_points(mmdb::Residue *residue_this_p,
                                                    mmdb::Residue *residue_next_p,
                                                    int n_steps, float grid_box_radius,
                                                    const std::string &gen_pts_file_name) const {

   std::pair<clipper::Coord_orth, std::vector<clipper::Coord_orth> > cb_pos_and_axes =
      get_residue_axes(residue_this_p);
   const std::vector<clipper::Coord_orth> &axes = cb_pos_and_axes.second;
   if (! axes.empty()) {
      // sample_masked density around 8A around CB
      clipper::Coord_orth cb_pt = cb_pos_and_axes.first;

      // make a block
      mode_t mode = GEN_USABLE_POINTS;
      clipper::Xmap<float> dummy;
      density_box_t block = sample_map(residue_this_p, residue_next_p, mode, cb_pt,
                                       axes, dummy, gen_pts_file_name);
   }
}

std::pair<std::string, std::string>
coot::side_chain_densities::map_key_to_residue_and_rotamer_names(const std::string &key) const {

   std::string::size_type pos = key.find_last_of(":");
   std::string residue_name = key.substr(0, pos);
   std::string rotamer_name = key.substr(pos+1, key.length());

   return std::pair<std::string, std::string>(residue_name, rotamer_name);
}

void
coot::side_chain_densities::fill_residue_blocks(mmdb::Manager *mol, const std::string &chain_id,
                                                int resno_start, int resno_end,
                                                const clipper::Xmap<float> &xmap) {

   const std::vector<mmdb::Residue *> residues = make_a_run_of_residues(mol, chain_id, resno_start, resno_end);
   fill_residue_blocks(residues, xmap);

}

// we can do better normalization of the grids for the user/test structure if
// we do them all at the same time.
void
coot::side_chain_densities::fill_residue_blocks(const std::vector<mmdb::Residue *> &residues,
                                                const clipper::Xmap<float> &xmap) {

   // call this function only once when multithreaded. We need only fill the residue blocks once.
   // I mean, don't wrap this function with multi-threads.

   if (density_block_map_cache.size() > 0) {
      return; // already done
   } else {
      auto tp_0 = std::chrono::high_resolution_clock::now();
      for (std::size_t i=0; i<residues.size(); i++) {
         mmdb::Residue *residue_p = residues[i];
         mode_t mode = SAMPLE_FOR_RESIDUE;
         std::pair<clipper::Coord_orth, std::vector<clipper::Coord_orth> > cb_pos_and_axes = get_residue_axes(residue_p);
         const clipper::Coord_orth &cb_pt = cb_pos_and_axes.first;
         const std::vector<clipper::Coord_orth> &axes = cb_pos_and_axes.second;
         density_box_t block = sample_map(residue_p, 0, mode, cb_pt, axes, xmap);
         if (false)
            std::cout << "debug:: fill_residue_blocks(): sample map for " << coot::residue_spec_t(residue_p)
                      << " mean " << block.mean << " var " << block.var << std::endl;
         block.normalize_using_ca_stats();
         density_block_map_cache[residue_p] = block;
      }
      add_mean_and_variance_to_individual_density_blocks();
      auto tp_1 = std::chrono::high_resolution_clock::now();
      auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
      std::cout << "TIMINGS:: fill_residue_blocks() " << d10 << " milliseconds" << std::endl;
   }

}

void
coot::side_chain_densities::normalize_density_blocks() {

   std::map<mmdb::Residue *, density_box_t>::const_iterator it;
   unsigned int n_grid_pts = 0;
   double sum = 0;

   stats::single s_all;
   for(it=density_block_map_cache.begin(); it!=density_block_map_cache.end(); ++it) {
      const density_box_t &block = it->second;
      if (! block.empty()) {
         int nnn = block.nnn();
         for (int i=0; i<nnn; i++) {
            if (block[i] > 0.0) {
               sum += block[i];
               n_grid_pts++;
            }
            s_all.add(block[i]);
         }
      }
   }
   if (n_grid_pts > 0) {
      double av = sum/static_cast<double>(n_grid_pts);
      double sc = mn_scale_for_normalized_density/av;
      std::map<mmdb::Residue *, density_box_t>::iterator it_inner;
      for(it_inner=density_block_map_cache.begin(); it_inner!=density_block_map_cache.end(); ++it_inner) {
         density_box_t &block = it_inner->second; // because not const
         block.scale_by(sc);
      }
   }
}

void
coot::side_chain_densities::add_mean_and_variance_to_individual_density_blocks() {

   stats::single s;
   stats::single s_positive;
   std::map<mmdb::Residue *, density_box_t>::iterator it;
   for(it=density_block_map_cache.begin(); it!=density_block_map_cache.end(); ++it) {
      density_box_t &block = it->second;
      if (! block.empty()) {
         int nnn = block.nnn();
         for (int i=0; i<nnn; i++) {
            const float &bi = block[i];
            if (bi > -1000.0)
               s.add(bi);
            if (bi > 0.0)
               s_positive.add(bi);
         }
         if (false)
            std::cout << "debug:: in add_mean_and_variance_to_individual_density_blocks() "
                      << residue_spec_t(it->first) << " mean: " << s.mean() << " variance: " << s.variance()
                      << std::endl;
         block.set_stats(s.mean(), s.variance(), s_positive.mean());
      }
   }
}

coot::density_box_t
coot::side_chain_densities::get_block(mmdb::Residue *residue_p) const {

   std::map<mmdb::Residue *, density_box_t>::const_iterator it;
   it = density_block_map_cache.find(residue_p); // this cannot (must not) fail - so make
                                                 // sure that fill_residue_blocks is called
                                                 // before this function is called.
   if (it == density_block_map_cache.end()) {
      std::cout << "ERROR:: in get_block(): Hideous failure!" << std::endl;
   }
   return it->second;
}



// the given residue needs to have a CB - caller should check and make one if
// the model doesn't have one
// limit_to_correct_rotamers_only is optional arg, default false, and is for debugging the llr for
// correct solutions - which is bad/low and why?
// verbose_output_mode is optional arg, default true
//
std::map<std::string, std::pair<std::string, double> >
coot::side_chain_densities::get_rotamer_likelihoods(mmdb::Residue *residue_p,
                                                    const clipper::Xmap<float> &xmap,
                                                    bool limit_to_correct_rotamers_only,
                                                    bool verbose_output_mode) {

   auto print_probability_map = [residue_p] (const std::map<std::string, std::pair<std::string, double> > &probability_map) {
      std::map<std::string, std::pair<std::string, double> >::const_iterator it;
      std::cout << "            " << coot::residue_spec_t(residue_p) << " " << residue_p->GetResName() << std::endl;
      for (it=probability_map.begin(); it!=probability_map.end(); ++it) {
         const std::string &key = it->first;
         const std::pair<std::string, double> &value = it->second;
         std::cout << "debug:: get_rotamer_likelihoods(): key " << std::setw(8) << key
                   << " rot " << value.first << " ll " << value.second << "\n";
      }
   };

   // To see what's going on with residue scoring.
   // verbose_output_mode = true;

   // fill_residue_blocks() has been called before we get here

   // key is residue name, value is rotamer name and score pair
   std::map<std::string, std::pair<std::string, double> > bs; // return this, best_score_for_res_type

   if (density_block_map_cache.size() == 0) {
      std::cout << "ERROR:: Cache is empty - fill it first" << std::endl;
      return bs;
   }

   // are axes needed?
   std::pair<clipper::Coord_orth, std::vector<clipper::Coord_orth> > cb_pos_and_axes =
      get_residue_axes(residue_p);
   const std::vector<clipper::Coord_orth> &axes = cb_pos_and_axes.second;
   if (! axes.empty()) {
      // sample_masked density around CB

      // retrieve the block for this residue (filled in fill_residue_blocks())

      density_box_t block = get_block(residue_p);

      if (verbose_output_mode)
         std::cout << "debug:: in get_rotamer_likelihoods() residue: "
                   << residue_spec_t(residue_p) << " block has mean " << block.mean
                   << " and var " << block.var << " sd " << sqrt(block.var) << std::endl;

      if (block.empty()) {

         std::cout << "WARNING:: failed to get a density block for "
                   << residue_spec_t(residue_p) << std::endl;

      } else {

         if (block.is_weird) {
            // do nothing
            std::cout << "WARNING:: weird density block for "
                      << residue_spec_t(residue_p) << std::endl;
         } else {

            // Happy Path

            // compare the block of density to every rotamer (of every residue)

            std::string res_name = residue_p->GetResName();
            std::string rot_name = get_rotamer_name(residue_p);

            if (verbose_output_mode)
               std::cout << "debug:: in get_rotamer_likelihoods(): residue "
                         << residue_spec_t(residue_p) << " " << res_name << " " << rot_name << std::endl;

            std::pair<bool, std::vector<std::pair<std::string, std::string> > > rotamer_limits;
            rotamer_limits.first = false;

            // just test the rotamer that this residue is, combine with
            // outputting the likelihoods for the grids in get_log_likelihood_ratio()
            // in the "check the engine" ouput.

            if (limit_to_correct_rotamers_only)
               rotamer_limits.first = true; // testing the likelihoods of the correct rotamers

            rotamer_limits.second.push_back(std::pair<std::string, std::string> (res_name, rot_name));

            std::map<std::string, std::pair<std::string, double> > probability_map =
               compare_block_vs_all_rotamers(block, residue_p, data_dir, rotamer_limits, xmap);

            if (false)
               print_probability_map(probability_map);

            // find the max value so that I can mark it and close others for screen output
            double best_score = -999999999999999.9;
            std::string best_rotamer;
            std::map<std::string, std::pair<std::string, double> >::const_iterator it;
            for (it=probability_map.begin(); it!=probability_map.end(); ++it) {
               if (it->second.second > best_score) {
                  best_score = it->second.second;
                  best_rotamer = it->second.first;
               }
            }

            std::map<std::string, std::pair<std::string, double> > best_score_for_res_type;
            // std::cout << "here with probability_map size " << probability_map.size() << std::endl;
            for (it=probability_map.begin(); it!=probability_map.end(); ++it) {
               const std::string &key = it->first; // resname:rot_name
               // std::cout << "key: " << key << std::endl;
               std::pair<std::string, std::string> rrp = map_key_to_residue_and_rotamer_names(key);
               const std::string &res_name_l = rrp.first;
               // const std::string &rot_name = rrp.second;
               const double &score = it->second.second;
               const std::string &table_rot_name = rrp.second;
               if (best_score_for_res_type.find(res_name_l) == best_score_for_res_type.end()) {
                  best_score_for_res_type[res_name_l] = std::make_pair(table_rot_name, score);
                  if (false)
                     std::cout << "DEBUG:: first  score for res-type " << std::setw(10) << std::left << key << " "
                               << std::fixed << std::right << std::setw(7) << std::setprecision(2)
                               << table_rot_name << " " << score << std::endl;
               } else {
                  if (score > best_score_for_res_type[res_name_l].second) {
                     best_score_for_res_type[res_name_l] = std::make_pair(table_rot_name, score);
                     if (false)
                        std::cout << "DEBUG:: better score for res-type " << std::setw(10) << std::left << key << " "
                                  << std::fixed << std::right << std::setw(7) << std::setprecision(2)
                                  << table_rot_name << " " << score << std::endl;
                  }
               }
            }
            bs = best_score_for_res_type;

            if (verbose_output_mode) {
               for (it=best_score_for_res_type.begin(); it!=best_score_for_res_type.end(); ++it) {
                  const std::string &res_type = it->first;
                  const double &score = it->second.second;
                  std::string m;
                  if (best_score < 0) {
                     if (score > 1.1 * best_score) m = " ooo";
                  } else {
                     if (score > 0)
                        if (score > 0.5 * best_score)
                           m = " ooo";
                  }
                  if (score == best_score) m = " ***";
                  std::cout << "   " << res_type << " "
                            << std::fixed << std::right << std::setw(8) << std::setprecision(2)
                            << score << m << std::endl;
               }
            }
         }
      }
   }

   return bs;
}


// key is residue type and value is rotamer name and score pair
//
std::map<std::string, std::pair<std::string, double> >
coot::side_chain_densities::compare_block_vs_all_rotamers(density_box_t block,
                                                          mmdb::Residue *residue_p, // for debugging
                                                          const std::string &data_dir,
                                                          const std::pair<bool, std::vector<std::pair<std::string, std::string> > > &rotamer_limits,
                                                          const clipper::Xmap<float> &xmap) {

   std::map<std::string, std::pair<std::string, double> > probability_map;

   std::string glob_pattern = "*";
   std::vector<std::string> dirs = coot::util::glob_files(data_dir, glob_pattern);
   // std::cout << "found " << dirs.size() << " files in " << data_dir << std::endl;

   for (std::size_t idir=0; idir<dirs.size(); idir++) {
      const std::string &res_dir = dirs[idir];

      std::vector<std::string> rot_dirs = coot::util::glob_files(res_dir, glob_pattern);
      for (std::size_t jdir=0; jdir<rot_dirs.size(); jdir++) {
         const std::string &rot_dir = rot_dirs[jdir];

         std::string res = util::file_name_non_directory(res_dir);
         std::string rot = util::file_name_non_directory(rot_dir);
         std::string key = res + ":" + rot;

         bool do_it = false;

         if (false)
            std::cout << "debug:: in compare_block_vs_all_rotamers() rotamer limits first "
                      << rotamer_limits.first << std::endl;
         if (! rotamer_limits.first) {
            // so, don't apply the limits
            do_it = true;
         } else {
            if (rotamer_limits.first) {
               // testing path
               if (! rotamer_limits.second.empty()) {
                  bool found = false;
                  for (std::size_t i=0; i<rotamer_limits.second.size(); i++) {
                     std::string limit_key = rotamer_limits.second[i].first + ":" + rotamer_limits.second[i].second;
                     if (limit_key == key) {
                        found = true;
                        break;
                     }
                  }
                  if (found)
                     do_it = true;
               }
            }
         }

         if (do_it)
            if (rot == "none")
               if (res != "GLY")
                  if (res != "ALA")
                     do_it = false;

         if (do_it) {

            if (false) // interesting but TMI.
               std::cout << "debug:: in compare_block_vs_all_rotamers(): rot_dir: " << rot_dir << std::endl;

            // std::cout << "debug:: in compare_block_vs_all_rotamers(): block var " << block.var << std::endl;
            std::pair<bool, double> p = compare_block_vs_rotamer(block, residue_p, rot_dir, xmap);
            if (p.first) {
               // std::cout << "debug:: in compare_block_vs_rotamer() pr: " << key << " " << p.second << std::endl;
               probability_map[key] = std::make_pair(res, p.second);
            }
         } else {
            // too noisy
            if (false)
               std::cout << "debug:: in compare_block_vs_rotamer() do_it was false " << rot_dir << std::endl;
         }
      }
   }
   return probability_map;
}

bool
coot::side_chain_densities::get_test_map_is_above_model_mean(const unsigned int &grid_idx,
                                                             const density_box_t &block,
                                                             const double &mean) const {
   const double &x = block[grid_idx];
   return (x > mean);

}

double
coot::side_chain_densities::get_log_likelihood(const unsigned int &grid_idx,
                                               const density_box_t &block,
                                               const double &mean,
                                               const double &variance,
                                               const double &skew) const {

   double x = block[grid_idx];
   double a = x - mean;
   double c_part = log(sqrt(1.0/(2.0 * M_PI * variance)));
   double e_part = -0.5*a*a/variance;

   if (false)
      std::cout << "get_log_likelihood() x " << x
                << " grid-idx " << grid_idx << " nz " << a/sqrt(variance) << std::endl;

   return c_part + e_part;
}

double
coot::side_chain_densities::get_grid_point_distance_from_grid_centre(const unsigned int &idx,
                                                                     const double &step_size) const {

   std::tuple<int, int, int> grid_coord = grid_index_to_grid(idx);
   // this is not the same axes system as the grid points - but they are on the same scale
   clipper::Coord_orth pt_in_grid(std::get<0>(grid_coord) * step_size,
                                  std::get<1>(grid_coord) * step_size,
                                  std::get<2>(grid_coord) * step_size);
   double l = std::sqrt(pt_in_grid.lengthsq());
   return l;
}

// log_likelihood ratio vs the gaussian sphere null hypothesis
double
coot::side_chain_densities::get_log_likelihood_ratio(const unsigned int &grid_idx,
                                                     const density_box_t &block,
                                                     const std::string &rotamer_dir,
                                                     const double &step_size,
                                                     const double &mean,
                                                     const double &variance_in,
                                                     const double &skew) const {

   double density_val = block[grid_idx];
   if (density_val > mn_density_block_sample_x_max)
      density_val = mn_density_block_sample_x_max;

   double variance = variance_in; // test/hack
   variance = 0.11; // observed variances are not useful.

   double var_scale = variance/block.var;
   double sd_scale = sqrt(var_scale);
   double mean_offset = mean - block.mean;

   if (false) {
      std::cout << "debug:: variance " << variance << " block.var " << block.var
                << " sd_scale " << sd_scale << std::endl;
      std::cout << "debug:: scaling density_val " << density_val << " with " << sd_scale
                << " and mean offset " << mean_offset << std::endl;
   }
   // double x = density_val * sd_scale - mean_offset;
   double x = density_val * 1.0; // weird scaling makes LLR better.

   double z = x - mean;
   double c_part = log(sqrt(1.0/(2.0 * M_PI * variance)));
   double e_part = -0.5*z*z/variance;

   // null hypothesis

   double nhs = null_hypothesis_scale;
   nhs = 2.0; // why is this 2.0?

   // distance between grid point and the CB
   double d = get_grid_point_distance_from_grid_centre(grid_idx, step_size);
   double x0_null_hypothesis = d;
   double c_part_null_normal = 1.0/(sqrt(2.0 * M_PI * null_hypothesis_sigma * null_hypothesis_sigma));
   double z_null_normal = d;
   double e_part_null_normal = exp(-((z_null_normal*z_null_normal)/(2.0 * null_hypothesis_sigma * null_hypothesis_sigma)));
   double x0_fake_density = nhs * c_part_null_normal * e_part_null_normal;

   double z_null = x0_fake_density - mean; // z value for the null hypothesis density value
   double e_part_normal = -0.5 * z_null * z_null / variance;

   if (false) {
      std::cout << "get_log_likelihood_ratio() x " << x
                << " grid-idx " << grid_idx << " nz " << z/sqrt(variance) << std::endl;
      std::cout << "in get_log_likelihood_ratio() null-hyp scale sigma " << null_hypothesis_scale
                << " " << null_hypothesis_sigma << std::endl;
      std::cout << "in get_log_likelihood_ratio() d " << d << std::endl;
      std::cout << "in get_log_likelihood_ratio() A " << e_part << " " << e_part_normal << std::endl;
      std::cout << "in get_log_likelihood_ratio() B " << z_null << " " << x0_fake_density << std::endl;
      std::cout << "in get_log_likelihood_ratio() C " << c_part_null_normal << " " << e_part_null_normal
                << std::endl;
   }
   double w = 1.0;
   // w = 2.3/(d + 1.0);
   // w = 1.0 - d * 0.166;
   // w = w * w * w;
   // w = 1.0;
   double diff = w * (e_part - e_part_normal);

   // std::cout << "for " << rotamer_dir << " diff " << diff << std::endl;

   // remove this hideous baddies: Magic number - needs optimizing
   double mn_log_likelihood_ratio_difference_max = 18.0;
   if (diff < mn_log_likelihood_ratio_difference_min)
      diff = mn_log_likelihood_ratio_difference_min;
   if (diff > mn_log_likelihood_ratio_difference_max)
      diff = mn_log_likelihood_ratio_difference_max;

   if (false) // debug/check the engine
       /*
                << " e_part: " << std::setw(10) << e_part
                << " e_part_normal: " << std::setw(8) << e_part_normal
      */
      std::cout << "engine: idx: " << grid_idx
                << " with density_val " << std::setw(8) << std::right << std::setprecision(5) << density_val
                << " dv-mean " << std::right << std::setw(8) << std::setprecision(5) << z
                << " x0_fake_density " << std::setw(8) << x0_fake_density
                << " mean " << std::setw(5) << mean << " sigma " << sqrt(variance)
                << " return " << std::fixed << std::right << std::setprecision(6) << diff
                << "\n";

   return diff;
}

std::tuple<int, int, int>
coot::side_chain_densities::grid_index_to_grid(int grid_idx) const {

   int n_z = 0;
   int n_y = 0;
   int n_x = 0;
   int n = 2 * n_steps + 1;

   while (grid_idx >= (n * n)) {
      grid_idx -= n * n;
      n_x++;
   }
   while (grid_idx >= n) {
      grid_idx -= n;
      n_y++;
   }
   n_z = grid_idx;

   std::tuple<int, int, int> t(n_x - n_steps , n_y - n_steps, n_z - n_steps);
   return t;

}

// Manhattan test - not very useful
bool
coot::side_chain_densities::in_sphere(int grid_idx, const int &n_steps) const {

   bool inside = true;

   std::tuple<int, int, int> t = grid_index_to_grid(grid_idx);
   int n_x = std::get<0>(t);
   int n_y = std::get<1>(t);
   int n_z = std::get<2>(t);
   int delta = abs(n_x) + abs(n_y) + abs(n_z);
   if (delta > n_steps)
      inside = false;

   return inside;

}

bool
coot::side_chain_densities::in_sphere(const clipper::Coord_orth &pt,
                                      const clipper::Coord_orth &cb,
                                      const double &d_max) const {
   return ((pt-cb).lengthsq() < (d_max * d_max));
}

// We fill the rotamer grid cache, so not const.
// Maybe there should be 2 separate functions
//
std::pair<bool, double>
coot::side_chain_densities::compare_block_vs_rotamer(density_box_t block,
                                                     mmdb::Residue *residue_p,
                                                     const std::string &rotamer_dir,
                                                     const clipper::Xmap<float> &xmap) {

   // This function (-is-) was slow - conversion of strings to numbers
   // So, to fix that, in the calling function (compare_block_vs_all_rotamers()),
   // read in all the stats files - write and read the stats file as binaries
   // if it's still slow, so that this function is not needed.

   // I want to match the density mean and sigma of this blocks mean and sigma with the
   // overall stats with a scale (to match sigmas) and offset (to match means)

   if (false)
      std::cout << "debug:: compare_block_vs_rotamer() residue "
                << coot::residue_spec_t(residue_p) << " " << rotamer_dir
                << std::endl;

   // bool do_debug_scoring = false; // count the number of times the the sample map is more than
                                  // the template stats mean. Should be about 50-50?
                                  // Turns out not. 100-150 seems to work.

   bool success = false; // initially
   double sum_log_likelihood = 0.0;
   double step_size = grid_box_radius/static_cast<float>(n_steps);

   get_results_addition_lock();

   // std::cout << "------- calling get_log_likelihood_ratio() for rotamer_dir " << rotamer_dir << std::endl;

   auto fill_rotamer_dir_grid_stats_map_cache = [rotamer_dir] (std::map<std::string, std::map<unsigned int, std::tuple<double, double, double> > > &rotamer_dir_grid_stats_map_cache) {
                                                   std::string glob_pattern = "stats.table";
                                                   std::vector<std::string> tables = coot::util::glob_files(rotamer_dir, glob_pattern);
                                                   if (tables.size() == 1) {
                                                      std::map<unsigned int, std::tuple<double, double, double> > stats_map;
                                                      std::string stats_table_file_name = tables[0];

                                                      // std::cout << "stats_table_file_name: " << stats_table_file_name << std::endl;
                                                      std::ifstream f(stats_table_file_name.c_str());
                                                      if (f) {
                                                         std::string line;
                                                         while (std::getline(f, line)) {
                                                            std::vector<std::string> words = coot::util::split_string_no_blanks(line);
                                                            if (words.size() == 4) { // 5 with kurtosis
                                                               unsigned int grid_idx = util::string_to_int(words[0]);
                                                               double mean = util::string_to_double(words[1]);
                                                               double var  = util::string_to_double(words[2]);
                                                               double skew = util::string_to_double(words[3]);
                                                               std::tuple<double, double, double> t(mean, var, skew);
                                                               stats_map[grid_idx] = t;
                                                            }
                                                         }
                                                         rotamer_dir_grid_stats_map_cache[rotamer_dir] = stats_map;
                                                      }
                                                   }
                                                };

   std::map<std::string, std::map<unsigned int, std::tuple<double, double, double> > >::const_iterator it;
   it = rotamer_dir_grid_stats_map_cache.find(rotamer_dir);
   if (it == rotamer_dir_grid_stats_map_cache.end())
      fill_rotamer_dir_grid_stats_map_cache(rotamer_dir_grid_stats_map_cache);
   if (it == rotamer_dir_grid_stats_map_cache.end())
      it = rotamer_dir_grid_stats_map_cache.find(rotamer_dir);
   if (it != rotamer_dir_grid_stats_map_cache.end()) {
      success = true;
      const std::map<unsigned int, std::tuple<double, double, double> > &stats_map = it->second;
      std::map<unsigned int, std::tuple<double, double, double> >::const_iterator it_inner;
      for (it_inner=stats_map.begin(); it_inner!=stats_map.end(); ++it_inner) {
         const unsigned int &grid_idx = it_inner->first;
         const std::tuple<double, double, double> &m_v_s = it_inner->second;
         const double &mean = std::get<0>(m_v_s);
         const double &var  = std::get<1>(m_v_s);
         const double &skew = std::get<2>(m_v_s);
         double llr = get_log_likelihood_ratio(grid_idx, block, rotamer_dir, step_size, mean, var, skew);
         sum_log_likelihood += llr;
      }
   }

   release_results_addition_lock();

   return std::pair<bool, double>(success, sum_log_likelihood);
}

void
coot::side_chain_densities::check_useable_grid_points(mmdb::Residue *residue_p,
                                                      const std::string &useable_grid_points_mapped_to_residue_file_name) const {

   int n_per_side = n_steps * 2 + 1;
   float step_size = grid_box_radius/static_cast<float>(n_steps);
   std::pair<clipper::Coord_orth, std::vector<clipper::Coord_orth> > cb_pos_and_axes =
      get_residue_axes(residue_p);
   const std::vector<clipper::Coord_orth> &axes = cb_pos_and_axes.second;
   const clipper::Coord_orth &cb_pt = cb_pos_and_axes.first;

   std::ofstream f(useable_grid_points_mapped_to_residue_file_name.c_str());
   if (f) {
      if (! axes.empty()) {
         for (int ix= -n_steps; ix<=n_steps; ix++) {
            for (int iy= -n_steps; iy<=n_steps; iy++) {
               for (int iz= -n_steps; iz<=n_steps; iz++) {
                  int idx =
                     (ix + n_steps) * n_per_side * n_per_side +
                     (iy + n_steps) * n_per_side +
                     (iz + n_steps);
                  if (useable_grid_points.find(idx) != useable_grid_points.end()) {
                     clipper::Coord_orth pt_in_grid = make_pt_in_grid(ix, iy, iz, step_size, axes);
                     clipper::Coord_orth pt_grid_point = cb_pt + pt_in_grid;
                     if (useable_grid_points.find(idx) != useable_grid_points.end()) {
                        f << "check-useable-grid-points " << idx << " "
                          << pt_grid_point.x() << " "
                          << pt_grid_point.y() << " "
                          << pt_grid_point.z() << "\n";
                     }
                  }
               }
            }
         }
      }
   }
   f.close();
}


void
coot::side_chain_densities::check_stats(mmdb::Residue *residue_p,
                                        const std::string &res_name,
                                        const std::string &rot_name) const {

   if (useable_grid_points.size() == 0) {
      std::cout << "ERROR:: useable_grid_points size is 0 " << std::endl;
      return;
   }

   int n_per_side = n_steps * 2 + 1;
   float step_size = grid_box_radius/static_cast<float>(n_steps);
   std::pair<clipper::Coord_orth, std::vector<clipper::Coord_orth> > cb_pos_and_axes =
      get_residue_axes(residue_p);
   const std::vector<clipper::Coord_orth> &axes = cb_pos_and_axes.second;
   const clipper::Coord_orth &cb_pt = cb_pos_and_axes.first;
   if (! axes.empty()) {
      std::string dir = "side-chain-data"; // here use a class data item (which is set on construction)
      std::string rot_dir = dir + "/" + res_name + "/" + rot_name;
      std::string file_name = rot_dir + "/" + "stats.table";
      std::ifstream f(file_name.c_str());
      if (f) {
         std::cout << "DEBUG:: check_stats() reading " << file_name << std::endl;
         std::string line;
         std::map<unsigned int, std::tuple<double, double, double> > grid_map;
         while (std::getline(f, line)) {
            std::vector<std::string> words = coot::util::split_string_no_blanks(line);
            if (words.size() == 4) { // 5 with kurtosis
               unsigned int idx = util::string_to_int(words[0]);
               double mean = util::string_to_double(words[1]);
               double var  = util::string_to_double(words[2]);
               double skew = util::string_to_double(words[3]);
               std::tuple<double, double, double> t(mean, var, skew);
               grid_map[idx] = t;
            }
         }

         for (int ix= -n_steps; ix<=n_steps; ix++) {
            for (int iy= -n_steps; iy<=n_steps; iy++) {
               for (int iz= -n_steps; iz<=n_steps; iz++) {
                  int idx =
                     (ix + n_steps) * n_per_side * n_per_side +
                     (iy + n_steps) * n_per_side +
                     (iz + n_steps);
                  if (useable_grid_points.find(idx) != useable_grid_points.end()) {
                     clipper::Coord_orth pt_in_grid = make_pt_in_grid(ix, iy, iz, step_size, axes);
                     clipper::Coord_orth pt_grid_point = cb_pt + pt_in_grid;

                     std::cout << "check-stats "
                               << idx << " "
                               << pt_grid_point.x() << " "
                               << pt_grid_point.y() << " "
                               << pt_grid_point.z() << " "
                               << "mean "   << std::get<0>(grid_map[idx])
                               << " stddev " << std::get<1>(grid_map[idx])
                               << std::endl;
                  }
               }
            }
         }
      } else {
         std::cout << "WARNING:: check_stats() file not found: " << file_name << std::endl;
      }
   } else {
      std::cout << "WARNING:: check_stats() empty axes" << std::endl;
   }
}


// static
void
coot::side_chain_densities::combine_directory(const std::string &rot_dir, int n_steps,
                                              double mn_unreliable_minimum_counts,
                                              double mn_unreliable_minimum_counts_for_low_variance,
                                              double mn_unreliable_minimum_variance,
                                              double mn_use_this_variance_for_unreliable) {

   // work out the mean and standard deviation of the grid point
   // for this rotamer (of this residue)

   // first check that the files are all standard format

   unsigned int n_per_side = 2 * n_steps + 1;
   unsigned int n_target = n_per_side * n_per_side * n_per_side;

   std::string glob_pattern = "*.tab";
   std::vector<std::string> files = coot::util::glob_files(rot_dir, glob_pattern);

   // first check that every file has the same number of numbers
   std::map<std::string, unsigned int> number_count_map;
   for (std::size_t i=0; i<files.size(); i++) {
      const std::string &fn = files[i];
      std::ifstream f(fn.c_str());
      if (f) {
         unsigned int word_count = 0;
         std::string line;
         while (std::getline(f, line)) {
            std::vector<std::string> words = coot::util::split_string_no_blanks(line);
            word_count += words.size();
         }
         number_count_map[fn] = word_count;
      } else {
         std::cout << "failed to open " << fn << std::endl;
      }
   }

   std::map<std::string, unsigned int>::const_iterator it;
   for (it=number_count_map.begin(); it!=number_count_map.end(); ++it) {
      // std::cout << "counts " << it->first << " " << it->second << std::endl;
      if (it->second != n_target) {
         std::cout << "combine_directory() fail " << rot_dir << " " << it->first
                   << " " << it->second << " c.f. " << n_target << std::endl;
         return;
      }
   }


   // OK, files were good

   // std::cout << "debug:: files were good for " << rot_dir << std::endl;

   // so that we can organize the memory a bit easier - on the stack
   const int nps = n_steps * 2 + 1;
   const int nnn = nps * nps * nps;
   std::vector<float> *x = new std::vector<float>[nnn]; // deleted
   for (int i=0; i<nnn; i++) x[i].reserve(40);

   for (std::size_t i=0; i<files.size(); i++) {
      const std::string &fn = files[i];
      std::ifstream ff(fn.c_str());
      if (ff) {
         unsigned int word_count = 0;
         std::string line;
         while (std::getline(ff, line)) {
            std::vector<std::string> words = util::split_string_no_blanks(line);
            for (std::size_t ii=0; ii<words.size(); ii++) {
               const std::string &w = words[ii];
               try {
                  float f = util::string_to_float(w);
                  if (f > -1000.0) { // mask value
                     x[word_count].push_back(f);
                  }
               }
               catch (const std::runtime_error &rte) {
                  std::cout << "ERROR:: rte: combine_directory() " << fn << " " << rte.what() << std::endl;
               }
               word_count++;
            }
         }
         // is fine
         // std::cout << "ending word count "<< word_count << " " << fn << std::endl;
      }
   }

   if (true) {

      stats::single s_all;

      for (int i=0; i<nnn; i++) {

         // stats
         {
            if (x[i].size() > 1) {
               stats::single s;
               for (std::size_t j=0; j<x[i].size(); j++) {
                  const float &val = x[i][j];
                  s.add(val);
                  s_all.add(val);
               }

               std::string fn_stats =  + "stats.table";
               fn_stats = util::append_dir_file(rot_dir, fn_stats);
               std::ios_base::openmode mode = std::ios_base::app;
               std::ofstream f_stats(fn_stats.c_str(), mode);
               double mean = s.mean();
               double var  = s.variance();
               double skew = s.skew();
               std::pair<double, double> mi = s.median_and_iqr();
               double median = mi.first;
               double irq    = mi.second;
               double lim_low  = median - 1.75 * irq;
               double lim_high = median + 1.75 * irq;

               // remove values that are more than 3 iqrs from the mean
               stats::single s_filtered;
               for (std::size_t j=0; j<x[i].size(); j++) {
                  const float &val = x[i][j];
                  if (val < lim_high)
                     if (val > lim_low)
                        s_filtered.add(val);
               }

               bool unreliable = false;
               unsigned int n_filtered = s.size() - s_filtered.size();
               std::cout << "debug:: n_filtered " << n_filtered << std::endl;
               if (s_filtered.size() < mn_unreliable_minimum_counts) // magic number 
                  unreliable = true;
               if (s_filtered.size() < mn_unreliable_minimum_counts_for_low_variance &&  // magic numbers
                   var < mn_unreliable_minimum_counts_for_low_variance)
                  unreliable = true;

               if (unreliable)
                  var = mn_use_this_variance_for_unreliable; // magic number

               if (f_stats) {
                  f_stats << i << " " << s_filtered.mean() << " " << s_filtered.variance()
                          << " " << s_filtered.skew() << "\n";
               } else {
                  std::cout << "failed to open " << fn_stats << std::endl;
               }
               f_stats.close();

            } else {

               if (x[i].size() == 1) {
                  std::cout << "only 1 point " << rot_dir << " " << i << std::endl;
                  double mean = x[i][0];
                  double var  = 1.0;
                  double skew = 0.0;
                  std::string fn_stats =  + "stats.table";
                  fn_stats = util::append_dir_file(rot_dir, fn_stats);
                  std::ios_base::openmode mode = std::ios_base::app;
                  std::ofstream f_stats(fn_stats.c_str(), mode);
                  f_stats << i << " " << mean << " " << var << " " << skew << "\n";

               } else {
                  // std::cout << "No grid point data for " << rot_dir << " grid point " << i << std::endl;
               }
            }
         }

         // file output
         bool file_output = true; // maybe not... (useful for R, not otherwise and makes many files)
         file_output = false; // Not at the moment, then.
         if (file_output) {
            std::string fn = "grid-point-" + util::int_to_string(i) + ".data";
            fn = util::append_dir_file(rot_dir, fn);
            if (x[i].size() > 0) {
               std::ofstream f(fn.c_str());
               if (f) {
                  for (std::size_t j=0; j<x[i].size(); j++)
                     f << x[i][j] << "\n";
                  f.close();
               } else {
                  std::cout << "Failed to open " << fn << " for writing " << std::endl;
               }
            } else {
               // std::cout << "No grid point data for " << rot_dir << " grid point " << i << std::endl;
            }
         }

         bool screen_output = false;
         if (screen_output) {
            std::cout << "sample point of nnn: " << i << " " << x[i].size() << " values" << std::endl;
            for (std::size_t j=0; j<x[i].size(); j++) {
               std::cout << x[i][j] << " ";
            }
            std::cout << std::endl;
         }
      }

      // all grid points, all sample data:
      std::string fn_summary_stats = "summary-stats";
      fn_summary_stats = util::append_dir_file(rot_dir, fn_summary_stats);
      // std::ios_base::openmode mode = std::ios_base::fwrite;
      std::ofstream f_stats(fn_summary_stats.c_str());
      double mean_all = s_all.mean();
      double var_all = s_all.variance();

      f_stats << s_all.size() << " " << mean_all << " " << var_all << std::endl;
      f_stats.close();

   }

   delete [] x;

}


void
coot::side_chain_densities::set_magic_number(const std::string &mn_name, double val) {

   get_results_addition_lock();
   if (mn_name == "mn_log_likelihood_ratio_difference_min") mn_log_likelihood_ratio_difference_min = val;
   if (mn_name == "mn_scale_for_normalized_density") mn_scale_for_normalized_density = val;
   if (mn_name == "mn_density_block_sample_x_max") mn_density_block_sample_x_max = val;
   release_results_addition_lock();

}

#include "coot-utils/fragment-container.hh"
#include "utils/split-indices.hh"

// This function seems to crash if I put it into side-chain-densities (which gets added to the
// coot-ligand library)
//

std::vector<std::pair<coot::fragment_container_t::fragment_range_t, std::vector<coot::side_chain_densities::results_t> > >
coot::get_fragment_sequence_scores(mmdb::Manager *mol,
                                   const coot::fasta_multi &fam,
                                   const clipper::Xmap<float> &xmap) {


   // score blocks of sequences (of which there are, say, 21000)
   auto proc_threads = [] (const std::pair<unsigned int, unsigned int> &start_stop_pair,
                           const coot::fasta_multi &fam,
                           const std::vector<mmdb::Residue *> &a_run_of_residues,
                           const clipper::Xmap<float> &xmap,
                           coot::side_chain_densities &scd) { // fill scd

                          for(unsigned int idx=start_stop_pair.first; idx!=start_stop_pair.second; ++idx) {
                             // std::cout << "proc_threads calls scd.test_sequence() with idx" << idx << std::endl;
                             scd.test_sequence(a_run_of_residues, xmap, fam[idx].name, fam[idx].sequence);
                          }
                       };

   std::vector<std::pair<coot::fragment_container_t::fragment_range_t, std::vector<coot::side_chain_densities::results_t> > > results_vec;

   unsigned int n_sequences = fam.size();
   // "analysis" constructor
   // coot::side_chain_densities scd(n_steps, grid_box_radius, useable_grid_points_file_name);
   // scd.set_data_dir("side-chain-data");

   coot::fragment_container_t fc = make_fragments(mol);

   std::cout << "get_fragment_sequence_scores() debug fragments" << std::endl;
   fc.print_fragments();

   std::cout << "INFO:: number of sequences in sequence file: " << n_sequences << std::endl;

   for (const auto &range : fc.ranges) {
      std::cout << "::: new-range" << std::endl;
      auto tp_0 = std::chrono::high_resolution_clock::now();
      coot::side_chain_densities scd;
      std::vector<coot::side_chain_densities::results_t> results;
      std::pair<std::string, std::vector<mmdb::Residue *> > a_run_of_residues =
         scd.setup_test_sequence(mol, range.chain_id, range.start_res.res_no, range.end_res.res_no, xmap);
      auto tp_1 = std::chrono::high_resolution_clock::now();
      if (! a_run_of_residues.first.empty()) {
         std::cout << "WARNING:: Failed to make a run of residue - due to missing atoms" << std::endl;
         std::cout << a_run_of_residues.first << std::endl;
      } else {
         // happy path
         scd.setup_likelihood_of_each_rotamer_at_every_residue(a_run_of_residues.second, xmap);

         auto tp_2 = std::chrono::high_resolution_clock::now();

#if 1 // threaded version
         unsigned int n_threads = get_max_number_of_threads();
         std::vector<std::pair<unsigned int, unsigned int> > seq_index_vector =
            coot::atom_index_ranges(n_sequences, n_threads);
         std::vector<std::thread> threads;

         for (unsigned int i=0; i<seq_index_vector.size(); i++) {
            std::pair<unsigned int, unsigned int> index_pair = seq_index_vector[i];
            threads.push_back(std::thread(proc_threads, index_pair, fam, a_run_of_residues.second, xmap, std::ref(scd)));
         }

         for (unsigned int i=0; i<seq_index_vector.size(); i++)
            threads[i].join();
#endif

#if 0 // the single threaded way

         for (unsigned int idx=0; idx<n_sequences; idx++) {
            std::string sequence = fam[idx].sequence;
            // std::cout << "Input Sequence:\n" << sequence << std::endl;
            const std::string &name = fam[idx].name;
            scd.test_sequence(a_run_of_residues, xmap, name, sequence);
         }
#endif

         std::map<std::string, std::vector<coot::side_chain_densities::results_t> >::const_iterator it;
         bool print_results = false;
         if (print_results) {
            for (it=scd.results_container.begin(); it!=scd.results_container.end(); ++it) {
               const std::string &sequence = it->first;
               std::cout << "sequence: " << sequence << std::endl;
               for (const auto &r : it->second) {
                  std::cout << "   " << r.offset << " " << r.sequence << " " << r.sum_score << std::endl;
               }
            }
         }

         auto tp_3 = std::chrono::high_resolution_clock::now();
         // transfer scd.results_container to returned results
         if (false) {
            for (it=scd.results_container.begin(); it!=scd.results_container.end(); ++it) {
               const std::vector<side_chain_densities::results_t> &v(it->second);
               if (! v.empty())
                  results.insert(results.end(), v.begin(), v.end());
            }
         } else {
            auto tp_3i = std::chrono::high_resolution_clock::now();
            double sum_data = 0.0;
            double sum_data_sqrd = 0.0;
            unsigned int n = 0;
            for (it=scd.results_container.begin(); it!=scd.results_container.end(); ++it) {
               const std::vector<side_chain_densities::results_t> &v(it->second);
               for (unsigned int i=0; i<v.size(); i++) {
                  sum_data += v[i].sum_score;
                  sum_data_sqrd += v[i].sum_score * v[i].sum_score;
               }
               n += v.size();
            }
            if (n > 1000) {
               results.reserve(10000);
               double mean = sum_data/static_cast<double>(n);
               double var = sum_data_sqrd/static_cast<double>(n) - mean * mean;
               if (var < 0) var = 0;
               double sd = std::sqrt(var);
               std::cout << "Mean: " << mean << " sd " << sd << std::endl;
               double lim_good_enough = mean + 1.5 * sd;
               for (it=scd.results_container.begin(); it!=scd.results_container.end(); ++it) {
                  const std::vector<side_chain_densities::results_t> &v(it->second);
                  for (unsigned int i=0; i<v.size(); i++) {
                     if (v[i].sum_score > lim_good_enough)
                        results.push_back(v[i]);
                  }
               }
            } else {
               for (it=scd.results_container.begin(); it!=scd.results_container.end(); ++it) {
                  const std::vector<side_chain_densities::results_t> &v(it->second);
                  if (! v.empty())
                     results.insert(results.end(), v.begin(), v.end());
               }
            }

            auto tp_3j = std::chrono::high_resolution_clock::now();
            auto dij = std::chrono::duration_cast<std::chrono::milliseconds>(tp_3j - tp_3i).count();
            std::cout << "TIMINGS:: get_fragment_sequence_scores() calc mean results " << dij
                      << " milliseconds" << std::endl;
         }

         auto tp_4 = std::chrono::high_resolution_clock::now();
         auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
         auto d21 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_2 - tp_1).count();
         auto d32 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_3 - tp_2).count();
         auto d43 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_4 - tp_3).count();
         std::cout << "TIMINGS:: get_fragment_sequence_scores() setup model: "
                   << d10 << " setup likelihoods: " << d21 << " proc_theads: " << d32 << " consolidate: " << d43
                   << " milliseconds" << std::endl;
      }
      results_vec.push_back(std::make_pair(range, results));
   }

   return results_vec; // results for each range/fragment
}



// do TRPs score higher than a VAL? Is that bad?
void
coot::side_chain_densities::check_vs_flat_density() const {

   std::string pdb_file_name = "tutorial-modern.pdb";
   std::string mtz_file_name = "rnasa-1.8-all_refmac1.mtz";

   std::string res_type = "VAL";
   std::string rot_name = "p";

   
}
