
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-map-utils.hh"
#include "utils/split-indices.hh"

#include "new-residue-by-3-phi-psi.hh"


coot::new_residue_by_3_phi_psi::new_residue_by_3_phi_psi(const std::string &terminus_type_in, mmdb::Residue *residue_p_in, const std::string &chain_id_in) {

   thread_pool_p = 0;
   n_threads = 0;

   chain_id = chain_id_in;
   terminus_type = terminus_type_in;
   residue_p = residue_p_in;

   init_phi_psi_plot();

}

void
coot::new_residue_by_3_phi_psi::init_phi_psi_plot() {

   rama.init(clipper::Ramachandran::All2);

   rama_pro.init(clipper::Ramachandran::Pro2);

   rama_max = 0.0;
   rama_pro_max = 0.0;

   upstream_neighbour_residue_p = 0;
   downstream_neighbour_residue_p = 0;

   for (float phi=0.0; phi<360.0; phi+=3.0) {
      for (float psi=0.0; psi<360.0; psi+=3.0) {
	 float v = rama.probability(clipper::Util::d2rad(phi),
				    clipper::Util::d2rad(psi));
	 if (v > rama_max)
	    rama_max = v;
      }
   }
   for (float phi=0.0; phi<360.0; phi+=3.0) {
      for (float psi=0.0; psi<360.0; psi+=3.0) {
	 float v = rama_pro.probability(clipper::Util::d2rad(phi),
					clipper::Util::d2rad(psi));
	 if (v > rama_pro_max)
	    rama_pro_max = v;
      }
   }
}

void
coot::new_residue_by_3_phi_psi::add_thread_pool(ctpl::thread_pool  *thread_pool_p_in, unsigned int n_threads_in) {
   thread_pool_p = thread_pool_p_in;
   n_threads = n_threads_in;
}

// If there are 3 positions on return their order is N, C, CA
//
coot::new_residue_by_3_phi_psi::connecting_atoms_t
coot::new_residue_by_3_phi_psi::get_connecting_residue_atoms() const {

   connecting_atoms_t atoms_in_residue;

   // Note, it seems that if there are 2 atoms in a residue with the
   // selected name (which happens if there is an altconf), then
   // GetAtom() returns NULL.  Ugh.
   // 20180403-PE But is that something we care about in this class?

   // std::vector<clipper::Coord_orth> pos; old

   mmdb::Atom *N_at = 0;
   mmdb::Atom *C_at = 0;
   mmdb::Atom *CA_at = 0;

   mmdb::PPAtom residue_atoms = 0;
   int nResidueAtoms;
   residue_p->GetAtomTable(residue_atoms, nResidueAtoms);
   for (int i=0; i<nResidueAtoms; i++) {
      mmdb::Atom *at = residue_atoms[i];
      std::string atom_name(at->name);
      if (atom_name == " N  ") N_at = at;
      if (atom_name == " C  ") C_at = at;
      if (atom_name == " CA ") CA_at = at;
   }

   if (N_at) {
      if (C_at) {
	 if (CA_at) {
	    clipper::Coord_orth  N_at_pos = co(N_at);
	    clipper::Coord_orth CA_at_pos = co(CA_at);
	    clipper::Coord_orth  C_at_pos = co(C_at);
	    atoms_in_residue = connecting_atoms_t(N_at_pos, CA_at_pos, C_at_pos);
	 }
      }
   }

   if (upstream_neighbour_residue_p) {
      mmdb::PPAtom residue_atoms_up = 0;
      int nResidueAtoms_up;
      upstream_neighbour_residue_p->GetAtomTable(residue_atoms_up, nResidueAtoms_up);
      for (int i=0; i<nResidueAtoms_up; i++) {
	 mmdb::Atom *at = residue_atoms_up[i];
	 std::string atom_name(at->GetAtomName());
	 if (atom_name == " C  ") {  // PDBv3 FIXME
	    clipper::Coord_orth pos = co(at);
	    atoms_in_residue.set_upstream_C(pos);
	    break;
	 }
      }
   }

   if (downstream_neighbour_residue_p) {
      mmdb::PPAtom residue_atoms_down = 0;
      int nResidueAtoms_down;
      downstream_neighbour_residue_p->GetAtomTable(residue_atoms_down, nResidueAtoms_down);
      for (int i=0; i<nResidueAtoms_down; i++) {
	 mmdb::Atom *at = residue_atoms_down[i];
	 std::string atom_name(at->GetAtomName());
	 if (atom_name == " N  ") {  // PDBv3 FIXME
	    clipper::Coord_orth pos = co(at);
	    atoms_in_residue.set_downstream_N(pos);
	    break;
	 }
      }
   }
   return atoms_in_residue;
}

#include<random>
std::random_device rd;

// between -1 and 1
float coot::get_random_float_rd() {

   // this function is slow when used with threads.

   int num = rd();
   float sf = 4.66e-10;
   float f = sf * static_cast<float>(num);
   return f;
}


float coot::get_random_float_mt(dsfmt_t *dsfmt) {
   double d = dsfmt_genrand_close_open(dsfmt);
   // std::cout << "in get_random_float_mt() d " << d << " from dsfmt " << dsfmt << std::endl;
   float f = static_cast<float>(d);
   return f; // between 0 and 1
}

float coot::get_random_float() {
   return get_random_float_rd();
}



// Xoroshiro128+
uint64_t shuffle_table[4];

void init_xoroshiro128plus() {
   for(unsigned int i=0; i<4; i++) shuffle_table[i] = i;
}

uint64_t xoroshiro128plus_next() {
    uint64_t s1 = shuffle_table[0];
    uint64_t s0 = shuffle_table[1];
    uint64_t result = s0 + s1;
    shuffle_table[0] = s0;
    s1 ^= s1 << 23;
    shuffle_table[1] = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5);
    return result;
}

float get_random_float() {

   uint64_t ii = xoroshiro128plus_next();
   float sf = 4.66e-10;
   float f = sf * static_cast<float>(ii);
   return f;

}

// static beacuse used in threads.
coot::phi_psi_t
coot::new_residue_by_3_phi_psi::get_phi_psi_by_random(const clipper::Ramachandran &rama_local,
                                                      const float &rama_max_local,
                                                      bool is_pro_rama,
                                                      dsfmt_t *dsfmt) {
   float phi = 0.0f;
   float psi = 0.0f;
   // bool force_positive_phi = false; // not sure if this will be needed.
   for (;;) {
      phi = 360.0f * fabs(get_random_float_mt(dsfmt));
      psi = 360.0f * fabs(get_random_float_mt(dsfmt));

      float r = rama_max_local * get_random_float_mt(dsfmt);
      float prob = rama_local.probability(clipper::Util::d2rad(phi),
					  clipper::Util::d2rad(psi));

      // std::cout << "            compare " << prob << " " << r << std::endl;

      if (prob > r)
         break;
   }
   float u1 = get_random_float_mt(dsfmt);
   // std::cout << "debug get_phi_psi_by_random(): u1 " << u1 << std::endl;
   float u2 = 2.0f *  u1 - 1.0; // -1 -> 1.
   float tau = 111.2 + u2 * 6.2; // 6.2 needs investigation
   return phi_psi_t(phi, psi, tau);
}

// phi and psi in radians
// When we are building backwards, we have psi for the selected residue, but not phi.
// To be clear, we only have psi, when we have also have the downstream (higher residue number) residue
// as well.
// static
double
coot::new_residue_by_3_phi_psi::get_phi_by_random_given_psi(double psi,
                                                            const clipper::Ramachandran &rama,
                                                            dsfmt_t *dsfmt) {

   double phi;

   // find rand_max_local (for probability of phi given psi
   //
   double step = M_PI / 36.0;
   double conditional_pr_rama_max = 0.0;
   for (unsigned int i=0; i<72; i++) {
      double phi_i = (i+0.5)*step;
      double pr = rama.probability(phi_i, psi);
      if (pr > conditional_pr_rama_max)
	 conditional_pr_rama_max = pr;
   }
   
   for (;;) {
      phi = 2.0 * M_PI * fabsf(get_random_float_mt(dsfmt));
      double r = conditional_pr_rama_max * fabsf(get_random_float_mt(dsfmt));
      double prob = rama.probability(phi, psi);
      if (prob > r) {
	 break;
      }
   }
   return phi;
}

// phi and psi in radians
// When we are building forward, we have phi for the selected residue, but not psi.
// To be clear, we only have phi, when we have also have the upstream (lower residue number) residue
// as well.
// static
double
coot::new_residue_by_3_phi_psi::get_psi_by_random_given_phi(double phi, const clipper::Ramachandran &rama,
                                                            dsfmt_t *dsfmt) {

   double psi;

   // find rand_max_local (for probability of phi given psi
   //
   double step = M_PI / 36.0;
   double conditional_pr_rama_max = 0.0;
   for (unsigned int i=0; i<72; i++) {
      double psi_i = (i+0.5)*step;
      double pr = rama.probability(phi, psi_i);
      if (pr > conditional_pr_rama_max)
	 conditional_pr_rama_max = pr;
   }

   if (conditional_pr_rama_max < 0.0001) {
      // something went wrong, hack a return value
      // so that we don't stay in the below loop forever
      //
      psi = 2.0 * M_PI * fabsf(get_random_float_mt(dsfmt));
   } else{

      for (;;) {
	 psi = 2.0 * M_PI * fabsf(get_random_float_mt(dsfmt));
	 double r = conditional_pr_rama_max * fabsf(get_random_float_mt(dsfmt));
	 double prob = rama.probability(phi, psi);
	 if (prob > r) {
	    break;
	 }
      }
   }
   return psi;
}

// static
float
coot::new_residue_by_3_phi_psi::score_fragment_basic(const minimol::fragment &frag,
                                                     const connecting_atoms_t &current_res_pos,
                                                     const clipper::Xmap<float> &xmap) {

   float score = 0.0;
   float w_sum = 0.0;

   for(int ires=frag.min_res_no(); ires<=frag.max_residue_number(); ires++) {
      for (unsigned int iat=0; iat<frag[ires].atoms.size(); iat++) {
         const clipper::Coord_orth atom_pos(frag[ires][iat].pos);
         float w = frag[ires][iat].occupancy;
         float f = util::density_at_point(xmap, atom_pos) * w;
         score += f; // or f * f * f  so that fragments that have atoms in low or negative density score hideously
                     // but we also need a roll-off - a f^5 say that stops fitting into sulphurs (for x-ray fitting)
         w_sum += w;
      }
   }
   return score/w_sum;
}

#include <fstream>

// static
float
coot::new_residue_by_3_phi_psi::score_fragment_using_peptide_fingerprint(const minimol::fragment &frag,
                                                                         const connecting_atoms_t &current_res_pos,
                                                                         const clipper::Xmap<float> &xmap,
                                                                         int res_no_base, int i_trial) {

   float score = 0.0;
   float w_sum = 0.0;

   auto roll_off = [] (float f) {
                      return f;
                   };


   // same again for the "dont_add" protection
   std::atomic<bool> print_lock(false);

   auto get_print_lock = [&print_lock] () {
                            bool unlocked = false;
                            while (! print_lock.compare_exchange_weak(unlocked, true)) {
                               std::this_thread::sleep_for(std::chrono::milliseconds(10)); // or faster?  (was 100ns)
                               unlocked = false;
                            }
                         };

   auto release_print_lock = [&print_lock] () {
                                print_lock = false;
                             };
   
   auto fingerprint_score = [&xmap, roll_off, i_trial, res_no_base] (const clipper::Coord_orth &ca_pos_1,
                                               const clipper::Coord_orth &ca_pos_2,
                                               const clipper::Coord_orth &o_pos_1,
                                               int ires) { // pass ires for debugging (writing out positions)

                               // these need to be optimized somehow - sounds fun
                               float scale_CO       =  0.2; // we are already counting the O position in the normal scoring
                               float scale_CO_low   = -0.8;
                               float scale_CO_anti  = -0.3;
                               float scale_N        =  0.2; // and the N position
                               float scale_N_low    = -1.0;
                               float scale_perp     = -0.9;

                               float s = 0.0;
                               // clipper::Coord_orth arb(0,0,1);
                               clipper::Coord_orth diff_p(ca_pos_2 - ca_pos_1);
                               clipper::Coord_orth diff_p_unit(diff_p.unit());

                               clipper::Coord_orth O_shifted = o_pos_1 - ca_pos_1;
                               double dp = clipper::Coord_orth::dot(diff_p_unit, O_shifted);
                               clipper::Coord_orth pt_on_line(ca_pos_1 + dp * diff_p_unit);
                               clipper::Coord_orth O_perp_vec_unit((o_pos_1-pt_on_line).unit());

                               // clipper::Coord_orth perp(clipper::Coord_orth::cross(arb, diff_p));
                               // clipper::Coord_orth perp_unit(perp.unit());

                               clipper::Coord_orth perp_unit = O_perp_vec_unit; // by another name

                               clipper::Coord_orth double_perp(clipper::Coord_orth::cross(diff_p, perp_unit));
                               clipper::Coord_orth double_perp_unit(double_perp.unit());

                               double along_CA_CA_pt_O = 1.53; // the C is lower down than the O.
                               double along_CA_CA_pt_for_perp = 2.33;

                               double along_CA_CA_pt_N = 2.5;
                               double ideal_peptide_length = 3.81;

                               // we don't want the peptide to be scrunched up on one side of a
                               // "long" peptide... let the atom positions expand along a long peptide.
                               //
                               double diff_p_len = sqrt(diff_p.lengthsq());
                               double f_ca_ca_o = along_CA_CA_pt_O * diff_p_len/ideal_peptide_length;
                               // double f_ca_ca_c = along_CA_CA_pt_C * diff_p_len/ideal_peptide_length;
                               double f_ca_ca_n = along_CA_CA_pt_N * diff_p_len/ideal_peptide_length;
                               double f_ca_ca_pt_for_perp = along_CA_CA_pt_for_perp * diff_p_len/ideal_peptide_length;

                               // there is good density 1.9A away from the mid-line in the direction of the CO
                               // (at the O). 
                               // there is little density 3.2A away from the mid-line in the direction of the CO
                               clipper::Coord_orth rel_line_pt_O(      diff_p_unit * f_ca_ca_o + perp_unit * 1.89);
                               clipper::Coord_orth rel_line_pt_O_low(  diff_p_unit * f_ca_ca_o + perp_unit * 3.2); // was 3.7
                               clipper::Coord_orth rel_line_pt_CO_anti(diff_p_unit * f_ca_ca_o * 0.9 - perp_unit * 0.6);
                               clipper::Coord_orth rel_line_pt_N(      diff_p_unit * f_ca_ca_n - perp_unit * 0.3);
                               clipper::Coord_orth rel_line_pt_N_low(  diff_p_unit * f_ca_ca_n - perp_unit * 1.45);
                               clipper::Coord_orth rel_line_pt_perp1(diff_p_unit * f_ca_ca_pt_for_perp  + double_perp_unit * 1.85);
                               clipper::Coord_orth rel_line_pt_perp2(diff_p_unit * f_ca_ca_pt_for_perp  - double_perp_unit * 1.72);

                               // we don't need to do the spin search here of course, but we need that the
                               // coordinates system is the same as the spin-search. i,e  perp is perpendicular to
                               // CA1-CA-2 in the plane of the peptide, toward the CO

                               // Let's copy across the variable names for consistency with spin search.

                               clipper::Coord_orth p_CO      = ca_pos_1 + rel_line_pt_O;
                               clipper::Coord_orth p_CO_low  = ca_pos_1 + rel_line_pt_O_low;
                               clipper::Coord_orth p_CO_anti = ca_pos_1 + rel_line_pt_CO_anti;
                               clipper::Coord_orth p_N       = ca_pos_1 + rel_line_pt_N;
                               clipper::Coord_orth p_N_low   = ca_pos_1 + rel_line_pt_N_low;
                               clipper::Coord_orth p_2       = ca_pos_1 + rel_line_pt_perp1;
                               clipper::Coord_orth p_3       = ca_pos_1 + rel_line_pt_perp2;
                               float rho_CO      = coot::util::density_at_point(xmap, p_CO);
                               float rho_CO_low  = coot::util::density_at_point(xmap, p_CO_low);
                               float rho_CO_anti = coot::util::density_at_point(xmap, p_CO_anti);
                               float rho_N       = coot::util::density_at_point(xmap, p_N);
                               float rho_N_low   = coot::util::density_at_point(xmap, p_N_low);
                               float rho_perp_1  = coot::util::density_at_point(xmap, p_2);
                               float rho_perp_2  = coot::util::density_at_point(xmap, p_3);
                               float this_score =
                                  scale_CO      * roll_off(rho_CO)      +
                                  scale_CO_low  * roll_off(rho_CO_low)  + 
                                  scale_CO_anti * roll_off(rho_CO_anti) +
                                  scale_N       * roll_off(rho_N)       +
                                  scale_N_low   * roll_off(rho_N_low)   +
                                  scale_perp    * roll_off(rho_perp_1)  +
                                  scale_perp    * roll_off(rho_perp_2);
                               s = this_score;

                               if (true) {
                                  std::string fn = "fp/peptide-fingerprint-" + std::to_string(res_no_base) + "/fp-" + std::to_string(ires) + "-" + std::to_string(i_trial) + ".points";
                                  std::ofstream f(fn.c_str());
                                  if (f) {
                                     f << " CO      " << p_CO.x()      << " " << p_CO.y()      << " " << p_CO.z()      << "\n";
                                     f << " CO_low  " << p_CO_low.x()  << " " << p_CO_low.y()  << " " << p_CO_low.z()  << "\n";
                                     f << " CO_anti " << p_CO_anti.x() << " " << p_CO_anti.y() << " " << p_CO_anti.z() << "\n";
                                     f << " N       " << p_N.x()       << " " << p_N.y()       << " " <<       p_N.z() << "\n";
                                     f << " N_low   " << p_N_low.x()   << " " << p_N_low.y()   << " " <<   p_N_low.z() << "\n";
                                     f << " p_2     " << p_2.x()       << " " << p_2.y()       << " " <<       p_2.z() << "\n";
                                     f << " p_3     " << p_3.x()       << " " << p_3.y()       << " " <<       p_3.z() << "\n";
                                     f << " score " << this_score << "\n";
                                     f.close();
                                  }
                               }
                               return s;
                            };

   for(int ires=frag.min_res_no(); ires<=frag.max_residue_number(); ires++) {
      for (unsigned int iat=0; iat<frag[ires].atoms.size(); iat++) {
         const clipper::Coord_orth atom_pos(frag[ires][iat].pos);
         float w = frag[ires][iat].occupancy;
         float f = util::density_at_point(xmap, atom_pos) * w;
         score += f; // or f * f * f  so that fragments that have atoms in low or negative density score hideously
                     // but we also need a roll-off - a f^5 say that stops fitting into sulphurs (for x-ray fitting)
         w_sum += w;
      }
   }

   for(int ires=frag.min_res_no(); ires<frag.max_residue_number(); ires++) {
      const coot::minimol::residue &res_1 = frag[ires  ];
      const coot::minimol::residue &res_2 = frag[ires+1];
      if (! res_1.is_empty() && ! res_2.is_empty()) {
         auto at_1 = res_1.get_atom(" CA ");
         auto at_2 = res_1.get_atom(" O  ");
         auto at_3 = res_2.get_atom(" CA ");
         get_print_lock();
         if (at_1.first && at_2.first && at_3.first) {
            // std::cout << "debug:: fingerprint_score(): atom at_1 " << at_1.second << std::endl;
            // std::cout << "debug:: type name(args) const;ingerprint_score(): atom at_2 " << at_2.second << std::endl;
            // std::cout << "debug:: fingerprint_score(): atom at_3 " << at_3.second << std::endl;
            const clipper::Coord_orth &pt_1 = at_1.second.pos;
            const clipper::Coord_orth &pt_2 = at_2.second.pos;
            const clipper::Coord_orth &pt_3 = at_3.second.pos;
            float fps = fingerprint_score(pt_1, pt_3, pt_2, ires); // CA, CA, O
            // std::cout << "fingerprint_score " << fps << " cf score " << score << "\n";
            score += fps;
            w_sum += at_1.second.occupancy * 1.0; // or so - because multiple fingerprint positions
         } else {
            std::cout << "Failed to extract CA or O atom from residue - heyho " << ires << std::endl;
         }
         release_print_lock();
      }
   }

   return score/w_sum;
}

// a forward-built residue is made from psi of the previous residue and a phi from "this" one.
// phi_this, psi_prev are in degrees (what a mess)
//
coot::minimol::residue
coot::new_residue_by_3_phi_psi::construct_next_res_from_rama_angles(float phi_this, float psi_prev,
                                                                    float tau, int seqno,
                                                                    const connecting_atoms_t &current_res_pos, float occupancy,
                                                                    dsfmt_t *dsfmt) {

   double jitter_scale = 4.0;

   const clipper::Coord_orth &previous_n  = current_res_pos.N_pos;
   const clipper::Coord_orth &previous_ca = current_res_pos.CA_pos;
   const clipper::Coord_orth &previous_c  = current_res_pos.C_pos;

   // +/- 10 degrees (but in radians)
   double omega_jitter =  20.0 * (M_PI/180.0) * 0.5 * get_random_float_mt(dsfmt) * jitter_scale;
   double O_torsion    = 720.0 * (M_PI/180.0) * 0.5 * get_random_float_mt(dsfmt);

   double r1 = get_random_float_mt(dsfmt) * 2.0 - 1.0;
   double r2 = get_random_float_mt(dsfmt) * 2.0 - 1.0;
   double r3 = get_random_float_mt(dsfmt) * 2.0 - 1.0;

   clipper::Coord_orth jittered_delta(jitter_scale * 0.1 * clipper::Coord_orth(r1,r2,r3));

   double angle, torsion;

   coot::minimol::residue mres(seqno);
   mres.name = "ALA";

   // N  of next is placed by psi_previous
   // CA of next is placed by omega
   // C  of next is placed by phi_next
   //
   // So, for clarity, say we are building on residue 54 forwards
   // 55 N  is placed by psi of residue 54
   // 55 CA is places by omega
   // 55 C  is placed by phi of residue 55
   // 55 O  is placed by psi of residue 55 (because it places the N of 56)

   angle = clipper::Util::d2rad(116.200); // Ca-C-N
   torsion = clipper::Util::d2rad(psi_prev);
   clipper::Coord_orth n_pos(previous_n, previous_ca, previous_c, 1.329, angle, torsion); // C-N bond

   n_pos += jittered_delta;

   angle = clipper::Util::d2rad(121.700); // C-N-Ca
   torsion = clipper::Util::d2rad(180.0) + omega_jitter;
   clipper::Coord_orth ca_pos(previous_ca, previous_c, n_pos, 1.458, angle, torsion); // N-CA bond

   angle = clipper::Util::d2rad(tau); // N-CA-C  // was 111.200
   torsion = clipper::Util::d2rad(phi_this);
   clipper::Coord_orth c_pos(previous_c, n_pos, ca_pos, 1.525, angle, torsion); // CA-C bond

   angle = clipper::Util::d2rad(120.800);  // CA-C-O
   torsion = O_torsion; // Random.  It depends on the *next* residue, not the prevous one.
                        // i.e. should be: ca_next, n_next, c_pos with torsion 0.0
   clipper::Coord_orth o_pos(n_pos, ca_pos, c_pos, 1.231, angle, torsion); // C-O bond

   float b_factor = 20.0;
   mres.addatom(coot::minimol::atom(" N  ", " N", n_pos,  "", occupancy, b_factor));
   mres.addatom(coot::minimol::atom(" C  ", " C", c_pos,  "", occupancy, b_factor));
   mres.addatom(coot::minimol::atom(" CA ", " C", ca_pos, "", occupancy, b_factor));
   mres.addatom(coot::minimol::atom(" O  ", " O", o_pos,  "", occupancy, b_factor));

   return mres;
}

// phi is phi of the residue to which we are joining.
//
// static
coot::minimol::residue
coot::new_residue_by_3_phi_psi::construct_prev_res_from_rama_angles(float phi, float psi, float tau,
                                                                    int seqno, const connecting_atoms_t &current_res_pos, float occupancy,
                                                                    dsfmt_t *dsfmt) {

   double jitter_scale = 4.0;

   coot::minimol::residue mres(seqno);
   mres.name = "ALA";

   double angle, torsion;

   double r1 = get_random_float_mt(dsfmt) * 2.0 - 1.0;
   double r2 = get_random_float_mt(dsfmt) * 2.0 - 1.0;
   double r3 = get_random_float_mt(dsfmt) * 2.0 - 1.0;

   clipper::Coord_orth jittered_delta(jitter_scale * 0.1 * clipper::Coord_orth(r1,r2,r3));

   // C
   angle = clipper::Util::d2rad(121.700); // C-N-Ca
   torsion = clipper::Util::d2rad(phi);
   clipper::Coord_orth c_pos(current_res_pos.C_pos, current_res_pos.CA_pos, current_res_pos.N_pos, 1.329, angle, torsion); // C-N bond

   c_pos += jittered_delta;

   // Ca
   angle = clipper::Util::d2rad(116.200); // Ca-C-N
   torsion = clipper::Util::d2rad(180.0);
   clipper::Coord_orth ca_pos(current_res_pos.CA_pos, current_res_pos.N_pos, c_pos, 1.525, angle, torsion); // Ca-C bond

   // N
   angle = clipper::Util::d2rad(tau); // N-Ca-C
   torsion = clipper::Util::d2rad(psi);
   clipper::Coord_orth n_pos(current_res_pos.N_pos, c_pos, ca_pos, 1.458, angle, torsion); // Ca-N

   // O
   angle = clipper::Util::d2rad(120.800);  // Ca-C-O
   torsion = clipper::Util::d2rad(0.0);
   clipper::Coord_orth o_pos(current_res_pos.CA_pos, current_res_pos.N_pos, c_pos, 1.231, angle, torsion); // C=0

   // CB
   angle = clipper::Util::d2rad(tau);
   torsion = clipper::Util::d2rad(psi+123.4);
   clipper::Coord_orth cb_pos(current_res_pos.N_pos, c_pos, ca_pos, 1.52, angle, torsion); // Ca-Cb


   // ---------------------------

   float b_factor = 20.0;
   mres.addatom(minimol::atom(" N  ", " N", n_pos,  "", occupancy, b_factor));
   mres.addatom(minimol::atom(" C  ", " C", c_pos,  "", occupancy, b_factor));
   mres.addatom(minimol::atom(" CA ", " C", ca_pos, "", occupancy, b_factor));
   mres.addatom(minimol::atom(" O  ", " O", o_pos,  "", occupancy, b_factor));

   return mres;
}



coot::minimol::fragment
coot::new_residue_by_3_phi_psi::make_3_res_joining_frag_forward(const std::string &chain_id,
                                                                const connecting_atoms_t &current_res_pos_in,
                                                                const double &psi_conditional_deg,
                                                                const phi_psi_t &pp_1,
                                                                const phi_psi_t &pp_2,
                                                                const phi_psi_t &pp_3,
                                                                int seq_num,
                                                                dsfmt_t *dsfmt) {

   auto make_CB_ideal_pos = [] (minimol::residue &res) {
                               std::pair<bool, coot::minimol::atom> CB = res.get_atom(" CB ");
                               bool generated = false;
                               clipper::Coord_orth pos;
                               if (! CB.first) {
                                  auto CA = res.get_atom(" CA ");
                                  auto C  = res.get_atom(" C  ");
                                  auto N  = res.get_atom(" N  ");
                                  if (CA.first) {
                                     if (C.first) {
                                        if (N.first) {
                                           clipper::Coord_orth C_to_N = N.second.pos - C.second.pos;
                                           clipper::Coord_orth C_to_N_mid_point(0.5 * (N.second.pos + C.second.pos));
                                           clipper::Coord_orth CA_to_CN_mid_point = C_to_N_mid_point - CA.second.pos;
                                           clipper::Coord_orth CA_to_CN_mid_point_uv(CA_to_CN_mid_point.unit());
                                           clipper::Coord_orth perp(clipper::Coord_orth::cross(C_to_N, CA_to_CN_mid_point));
                                           clipper::Coord_orth perp_uv(perp.unit());
                                           // guess and fiddle these - good enough (function copied from res-trace.cc)
                                           clipper::Coord_orth CB_pos(CA.second.pos + 1.21 * perp_uv - 0.95 * CA_to_CN_mid_point_uv);
                                           pos = CB_pos;
                                           generated = true;
                                        }
                                     }
                                  }
                               }
                               return std::make_pair(generated, pos);
                            };


   coot::minimol::fragment frag(chain_id);
   // we need edit the postions of the reference atoms to add a bit of (useful) jitter
   connecting_atoms_t current_res_pos(current_res_pos_in);
   double rand_lim = 0.1;
   clipper::Coord_orth &current_n  = current_res_pos.N_pos;
   clipper::Coord_orth &current_ca = current_res_pos.CA_pos;

   // Am I double jittering here!?

   double r1 = get_random_float_mt(dsfmt) * 2.0 - 1.0;
   double r2 = get_random_float_mt(dsfmt) * 2.0 - 1.0;
   double r3 = get_random_float_mt(dsfmt) * 2.0 - 1.0;
   double r4 = get_random_float_mt(dsfmt) * 2.0 - 1.0;
   double r5 = get_random_float_mt(dsfmt) * 2.0 - 1.0;
   double r6 = get_random_float_mt(dsfmt) * 2.0 - 1.0;

   current_n  += clipper::Coord_orth(rand_lim * 0.5 * r1, rand_lim * 0.5 * r2, rand_lim * 0.5 * r3);
   current_ca += clipper::Coord_orth(rand_lim * 0.5 * r4, rand_lim * 0.5 * r5, rand_lim * 0.5 * r6);

   if (true) {
      minimol::residue res_1 = construct_next_res_from_rama_angles(pp_1.phi, psi_conditional_deg, pp_1.tau, seq_num + 1, current_res_pos, 1.0, dsfmt);

      connecting_atoms_t just_built_res_1(res_1[" N  "].pos, res_1[" CA "].pos, res_1[" C  "].pos);
      just_built_res_1.set_upstream_C(current_res_pos.C_pos);
      minimol::residue res_2 = construct_next_res_from_rama_angles(pp_2.phi, pp_1.psi, pp_2.tau, seq_num + 2, just_built_res_1, 0.8, dsfmt);
      connecting_atoms_t just_built_res_2(res_2[" N  "].pos, res_2[" CA "].pos, res_2[" C  "].pos);
      minimol::residue res_3 = construct_next_res_from_rama_angles(pp_3.phi, pp_3.psi, pp_3.tau, seq_num + 3, just_built_res_2, 0.5, dsfmt);

      // now set set the occupancy of res2 to 0.5 or so, because we care more that the
      // first residue is in density
      // ... or do we?
      // for (unsigned int iat=0; iat<res2.atoms.size(); iat++)
      // res2.atoms[iat].occupancy = 0.5;

      // now let's correct the position of the O atom of residue 1,
      // now that we have residue 2. And likewise, the position of the O
      // in residue 2 now that we have residue 3.
      //
      double angle   = clipper::Util::d2rad(123.0); // N-C-O
      double torsion = clipper::Util::d2rad(0.0);
      clipper::Coord_orth o_pos_1(res_2[" CA "].pos, res_2[" N  "].pos, res_1[" C  "].pos, 1.231, angle, torsion);
      clipper::Coord_orth o_pos_2(res_3[" CA "].pos, res_3[" N  "].pos, res_2[" C  "].pos, 1.231, angle, torsion);
      for (unsigned int iat=0; iat<res_1.atoms.size(); iat++)
	 if (res_1[iat].name == " O  ")
	    res_1[iat].pos = o_pos_1;
      for (unsigned int iat=0; iat<res_2.atoms.size(); iat++)
	 if (res_2[iat].name == " O  ")
	    res_2[iat].pos = o_pos_2;

      try {
         const std::pair<bool, clipper::Coord_orth> &CB_res_1 = make_CB_ideal_pos(res_1);
         const std::pair<bool, clipper::Coord_orth> &CB_res_2 = make_CB_ideal_pos(res_2);
         const std::pair<bool, clipper::Coord_orth> &CB_res_3 = make_CB_ideal_pos(res_3);
         if (CB_res_1.first) res_1.addatom(minimol::atom(" CB ", " C", CB_res_1.second, "", 0.5f, 20.0f));
         if (CB_res_2.first) res_2.addatom(minimol::atom(" CB ", " C", CB_res_2.second, "", 0.4f, 20.0f));
         if (CB_res_3.first) res_3.addatom(minimol::atom(" CB ", " C", CB_res_3.second, "", 0.3f, 20.0f));

	 frag.addresidue(res_1, 0);
	 frag.addresidue(res_2, 0);
	 frag.addresidue(res_3, 0);

         if (false) {
            std::string fn = "make_3_res_joining_frag_forward-" + chain_id + "-" + std::to_string(seq_num) + ".pdb";
            frag.write_file(fn);
         }
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "ERROR:: make_3_res_joining_frag_new_building_forwards() "
		   << rte.what() << std::endl;
      }
   }
   return frag;
}


coot::minimol::fragment
coot::new_residue_by_3_phi_psi::make_3_res_joining_frag_backward(const std::string &chain_id,
                                                                 const connecting_atoms_t &current_res_pos,
                                                                 const double &phi_conditional_deg,
                                                                 const phi_psi_t &pp_1,
                                                                 const phi_psi_t &pp_2,
                                                                 const phi_psi_t &pp_3,
                                                                 int seq_num,
                                                                 dsfmt_t *dsfmt) {
   coot::minimol::fragment frag(chain_id);

   clipper::Coord_orth current_n  = current_res_pos.N_pos;
   clipper::Coord_orth current_ca = current_res_pos.CA_pos;

   if (true) {

      // add a bit of jitter
      double rand_lim = 0.1;

   // Am I double jittering here!?

      double r1 = get_random_float_mt(dsfmt) * 2.0 - 1.0;
      double r2 = get_random_float_mt(dsfmt) * 2.0 - 1.0;
      double r3 = get_random_float_mt(dsfmt) * 2.0 - 1.0;
      double r4 = get_random_float_mt(dsfmt) * 2.0 - 1.0;
      double r5 = get_random_float_mt(dsfmt) * 2.0 - 1.0;
      double r6 = get_random_float_mt(dsfmt) * 2.0 - 1.0;

      current_n  += clipper::Coord_orth(rand_lim * 0.5 * r1, rand_lim * 0.5 * r2, rand_lim * 0.5 * r3);
      current_ca += clipper::Coord_orth(rand_lim * 0.5 * r4, rand_lim * 0.5 * r5, rand_lim * 0.5 * r6);

      minimol::residue res1 = construct_prev_res_from_rama_angles(phi_conditional_deg, pp_1.psi, pp_1.tau, seq_num - 1, current_res_pos, 1.0, dsfmt);

      connecting_atoms_t just_built_res_1(res1[" N  "].pos, res1[" CA "].pos, res1[" C  "].pos);
      just_built_res_1.set_downstream_N(current_res_pos.N_pos);
      minimol::residue res2 = construct_prev_res_from_rama_angles(pp_1.phi, pp_2.psi, pp_2.tau, seq_num - 2, just_built_res_1, 0.8, dsfmt);
      connecting_atoms_t just_built_res_2(res2[" N  "].pos, res2[" CA "].pos, res2[" C  "].pos);
      minimol::residue res3 = construct_prev_res_from_rama_angles(pp_2.phi, pp_3.psi, pp_3.tau, seq_num - 3, just_built_res_2, 0.5, dsfmt);

      // for (unsigned int iat=0; iat<res2.atoms.size(); iat++)
      // res2.atoms[iat].occupancy = 0.5;

      try {
	 frag.addresidue(res3, 0);
	 frag.addresidue(res2, 0);
	 frag.addresidue(res1, 0);
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "ERROR:: make_2_res_joining_frag_new_building_forwards() "
		   << rte.what() << std::endl;
      }

   }
   return frag;
}

coot::minimol::fragment
coot::new_residue_by_3_phi_psi::best_fit_phi_psi(unsigned int n_trials, const clipper::Xmap<float> &xmap,
                                                 float min_density_level_for_connecting_atom) const {

   auto run_forward_trials = [] (int thread_id, std::pair<unsigned int, unsigned int> trial_start_stop, 
                                 double phi_current, const clipper::Ramachandran &rama, float rama_max,
                                 const std::string &chain_id, const connecting_atoms_t &current_res_pos, int seq_num,
                                 const clipper::Xmap<float> *xmap, float min_density_level_for_connecting_atom,
                                 std::pair<minimol::fragment, float> &best_frag_result,
                                 std::atomic<unsigned int> &count) {

                                float best_score = -9999.9;
                                dsfmt_t dsfmt;
                                uint32_t seed = 1;
                                dsfmt_gv_init_gen_rand(seed);
                                dsfmt_init_gen_rand(&dsfmt, seed);

                                for (unsigned int i_trial=trial_start_stop.first; i_trial<trial_start_stop.second; i_trial++) {

                                   //std::cout << "debug:: in run_forward_trials() i_trial " << i_trial << std::endl;

                                   double psi_conditional = get_psi_by_random_given_phi(phi_current, rama, &dsfmt); // in radians

                                   phi_psi_t pp_1 = get_phi_psi_by_random(rama, rama_max, false, &dsfmt);
                                   phi_psi_t pp_2 = get_phi_psi_by_random(rama, rama_max, false, &dsfmt);
                                   phi_psi_t pp_3 = get_phi_psi_by_random(rama, rama_max, false, &dsfmt);

                                   minimol::fragment frag = make_3_res_joining_frag_forward(chain_id, current_res_pos, clipper::Util::rad2d(psi_conditional),
                                                                                            pp_1, pp_2, pp_3, seq_num, &dsfmt);
                                   float score = score_fragment_using_peptide_fingerprint(frag, current_res_pos, *xmap, seq_num, i_trial); // pass i_trial for debugging
                                   float score_basic = score_fragment_basic(frag, current_res_pos, *xmap);
                                   if (score > best_score) {

                                      // std::cout << "residue " << seq_num << " i_trial " << i_trial << " score_basic " << score_basic
                                      // << " score-with-fp " << score << std::endl;
                                      best_score = score;
                                      best_frag_result.first = frag;
                                      best_frag_result.second = score;
                                      if (false) {
                                         std::string fn = "run_forward_trials_" + std::to_string(seq_num) + "_" + std::to_string(i_trial) + ".pdb";
                                         frag.write_file(fn);
                                      }
                                   }
                                }
                                count++;
                             };

   auto run_backward_trials = [] (int thread_id, std::pair<unsigned int, unsigned int> trial_start_stop, 
                                  double psi_current, const clipper::Ramachandran &rama, float rama_max,
                                  const std::string &chain_id, const connecting_atoms_t &current_res_pos, int seq_num,
                                  const clipper::Xmap<float> *xmap, float min_density_level_for_connecting_atom,
                                  std::pair<minimol::fragment, float> &best_frag_result,
                                  std::atomic<unsigned int> &count) {

                                 dsfmt_t dsfmt;
                                 uint32_t seed = 1;
                                 dsfmt_gv_init_gen_rand(seed);
                                 dsfmt_init_gen_rand(&dsfmt, seed);
                                 float best_score = -9999.9;
                                 for (unsigned int i_trial=trial_start_stop.first; i_trial<trial_start_stop.second; i_trial++) {

                                    double phi_conditional = get_phi_by_random_given_psi(psi_current, rama, &dsfmt); // in radians

                                    phi_psi_t pp_1 = get_phi_psi_by_random(rama, rama_max, false, &dsfmt);
                                    phi_psi_t pp_2 = get_phi_psi_by_random(rama, rama_max, false, &dsfmt);
                                    phi_psi_t pp_3 = get_phi_psi_by_random(rama, rama_max, false, &dsfmt);

                                    minimol::fragment frag = make_3_res_joining_frag_backward(chain_id, current_res_pos,
                                                                                              clipper::Util::rad2d(phi_conditional),
                                                                                              pp_1, pp_2, pp_3, seq_num, &dsfmt);
                                    // inital previous-residue C density test:
                                    // don't bother with scoring the triple peptide if this is not at least 1 rmsd (min density level for connecting atom)
                                    // We can do this *after* make_3_res_joining_frag_backward if make_3_res_joining_frag_backward() is fast
                                    // otherwise we put this test *inside* make_3_res_joining_frag_backward()
                                    //
                                    // Put this in the forward build too
                                    //
                                    std::pair<bool, coot::minimol::atom> pos_pair = frag[seq_num-1].get_atom(" C  ");
                                    if (pos_pair.first) {
                                       float f = coot::util::density_at_point(*xmap, pos_pair.second.pos);
                                       if (f < min_density_level_for_connecting_atom)
                                          continue;
                                    } else {
                                       std::cout << "Hideous failure in run_backward_trials() " << std::endl;
                                    }
                                    
                                    float score = score_fragment_using_peptide_fingerprint(frag, current_res_pos, *xmap, seq_num, i_trial); // pass i_trial for debugging
                                    if (score > best_score) {
                                       best_score = score;
                                       best_frag_result.first = frag;
                                       best_frag_result.second = score;
                                       if (false) {
                                          std::string fn = "run_backward_trials_" + chain_id + "-" + std::to_string(seq_num) + "_" +
                                             std::to_string(i_trial) + ".pdb";
                                          frag.write_file(fn);
                                      }
                                    }
                                 }
                                 count++;
                              };

   coot::minimol::fragment best_frag;

   // not to self: try to restore the version without threading for speed comparison - because this current implementationn
   // seems slow.  Mabye grid scan phi and phis, rather than random selection?

   if (terminus_type == "C") {
      connecting_atoms_t current_res_pos = get_connecting_residue_atoms();
      std::pair<bool,double> phi_current = current_res_pos.get_phi();
      // forwards,  we have a phi and need to generate a psi to place the N

      if (false)
         std::cout << "debug:: best_fit_phi_psi(): C extension current_phi: " << coot::residue_spec_t(residue_p) << " phi: "
                   << phi_current.first << " " << phi_current.second << std::endl;

      std::atomic<unsigned int> count(0);

      int seq_num = residue_p->GetSeqNum();
      if (phi_current.first) {
      } else {
         std::cout << "INFO:: Missing phi " << residue_spec_t(residue_p) << " inventing -120" << std::endl;
         phi_current.second = clipper::Util::d2rad(-120.0);
      }
      float score_for_best_frag = -9999.9;
      std::vector<std::pair<unsigned int, unsigned int> > ranges = atom_index_ranges(n_trials, n_threads);
      std::vector<std::pair<minimol::fragment, float> > best_frag_vec(ranges.size()); // best frag for that thread/trial-range
      for (unsigned int ir=0; ir<ranges.size(); ir++) {
         thread_pool_p->push(run_forward_trials, ranges[ir], phi_current.second, std::cref(rama), rama_max, std::cref(chain_id),
                             std::cref(current_res_pos), seq_num, &xmap, min_density_level_for_connecting_atom,
                             std::ref(best_frag_vec[ir]), std::ref(count));
      }
      while (count != ranges.size()) {
         // std::cout << "waiting for trial sets: done " << count << " of " << ranges.size() << " ranges " << std::endl;
         std::this_thread::sleep_for(std::chrono::milliseconds(1));
      }
      for (unsigned int ir=0; ir<ranges.size(); ir++) {
         if (best_frag_vec[ir].second > score_for_best_frag) {
            best_frag = best_frag_vec[ir].first;
            score_for_best_frag = best_frag_vec[ir].second;
         }
      }
   }

   if (terminus_type == "N") {
      // backwards we have a psi and need to generate a phi to place the C

      connecting_atoms_t current_res_pos = get_connecting_residue_atoms();
      std::pair<bool,double> psi_current = current_res_pos.get_psi();

      int seq_num = residue_p->GetSeqNum();
      std::atomic<unsigned int> count(0);

      if (psi_current.first) {
      } else {
         std::cout << "INFO:: Missing psi " << residue_spec_t(residue_p) << " inventing +120" << std::endl;
         psi_current.second = clipper::Util::d2rad(120.0);
      }
            
      float score_for_best_frag = -9999.9;
      std::vector<std::pair<unsigned int, unsigned int> > ranges = atom_index_ranges(n_trials, n_threads);
      std::vector<std::pair<minimol::fragment, float> > best_frag_vec(ranges.size()); // best frag for that thread/trial-range
      for (unsigned int ir=0; ir<ranges.size(); ir++) {
         thread_pool_p->push(run_backward_trials, ranges[ir], psi_current.second, std::cref(rama), rama_max, std::cref(chain_id),
                             std::cref(current_res_pos), seq_num, &xmap,
                             min_density_level_for_connecting_atom,
                             std::ref(best_frag_vec[ir]), std::ref(count));
      }
      while (count != ranges.size()) {
         // std::cout << "waiting for trial sets: done " << count << " of " << ranges.size() << " ranges " << std::endl;
         std::this_thread::sleep_for(std::chrono::microseconds(10));
      }
      for (unsigned int ir=0; ir<ranges.size(); ir++) {
         if (best_frag_vec[ir].second > score_for_best_frag) {
            best_frag = best_frag_vec[ir].first;
            score_for_best_frag = best_frag_vec[ir].second;
         }
      }
   }

   return best_frag;
}
