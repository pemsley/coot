// Perform structure factor calculation
/* Copyright 2003-2004 Kevin Cowtan & University of York all rights reserved */
/* Cootenized by Paul Emsley */


#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/core/clipper_util.h>

#include "sfcalc-genmap.hh"

// calculate structure factors from the given model and data
// and update the map xmap_p
void coot::util::sfcalc_genmap(mmdb::Manager *mol,
                               const clipper::HKL_data<clipper::data32::F_sigF> &fobs_in,
                               const clipper::HKL_data<clipper::data32::Flag> &free,
                               clipper::Xmap<float> *xmap_p) {

   if (fobs_in.num_obs() == 0) {
      std::cout << "sfcalc_genmap(): No Fobs reflections\n";
      return;
   }
   enum ANISO { NONE, FOBS, FCAL };
   bool bulk  = true;
   ANISO aniso = NONE;
   int freeflag = 0;
   int n_refln = 1000;
   int n_param = 20;
   int verbose = 0;
   typedef clipper::HKL_data_base::HKL_reference_index HRI;

   // modifyable fobs (because we want to anisotropically scale Fobs)
   clipper::HKL_data<clipper::data32::F_sigF> fobs(fobs_in);

   // get a list of all the atoms
   clipper::mmdb::CAtom **atom_sel = 0;
   int nsel = 0;
   int hndl = mol->NewSelection();  // d
   mol->SelectAtoms(hndl, 0, 0, ::mmdb::SKEY_NEW);
   mol->GetSelIndex(hndl, atom_sel, nsel);
   clipper::MMDBAtom_list atoms(atom_sel, nsel);

   // clipper::MTZcrystal cxtl;
   clipper::HKL_info hkls;
   hkls.init(fobs.spacegroup(),fobs.cell(),fobs.hkl_sampling(), true); // init this correctly - how?
   double bulkfrc, bulkscl;

   // calculate structure factors
   //
   // clipper::HKL_data<clipper::data32::F_phi> fc(hkls, cxtl);
   clipper::HKL_data<clipper::data32::F_phi> fc(fobs.spacegroup(),fobs.cell(),fobs.hkl_sampling());
   if (bulk) {
      clipper::SFcalc_obs_bulk<float> sfcb;
      sfcb(fc, fobs, atoms);
      bulkfrc = sfcb.bulk_frac();
      bulkscl = sfcb.bulk_scale();
   } else {
      clipper::SFcalc_aniso_fft<float> sfc;
      sfc(fc, atoms);
      bulkfrc = bulkscl = 0.0;
   }

  // do anisotropic scaling
   if (aniso != NONE)  {
      clipper::SFscale_aniso<float>::TYPE F = clipper::SFscale_aniso<float>::F;
      clipper::SFscale_aniso<float> sfscl;
      if ( aniso == FOBS ) sfscl(fobs, fc);  // scale Fobs
      if ( aniso == FCAL ) sfscl(fc, fobs);  // scale Fcal
      std::cout << "\nAnisotropic scaling:\n" << sfscl.u_aniso_orth(F).format() << "\n";
   }

   // now do sigmaa calc
   // clipper::HKL_data<clipper::data32::F_phi> fb(hkls, cxtl);
   // clipper::HKL_data<clipper::data32::F_phi> fd( hkls, cxtl);
   // clipper::HKL_data<clipper::data32::Phi_fom> phiw( hkls, cxtl );
   // clipper::HKL_data<clipper::datatypes::Flag>    flag( hkls, cxtl );
   clipper::HKL_data<clipper::data32::F_phi> f_best(fobs.spacegroup(),fobs.cell(),fobs.hkl_sampling());
   clipper::HKL_data<clipper::data32::F_phi> f_diff(fobs.spacegroup(),fobs.cell(),fobs.hkl_sampling());
   clipper::HKL_data<clipper::data32::Phi_fom> phiw(fobs.spacegroup(),fobs.cell(),fobs.hkl_sampling());
   clipper::HKL_data<clipper::data32::Flag>    flag(fobs.spacegroup(),fobs.cell(),fobs.hkl_sampling());

   for ( HRI ih = flag.first(); !ih.last(); ih.next() )
      if (!fobs[ih].missing() && (free[ih].missing()||free[ih].flag()==freeflag))
         flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
      else
         flag[ih].flag() = clipper::SFweight_spline<float>::NONE;

   // do sigmaa calc
   clipper::SFweight_spline<float> sfw(n_refln, n_param);
   sfw(f_best, f_diff, phiw, fobs, fc, flag);

   if (true) {
      unsigned int n_nans_fobs   = 0;
      unsigned int n_nans_fc     = 0;
      unsigned int n_nans_f_diff = 0;
      for (HRI ih = fobs.first(); !ih.last(); ih.next()) {
         if (!fobs[ih].missing() && (free[ih].missing()||free[ih].flag()==freeflag)) {
            if (clipper::Util::isnan(fobs[ih].f())) n_nans_fobs++;
            if (clipper::Util::isnan(fc[ih].f()))   n_nans_fc++;
            if (clipper::Util::isnan(f_diff[ih].f())) n_nans_f_diff++;
         }
      }
      // std::cout << "DEBUG:: the nan count: " << n_nans_fobs << " " << n_nans_fc << " " << n_nans_f_diff << std::endl;
   }

   // calc abcd (needed?)
   clipper::HKL_data<clipper::data32::ABCD> abcd(hkls);
   abcd.compute(phiw, clipper::data32::Compute_abcd_from_phifom());

   mol->DeleteSelection(hndl);

   // now calc R and R-free
   std::vector<double> params(n_param, 1.0);
   clipper::BasisFn_spline basisfn(fobs, n_param, 1.0);
   clipper::TargetFn_scaleF1F2<clipper::data32::F_phi,clipper::data32::F_sigF> targetfn(fc, fobs);
   clipper::ResolutionFn rfn(hkls, basisfn, targetfn, params);
   double r1w, f1w, r1f, f1f, Fo, Fc;
   r1w = f1w = r1f = f1f = 0.0;
   for (HRI ih = fobs.first(); !ih.last(); ih.next()) {
      if (!fobs[ih].missing()) {
         Fo = fobs[ih].f();
         Fc = sqrt( rfn.f(ih) ) * fc[ih].f();
         if (free[ih].flag() == freeflag) {
            r1f += fabs(Fo - Fc);
            f1f += Fo;
         } else {
            r1w += fabs(Fo - Fc);
            f1w += Fo;
         }
      }
   }
   r1f /= clipper::Util::max(f1f, 0.1);
   r1w /= clipper::Util::max(f1w, 0.1);
   if (false)
      std::cout << "\n R-factor      : " << r1w << "\n Free R-factor : " << r1f << "\n";

   // now make a map

   // test
   // xmap_p->fft_from(fc);
   xmap_p->fft_from(f_diff);

   // DIAGNOSTIC OUTPUT
   if (false) {
      std::cout << "\n Bulk Correction Volume: " << bulkfrc;
      std::cout << "\n Bulk Correction Factor: " << bulkscl << "\n";
      std::cout << "\nNumber of spline params: " << sfw.params_scale().size() << "\n";
      clipper::BasisFn_spline basisfn( hkls, sfw.params_scale().size(), 1.0 );
      printf("\n $TABLE: Sigmaa statistics :\n $GRAPHS:scale vs resolution:N:1,2:\n        :lack of closure vs resolution:N:1,3:\n $$\n 1/resol^2   scale   lack_of_closure $$\n $$\n");
      for (int i = 0; i <= 20; i++) {
         double s = hkls.resolution().invresolsq_limit()*double(i)/20.0;
         printf("%6.3f %12.3f %12.3f\n", s, basisfn.f_s(s,sfw.params_scale()),
                                         basisfn.f_s(s,sfw.params_error()));
      }
      printf(" $$\n");
   }
}

#include <chrono>

coot::util::sfcalc_genmap_stats_t
coot::util::sfcalc_genmaps_using_bulk_solvent(mmdb::Manager *mol,
                                              const clipper::HKL_data<clipper::data32::F_sigF> &fobs_in,
                                              const clipper::HKL_data<clipper::data32::Flag> &free,
                                              const clipper::Cell &cell_for_fobs,
                                              clipper::Xmap<float> *xmap_2fofc_p,
                                              clipper::Xmap<float> *xmap_fofc_p) {

   auto tp_0 = std::chrono::high_resolution_clock::now();

   sfcalc_genmap_stats_t sfcgs;

   if (fobs_in.num_obs() == 0) {
      std::cout << "sfcalc_genmap(): No Fobs reflections\n";
      return sfcgs;
   }

   if (false) {

      // sanity check
      const clipper::HKL_info &hkls_check = fobs_in.base_hkl_info();
      const clipper::Spacegroup &spgr_check = hkls_check.spacegroup();

      std::cout << "DEBUG:: Sanity check A in sfcalc_genmaps_using_bulk_solvent(): HKL_info: "
                << "cell_for_fobs: " << cell_for_fobs.format() << " "
                << "cell of fobs: " << hkls_check.cell().format() << " "
                << "spacegroup: " << spgr_check.symbol_xhm() << " "
                << "resolution: " << hkls_check.resolution().limit() << " "
                << "invsqreslim: " << hkls_check.resolution().invresolsq_limit() << " "
                << std::endl;
   }


   enum ANISO { NONE, FOBS, FCAL };
   bool bulk  = true;
   ANISO aniso = NONE;
   int freeflag = 0;
   int n_refln = 1000;
   int n_param = 20;

   typedef clipper::HKL_data_base::HKL_reference_index HRI;

   // modifyable fobs (because we want to anisotropically scale Fobs)
   clipper::HKL_data<clipper::data32::F_sigF> fobs(fobs_in);

   auto tp_1 = std::chrono::high_resolution_clock::now();

   // get a list of all the atoms
   clipper::mmdb::CAtom **atom_sel = 0;
   int nsel = 0;
   int hndl = mol->NewSelection();  // d
   mol->SelectAtoms(hndl, 0, 0, ::mmdb::SKEY_NEW);
   mol->GetSelIndex(hndl, atom_sel, nsel);
   clipper::MMDBAtom_list atoms(atom_sel, nsel);

   auto tp_2 = std::chrono::high_resolution_clock::now();
   // clipper::MTZcrystal cxtl;
   clipper::HKL_info hkls;
   hkls.init(fobs.spacegroup(), cell_for_fobs, fobs.hkl_sampling(), true); // init this correctly - how?
   double bulkfrc, bulkscl;

   // calculate structure factors
   //
   // clipper::HKL_data<clipper::data32::F_phi> fc(hkls, cxtl);
   clipper::HKL_data<clipper::data32::F_phi> fc(fobs.spacegroup(), cell_for_fobs, fobs.hkl_sampling());

   auto tp_3 = std::chrono::high_resolution_clock::now();
   if (bulk) {
      clipper::SFcalc_obs_bulk<float> sfcb;
      sfcb(fc, fobs, atoms);
      bulkfrc = sfcb.bulk_frac();
      bulkscl = sfcb.bulk_scale();
   } else {
      clipper::SFcalc_aniso_fft<float> sfc;
      sfc(fc, atoms);
      bulkfrc = bulkscl = 0.0;
   }
   auto tp_4 = std::chrono::high_resolution_clock::now();

  // do anisotropic scaling
   if (aniso != NONE)  {
      clipper::SFscale_aniso<float>::TYPE F = clipper::SFscale_aniso<float>::F;
      clipper::SFscale_aniso<float> sfscl;
      if ( aniso == FOBS ) sfscl(fobs, fc);  // scale Fobs
      if ( aniso == FCAL ) sfscl(fc, fobs);  // scale Fcal
      std::cout << "\nAnisotropic scaling:\n" << sfscl.u_aniso_orth(F).format() << "\n";
   }

   // now do sigmaa calc
   // clipper::HKL_data<clipper::data32::F_phi> fb(hkls, cxtl);
   // clipper::HKL_data<clipper::data32::F_phi> fd( hkls, cxtl);
   // clipper::HKL_data<clipper::data32::Phi_fom> phiw( hkls, cxtl );
   // clipper::HKL_data<clipper::datatypes::Flag>    flag( hkls, cxtl );
   clipper::HKL_data<clipper::data32::F_phi> f_best(fobs.spacegroup(), cell_for_fobs, fobs.hkl_sampling());
   clipper::HKL_data<clipper::data32::F_phi> f_diff(fobs.spacegroup(), cell_for_fobs, fobs.hkl_sampling());
   clipper::HKL_data<clipper::data32::F_phi> f_2mfodfc(fobs.spacegroup(), cell_for_fobs, fobs.hkl_sampling());
   clipper::HKL_data<clipper::data32::Phi_fom> phiw(fobs.spacegroup(), cell_for_fobs, fobs.hkl_sampling());
   clipper::HKL_data<clipper::data32::Flag>    flag(fobs.spacegroup(), cell_for_fobs, fobs.hkl_sampling());

   auto tp_5 = std::chrono::high_resolution_clock::now();

   for ( HRI ih = flag.first(); !ih.last(); ih.next() )
      if (!fobs[ih].missing() && (free[ih].missing()||free[ih].flag()==freeflag))
         flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
      else
         flag[ih].flag() = clipper::SFweight_spline<float>::NONE;

   auto tp_6 = std::chrono::high_resolution_clock::now();

   // do sigmaa calc
   clipper::SFweight_spline<float> sfw(n_refln, n_param);
   sfw(f_best, f_diff, phiw, fobs, fc, flag);

   if (true) {
      unsigned int n_nans_fobs   = 0;
      unsigned int n_nans_fc     = 0;
      unsigned int n_nans_f_diff = 0;
      for (HRI ih = fobs.first(); !ih.last(); ih.next()) {
         if (!fobs[ih].missing() && (free[ih].missing()||free[ih].flag()==freeflag)) {
            if (clipper::Util::isnan(fobs[ih].f())) n_nans_fobs++;
            if (clipper::Util::isnan(fc[ih].f()))   n_nans_fc++;
            if (clipper::Util::isnan(f_diff[ih].f())) n_nans_f_diff++;
         }
      }
      // std::cout << "DEBUG:: the nan count: " << n_nans_fobs << " " << n_nans_fc << " " << n_nans_f_diff << std::endl;
   }

   auto tp_7 = std::chrono::high_resolution_clock::now();

   // calc abcd (needed?)
   // clipper::HKL_data<clipper::data32::ABCD> abcd(hkls);
   // abcd.compute(phiw, clipper::data32::Compute_abcd_from_phifom());

   mol->DeleteSelection(hndl);

   // now calc R and R-free
   std::vector<double> params(n_param, 1.0);
   clipper::BasisFn_spline basisfn(fobs, n_param, 1.0);
   clipper::TargetFn_scaleF1F2<clipper::data32::F_phi,clipper::data32::F_sigF> targetfn(fc, fobs);
   clipper::ResolutionFn rfn(hkls, basisfn, targetfn, params);
   double r1w, f1w, r1f, f1f, Fo, Fc;
   r1w = f1w = r1f = f1f = 0.0;
   for (HRI ih = fobs.first(); !ih.last(); ih.next()) {
      if (!fobs[ih].missing()) {
         Fo = fobs[ih].f();
         Fc = sqrt( rfn.f(ih) ) * fc[ih].f();
         if (free[ih].flag() == freeflag) {
            r1f += fabs(Fo - Fc);
            f1f += Fo;
         } else {
            r1w += fabs(Fo - Fc);
            f1w += Fo;
         }
      }
   }
   auto tp_8 = std::chrono::high_resolution_clock::now();
   r1f /= clipper::Util::max(f1f, 0.1);
   r1w /= clipper::Util::max(f1w, 0.1);
   sfcalc_genmap_stats_t::loc_table_t loct;
   for (int i = 0; i <= 20; i++) {
      double s = hkls.resolution().invresolsq_limit()*double(i)/20.0;
      // printf("%6.3f %12.3f %12.3f\n", s, basisfn.f_s(s,sfw.params_scale()), basisfn.f_s(s,sfw.params_error()));
      sfcalc_genmap_stats_t::loc_table_t::loc_table_item_t item(s, basisfn.f_s(s,sfw.params_scale()), basisfn.f_s(s,sfw.params_error()));
      loct.add(item);
   }

   auto tp_9 = std::chrono::high_resolution_clock::now();

   if (false)
      std::cout << "debug:: capturing rfactors " << r1w << " and " << r1f << std::endl;
   sfcgs = sfcalc_genmap_stats_t(r1w, r1f, bulkfrc, bulkscl, sfw.params_scale().size(), loct);
   if (false)
      std::cout << "\n R-factor      : " << r1w << "\n Free R-factor : " << r1f << "\n";
   auto tp_10 = std::chrono::high_resolution_clock::now();

   // now make a map

   // test
   // xmap_p->fft_from(fc);
   xmap_2fofc_p->fft_from(f_best);
   auto tp_11 = std::chrono::high_resolution_clock::now();
   xmap_fofc_p->fft_from(f_diff);
   auto tp_12 = std::chrono::high_resolution_clock::now();

   // DIAGNOSTIC OUTPUT
   if (false) {
      std::cout << "\n Bulk Correction Volume: " << bulkfrc;
      std::cout << "\n Bulk Correction Factor: " << bulkscl << "\n";
      std::cout << "\nNumber of spline params: " << sfw.params_scale().size() << "\n";
      clipper::BasisFn_spline basisfn_table( hkls, sfw.params_scale().size(), 1.0 );
      printf("\n $TABLE: Sigmaa statistics :\n $GRAPHS:scale vs resolution:N:1,2:\n        :lack of closure vs resolution:N:1,3:\n $$\n 1/resol^2   scale   lack_of_closure $$\n $$\n");
      for (int i = 0; i <= 20; i++) {
         double s = hkls.resolution().invresolsq_limit()*double(i)/20.0;
         printf("%6.3f %12.3f %12.3f\n", s, basisfn_table.f_s(s,sfw.params_scale()), basisfn_table.f_s(s,sfw.params_error()));
      }
      printf(" $$\n");
   }

   bool show_timings = true;
   if (show_timings) {
      auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
      auto d21 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_2 - tp_1).count();
      auto d32 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_3 - tp_2).count();
      auto d43 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_4 - tp_3).count();
      auto d54 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_5 - tp_4).count();
      auto d65 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_6 - tp_5).count();
      auto d76 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_7 - tp_6).count();
      auto d87 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_8 - tp_7).count();
      auto d98 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_9 - tp_8).count();
      auto d109 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_10 - tp_9).count();
      auto d1110 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_11 - tp_10).count();
      auto d1211 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_12 - tp_11).count();

      std::cout << "sfs calc timings: d10 " << d10 << " d21 " << d21 << " d32 " << d32 << " d43 " << d43 << " d54 " << d54
                << " d65 " << d65 << " d76 " << d76 << " " << d87 << " d98 " << d98 << " d10-9 " << d109
                << " d11-10 " << d1110 << " d12-11 " << d1211 << "\n";
   }

   return sfcgs;

}
