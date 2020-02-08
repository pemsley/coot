// Perform structure factor calculation
/* Copyright 2003-2004 Kevin Cowtan & University of York all rights reserved */
/* Cootenized by Paul Emsley */


#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-mmdb.h>

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
   std::cout << "\n R-factor      : " << r1w << "\n Free R-factor : " << r1f << "\n";

   // now make a map

   // test
   // xmap_p->fft_from(fc);
   xmap_p->fft_from(f_diff);

   // DIAGNOSTIC OUTPUT
   if (true) {
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
