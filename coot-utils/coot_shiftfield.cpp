// Clipper app to perform shift field refinement
/* Copyright 2018 Kevin Cowtan & University of York all rights reserved */


#include "coot_shiftfield.h"
#include "shiftfield.h"

void
test_import_minimol(clipper::MMDBfile* mfile, clipper::MiniMol& minimol, const int hnd ) {

   // clear the model
   minimol = clipper::MiniMol( mfile->spacegroup(), mfile->cell() );

   // make atom, residue, chain selections
   int h_atm = mfile->NewSelection();
   int h_res = mfile->NewSelection();
   int h_chn = mfile->NewSelection();
   int h_mod = mfile->NewSelection();
   if ( hnd < 0 ) {
      mfile->SelectAtoms( h_atm, 0, 0, ::mmdb::SKEY_NEW );
   } else {
      mfile->Select( h_atm, ::mmdb::STYPE_ATOM, hnd, ::mmdb::SKEY_NEW );
   }
   mfile->Select( h_res, ::mmdb::STYPE_RESIDUE, h_atm, ::mmdb::SKEY_NEW );
   mfile->Select( h_chn, ::mmdb::STYPE_CHAIN,   h_atm, ::mmdb::SKEY_NEW );
   mfile->Select( h_mod, ::mmdb::STYPE_MODEL,   h_atm, ::mmdb::SKEY_NEW );

   std::cout << "debug:: h_atm " << h_atm << " h_res " << h_res << " h_chn " << h_chn << " h_mod " << h_mod << std::endl;

   // now import objects
   char txt[256];
   clipper::MModel& mol = minimol.model();
   mmdb::Model *p_mod = mfile->GetModel(1);
   if ( p_mod != NULL ) {
      for ( int c = 0; c < p_mod->GetNumberOfChains(); c++ ) {
         mmdb::Chain *p_chn = p_mod->GetChain(c);
         if ( p_chn != NULL ) if ( p_chn->isInSelection( h_chn ) ) {
               // import the chain
               clipper::MPolymer pol;
               for ( int r = 0; r < p_chn->GetNumberOfResidues(); r++ ) {
                  mmdb::Residue *p_res = p_chn->GetResidue(r);
                  if ( p_chn != NULL ) if ( p_res->isInSelection( h_res ) ) {
                        // import the residue
                        clipper::MMonomer mon;
                        for ( int a = 0; a < p_res->GetNumberOfAtoms(); a++ ) {
                           mmdb::Atom *p_atm = p_res->GetAtom(a);
                           if ( p_atm != NULL ) if ( p_atm->isInSelection( h_atm ) )
                                                   if ( !p_atm->Ter ) {
                                                      // import the atom
                                                      clipper::MAtom atm( clipper::Atom::null() );
                                                      atm.set_name( p_atm->GetAtomName(), p_atm->altLoc );
                                                      atm.set_element( p_atm->element );
                                                      if ( p_atm->WhatIsSet & ::mmdb::ASET_Coordinates )
                                                         atm.set_coord_orth(clipper::Coord_orth( p_atm->x, p_atm->y, p_atm->z ) );
                                                      if ( p_atm->WhatIsSet & ::mmdb::ASET_Occupancy )
                                                         atm.set_occupancy( p_atm->occupancy );
                                                      if ( p_atm->WhatIsSet & ::mmdb::ASET_tempFactor )
                                                         atm.set_u_iso( clipper::Util::b2u( p_atm->tempFactor ) );
                                                      if ( p_atm->WhatIsSet & ::mmdb::ASET_Anis_tFac )
                                                         atm.set_u_aniso_orth(
                                                                              clipper::U_aniso_orth( p_atm->u11, p_atm->u22, p_atm->u33,
                                                                                                     p_atm->u12, p_atm->u13, p_atm->u23 ) );
                                                      p_atm->GetAtomID( txt );
                                                      atm.set_property("CID",clipper::Property<clipper::String>(clipper::String(txt)));
                                                      if ( p_atm->altLoc[0] != '\0' )
                                                         atm.set_property("AltConf",
                                                                          clipper::Property<clipper::String>(clipper::String(p_atm->altLoc)));
                                                      mon.insert( atm );  // store the atom
                                                   }
                        }
                        mon.set_seqnum( p_res->GetSeqNum(), clipper::String(p_res->GetInsCode()) );
                        mon.set_type( p_res->GetResName() );
                        p_res->GetResidueID( txt );
                        mon.set_property("CID",clipper::Property<clipper::String>(clipper::String(txt)));
                        pol.insert( mon );  // store the residue
                     }
               }
               pol.set_id( p_chn->GetChainID() );
               p_chn->GetChainID( txt );
               pol.set_property("CID",clipper::Property<clipper::String>(clipper::String(txt)));
               mol.insert( pol );  // store the chain
            }
      }
      p_mod->GetModelID( txt );
      mol.set_property("CID",clipper::Property<clipper::String>(clipper::String(txt)));
   }

   // clean up
   mfile->DeleteSelection( h_atm );
   mfile->DeleteSelection( h_res );
   mfile->DeleteSelection( h_chn );
   mfile->DeleteSelection( h_mod );
}



// update the temperature-factors of the atoms in mol
void
coot::shift_field_b_factor_refinement(const clipper::HKL_data< clipper::datatypes::F_sigF<float> > &fo0,
                                      const clipper::HKL_data<clipper::data32::Flag> &free,
                                      mmdb::Manager *mol, int ncyc)
{
  double radscl = 3.0;
  int freeflag = 0;
  int filter = 2;
  int n_refln = 1000;
  int n_param = 10;

  typedef clipper::HKL_data_base::HKL_reference_index HRI;
  clipper::Spacegroup spgr0 = fo0.spacegroup();
  clipper::Cell       cell  = fo0.cell();
  clipper::Resolution reso  = fo0.resolution();
  double rcyc = reso.limit();

  // get a list of all the atoms
  clipper::MMDBfile* mfile = static_cast<clipper::MMDBfile*>(mol);
  clipper::MiniMol mmol;

  // Either of these cause a crash when the molecule is imported using gemmi
  // test_import_minimol(mfile, mmol, -1);
  mfile->import_minimol( mmol );

  for ( int cyc = 0; cyc < ncyc; cyc++ ) {
    double radcyc = radscl * rcyc;
    std::cout << std::endl << "Shiftfield Cycle: " << cyc+1 << "  Resolution: " << rcyc << "  Radius: " << radcyc << std::endl;

    // truncate resolution
    clipper::Resolution rescyc( rcyc );
    clipper::HKL_sampling hklsam( cell, rescyc );
    clipper::HKL_data<clipper::data32::F_sigF> fo_all( spgr0, cell, hklsam ), fo( spgr0, cell, hklsam );
    for ( HRI ih = fo.first(); !ih.last(); ih.next() ) fo_all[ih] = fo0[ih.hkl()];
    for ( HRI ih = fo.first(); !ih.last(); ih.next() ) if (free[ih.hkl()].flag() != freeflag) fo[ih] = fo0[ih.hkl()];

    // calculate structure factors
    clipper::Atom_list atoms = mmol.atom_list();
    clipper::HKL_data<clipper::data32::F_phi> fc( fo );
    clipper::SFcalc_obs_bulk<float> sfcb;
    sfcb( fc, fo, atoms );
    //std::cerr << "Done sfcalc" << std::endl;

    // now calc R and R-free
    std::vector<double> params( 2, 1.0 );
    clipper::BasisFn_log_gaussian basisfn; // we're using a Gaussian because using a spline can result in negative square roots {see below}
    clipper::TargetFn_scaleLogF1F2<clipper::data32::F_phi,clipper::data32::F_sigF> targetfn( fc, fo );

    clipper::ResolutionFn rfn( fo.hkl_info(), basisfn, targetfn, params);
    double r1w, f1w, r1f, f1f, Fo, Fc;
    r1w = f1w = r1f = f1f = 0.0;
    for ( HRI ih = fo_all.first(); !ih.last(); ih.next() )
      if ( !fo_all[ih].missing() ) {
        float sfac = exp(0.5*basisfn.f(ih.hkl(),cell,rfn.params()));
        Fo = fo_all[ih].f();
        Fc = sfac * fc[ih].f();

        if ( free[ih].flag() == freeflag ) {
                r1f += fabs( Fo - Fc );
                f1f += Fo;
        } else {
                r1w += fabs( Fo - Fc );
                f1w += Fo;
        }
      }

    r1w /= clipper::Util::max( f1w, 0.1 );
    r1f /= clipper::Util::max( f1f, 0.1 );

    std::cout << "R-factor      : " << r1w << std::endl
              << "Free R-factor : " << r1f << std::endl;

    // now do sigmaa calc
    clipper::HKL_data<clipper::data32::F_phi>   fb( fo ), fd( fo );
    clipper::HKL_data<clipper::data32::Phi_fom> phiw( fo );
    clipper::HKL_data<clipper::data32::Flag>    flag( fo );
    for ( HRI ih = flag.first(); !ih.last(); ih.next() )
      flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;

    // do sigmaa calc
    clipper::SFweight_spline<float> sfw( n_refln, n_param );
    sfw( fb, fd, phiw, fo, fc, flag );
    //std::cerr << "Done sigmaa" << std::endl;

    // expand to P1
    clipper::Spacegroup spgr1 = clipper::Spacegroup::p1();
    clipper::HKL_data<clipper::data32::F_phi> fphi( spgr1, cell, hklsam );
    clipper::HKL_data<clipper::data32::F_phi> dphi( fphi );
    for ( HRI ih = fphi.first(); !ih.last(); ih.next() ) {
      fphi[ih] = fc[ih.hkl()];
      dphi[ih] = fd[ih.hkl()];
    }
    std::cout << "Reflections: " << fo0.base_hkl_info().num_reflections() << " P1: " << fo.base_hkl_info().num_reflections() << std::endl;

    //for ( HRI ih = fphi.first(); !ih.last(); ih.next() ) std::cout << ih.hkl().format() << fphi[ih].f() << " " << fphi[ih].phi() << std::endl;
    //for ( HRI ih = dphi.first(); !ih.last(); ih.next() ) std::cout << ih.hkl().format() << dphi[ih].f() << " " << dphi[ih].phi() << std::endl;

    // apply U value
    //fphi.compute( fphi, clipper::data32::Compute_scale_u_iso_fphi(1.0,-uvalue) );

    // make grid if necessary
    clipper::Grid_sampling grid( spgr1, cell, rescyc );

    // make xmaps
    clipper::Xmap<float> cmap( spgr1, cell, grid );
    clipper::Xmap<float> dmap( spgr1, cell, grid );
    clipper::Xmap<float> mmap( spgr1, cell, grid );
    clipper::Xmap<float> x1map( spgr1, cell, grid );

    cmap.fft_from( fphi );
    dmap.fft_from( dphi );
    mmap = 1;

    // isotropic u refinement
    std::cout << "REFINE U ISO" << std::endl;

    // make shift field
    Shift_field_refine::shift_field_u_iso( cmap, dmap, mmap, x1map, radcyc, filter );

    std::cout << "DONE REFINE U ISO" << std::endl;

    // read pdb and update
    for ( int p = 0; p < mmol.size(); p++ )
      for ( int m = 0; m < mmol[p].size(); m++ )
        for ( int a = 0; a < mmol[p][m].size(); a++ ) {
          const clipper::Coord_frac cf = mmol[p][m][a].coord_orth().coord_frac(cell);
          const float du = 1.0*x1map.interp<clipper::Interp_cubic>( cf );
          mmol[p][m][a].set_u_iso( mmol[p][m][a].u_iso() - du );
          // std::cout << p << " " << m << " " << a << " " << du << std::endl;
        }
  }

  mfile->export_minimol( mmol );
}


// update the coordinates of the atoms in mol
void
coot::shift_field_xyz_refinement(const clipper::HKL_data< clipper::datatypes::F_sigF<float> > &fo0,
                                 const clipper::HKL_data<clipper::data32::Flag> &free,
                                 mmdb::Manager *mol,
                                 float resolution) {

}


