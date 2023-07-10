// Clipper sequins
/* Copyright 2007 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include "simulate-lib.h"
#include "buccaneer-sequence.h"

extern "C" {
#include <stdlib.h>
}
//#include "buccaneer-correct.h"
//#include "buccaneer-build.h"


int main( int argc, char** argv )
{
  CCP4Program prog( "csequins", "0.1.1", "$Date: 2007/09/11" );

  std::cout << "\nCopyright 2007 Kevin Cowtan and University of York. All rights reserved.\n\n";
  std::cout << "Please reference:\n Cowtan K. (2006) Acta Cryst. D62, 1002-1011.\n\n";

  // defaults
  clipper::String title;
  clipper::String ipmtz_ref = "NONE";
  clipper::String ipmtz_wrk = "NONE";
  clipper::String ippdb_ref = "NONE";
  clipper::String ippdb_wrk = "NONE";
  clipper::String ipcol_ref_fo = "/*/*/FP";
  clipper::String ipcol_ref_hl = "/*/*/FC";
  clipper::String ipcol_wrk_fo = "NONE";
  clipper::String ipcol_wrk_hl = "NONE";
  clipper::String opmtz = "NONE";
  clipper::String oppdb = "NONE";
  clipper::String chnid = "NONE";
  double res_in = 1.0;         // Resolution limit
  bool scomit = false;
  bool correl = false;
  double side_tgt_rad = 5.5;
  int n_refln = 1000;
  int n_param = 20;
  int verbose = 0;

  // command input
  CCP4CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if        ( args[arg] == "-title" ) {
      if ( ++arg < args.size() ) title = args[arg];
    } else if ( args[arg] == "-mtzin-ref" ) {
      if ( ++arg < args.size() ) ipmtz_ref = args[arg];
    } else if ( args[arg] == "-mtzin-wrk" ) {
      if ( ++arg < args.size() ) ipmtz_wrk = args[arg];
    } else if ( args[arg] == "-pdbin-ref" ) {
      if ( ++arg < args.size() ) ippdb_ref = args[arg];
    } else if ( args[arg] == "-pdbin-wrk" ) {
      if ( ++arg < args.size() ) ippdb_wrk = args[arg];
    } else if ( args[arg] == "-colin-ref-fo" ) {
      if ( ++arg < args.size() ) ipcol_ref_fo = args[arg];
    } else if ( args[arg] == "-colin-ref-hl" ) {
      if ( ++arg < args.size() ) ipcol_ref_hl = args[arg];
    } else if ( args[arg] == "-colin-wrk-fo" ) {
      if ( ++arg < args.size() ) ipcol_wrk_fo = args[arg];
    } else if ( args[arg] == "-colin-wrk-hl" ) {
      if ( ++arg < args.size() ) ipcol_wrk_hl = args[arg];
    } else if ( args[arg] == "-mtzout-wrk" ) {
      if ( ++arg < args.size() ) opmtz = args[arg];
    } else if ( args[arg] == "-pdbout-wrk" ) {
      if ( ++arg < args.size() ) oppdb = args[arg];
    } else if ( args[arg] == "-resolution" ) {
      if ( ++arg < args.size() ) res_in = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-side-chain-likelihood-radius" ) {
      if ( ++arg < args.size() ) side_tgt_rad = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-chain-id" ) {
      if ( ++arg < args.size() ) chnid = args[arg];
    } else if ( args[arg] == "-side-chain-omit" ) {
      scomit = true;
    } else if ( args[arg] == "-correlation-mode" ) {
      correl = true;
    } else if ( args[arg] == "-verbose" ) {
      if ( ++arg < args.size() ) verbose = clipper::String(args[arg]).i();
    } else {
      std::cout << "\nUnrecognized:\t" << args[arg] << std::endl;
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "\nUsage: csequins\n\t-mtzin-ref <filename>\t\tCOMPULSORY\n\t-pdbin-ref <filename>\t\tCOMPULSORY\n\t-mtzin-wrk <filename>\t\tCOMPULSORY\n\t-pdbin-wrk <filename>\t\tCOMPULSORY\n\t-mtzout-wrk <filename>\n\t-pdbout-wrk <filename>\n\t-colin-ref-fo <colpath>\n\t-colin-ref-hl <colpath>\n\t-colin-wrk-fo <colpath>\t\tCOMPULSORY\n\t-colin-wrk-hl <colpath>\n\t-resolution <resolution/A>\n\t-side-chain-likelihood-radius <radius/A>\n\t-chain-id <id>\n\t-side-chain-omit\n\t-correlation-mode\nValidate sequence against structure factors.\nAn input pdb and mtz are required for the reference structure, and \nan input mtz file for the work structure.\nExperimental phases are used if provided, otherwise phases are calculated\nfrom the model omitting the side chains.\n";
    exit(1);
  }

  // other initialisations
  clipper::Resolution resol;
  clipper::CCP4MTZfile mtzfile;
  using clipper::data32::F_sigF;  using clipper::data32::F_phi;
  using clipper::data32::Phi_fom; using clipper::data32::Flag;
  using clipper::data32::ABCD;
  const int mmdbflags = ::mmdb::MMDBF_IgnoreBlankLines | ::mmdb::MMDBF_IgnoreDuplSeqNum | ::mmdb::MMDBF_IgnoreNonCoorPDBErrors | ::mmdb::MMDBF_IgnoreRemarks;
 
  // Get resolution for calculation
  mtzfile.open_read( ipmtz_ref );
  double res_ref = clipper::Util::max( mtzfile.resolution().limit(), res_in );
  mtzfile.close_read();
  mtzfile.open_read( ipmtz_wrk );
  double res_wrk = clipper::Util::max( mtzfile.resolution().limit(), res_in );
  mtzfile.close_read();
  resol = clipper::Resolution( clipper::Util::max( res_ref, res_wrk ) );
  if ( res_ref > res_wrk ) std::cout << "\nWARNING: resolution of work structure truncated to reference:\n Ref: " << res_ref << " Wrk: " << res_wrk << std::endl;

  // Get reference reflection data
  clipper::HKL_info hkls_ref;
  mtzfile.open_read( ipmtz_ref );
  hkls_ref.init( mtzfile.spacegroup(), mtzfile.cell(), resol, true );
  clipper::HKL_data<F_sigF> ref_f( hkls_ref );
  clipper::HKL_data<ABCD> ref_hl( hkls_ref );
  mtzfile.import_hkl_data( ref_f, ipcol_ref_fo );
  mtzfile.import_hkl_data( ref_hl, ipcol_ref_hl );
  mtzfile.close_read();

  // Get work reflection data
  clipper::HKL_info hkls_wrk;
  mtzfile.open_read( ipmtz_wrk );
  hkls_wrk.init( mtzfile.spacegroup(), mtzfile.cell(), resol, true );
  clipper::HKL_data<F_sigF> wrk_f ( hkls_wrk );
  clipper::HKL_data<F_phi>  wrk_fc( hkls_wrk );
  clipper::HKL_data<ABCD>   wrk_hl( hkls_wrk );
  mtzfile.import_hkl_data( wrk_f, ipcol_wrk_fo );
  if ( ipcol_wrk_hl != "NONE" ) mtzfile.import_hkl_data( wrk_hl, ipcol_wrk_hl );
  mtzfile.close_read();

  // Get reference model
  clipper::MiniMol mol_ref;
  clipper::MMDBfile mmdb_ref;
  mmdb_ref.SetFlag( mmdbflags );
  mmdb_ref.read_file( ippdb_ref );
  mmdb_ref.import_minimol( mol_ref );

  // Get work model (optional)
  clipper::MiniMol mol_wrk, mol_tmp;
  mol_wrk.init( hkls_wrk.spacegroup(), hkls_wrk.cell() );
  clipper::MMDBfile mmdb_wrk;
  mmdb_wrk.SetFlag( mmdbflags );
  mmdb_wrk.read_file( ippdb_wrk );
  mmdb_wrk.import_minimol( mol_tmp );
  mol_wrk.copy( mol_tmp, clipper::MM::COPY_MPC );

  // get map coefficients
  // if no HL coeffs, generate them
  if ( ipcol_wrk_hl == "NONE" ) {
    // optionally delete side chain atoms
    mol_tmp = mol_wrk;
    if ( scomit ) {
      for ( int c = 0; c < mol_tmp.size(); c++ ) {
	for ( int r = 0; r < mol_tmp[c].size(); r++ ) {
	  if ( ProteinTools::residue_index_3( mol_tmp[c][r].type() ) >= 0 ) {
	    clipper::MMonomer mm;
	    mm.set_id( mol_tmp[c][r].id() );
	    mm.set_type( mol_tmp[c][r].type() );
	    for ( int a = 0; a < mol_tmp[c][r].size(); a++ )
	      if ( mol_tmp[c][r][a].name() == " N  " ||
		   mol_tmp[c][r][a].name() == " CA " || 
		   mol_tmp[c][r][a].name() == " C  " || 
		   mol_tmp[c][r][a].name() == " O  " )
		mm.insert( mol_tmp[c][r][a] );
	    mol_tmp[c][r] = mm;
	  }
	}
      }
    }
    clipper::Atom_list atoms = mol_tmp.atom_list();
    // calculate structure factors
    clipper::HKL_data<F_phi> fc( hkls_wrk );
    clipper::SFcalc_obs_bulk<float> sfcb;
    sfcb( fc, wrk_f, atoms );
    // do anisotropic scaling
    clipper::SFscale_aniso<float> sfscl;
    sfscl( wrk_f, fc );
    // do sigmaa calc
    clipper::HKL_data<F_phi>   fd( hkls_wrk );
    clipper::HKL_data<Phi_fom> phiw( hkls_wrk );
    clipper::HKL_data<Flag>    flag( hkls_wrk );
    typedef clipper::HKL_data_base::HKL_reference_index HRI;
    for ( HRI ih = flag.first(); !ih.last(); ih.next() )
      if ( !wrk_f[ih].missing() )
	flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
      else
	flag[ih].flag() = clipper::SFweight_spline<float>::NONE;
    // do sigmaa calc
    clipper::SFweight_spline<float> sfw( n_refln, n_param );
    sfw( wrk_fc, fd, phiw, wrk_f, fc, flag );
    // calc abcd
    wrk_hl.compute( phiw, clipper::data32::Compute_abcd_from_phifom() );
    // optional output
    if ( opmtz != "NONE" ) {
      mtzfile.open_append( ipmtz_wrk, opmtz );
      mtzfile.export_hkl_data( wrk_fc, "/*/*/sideomit" );
      mtzfile.export_hkl_data( wrk_hl, "/*/*/sideomit" );
      mtzfile.close_append();
    }
    if ( oppdb != "NONE" ) {
      clipper::MMDBfile mmdb;
      mmdb.export_minimol( mol_tmp );
      mmdb.write_file( oppdb );
    }
  } else {
    // get map coefficients from input HL coeffs
    clipper::HKL_data<Phi_fom> phiw( hkls_wrk );
    phiw.compute( wrk_hl, clipper::data32::Compute_phifom_from_abcd() );
    wrk_fc.compute( wrk_f, phiw, clipper::data32::Compute_fphi_from_fsigf_phifom() );
  }

  // DO INITIAL MAP SIMULATION
  clipper::HKL_data<F_sigF> sim_f( hkls_ref );
  clipper::HKL_data<ABCD> sim_hl( hkls_ref );
  MapSimulate mapsim( 100, 20 );
  mapsim( sim_f, sim_hl, ref_f, ref_hl, wrk_f, wrk_hl );

  // make llk target objects
  LLK_map_target::TYPE tgttyp =
    correl ? LLK_map_target::CORREL : LLK_map_target::NORMAL;
  std::vector<LLK_map_target> llkcls( 20 );
  for ( int t = 0; t < 20; t++ ) llkcls[t].init( side_tgt_rad, 0.5, tgttyp );

  // STAGE 1: Calculate target from reference data

  {
    // reference map
    clipper::HKL_data<F_phi> fphi( hkls_ref );
    clipper::HKL_data<Phi_fom> phiw( hkls_ref );
    phiw.compute( sim_hl, clipper::data32::Compute_phifom_from_abcd() );
    fphi.compute( sim_f, phiw, clipper::data32::Compute_fphi_from_fsigf_phifom() );
    clipper::Grid_sampling grid( hkls_ref.spacegroup(), hkls_ref.cell(), hkls_ref.resolution() );
    clipper::Xmap<float> xref( hkls_ref.spacegroup(), hkls_ref.cell(), grid );
    xref.fft_from( fphi );

    // prepare llk targets
    typedef clipper::MMonomer MM;
    for ( int chn = 0; chn < mol_ref.size(); chn++ )
      for ( int res = 1; res < mol_ref[chn].size()-1; res++ ) {
	const clipper::MMonomer& mm1 = mol_ref[chn][res-1];
	const clipper::MMonomer& mm2 = mol_ref[chn][res  ];
	const clipper::MMonomer& mm3 = mol_ref[chn][res+1];
	bool b1 = MM::protein_peptide_bond( mm1, mm2 );
	bool b2 = MM::protein_peptide_bond( mm2, mm3 );
	int index_n  = mm2.lookup( " N  ", clipper::MM::ANY );
	int index_ca = mm2.lookup( " CA ", clipper::MM::ANY );
	int index_c  = mm2.lookup( " C  ", clipper::MM::ANY );
	if ( b1 && b2 && index_ca >= 0 && index_c >= 0 && index_n >= 0 ) {
	  Ca_group ca( mm2[index_n ].coord_orth(),
		       mm2[index_ca].coord_orth(),
		       mm2[index_c ].coord_orth() );
	  // side chain target:
	  int type = ProteinTools::residue_index_3( mm2.type() );
	  if ( type >= 0 )
	    llkcls[type].accumulate( xref, ca.rtop_beta_carbon() );
	}
      }

    for ( int t = 0; t < llkcls.size(); t++ ) llkcls[t].prep_llk();
  }


  // STAGE 2: Apply target to work data

  {
    // calc work map
    clipper::Grid_sampling grid( hkls_wrk.spacegroup(), hkls_wrk.cell(), hkls_wrk.resolution() );
    clipper::Xmap<float> xwrk( hkls_wrk.spacegroup(), hkls_wrk.cell(), grid );
    xwrk.fft_from( wrk_fc );

    // extract the necessary bits of the likelihood targets
    std::vector<LLK_map_target::Sampled> llksmp( llkcls.size() );
    for ( int t = 0; t < llkcls.size(); t++ )
      llksmp[t] = llkcls[t].sampled();

    // sequence
    for ( int c = 0; c < mol_wrk.size(); c++ ) {
      // extract input sequence
      clipper::String seq = "";
      clipper::MPolymer mp;
      for ( int r = 0; r < mol_wrk[c].size(); r++ ) {
	int t = ProteinTools::residue_index_3( mol_wrk[c][r].type() );
	if ( mol_wrk[c][r].type().length() == 3 && t >= 0 ) { 
	  seq += ProteinTools::residue_code_1( t );
	  mp.insert( mol_wrk[c][r] );
	}
      }
      if ( seq.length() > 1 && ( chnid==mol_wrk[c].id() || chnid=="NONE" ) ) {
	clipper::MMoleculeSequence seq_mol;
	clipper::MPolymerSequence  seq_chn;
	seq_chn.set_id( "TEST" );
	seq_chn.set_sequence( seq );
	seq_mol.insert( seq_chn );
	// reapply it
	Ca_sequence::prepare_scores( mp, xwrk, llksmp );
	Score_list<clipper::String> result = 
	  Ca_sequence::sequence_chain( mp, seq_mol );
	clipper::String s0 = ProteinTools::chain_sequence( mp );
        clipper::String s1 = "";
        double prob = 0.0;
	if ( result.size() > 1 ) {
          s1 = result[0];
	  prob = Ca_sequence::phi_approx(result.score(1)-result.score(0));
        }
	// basic output
	std::cout << std::endl << "------------------------------------------------------------------------" << std::endl;
	std::cout << std::endl << "CHAIN: " << mol_wrk[c].id() << std::endl;
	std::cout << "Original: " << s0 << std::endl;
	std::cout << "Modified: " << s1 << std::endl;
	std::cout << "Confidence: " << prob << std::endl;
	// fancy output
	if ( prob < 0.95 ) {
	  std::cout << std::endl << "The data is insufficient to either confirm or deny the sequence for this chain." << std::endl;
	} else {
	  if ( s0.length() == s1.length() ) {
	    int n, nmtch, nmiss, nshft;
	    n =  s0.length();
	    nmtch = nmiss = nshft = 0;
	    for ( int i = 0; i < n; i++ ) {
	      if ( s0[i] == s1[i] ) nmtch++;
	      else { 
		if ( isupper( s1[i] ) ) nmiss++;
		if ( s1[i] == '+' || s1[i] == '-' ) nshft++;
	      }
	    }
	    clipper::String conf;
	    if ( prob > 0.999 ) conf = "strongly";
	    else                conf = "weakly";
	    if ( nmtch > (3*n)/4 && nmiss == 0 )
	      std::cout << std::endl << "The supplied data "+conf+" supports the assigned sequence." << std::endl;
	    else if ( nmtch > n/4 && nmiss > 6 && nshft > 0 )
	      std::cout << std::endl << "WARNING: The supplied data suggests a partial sequence register error." << std::endl;
	    else if ( nmtch < n/4 && nmiss > n/4 )
	      std::cout << std::endl << "WARNING: The supplied data suggests a tracing error or major sequence register error." << std::endl;
	    else
	      std::cout << std::endl << "Manual analysis is required." << std::endl;
	  } else {
	    std::cout << std::endl << "Manual analysis is required." << std::endl;
	  }
	}
      }
    }
    std::cout << std::endl << "------------------------------------------------------------------------" << std::endl;
  }
}

