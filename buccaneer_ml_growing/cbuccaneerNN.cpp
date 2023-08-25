// Clipper buccaneer
/* Copyright 2002-2008 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include "simulate-lib.h"
#include "buccaneer-prep.h"
#include "buccaneer-find.h"
#include "buccaneer-grow.h"
#include "buccaneer-join.h"
#include "buccaneer-link.h"
#include "buccaneer-sequence.h"
#include "buccaneer-correct.h"
#include "buccaneer-filter.h"
#include "buccaneer-ncsbuild.h"
#include "buccaneer-prune.h"
#include "buccaneer-build.h"
#include "buccaneer-merge.h"
#include "buccaneer-known.h"
#include "buccaneer-tidy.h"
#include "buccaneer-util.h"
#include "buccaneer-NN-Select.h"
#include "buccaneer-NN-build.h"
#include "buccaneer-NN-build-cyc.h"

extern "C" {
  #include <stdlib.h>
}


int main( int argc, char** argv )
{
  CCP4Program prog( "cbuccaneer", "1.7.0", "$Date: 2023/06/12" );
  prog.set_termination_message( "Failed" );

  std::cout << std::endl << "Copyright 2002-2020 Kevin Cowtan and University of York." << std::endl << std::endl;
  prog.summary_beg();
  std::cout << "$TEXT:Reference: $$ Please reference $$" << std::endl << std::endl << " 'Fitting molecular fragments into electron density'" << std::endl << " Cowtan K. (2008) Acta Cryst. D64, 83-89." << std::endl << std::endl << "$$" << std::endl;
  std::cout << "$TEXT:Reference: $$ Please reference $$" << std::endl << std::endl << " 'The Buccaneer software for automated model building'" << std::endl << " Cowtan K. (2006) Acta Cryst. D62, 1002-1011." << std::endl << std::endl << "$$" << std::endl;
  prog.summary_end();

  // defaults
  clipper::String title;
  clipper::String ipmtz_ref = "NONE";
  clipper::String ipmtz_wrk = "NONE";
  clipper::String ippdb_ref = "NONE";
  clipper::String ippdb_wrk = "NONE";
  clipper::String ippdb_mr  = "NONE";
  clipper::String ippdb_seq = "NONE";
  clipper::String ipseq_wrk = "NONE";
  clipper::String ipcol_ref_fo = "/*/*/FP";
  clipper::String ipcol_ref_hl = "/*/*/FC";
  clipper::String ipcol_wrk_fo = "NONE";
  clipper::String ipcol_wrk_hl = "NONE";
  clipper::String ipcol_wrk_pw = "NONE";
  clipper::String ipcol_wrk_fc = "NONE";
  clipper::String ipcol_wrk_fr = "NONE";
  clipper::String oppdb = "buccaneer.pdb";
  clipper::String opmap = "NONE";
  clipper::String opxml = "NONE";
  clipper::String newresname = "UNK";
  clipper::String newrestype = "ALA";
  double res_in = 1.0;         // Resolution limit
  int ncyc = 3;
  int ncpu = 0;
  int nfrag  = 500;
  int nfragr = 20;
  int freerindex = 0;
  int modelindex = 0;
  int randomseed = -1;
  bool merge = false;  // multimodel merge
  bool find  = false;  // calculation steps
  bool grow  = false;
  bool join  = false;
  bool link  = false;
  bool seqnc = false;
  bool corct = false;
  bool filtr = false;
  bool ncsbd = false;
  bool prune = false;
  bool build = false;
  bool tidy  = false;
  bool fast   = false;  // further options
  bool semet  = false;
  bool doanis = false;
  bool optemp = false;
  bool fixpos = false;
  bool correl = false;
  bool nocorrel = false;
  clipper::MMDBManager::TYPE cifflag = clipper::MMDBManager::Default;
  double nprad = -1.0;
  double main_tgt_rad = 4.0;
  double side_tgt_rad = 5.5;
  double seq_rel = 0.95;
  double moffset = 0.0;
  double model_filter_sig = 3.0;
  double mr_filter_sig = 3.0;
  bool model_filter = false;
  bool mr_model = false;
  bool mr_model_filter = false;
  bool mr_model_seed = false;
  Ca_prep::Rama_flt rama_flt = Ca_prep::rama_flt_all;
  bool rama_grow = false;
  std::vector<std::pair<clipper::String,double> > known_ids;
  int verbose = 0;

  //NN
  double model_num=10;
  bool nnselect=false;
  clipper::String threshold_select_method = "NONE";
  bool nn_confirmation_cycles=true;
  int confirmation_cycles_num=-1;// -1 meaning here the number of confirmation cycles will be the number of remaining building cycles


  // command input
  CCP4CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    std::string key = args[arg];
    if        ( key == "-title" ) {
      if ( ++arg < args.size() ) title = args[arg];
    } else if ( key == "-mtzin-ref" ) {
      if ( ++arg < args.size() ) ipmtz_ref = args[arg];
    } else if ( key == "-pdbin-ref" ) {
      if ( ++arg < args.size() ) ippdb_ref = args[arg];
    } else if ( key == "-mtzin"        || key == "-mtzin-wrk" ) {
      if ( ++arg < args.size() ) ipmtz_wrk = args[arg];
    } else if ( key == "-seqin"        || key == "-seqin-wrk" ) {
      if ( ++arg < args.size() ) ipseq_wrk = args[arg];
    } else if ( key == "-pdbin"        || key == "-pdbin-wrk" ) {
      if ( ++arg < args.size() ) ippdb_wrk = args[arg];
    } else if ( key == "-pdbin-mr"     || key == "-pdbin-wrk-mr" ) {
      if ( ++arg < args.size() ) ippdb_mr  = args[arg];
    } else if ( key == "-pdbin-sequence-prior" || key == "-pdbin-wrk-sequence-prior" ) {
      if ( ++arg < args.size() ) ippdb_seq = args[arg];
    } else if ( key == "-pdbout"      || key == "-pdbout-wrk" ) {
      if ( ++arg < args.size() ) oppdb = args[arg];
    } else if ( key == "-mapout" ) {
      if ( ++arg < args.size() ) opmap  = args[arg];
    } else if ( key == "-xmlout" ) {
      if ( ++arg < args.size() ) opxml  = args[arg];
    } else if ( key == "-colin-ref-fo" ) {
      if ( ++arg < args.size() ) ipcol_ref_fo = args[arg];
    } else if ( key == "-colin-ref-hl" ) {
      if ( ++arg < args.size() ) ipcol_ref_hl = args[arg];
    } else if ( key == "-colin-fo"     || key == "-colin-wrk-fo" ) {
      if ( ++arg < args.size() ) ipcol_wrk_fo = args[arg];
    } else if ( key == "-colin-hl"     || key == "-colin-wrk-hl" ) {
      if ( ++arg < args.size() ) ipcol_wrk_hl = args[arg];
    } else if ( key == "-colin-phifom" || key == "-colin-wrk-phifom" ) {
      if ( ++arg < args.size() ) ipcol_wrk_pw = args[arg];
    } else if ( key == "-colin-fc"     || key == "-colin-wrk-fc" ) {
      if ( ++arg < args.size() ) ipcol_wrk_fc = args[arg];
    } else if ( key == "-colin-free"   || key == "-colin-wrk-free" ) {
      if ( ++arg < args.size() ) ipcol_wrk_fr = args[arg];
    } else if ( key == "-resolution" ) {
      if ( ++arg < args.size() ) res_in = clipper::String(args[arg]).f();
    } else if ( key == "-cycles" ) {
      if ( ++arg < args.size() ) ncyc  = clipper::String(args[arg]).i();
    } else if ( key == "-merge" ) {
      merge = true;
    } else if ( key == "-find" ) {
      find = true;
    } else if ( key == "-grow" ) {
      grow = true;
    } else if ( key == "-join" ) {
      join = true;
    } else if ( key == "-link" ) {
      link = true;
    } else if ( key == "-sequence" ) {
      seqnc = true;
    } else if ( key == "-correct" ) {
      corct = true;
    } else if ( key == "-filter" ) {
      filtr = true;
    } else if ( key == "-ncsbuild" ) {
      ncsbd = true;
    } else if ( key == "-prune" ) {
      prune = true;
    } else if ( key == "-rebuild" ) {
      build = true;
    } else if ( key == "-tidy" ) {
      tidy = true;
    } else if ( key == "-fast" ) {
      fast = true;
    } else if ( key == "-nn-select" ) { //neural network
      nnselect = true;
    } else if ( key == "-nn-model-num" ) {
      if ( ++arg < args.size() ) model_num  = clipper::String(args[arg]).f();
    } else if ( key == "-nn-threshold-select-method" ) {
      if ( ++arg < args.size() ) threshold_select_method  = clipper::String(args[arg]);
    } else if ( key == "-nn-no-confirmation-cycles" ) { //neural network
      nn_confirmation_cycles = false;
    } else if ( key == "-nn-confirmation-cycles-num" ) { //neural network
      if ( ++arg < args.size() ) confirmation_cycles_num = clipper::String(args[arg]).i();
    } else if ( key == "-build-semet" ) {
      semet = true;
    } else if ( key == "-anisotropy-correction" ) {
      doanis = true;
    } else if ( key == "-fix-position" ) {
      fixpos = true;
    } else if ( key == "-output-intermediate-models" ) {
      optemp = true;
    } else if ( key == "-fragments" ) {
      if ( ++arg < args.size() ) nfrag  = clipper::String(args[arg]).i();
    } else if ( key == "-fragments-per-100-residues" ) {
      if ( ++arg < args.size() ) nfragr = clipper::String(args[arg]).i();
    } else if ( key == "-free-flag" ) {
      if ( ++arg < args.size() ) freerindex = clipper::String(args[arg]).i();
    } else if ( key == "-model-index" ) {
      if ( ++arg < args.size() ) modelindex = clipper::String(args[arg]).i();
    } else if ( key == "-random-seed" ) {
      if ( ++arg < args.size() ) randomseed = clipper::String(args[arg]).i();
    } else if ( key == "-ramachandran-filter" ) {
      if ( ++arg < args.size() ) {
        if ( args[arg] == "all"      ) rama_flt = Ca_prep::rama_flt_all;
        if ( args[arg] == "helix"    ) rama_flt = Ca_prep::rama_flt_helix;
        if ( args[arg] == "strand"   ) rama_flt = Ca_prep::rama_flt_strand;
        if ( args[arg] == "nonhelix" ) rama_flt = Ca_prep::rama_flt_nonhelix;
      }
    } else if ( key == "-ramachandran-growing" ) {
      rama_grow = true;
    } else if ( key == "-known-structure" ) {
      if ( ++arg < args.size() )
        known_ids.push_back( KnownStructure::parse(args[arg] ) );
    } else if ( key == "-nonprotein-radius" ) {
      if ( ++arg < args.size() ) nprad = clipper::String(args[arg]).f();
    } else if ( key == "-main-chain-likelihood-radius" ) {
      if ( ++arg < args.size() ) main_tgt_rad = clipper::String(args[arg]).f();
    } else if ( key == "-side-chain-likelihood-radius" ) {
      if ( ++arg < args.size() ) side_tgt_rad = clipper::String(args[arg]).f();
    } else if ( key == "-sequence-reliability" ) {
      if ( ++arg < args.size() ) seq_rel = clipper::String(args[arg]).f();
    } else if ( key == "-new-residue-name" ) {
      if ( ++arg < args.size() ) newresname = args[arg];
    } else if ( key == "-new-residue-type" ) {
      if ( ++arg < args.size() ) newrestype = args[arg];
    } else if ( key == "-offset" ) {
      if ( ++arg < args.size() ) moffset = clipper::String(args[arg]).f();
    } else if ( key == "-correlation-mode" ) {
      correl = true;
    } else if ( key == "-no-correlation-mode" ) {
      nocorrel = true;
    } else if ( key == "-cif" ) {
      cifflag =  clipper::MMDBManager::CIF;
    } else if ( key == "-model-filter" ) {
      model_filter = true;
    } else if ( key == "-mr-model" ) {
      mr_model = true;
    } else if ( key == "-mr-model-filter" ) {
      mr_model_filter = mr_model = true;
    } else if ( key == "-mr-model-seed" ) {
      mr_model_seed = mr_model = true;
    } else if ( key == "-model-filter-sigma" ) {
      if (++arg<args.size()) model_filter_sig = clipper::String(args[arg]).f();
    } else if ( key == "-mr-model-filter-sigma" ) {
      if (++arg<args.size()) mr_filter_sig = clipper::String(args[arg]).f();
    } else if ( key == "-jobs" || key == "-j" ) {
      if ( ++arg < args.size() ) ncpu = clipper::String(args[arg]).i();
    } else if ( key == "-verbose" ) {
      if ( ++arg < args.size() ) verbose = clipper::String(args[arg]).i();
    } else {
      std::cout << "\nUnrecognized:\t" << args[arg] << std::endl;
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "\nUsage: cbuccaneer\n\t-mtzin-ref <filename>\n\t-pdbin-ref <filename>\n\t-mtzin <filename>\t\tCOMPULSORY\n\t-seqin <filename>\n\t-pdbin <filename>\n\t-pdbin-mr <filename>\n\t-pdbin-sequence-prior <filename>\n\t-pdbout <filename>\n\t-xmlout <filename>\n\t-colin-ref-fo <colpath>\n\t-colin-ref-hl <colpath>\n\t-colin-fo <colpath>\t\tCOMPULSORY\n\t-colin-hl <colpath> or -colin-phifom <colpath>\tCOMPULSORY\n\t-colin-fc <colpath>\n\t-colin-free <colpath>\n\t-resolution <resolution/A>\n\t-find\n\t-grow\n\t-join\n\t-link\n\t-sequence\n\t-correct\n\t-filter\n\t-ncsbuild\n\t-prune\n\t-rebuild\n\t-tidy\n\t-fast\n\t-anisotropy-correction\n\t-build-semet\n\t-fix-position\n\t-cycles <num_cycles>\n\t-fragments <max_fragments>\n\t-fragments-per-100-residues <num_fragments>\n\t-ramachandran-filter <type>\n\t-ramachandran-growing\n\t-known-structure <atomid[,radius]>\n\t-main-chain-likelihood-radius <radius/A>\n\t-side-chain-likelihood-radius <radius/A>\n\t-sequence-reliability <value>\n\t-new-residue-name <type>\n\t-new-residue-type <type>\n\t-correlation-mode\n\t-no-correlation-mode\n\t-mr-model\n\t-mr-model-seed\n\t-mr-model-filter\n\t-mr-model-filter-sigma <sigma>\n\t-model-filter\n\t-model-filter-sigma <sigma>\n\t-jobs <CPUs>\n\t-cif\t\t*will only output model in cif format\n\t-output-intermediate-models\n\t-verbose verbosity\n";
    std::cout << "\nAn input pdb and mtz are required for the reference structure, and an input mtz file for the work structure. Chains will be located and grown for the work structure and written to the output pdb file. \nThis involves 10 steps:\n finding, growing, joining, linking, sequencing, correcting, filtering, ncs building, pruning, and rebuilding. \nIf the optional input pdb file is provided for the work structure, then the input model is extended." << std::endl;
    exit(1);
  }

  // other initialisations
  if ( randomseed >= 0 ) srand( randomseed );
  using clipper::data32::Compute_abcd_from_phifom;
  using clipper::data32::Compute_phifom_from_abcd;
  using clipper::data32::Compute_fphi_from_fsigf_phifom;
  using clipper::data32::Compute_scale_u_iso_fsigf;
  using clipper::data32::Compute_scale_u_aniso_fphi;
  using clipper::data32::F_sigF;
  using clipper::data32::F_phi;
  using clipper::data32::Phi_fom;
  using clipper::data32::ABCD;
  using clipper::data32::Flag;
  clipper::Resolution resol;
  clipper::CCP4MTZfile mtzfile;
  mtzfile.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
  std::string msg;
  const int mmdbflags = ::mmdb::MMDBF_IgnoreBlankLines | ::mmdb::MMDBF_IgnoreDuplSeqNum | ::mmdb::MMDBF_IgnoreNonCoorPDBErrors | ::mmdb::MMDBF_IgnoreRemarks;
  ProteinTools proteintools = ProteinTools();
  Ca_prep::set_cpus( ncpu );
  Ca_find::set_cpus( ncpu );
  Ca_grow::set_cpus( ncpu );
  Ca_sequence::set_cpus( ncpu );
  Ca_sequence::set_semet( semet );
  Ca_find::TYPE findtype = fast ? Ca_find::SECSTRUC : Ca_find::LIKELIHOOD;
  if ( !(find||grow||join||link||seqnc||corct||filtr||ncsbd||prune||build||tidy) )
    find=grow=join=link=seqnc=corct=filtr=ncsbd=prune=build=tidy=true;
  if ( !( correl || nocorrel ) )
    correl = ( ippdb_wrk != "NONE" || ippdb_mr != "NONE" );
  if ( ipmtz_ref == "NONE" || ippdb_ref == "NONE" )
    BuccaneerUtil::set_reference( ipmtz_ref, ippdb_ref );
    
  BuccaneerLog log( title );
  // messages
  std::cout << std::endl << ((correl)?std::string("Correlation mode selected"):std::string("Correlation mode not selected")) << std::endl << std::endl;

  // Get resolution for calculation
  mtzfile.open_read( ipmtz_ref );
  double res_ref = std::max( mtzfile.resolution().limit(), res_in );
  mtzfile.close_read();
  mtzfile.open_read( ipmtz_wrk );
  double res_wrk = std::max( mtzfile.resolution().limit(), res_in );
  mtzfile.close_read();
  resol = clipper::Resolution( std::max( res_ref, res_wrk ) );
  if ( res_ref > res_wrk ) std::cout << std::endl << "WARNING: resolution of work structure truncated to reference:\n Ref: " << res_ref << " Wrk: " << res_wrk << std::endl;

  // Get reference reflection data
  clipper::HKL_info hkls_ref;
  mtzfile.open_read( ipmtz_ref );
  hkls_ref.init( mtzfile.spacegroup(), mtzfile.cell(), resol, true );
  clipper::HKL_data<F_sigF> ref_f( hkls_ref );
  clipper::HKL_data<ABCD>   ref_hl( hkls_ref );
  mtzfile.import_hkl_data( ref_f,  ipcol_ref_fo );
  mtzfile.import_hkl_data( ref_hl, ipcol_ref_hl );
  mtzfile.close_read();

  // Get work reflection data
  clipper::MTZcrystal cxtl;
  clipper::HKL_info hkls_wrk;
  mtzfile.set_verbose( (verbose>0) ? 3 : 2 );
  mtzfile.open_read( ipmtz_wrk );
  hkls_wrk.init( mtzfile.spacegroup(), mtzfile.cell(), resol, true );
  mtzfile.import_crystal( cxtl, ipcol_wrk_fo+".F_sigF.F" );
  clipper::HKL_data<F_sigF>  wrk_f ( hkls_wrk, cxtl );
  clipper::HKL_data<ABCD>    wrk_hl( hkls_wrk, cxtl );
  clipper::HKL_data<Phi_fom> wrk_pw( hkls_wrk, cxtl );
  clipper::HKL_data<F_phi>   wrk_fp( hkls_wrk, cxtl );
  clipper::HKL_data<Flag>    flag( hkls_wrk, cxtl );
  mtzfile.import_hkl_data( wrk_f , ipcol_wrk_fo );
  if ( ipcol_wrk_hl != "NONE" ) mtzfile.import_hkl_data( wrk_hl,ipcol_wrk_hl );
  if ( ipcol_wrk_pw != "NONE" ) mtzfile.import_hkl_data( wrk_pw,ipcol_wrk_pw );
  if ( ipcol_wrk_fc != "NONE" ) mtzfile.import_hkl_data( wrk_fp,ipcol_wrk_fc );
  if ( ipcol_wrk_fr != "NONE" ) mtzfile.import_hkl_data( flag,  ipcol_wrk_fr );
  mtzfile.close_read();

  // do anisotropy correction
  clipper::U_aniso_orth uaniso( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
  if ( doanis ) {
    // scale obs data
    typedef clipper::SFscale_aniso<float> SFscale;
    SFscale sfscl( 3.0, SFscale::SHARPEN );
    sfscl( wrk_f );
    uaniso = sfscl.u_aniso_orth( SFscale::F );
    // scale map coeffs
    Compute_scale_u_aniso_fphi compute_aniso( 1.0, -uaniso );
    if ( ipcol_wrk_fc != "NONE" ) wrk_fp.compute( wrk_fp, compute_aniso );
    // output
    std::cout << std::endl << "Applying anisotropy correction:"
              << std::endl << uaniso.format() << std::endl << std::endl;
  }

  // apply free flag
  clipper::HKL_data<F_sigF> wrk_fwrk = wrk_f;
  clipper::HKL_data<F_sigF>  wrk_f_copy=wrk_f;// to use in R-factor calculation as the following make changes on wrk_f

  //wrk_fwrk.mask( flag != freerindex );
  for ( clipper::HKL_data_base::HKL_reference_index ih = hkls_wrk.first(); !ih.last(); ih.next() ) if ( flag[ih].flag() == freerindex ) wrk_fwrk[ih] = F_sigF();  //ugly hack for broken SGI compilers
  // and fill in hl
  if ( ipcol_wrk_hl == "NONE" )
    wrk_hl.compute( wrk_pw, Compute_abcd_from_phifom() );

  // Get reference model
  clipper::Spacegroup cspg = hkls_wrk.spacegroup();
  clipper::MiniMol mol_ref, mol_wrk_in, mol_tmp;
  clipper::MiniMol mol_wrk(cspg,cxtl), mol_mr(cspg,cxtl), mol_seq(cspg,cxtl);
  clipper::MMDBfile mmdb_ref;
  mmdb_ref.SetFlag( mmdbflags );
  mmdb_ref.read_file( ippdb_ref );
  mmdb_ref.import_minimol( mol_ref );

  // Get work model (optional)
  BuccaneerUtil::read_model( mol_wrk, ippdb_wrk, verbose>5 );
  // Get MR model - to help rebuilding
  BuccaneerUtil::read_model( mol_mr,  ippdb_mr,  verbose>5 );
  // Get sequencing model - heavy atoms or MR (optional)
  BuccaneerUtil::read_model( mol_seq, ippdb_seq, verbose>5 );
  if ( mol_seq.size() > 0 ) Ca_sequence::set_prior_model( mol_seq );
  // store a copy of the input model
  mol_wrk_in = mol_wrk;
  // store known structure info
  KnownStructure knownstruc( mol_wrk_in, known_ids, nprad );
  if ( verbose >= 1 ) knownstruc.debug();

  // Get work sequence (optional)
  clipper::MMoleculeSequence seq_wrk;
  if ( ipseq_wrk != "NONE" ) {
    clipper::SEQfile seqf_wrk;
    seqf_wrk.read_file( ipseq_wrk );
    seqf_wrk.import_molecule_sequence( seq_wrk );
  }

  // check input files match mode
  if ( !find && mol_wrk.is_null() )
    clipper::Message::message(clipper::Message_fatal("Missing work model."));
  if ( seqnc && seq_wrk.is_null() )
    clipper::Message::message(clipper::Message_fatal("Missing work sequence."));

  // DO INITIAL MAP SIMULATION
  clipper::HKL_data<F_sigF> sim_f( hkls_ref );
  clipper::HKL_data<ABCD> sim_hl( hkls_ref );
  MapSimulate mapsim( 100, 20 );
  mapsim( sim_f, sim_hl, ref_f, ref_hl, wrk_f, wrk_hl );

  // make llk target objects
  LLK_map_target llktgt;
  std::vector<LLK_map_target> llkcls( 20 );

  // STAGE 1: Calculate target from reference data

  {
    // reference map
    clipper::HKL_data<F_phi>   ref_fp( hkls_ref );
    clipper::HKL_data<Phi_fom> ref_pw( hkls_ref );
    ref_pw.compute( sim_hl, Compute_phifom_from_abcd() );
    ref_fp.compute( sim_f, ref_pw, Compute_fphi_from_fsigf_phifom() );
    clipper::Grid_sampling grid( hkls_ref.spacegroup(), hkls_ref.cell(), hkls_ref.resolution() );
    clipper::Xmap<float> xref( hkls_ref.spacegroup(), hkls_ref.cell(), grid );
    xref.fft_from( ref_fp );

    // prepare llk targets
    Ca_prep caprep( main_tgt_rad, side_tgt_rad, rama_flt, correl, seqnc,
                    verbose>3 );
    caprep( llktgt, llkcls, mol_ref, xref );

    log.log( "PREP" );
  }


  // STAGE 2: Apply target to work data

  {
    // work map
    wrk_pw.compute( wrk_hl, Compute_phifom_from_abcd() );
    if ( ipcol_wrk_fc == "NONE" )
      wrk_fp.compute( wrk_fwrk, wrk_pw, Compute_fphi_from_fsigf_phifom() );
    clipper::Grid_sampling grid( cspg, cxtl, hkls_wrk.resolution() );
    clipper::Xmap<float>   xwrk( cspg, cxtl, grid );
    xwrk.fft_from( wrk_fp );

    // optionlly write work map
    if ( opmap != "NONE" ) {
      clipper::CCP4MAPfile mapfile;
      mapfile.open_write( opmap );
      mapfile.export_xmap( xwrk );
      mapfile.close_write();
    }

    // number of residues to find
    double vol = xwrk.cell().volume() / double(xwrk.spacegroup().num_symops());
    int nres   = int( vol / 320.0 );  // 320A^3/residue on average (inc solvent)
    nfrag = std::min( nfrag, (nfragr*nres)/100 );

    // offset the map density
    clipper::Map_stats stats( xwrk );
    clipper::Xmap_base::Map_reference_index ix;
    for ( ix = xwrk.first(); !ix.last(); ix.next() )
      xwrk[ix] += moffset * stats.std_dev();

    // generate distribution of llk target values for cutoffs
    llktgt.prep_llk_distribution( xwrk );

    // tidy input model
    ProteinTools::split_chains_at_gap( mol_wrk );

    // prepare search target
    Ca_find cafind( nfrag, resol.limit() );

    // merge multi-model results
    if ( merge ) {
      Ca_merge camerge( seq_rel );
      std::cout << " C-alphas before model merge:  " << mol_wrk.select("*/*/CA").atom_list().size() << std::endl;
      camerge( mol_wrk, xwrk, llkcls, seq_wrk );
      std::cout << " C-alphas after model merge:  " << mol_wrk.select("*/*/CA").atom_list().size() << std::endl;
      log.log( "MRG ", mol_wrk, verbose>9 );
    }

    // Filter input model
    if ( model_filter ) {
      std::cout << " C-alphas before model filter:  " << mol_wrk.select("*/*/CA").atom_list().size() << std::endl;
      Ca_filter::filter( mol_wrk, xwrk, model_filter_sig );
      std::cout << " C-alphas after model filter:  " << mol_wrk.select("*/*/CA").atom_list().size() << std::endl;
      log.log( "FLT ", mol_wrk, verbose>9 );
    }

    // Augment input model with MR model
    if ( mr_model ) {
      std::vector<int> result = Ca_merge::merge_mr( mol_wrk, mol_mr, mr_filter_sig, 3, mr_model_filter, mr_model_seed );
      std::cout << " MR residues input: " << result[0] << ", after filter: " << result[1] << ", after seeding: " << result[2] << std::endl;
    }

    // Trim input model to only protein residues
    ProteinTools::trim_to_protein( mol_wrk );

    // model building loop
    for ( int cyc = 0; cyc < ncyc; cyc++ ) {
      std::cout << std::endl << "Cycle: " << cyc+1 << std::endl << std::endl;
      clipper::String history = "";

      // find C-alphas by slow likelihood search
      if ( find ) {
        cafind( mol_wrk, knownstruc, xwrk, llktgt, findtype, modelindex );
        std::cout << " C-alphas after finding:    " << mol_wrk.select("*/*/CA").atom_list().size() << std::endl;
        log.log( "FIND", mol_wrk, verbose>9 );
      }

      // grow C-alphas
      if ( grow ) {
        Ca_grow cagrow( 25 );
        cagrow( mol_wrk, xwrk, llktgt, rama_grow );
        std::cout << " C-alphas after growing:    " << mol_wrk.select("*/*/CA").atom_list().size() << std::endl;
        log.log( "GROW", mol_wrk, verbose>9 );
      }

      // omit incorrect C-alphas using NN
      if (nnselect) {
        NNbuild *nn = new NNbuild(nfrag, resol, mol_wrk, knownstruc, xwrk, llktgt, findtype, modelindex, llkcls, seq_wrk, seq_rel, newrestype,  cyc,  ncyc, mol_mr,freerindex,mol_wrk_in,ipmtz_wrk,ipcol_wrk_fo,ipcol_wrk_fr,cafind,mtzfile.resolution().limit(),model_num,threshold_select_method,nn_confirmation_cycles,confirmation_cycles_num,rama_grow);
      } else { // begin NN else clause

      // join C-alphas
      if ( join ) {
        Ca_join cajoin( 2.0, 2.0 );
        cajoin( mol_wrk );
        std::cout << " C-alphas after joining:    " << mol_wrk.select("*/*/CA").atom_list().size() << std::endl;
        log.log( "JOIN", mol_wrk, verbose>9 );
      }

      // link C-alphas
      if ( link ) {
        Ca_link calnk( 10.0, 24 );
        calnk( mol_wrk, xwrk, llktgt );
        std::cout << " C-alphas linked:           " << calnk.num_linked() << std::endl;
        log.log( "LINK", mol_wrk, verbose>9 );
      }

      // assign sequences
      if ( seqnc ) {
        Ca_sequence caseq( seq_rel );
        caseq( mol_wrk, xwrk, llkcls, seq_wrk );
        std::cout << " C-alphas sequenced:        " << caseq.num_sequenced() << std::endl;
        history = caseq.format();
        log.log( "SEQU", mol_wrk, verbose>9 );
      }

      // correct insertions/deletions
      if ( corct ) {
        Ca_correct cacor( 12 );
        cacor( mol_wrk, xwrk, llkcls, seq_wrk );
        std::cout << " C-alphas corrected:        " << cacor.num_corrected() << std::endl;
        log.log( "CORR", mol_wrk, verbose>9 );
      }

      //  poor density
      if ( filtr ) {
        Ca_filter cafiltr( 1.0 );
        cafiltr( mol_wrk, xwrk );
        std::cout << " C-alphas after filtering:  " << mol_wrk.select("*/*/CA").atom_list().size() << std::endl;
        log.log( "FILT", mol_wrk, verbose>9 );
      }

      // ncsbuild C-alphas
      if ( ncsbd ) {
        Ca_ncsbuild cancsbuild( seq_rel, 1.0, 12 );
        cancsbuild( mol_wrk, xwrk, llkcls, seq_wrk );
        std::cout << " C-alphas after NCS build:  " << mol_wrk.select("*/*/CA").atom_list().size() << std::endl;
        log.log( "NCSB", mol_wrk, verbose>9 );
      }

      // prune C-alphas for clashes with other chains and known model
      if ( prune ) {
        Ca_prune caprune( 3.0 );
        caprune( mol_wrk, xwrk );
        std::cout << " C-alphas after pruning:    " << mol_wrk.select("*/*/CA").atom_list().size() << std::endl;
        log.log( "PRUN", mol_wrk, verbose>9 );
      }

      // build side chains/atoms
      if ( build ) {
        Ca_build cabuild( newrestype );
        cabuild( mol_wrk, xwrk );
        std::cout << " C-alphas after rebuilding: " << mol_wrk.select("*/*/CA").atom_list().size() << std::endl;
        log.log( "REBU", mol_wrk, verbose>9 );
      }

      // tidy output model
      knownstruc.prune( mol_wrk );
      ProteinTools::split_chains_at_gap( mol_wrk );
      ProteinTools::chain_number( mol_wrk );
      ProteinTools::chain_label( mol_wrk, cifflag );
      log.log( "TDYI", mol_wrk, verbose>9 );

      } // end NN else clause

      // user output
      std::cout << std::endl;
      if ( verbose > 7 ) std::cout << history << std::endl;
      msg = log.log( mol_wrk, mol_mr, seq_wrk );
      prog.summary_beg();
      std::cout << "Internal cycle " << clipper::String( cyc+1, 3 )
                << std::endl << msg << std::endl;
      prog.summary_end();

      // file output
      if ( opxml != "NONE" ) log.xml( opxml );
      if ( optemp ) {
        clipper::String c( cyc, 1 );
        clipper::MMDBfile mmdb;
        mmdb.export_minimol( mol_wrk );
        mmdb.write_file( oppdb.substr(0,oppdb.rfind(".")+1) + c + oppdb.substr(oppdb.rfind(".")), cifflag );
      }
    } // next cycle

    // move model
    if ( fixpos )
      ProteinTools::symm_match( mol_wrk, mol_wrk_in );

    // model tidy
    if ( tidy ) {
      // tidy chains and assign chain IDs
      ModelTidy mtidy( 1.0, 12, newrestype, verbose > 6 );
      bool chk = mtidy.tidy( mol_wrk, mol_mr, seq_wrk );
      if ( !chk ) std::cout << "ModelTidy error" << std::endl; // can't happen
      // rebuild side chains again due to possible change in sequence. FUTURE: flex side chains
      //Ca_build cabuild( newrestype );  // , True?
      //cabuild( mol_wrk, xwrk );
      log.log("TIDY", mol_wrk, verbose>9 );
    }

    // Assign default B-factors to missing values
    float default_u_iso = ProteinTools::main_chain_u_mean( mol_wrk_in );
    for ( int c = 0; c < mol_wrk.size(); c++ )
      for ( int r = 0; r < mol_wrk[c].size(); r++ )
        for ( int a = 0; a < mol_wrk[c][r].size(); a++ )
          if ( clipper::Util::is_nan( mol_wrk[c][r][a].u_iso() ) )
            mol_wrk[c][r][a].set_u_iso( default_u_iso );

    // adjust residue names
    if ( newresname != "NONE" )
      for ( int c = 0; c < mol_wrk.size(); c++ )
        for ( int r = 0; r < mol_wrk[c].size(); r++ )
          if ( ProteinTools::residue_index_3( mol_wrk[c][r].type() ) < 0 )
            mol_wrk[c][r].set_type( newresname );

    // add known structure from input model
    bool copy = knownstruc.copy_to( mol_wrk );
    if ( !copy ) std::cout << "$TEXT:Warning: $$ $$\nWARNING: chain ID clash between known-structure and pdbin-mr. Chains renamed.\n$$" << std::endl;
    
    // label unlabelled chains
    ProteinTools::chain_label( mol_wrk, cifflag );
    log.log("LABE",mol_wrk,verbose>9);
    
    // write answers
    clipper::MMDBfile mmdb;
    mmdb.export_minimol( mol_wrk );
    mmdb.write_file( oppdb, cifflag );
  }

  std::cout << "$TEXT:Result: $$ $$" << std::endl << msg << "$$" << std::endl;
  log.profile();
  prog.set_termination_message( "Normal termination" );
}
