/*
 * buccaneer-NN-build.cpp
 *
 *  Created on: 16 May 2021
 *      Author: Emad Alharbi
 */

#include "buccaneer-NN-build.h"

#include "buccaneer-NN-build-cyc.h"
#include "buccaneer-NN-Select.h"
#include "buccaneer-util.h"
#include "buccaneer-tidy.h"
#include <limits>
#include "buccaneer-known.h"

bool NNbuild::nn_init;
NNbuild::NNbuild(int nfrag, clipper::Resolution resol,
		clipper::MiniMol &mol_wrk, KnownStructure knownstruc,
		clipper::Xmap<float> xwrk, LLK_map_target llktgt,
		Ca_find::TYPE findtype, int modelindex,
		std::vector<LLK_map_target> llkcls, clipper::MMoleculeSequence seq_wrk,
		double seq_rel, clipper::String newrestype, int cyc, int ncyc,
		 clipper::MiniMol mol_mr ,int freeflag,clipper::MiniMol  mol_wrk_in ,clipper::String ipfile,clipper::String ipcolfo,clipper::String ipcolfree,Ca_find cafind,double reso,double model_num, clipper::String threshold_select_method, bool nn_confirmation_cycles, int confirmation_cycles_num, bool rama_grow) {
	clipper::String title;
	BuccaneerLog lognn(title);

	//NN init
	SelectNN::update_predictions=false;
	SelectNN::select=false;
	SelectNN::threshold=0;
    SelectNN::reso=reso;
    SelectNN::model_num=model_num;
    SelectNN::threshold_select_method =threshold_select_method;

    const char *nnweightsptr = getenv("CLIBD");


    SelectNN::model = clipper::String(nnweightsptr) + "/buccaneer_nn_weights.csv";
    if(SelectNN::threshold_select_method != "NONE")
       SelectNN::model_num=5000; // in case of selecting non-default method, we need to increase the number of model to be sure the determined number of thresholds is less than the number of models. However not all the number of models will be used.

    std::vector<double> thresholds;
    for(double thrs=0.0; thrs < SelectNN::model_num ; ++thrs)
    	thresholds.push_back(thrs);



	//model evaluation indicators
	int unq_num = -1;// in case the best model has zero sequenced residues
	double best_rwork=1;



	clipper::MiniMol mol_wrk_best(mol_wrk.spacegroup(), mol_wrk.cell());
	clipper::MiniMol mol_wrk_current_cyc(mol_wrk.spacegroup(), mol_wrk.cell());


	SelectNN::update_predictions = true;
	SelectNN::llktgt = llktgt;
	SelectNN::xwrk = xwrk;

	clipper::MiniMol base_mol_wrk_NN(mol_wrk.spacegroup(), mol_wrk.cell());


	clipper::String best_mol_log="";// to save the log and print out at the end of the cycle



	for (int thrs = 0; thrs < thresholds.size(); ++thrs) {

		if (thrs != 0)
		SelectNN::update_predictions = false; // to not run NN every time when test different thresholds
		SelectNN::threshold = thresholds[thrs];
		clipper::MiniMol mol_wrk_NN(mol_wrk.spacegroup(), mol_wrk.cell());
		mol_wrk_NN = mol_wrk;// take a copy of the model
		SelectNN::select = true;
		clipper::String current_mol_log="";



		//Now complete  the model building starting from joining step
		runcyc buildcyc(nfrag, resol, mol_wrk_NN, knownstruc, xwrk, llktgt,
				findtype, modelindex, llkcls, seq_wrk, seq_rel, newrestype, 0,
				2, 3, 10,current_mol_log,cafind,residues_container,true,rama_grow); // 0,2 will do only one building cycle





		mol_wrk_current_cyc = mol_wrk_NN; //save a copy before running confirmation cycles

		clipper::String current_mol_log_confirmation="";





		//Now do confirmation cycles starting from finding step
		if(nn_confirmation_cycles){
			int current_cycle=cyc;
			int cycles_num=ncyc;
			if(confirmation_cycles_num !=-1){
				 current_cycle=0;
				 cycles_num=confirmation_cycles_num+1; // +1 because in build-cyc will add 1 to current_cycle
			}
			runcyc buildcycconf(nfrag, resol, mol_wrk_NN, knownstruc, xwrk, llktgt,
					 	findtype, modelindex, llkcls, seq_wrk, seq_rel, newrestype, current_cycle,
						cycles_num, 1, 10,current_mol_log_confirmation,cafind,residues_container,true,rama_grow);
		}



		std::vector<double> evaluated_mol=lognn.evaluate(mol_wrk_NN, mol_mr, seq_wrk);







		std::vector<double> rwork_rfree=rfactors(mol_wrk_NN,freeflag,mol_wrk_in,mol_mr , ipfile, ipcolfo, ipcolfree);



		// evaluate the model
		if ( (evaluated_mol[0] > unq_num &&  rwork_rfree[0]<best_rwork)   ||  unq_num==-1){ // UnqNum==-1 in case R-work is nan for first built model



			unq_num = evaluated_mol[0];
			best_rwork=rwork_rfree[0];
			mol_wrk_best = mol_wrk_current_cyc; // Save the best model. The copy before the confirmation cycles
			best_mol_log=current_mol_log;
		}



		if(SelectNN::threshold==-1)break; // it will be set to -1 in buccaneer-NN-Select.cpp in case of the number of fragments higher than the NN input limit or when the thresholds reached the highest NN predicted probability

	}

	mol_wrk = mol_wrk_best;
	std::cout<<best_mol_log<<std::endl;
	residues_container.clear();

}

NNbuild::~NNbuild() {
	// TODO Auto-generated destructor stub
}


std::vector<double> NNbuild::rfactors(clipper::MiniMol mol1,int freeflag,clipper::MiniMol  mol_wrk_in,clipper::MiniMol mol_mr, clipper::String ipfile,clipper::String ipcolfo,clipper::String ipcolfree){

/*
 * using  clipper::HKL_data<F_sigF> fo from cbuccaneer leads to wrong R-factors compared to  csfcalc
 */


	// Assign default B-factors to missing values
			    float default_u_iso = ProteinTools::main_chain_u_mean( mol_wrk_in );
			    for ( int c = 0; c < mol1.size(); c++ )
			      for ( int r = 0; r < mol1[c].size(); r++ )
			        for ( int a = 0; a < mol1[c][r].size(); a++ )
			          if ( clipper::Util::is_nan( mol1[c][r][a].u_iso() ) )
			        	  mol1[c][r][a].set_u_iso( default_u_iso );



	  // defaults
	  enum ANISO { NONE, FOBS, FCAL };
	  clipper::String title;
	  clipper::String ippdb = "NONE";

	  clipper::String opfile = "sfcalc.mtz";
	  clipper::String opcol = "sfcalc";
	  bool bulk  = true;
	  ANISO aniso = NONE;

	  int n_refln = 1000;
	  int n_param = 20;
	  int verbose = 0;



	  // make data objects
	  clipper::CCP4MTZfile mtzin, mtzout;
	  clipper::MTZcrystal cxtl;
	  clipper::HKL_info hkls;
	  double bulkfrc, bulkscl;
	  typedef clipper::HKL_data_base::HKL_reference_index HRI;
	  using clipper::data32::F_sigF;  using clipper::data32::F_phi;
	  using clipper::data32::Phi_fom; using clipper::data32::Flag;
	  mtzin.set_column_label_mode( clipper::CCP4MTZfile::Legacy );

	  // open file
	  mtzin.open_read( ipfile );
	  mtzin.import_hkl_info( hkls );
	  mtzin.import_crystal( cxtl, ipcolfo+".F_sigF.F" );
	  clipper::HKL_data<F_sigF> fo( hkls, cxtl );
	  clipper::HKL_data<Flag> free( hkls, cxtl );
	  mtzin.import_hkl_data( fo, ipcolfo );
	  if ( ipcolfree != "NONE" ) mtzin.import_hkl_data( free, ipcolfree );
	  if ( opcol[0] != '/' ) opcol = mtzin.assigned_paths()[0].notail()+"/"+opcol;
	  mtzin.close_read();


	  clipper::Atom_list atoms = mol1.atom_list();

	  // calculate structure factors
	  clipper::HKL_data<F_phi> fc( hkls, cxtl );
	  if ( bulk ) {
	    clipper::SFcalc_obs_bulk<float> sfcb;
	    sfcb( fc, fo, atoms );
	    bulkfrc = sfcb.bulk_frac();
	    bulkscl = sfcb.bulk_scale();
	  } else {
	    clipper::SFcalc_aniso_fft<float> sfc;
	    sfc( fc, atoms );
	    bulkfrc = bulkscl = 0.0;
	  }

	  // do anisotropic scaling
	  if ( aniso != NONE )  {
	    clipper::SFscale_aniso<float>::TYPE F = clipper::SFscale_aniso<float>::F;
	    clipper::SFscale_aniso<float> sfscl;
	    if ( aniso == FOBS ) sfscl( fo, fc );  // scale Fobs
	    if ( aniso == FCAL ) sfscl( fc, fo );  // scale Fcal
	    std::cout << "\nAnisotropic scaling:\n"
		      << sfscl.u_aniso_orth(F).format() << "\n";
	  }

	  // now do sigmaa calc
	  clipper::HKL_data<F_phi>   fb( hkls, cxtl ), fd( hkls, cxtl );
	  clipper::HKL_data<Phi_fom> phiw( hkls, cxtl );
	  clipper::HKL_data<Flag>    flag( hkls, cxtl );
	  for ( HRI ih = flag.first(); !ih.last(); ih.next() )
	    if ( !fo[ih].missing() && (free[ih].missing()||free[ih].flag()==freeflag) )
	      flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
	    else
	      flag[ih].flag() = clipper::SFweight_spline<float>::NONE;

	  // do sigmaa calc
	  clipper::SFweight_spline<float> sfw( n_refln, n_param );
	  sfw( fb, fd, phiw, fo, fc, flag );

	  // calc abcd
	  clipper::HKL_data<clipper::data32::ABCD> abcd( hkls );
	  abcd.compute( phiw, clipper::data32::Compute_abcd_from_phifom() );

	   std::vector<clipper::ftype> params( 2, 0.0 );  // initial parameters
	    clipper::BasisFn_gaussian basisfn;
	    clipper::TargetFn_scaleF1F2<F_phi,F_sigF> targetfn( fc, fo );
	   // clipper::TargetFn_scaleF1F2<clipper::data32::F_sigF,clipper::data32::F_sigF>targetfn( fsig1, fsig2 );
	    clipper::ResolutionFn_nonlinear rfn( hkls, basisfn, targetfn, params );

	  // now calc R and R-free
	  //std::vector<double> params( n_param, 1.0 );
	 // clipper::BasisFn_spline basisfn( fo, n_param, 1.0 );
	  //clipper::TargetFn_scaleF1F2<F_phi,F_sigF> targetfn( fc, fo );
	  //clipper::ResolutionFn rfn( hkls, basisfn, targetfn, params );
	  double r1w, f1w, r1f, f1f, Fo, Fc;
	  r1w = f1w = r1f = f1f = 0.0;
	  for ( HRI ih = fo.first(); !ih.last(); ih.next() )
	    if ( !fo[ih].missing() ) {
	      Fo = fo[ih].f();
	      Fc = sqrt( rfn.f(ih) ) * fc[ih].f();
	      if ( free[ih].flag() == freeflag ) {
		r1f += fabs( Fo - Fc );
		f1f += Fo;
	      } else {
		r1w += fabs( Fo - Fc );
		f1w += Fo;
	      }
	    }
	  r1f /= clipper::Util::max( f1f, 0.1 );
	  r1w /= clipper::Util::max( f1w, 0.1 );
	  //std::cout << "\n R-factor      : " <<  r1w
	//	    << "\n Free R-factor : " <<  r1f  << "\n";
	  std::vector<double> rworkfree;
	  	  rworkfree.push_back (r1w );
	  	  	  rworkfree.push_back( r1f );


	  	  	  return rworkfree;


}


