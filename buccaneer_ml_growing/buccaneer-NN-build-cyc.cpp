/*
 * buccaneer-build-cyc.cpp
 *
 *  Created on: 12 Mar 2021
 *      Author: Emad Alharbi
 */

#include "buccaneer-NN-build-cyc.h"

#include "buccaneer-NN-Select.h"
#include <map>
runcyc::runcyc(int nfrag, clipper::Resolution resol, clipper::MiniMol &mol_wrk,
		KnownStructure knownstruc, clipper::Xmap<float> xwrk,
		LLK_map_target llktgt, Ca_find::TYPE findtype, int modelindex,
		std::vector<LLK_map_target> llkcls, clipper::MMoleculeSequence seq_wrk,
		double seq_rel, clipper::String newrestype, int cyc, int ncyc,
		int startstep, int endstep,clipper::String &mollog,Ca_find cafind,std::map<std::pair<clipper::String,clipper::String>, clipper::MMonomer> &seqscore, bool saveseqscore, bool rama_grow ) {


	// TODO Auto-generated constructor stub


	for (int c = cyc + 1; c < ncyc; ++c) {

		if (1 >= startstep && 1 <= endstep) {
			cafind(mol_wrk, knownstruc, xwrk, llktgt, findtype, modelindex);
			mollog+=" C-alphas after finding:    ";
			mollog+=std::to_string(mol_wrk.select("*/*/CA").atom_list().size());
			mollog+="\n";
		}



		if (2 >= startstep && 2 <= endstep) {
			Ca_grow cagrow(25);
			cagrow(mol_wrk, xwrk, llktgt, rama_grow);
			mollog+=" C-alphas after growing:    ";
			mollog+=std::to_string(mol_wrk.select("*/*/CA").atom_list().size());
			mollog+="\n";
		}



		if (3 >= startstep && 3 <= endstep) {
			Ca_join cajoin(2.0, 2.0);
			cajoin(mol_wrk);

			mollog+=" C-alphas after joining:    ";
			mollog+=std::to_string(mol_wrk.select("*/*/CA").atom_list().size());
			mollog+="\n";
		}


		SelectNN::select = false; // avoid use NN in ncsbuild and merge



		if (4 >= startstep && 4 <= endstep) {
			Ca_link calnk(10.0, 24);
			calnk(mol_wrk, xwrk, llktgt);
			mollog+=" C-alphas linked:           ";
			mollog+=std::to_string(calnk.num_linked());
			mollog+="\n";


		}

		if (5 >= startstep && 5 <= endstep) {
			if(saveseqscore==true)
			apply_sequence_score(mol_wrk, "SEQDAT",seqscore);

			Ca_sequence caseq(seq_rel);
			caseq(mol_wrk, xwrk, llkcls, seq_wrk);

			if(saveseqscore==true)
				save_sequence_score(mol_wrk, "SEQDAT",seqscore);

			mollog+=" C-alphas sequenced:        ";
			mollog+=std::to_string(caseq.num_sequenced());
			mollog+="\n";
		}

		if (6 >= startstep && 6 <= endstep) {
			Ca_correct cacor(12);
			cacor(mol_wrk, xwrk, llkcls, seq_wrk);

			mollog+=" C-alphas corrected:        ";
			mollog+=std::to_string(cacor.num_corrected());
			mollog+="\n";
		}

		if (7 >= startstep && 7 <= endstep) {
			Ca_filter cafiltr(1.0);
			cafiltr(mol_wrk, xwrk);
			mollog+=" C-alphas after filtering:  ";
			mollog+=std::to_string(mol_wrk.select("*/*/CA").atom_list().size());
			mollog+="\n";
		}

		if (8 >= startstep && 8 <= endstep) {
			Ca_ncsbuild cancsbuild(seq_rel, 1.0, 12);
			cancsbuild(mol_wrk, xwrk, llkcls, seq_wrk);
			mollog+=" C-alphas after NCS build:  ";
			mollog+=std::to_string(mol_wrk.select("*/*/CA").atom_list().size());
			mollog+="\n";
		}

		if (9 >= startstep && 9 <= endstep) {
			Ca_prune caprune(3.0);
			caprune(mol_wrk, xwrk);
			mollog+=" C-alphas after pruning:    ";
			mollog+=std::to_string(mol_wrk.select("*/*/CA").atom_list().size());
			mollog+="\n";
		}

		if (10 >= startstep && 10 <= endstep) {
			Ca_build cabuild(newrestype);
			cabuild(mol_wrk, xwrk);
			mollog+=" C-alphas after rebuilding: ";
			mollog+=std::to_string(mol_wrk.select("*/*/CA").atom_list().size());
			mollog+="\n";
		}

		knownstruc.prune(mol_wrk);
		ProteinTools::split_chains_at_gap(mol_wrk);
		ProteinTools::chain_number(mol_wrk);
		ProteinTools::chain_label(mol_wrk, clipper::MMDBManager::CIF);
	}

}

runcyc::~runcyc() {
	// TODO Auto-generated destructor stub
}
void runcyc::apply_sequence_score(clipper::MiniMol &Current_wrk, clipper::String propertylabel,std::map<std::pair<clipper::String,clipper::String>, clipper::MMonomer> &seqscore){

	// if the residue sequenced before, then we set the residue score that we have to avoid calculating the score again.
	int count_sequnced=0;

	for(int c=0; c < Current_wrk.size(); ++c){
		    	for(int r=0; r < Current_wrk[c].size(); ++r){
		    		Ca_group ca( Current_wrk[c][r] );
		    		clipper::String key = std::to_string(ca.coord_ca().x())+std::to_string(ca.coord_ca().y())+std::to_string(ca.coord_ca().z())+std::to_string(ca.coord_n().x())+std::to_string(ca.coord_n().y())+std::to_string(ca.coord_n().z())+std::to_string(ca.coord_c().x())+std::to_string(ca.coord_c().y())+std::to_string(ca.coord_c().z());


		    		if(seqscore.count( std::make_pair(key,propertylabel) )>0){//found
		    		if(seqscore[std::make_pair(key,propertylabel)].exists_property( propertylabel )){
		    		Current_wrk[c][r].delete_property(propertylabel);
		    		Current_wrk[c][r].set_property( propertylabel, seqscore[std::make_pair(key,propertylabel)].get_property( propertylabel ) );
		    		count_sequnced++;
		    		}
		    		}
		    	}
	}



}
void runcyc::save_sequence_score(clipper::MiniMol &Current_wrk, clipper::String propertylabel,std::map<std::pair<clipper::String,clipper::String>, clipper::MMonomer> &seqscore){

	// Save sequence score to speed up the multiple model building
 	for(int c=0; c < Current_wrk.size(); ++c){
		      	for(int r=0; r < Current_wrk[c].size(); ++r){
		      		Ca_group ca( Current_wrk[c][r] );
		      			    		clipper::String key = std::to_string(ca.coord_ca().x())+std::to_string(ca.coord_ca().y())+std::to_string(ca.coord_ca().z())+std::to_string(ca.coord_n().x())+std::to_string(ca.coord_n().y())+std::to_string(ca.coord_n().z())+std::to_string(ca.coord_c().x())+std::to_string(ca.coord_c().y())+std::to_string(ca.coord_c().z());







		      		if(seqscore.count( std::make_pair(key,propertylabel) )==0){// not found, then save it
		      			seqscore[std::make_pair(key,propertylabel)]=Current_wrk[c][r];
		      		}

		      	}
	  }



}
