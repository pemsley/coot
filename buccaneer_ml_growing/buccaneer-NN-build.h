/*
 * buccaneer-NN-build.h
 *
 *  Created on: 16 May 2021
 *      Author: Emad Alharbi
 */

#ifndef SRC_BUCCANEER_BUCCANEER_NN_BUILD_H_
#define SRC_BUCCANEER_BUCCANEER_NN_BUILD_H_
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include "simulate-lib.h"
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
#include <map>
class NNbuild {
public:
	NNbuild(int nfrag ,clipper::Resolution resol ,clipper::MiniMol &mol_wrk,KnownStructure knownstruc,clipper::Xmap<float>   xwrk, LLK_map_target llktgt,Ca_find::TYPE findtype,int modelindex,std::vector<LLK_map_target> llkcls,clipper::MMoleculeSequence seq_wrk,double seq_rel,clipper::String newrestype, int cyc, int ncyc,clipper::MiniMol mol_mr,int freeflag,clipper::MiniMol  mol_wrk_in,clipper::String ipfile,clipper::String ipcolfo,clipper::String ipcolfree,Ca_find cafind, double reso,double model_num, clipper::String threshold_select_method,bool nn_confirmation_cycles, int confirmation_cycles_num, bool rama_grow);
	virtual ~NNbuild();
	std::vector<double> rfactors(clipper::MiniMol mol1,int freeflag,clipper::MiniMol  mol_wrk_in,clipper::MiniMol mol_mr, clipper::String ipfile,clipper::String ipcolfo,clipper::String ipcolfree);
	static bool nn_init;
private:
	std::map<std::pair<clipper::String,clipper::String>, clipper::MMonomer> residues_container;

};

#endif /* SRC_BUCCANEER_BUCCANEER_NN_BUILD_H_ */
