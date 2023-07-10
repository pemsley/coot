/*
 * buccaneer-build-cyc.h
 *
 *  Created on: 12 Mar 2021
 *      Author: Emad Alharbi
 */

#ifndef SRC_BUCCANEER_BUCCANEER_NN_BUILD_CYC_H_
#define SRC_BUCCANEER_BUCCANEER_NN_BUILD_CYC_H_
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
class runcyc {
public:
	runcyc(int nfrag ,clipper::Resolution resol ,clipper::MiniMol &mol_wrk,KnownStructure knownstruc,clipper::Xmap<float>   xwrk, LLK_map_target llktgt,Ca_find::TYPE findtype,int modelindex,std::vector<LLK_map_target> llkcls,clipper::MMoleculeSequence seq_wrk,double seq_rel,clipper::String newrestype, int cyc , int ncyc, int startstep, int endstep,clipper::String &mollog,Ca_find cafind,
			std::map<std::pair<clipper::String,clipper::String>, clipper::MMonomer> &seqscore, bool saveseqscore, bool rama_grow);
	void save_sequence_score(clipper::MiniMol &Current_wrk, clipper::String propertylabel,std::map<std::pair<clipper::String,clipper::String>, clipper::MMonomer> &seqscore);
	void apply_sequence_score(clipper::MiniMol &Current_wrk, clipper::String propertylabel,std::map<std::pair<clipper::String,clipper::String>, clipper::MMonomer> &seqscore);
	virtual ~runcyc();
};

#endif /* SRC_BUCCANEER_BUCCANEER_NN_BUILD_CYC_H_ */
