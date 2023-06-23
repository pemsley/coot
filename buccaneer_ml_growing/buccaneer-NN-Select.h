/*
 * buccaneer-NN-Select.h
 *
 *  Created on: 8 Dec 2020
 *      Author: Emad Alharbi
 */

#ifndef SRC_BUCCANEER_BUCCANEER_NN_SELECT_H_
#define SRC_BUCCANEER_BUCCANEER_NN_SELECT_H_
#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper.h>
#include "buccaneer-prep.h"

#include "buccaneer-join.h"

#include "k2c_include.h"
#include <vector>
class SelectNN {
public:

	 SelectNN(){};
	virtual ~SelectNN();
	static void calculatethresholds(std::vector<int> frgs);
	static int countremovedfrg(std::vector<int> frgs, double ths);
	static  LLK_map_target llktgt;
	static  clipper::Xmap<float>   xwrk;
	static  double reso;
	static 	void predict(std::vector<Ca_join::Tri_residue> &fragments, clipper::MiniMol mol_wrk ,  double rjoin);
	static  std::string model;
	static double threshold;
	static bool select;
	static std::vector<float> predictions;
	static bool update_predictions;
	static std::vector<float> moving_thresholds;
	static double model_num;

	static void removeuselessthresholds(std::vector<int> frgs,std::vector<float> thresholds); // to remove duplicated thresholds
	static void fixedthresholds(std::vector<int> frgs);
	static void freedmandiaconisrule(std::vector<int> frgs);
	static void defaultthreshold(std::vector<int> frgs);

	static void minmaxprobability(std::vector<int> frgs,float &maxprobability,float &minprobability);
	static clipper::String threshold_select_method;

};

#endif /* SRC_BUCCANEER_BUCCANEER_NN_SELECT_H_ */


