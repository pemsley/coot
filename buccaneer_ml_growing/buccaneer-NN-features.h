/*
 * buccaneer-NN-features.h
 *
 *  Created on: 18 Sep 2020
 *      Author: Emad Alharbi
 */

#ifndef FEATURES_H_
#define FEATURES_H_
#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper.h>
#include "buccaneer-prep.h"
class Score{
public:
	double max;
		double mean;
		double total;
		double lowest;
		double sd;
	Score(double Max , double Mean, double Total , double Lowest, double SD){
		max=Max;
		mean=Mean;
		total=Total;
		lowest=Lowest;
		sd=SD;
	}


};
class Features {
public:
	Features();
	virtual ~Features();
	Score * ramachandran_score( clipper::MPolymer mol);
	Score * llk_score(clipper::MPolymer fragment , LLK_map_target llktgt , clipper::Xmap<float>   xwrk);
	Score * rmsd(clipper::MPolymer fragment,int nresf);
	Score * main_chain_densities_Score(clipper::MPolymer fragment ,  clipper::Xmap<float>   xwrk);

private:
	double standard_deviation(std::vector<float> num,double mean);


};

#endif /* FEATURES_H_ */
