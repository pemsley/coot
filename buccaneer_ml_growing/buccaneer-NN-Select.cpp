/*
 * buccaneer-NN-Select.cpp
 *
 *  Created on: 8 Dec 2020
 *      Author: Emad Alharbi
 */

#include "buccaneer-NN-Select.h"

#include "buccaneer-NN-features.h"
#include <algorithm>

#include "buccaneer-NN-model.h"
LLK_map_target SelectNN::llktgt;
clipper::Xmap<float> SelectNN::xwrk;
double SelectNN::reso;
std::string SelectNN::model;
double SelectNN::threshold;
bool SelectNN::select;
std::vector<float> SelectNN::predictions;
bool SelectNN::update_predictions;
std::vector<float> SelectNN::moving_thresholds;
double SelectNN::model_num;
clipper::String SelectNN::threshold_select_method;

void SelectNN::predict(std::vector<Ca_join::Tri_residue> &fragments,
		clipper::MiniMol mol_wrk, double rjoin) {
	// TODO Auto-generated constructor stub

	if (SelectNN::select == false) // if NN keyword not set
		return;

	if (SelectNN::update_predictions == false && SelectNN::threshold>SelectNN::moving_thresholds.size()-1){ // when number of the thresholds are less than the number of models that we want to build
		SelectNN::threshold=-1;

		return;
	}

	std::vector<int> NonZeroFrg;
	for (int f1 = 0; f1 < fragments.size(); ++f1) // save the index of a fragment that its flag not zero
		if (fragments[f1].flag != 0)
			NonZeroFrg.push_back(f1);

/*
	if (NonZeroFrg.size() > 2849) {
			std::cout
					<< "NN: Number of fragments higher than the maximum! (maximum 2849)"
					<< std::endl;
			SelectNN::threshold=-1;//set to -1 to stop
			return;
		} // can not predict more than 2849  fragments
*/


	if (SelectNN::update_predictions == true) {

		SelectNN::predictions.clear(); // clear if contains results from previous cycle

		// make a list of joins
		std::vector<Ca_join::Node> joins(fragments.size());
		for (int f1 = 0; f1 < fragments.size(); f1++)
			if (fragments[f1].flag != 0) {
				joins[f1].score = 1.0;
				clipper::Coord_frac cx1 = fragments[f1].res[1][1].coord_frac(
						mol_wrk.cell());
				clipper::Coord_frac cx2 = fragments[f1].res[2][1].coord_frac(
						mol_wrk.cell());
				for (int f2 = 0; f2 < fragments.size(); f2++)
					if (fragments[f2].flag != 0) {
						if (f1 != f2) {
							clipper::Coord_frac cy0 =
									fragments[f2].res[0][1].coord_frac(
											mol_wrk.cell());
							clipper::Coord_frac cy1 =
									fragments[f2].res[1][1].coord_frac(
											mol_wrk.cell());
							clipper::Coord_frac cy2 =
									fragments[f2].res[2][1].coord_frac(
											mol_wrk.cell());
							cy0 = cy0.symmetry_copy_near(mol_wrk.spacegroup(),
									mol_wrk.cell(), cx1);
							cy1 = cy1.symmetry_copy_near(mol_wrk.spacegroup(),
									mol_wrk.cell(), cx1);
							cy2 = cy2.symmetry_copy_near(mol_wrk.spacegroup(),
									mol_wrk.cell(), cx1);
							if ((cx1 - cy0).lengthsq(mol_wrk.cell()) < rjoin
									&& (cx2 - cy1).lengthsq(mol_wrk.cell())
											< rjoin) {
								if (fragments[f1].flag == 1
										&& fragments[f2].flag == 1)
									joins[f1].ptrs.push_back(f2);
								else if (f2 == f1 + 1)
									joins[f1].ptrs.push_back(f2);
							}
						}
					}
			}

		//visualize into 2d graph
		// for(int j=0 ; j < joins.size();++j){
		//	  for(int child=0 ; child < joins[j].ptrs.size();++child)
		//		  std::cout<<j<<","<<joins[j].ptrs[child]<<",1"<<std::endl;
		//  }

		// Find which Ca-alpha can be linked to which based on the list of joins
		std::vector<std::vector<float>> links;
		for (int j = 0; j < joins.size(); ++j) {

			std::vector<float> tmpdistance(joins.size(), 0);

			for (int child = 0; child < joins[j].ptrs.size(); ++child) {

				tmpdistance[joins[j].ptrs[child]] = 1;

			}
			links.push_back(tmpdistance);

		}

		const int tensor_rows = 2849; // number of max fragments in a model
		//const int tensor_rows = 50; // number of max fragments in a model
		const int tensor_cols = 14; // number of features

		//float t[NonZeroFrg.size() * tensor_cols]; // limit the input to the number of fragments to speed up the NN calculation
		float t[tensor_rows * tensor_cols]; // limit the input to the number of fragments to speed up the NN calculation

		int tensorcounter = 0;

		// Normalised resolution here to avoid duplicate  normalisation in the for loop
		double reso = SelectNN::reso;


		reso = (reso - 2.27850202) / 0.65511165;
		int r = 0;

		for (int f1 = 0; f1 < fragments.size(); ++f1) {

			if (fragments[f1].flag != 0) {

				bool startRes = false;
				bool endRes = true;
				int OnlyOneLink = 0;

				for (int dd = 0; dd < links[f1].size(); ++dd) {
					if (links[dd][f1] != 0) {
						OnlyOneLink++;

					}
					if (links[f1][dd] != 0) {
						endRes = false;

					}

				}
				if (OnlyOneLink == 0) {

					startRes = true;
				}

				bool midRes = false;
				if (startRes == false && endRes == false)
					midRes = true;


				clipper::MPolymer frg;
				Ca_chain chain;
				std::vector<Ca_chain> chains;
				Ca_group ca1(fragments[f1].res[0][0], fragments[f1].res[0][1],
						fragments[f1].res[0][2]);
				chain.push_back(ca1);
				Ca_group ca2(fragments[f1].res[1][0], fragments[f1].res[1][1],
						fragments[f1].res[1][2]);
				chain.push_back(ca2);
				Ca_group ca3(fragments[f1].res[2][0], fragments[f1].res[2][1],
						fragments[f1].res[2][2]);
				chain.push_back(ca3);
				chains.push_back(chain);
				clipper::MiniMol Current_wrk(mol_wrk.spacegroup(),
						mol_wrk.cell());

				ProteinTools::insert_ca_chains(Current_wrk, chains);

				Features *features = new Features();

				Score *RMSD = features->rmsd(Current_wrk[0], 3);

				Score *llk = features->llk_score(Current_wrk[0],
						SelectNN::llktgt, SelectNN::xwrk);

				Score *densitiyScore = features->main_chain_densities_Score(
						Current_wrk[0], SelectNN::xwrk);

				Score *rma = features->ramachandran_score(Current_wrk[0]);

				rma->max = (rma->max - 0.72345617) / 0.44728888;

				t[tensorcounter] = rma->max;

				tensorcounter++;

				t[tensorcounter] = reso;
				tensorcounter++;

				llk->mean = (llk->mean - -0.64929002) / 0.11796437;

				t[tensorcounter] = llk->mean;
				tensorcounter++;

				llk->max = (llk->max - -0.71953489) / 0.11392434;

				t[tensorcounter] = llk->max;
				tensorcounter++;

				densitiyScore->mean = (densitiyScore->mean - 0.48943167)
						/ 1.02339344;

				t[tensorcounter] = densitiyScore->mean;
				tensorcounter++;

				llk->lowest = (llk->lowest - -0.5777621) / 0.13567325;

				t[tensorcounter] = llk->lowest;
				tensorcounter++;

				rma->mean = (rma->mean - 0.25391793) / 0.43525121;

				t[tensorcounter] = rma->mean;
				tensorcounter++;

				llk->sd = (llk->sd - 0.06085652) / 0.03867459;

				t[tensorcounter] = llk->sd;
				tensorcounter++;

				densitiyScore->sd = (densitiyScore->sd - 0.08716381)
						/ 0.21006441;

				t[tensorcounter] = densitiyScore->sd;
				tensorcounter++;

				densitiyScore->lowest = (densitiyScore->lowest - 0.38655551)
						/ 0.85745071;

				t[tensorcounter] = densitiyScore->lowest;
				tensorcounter++;



				RMSD->max = (RMSD->max - 0.0036608) / 0.01043267;

				t[tensorcounter] = RMSD->max;
				tensorcounter++;




				double start = 0;
				if (startRes == true)
					start = 1;

				double end = 0;
				if (endRes == true)
					end = 1;

				double mid = 0;
				if (midRes == true)
					mid = 1;

				start = (start - 0.08658568) / 0.28122696;

				end = (end - 0.07891443) / 0.26960517;

				mid = (mid - 0.84111921) / 0.36556489;

				t[tensorcounter] = start;
				tensorcounter++;

				t[tensorcounter] = end;
				tensorcounter++;

				t[tensorcounter] = mid;
				tensorcounter++;

				r++;




			}
		if(r==tensor_rows || (f1+1==fragments.size() && r!=0)){ // r!=0 in case no more new fragments were not predicted

				// Set up NN weights and bias
				float *lstm_output_array;
				float *lstm_kernel_array;
				float *lstm_recurrent_kernel_array;
				float *lstm_bias_array;
				float *lstm_1_output_array;
				float *lstm_1_kernel_array;
				float *lstm_1_recurrent_kernel_array;
				float *lstm_1_bias_array;
				float *lstm_2_output_array;
				float *lstm_2_kernel_array;
				float *lstm_2_recurrent_kernel_array;
				float *lstm_2_bias_array;
				float *lstm_3_output_array;
				float *lstm_3_kernel_array;
				float *lstm_3_recurrent_kernel_array;
				float *lstm_3_bias_array;
				float *lstm_4_output_array;
				float *lstm_4_kernel_array;
				float *lstm_4_recurrent_kernel_array;
				float *lstm_4_bias_array;
				float *dense_kernel_array;
				float *dense_bias_array;

				nnmodel::nnmodel_initialize(&lstm_output_array, &lstm_kernel_array,
						&lstm_recurrent_kernel_array, &lstm_bias_array,
						&lstm_1_output_array, &lstm_1_kernel_array,
						&lstm_1_recurrent_kernel_array, &lstm_1_bias_array,
						&lstm_2_output_array, &lstm_2_kernel_array,
						&lstm_2_recurrent_kernel_array, &lstm_2_bias_array,
						&lstm_3_output_array, &lstm_3_kernel_array,
						&lstm_3_recurrent_kernel_array, &lstm_3_bias_array,
						&lstm_4_output_array, &lstm_4_kernel_array,
						&lstm_4_recurrent_kernel_array, &lstm_4_bias_array,
						&dense_kernel_array, &dense_bias_array,
						SelectNN::model.c_str());
		k2c_tensor nninput = { &t[0], 2, size_t(r * tensor_cols), {
				size_t(r), 14, 1, 1, 1 } };



		float nnputarr[tensor_rows] = { 0 };
		k2c_tensor nnoutput = { &nnputarr[0], 2, tensor_rows, { tensor_rows, 1, 1, 1, 1 } };

		// predict
		nnmodel *select = new nnmodel(&nninput, &nnoutput, lstm_output_array,
						lstm_kernel_array, lstm_recurrent_kernel_array, lstm_bias_array,
						lstm_1_output_array, lstm_1_kernel_array,
						lstm_1_recurrent_kernel_array, lstm_1_bias_array,
						lstm_2_output_array, lstm_2_kernel_array,
						lstm_2_recurrent_kernel_array, lstm_2_bias_array,
						lstm_3_output_array, lstm_3_kernel_array,
						lstm_3_recurrent_kernel_array, lstm_3_bias_array,
						lstm_4_output_array, lstm_4_kernel_array,
						lstm_4_recurrent_kernel_array, lstm_4_bias_array,
						dense_kernel_array, dense_bias_array, r);//r number of fragments to send to the NN



		//save predictions to use when test another threshold
		for (int i = 0; i < r; ++i) {
			SelectNN::predictions.push_back(nnoutput.array[i]);


		}


		//clear up the NN  weights and bias
				nnmodel::nnmodel_terminate(lstm_output_array, lstm_kernel_array,
						lstm_recurrent_kernel_array, lstm_bias_array,
						lstm_1_output_array, lstm_1_kernel_array,
						lstm_1_recurrent_kernel_array, lstm_1_bias_array,
						lstm_2_output_array, lstm_2_kernel_array,
						lstm_2_recurrent_kernel_array, lstm_2_bias_array,
						lstm_3_output_array, lstm_3_kernel_array,
						lstm_3_recurrent_kernel_array, lstm_3_bias_array,
						lstm_4_output_array, lstm_4_kernel_array,
						lstm_4_recurrent_kernel_array, lstm_4_bias_array,
						dense_kernel_array, dense_bias_array);
				free(select);

				tensorcounter = 0;
				r=0;

				for (int t_input=0; t_input < tensor_rows * tensor_cols ; t_input++)
					t[t_input]=0;


			}


				}

		SelectNN::calculatethresholds(NonZeroFrg);
	}




// Now set the fragments that have probability  lower than threshold to zero to not use in joining step
	int countremovedfrg=0;
	for (int i = 0; i < NonZeroFrg.size(); ++i) {

		float pre = SelectNN::predictions[i];


		if (pre < SelectNN::moving_thresholds[SelectNN::threshold] ) {
			fragments[NonZeroFrg[i]].flag = 0;
			countremovedfrg++;

		}




	}


}

void SelectNN::calculatethresholds(std::vector<int> frgs){
	moving_thresholds.clear();


	if(SelectNN::threshold_select_method == "NONE")
	 SelectNN::fixedthresholds(frgs);
	if(SelectNN::threshold_select_method == "freedmandiaconisrule")
	 SelectNN::freedmandiaconisrule(frgs);
	if(SelectNN::threshold_select_method == "defaultthreshold")// default threshold is 0.5
	 SelectNN::defaultthreshold(frgs);

}
void SelectNN::defaultthreshold(std::vector<int> frgs){

	float maxprobability=0;
	float minprobability=1;

  SelectNN::minmaxprobability(frgs,maxprobability,minprobability);


  moving_thresholds.push_back(minprobability);
  moving_thresholds.push_back(0.5);
}
void SelectNN::freedmandiaconisrule(std::vector<int> frgs){
		//IQR
		std::vector<float> thresholds;
		float maxprobability=0;
		float minprobability=1;

		SelectNN::minmaxprobability(frgs,maxprobability,minprobability);


		std::vector<float> predictions_copy(SelectNN::predictions);
		std::sort (predictions_copy.begin(), predictions_copy.end());

		 //https://stackoverflow.com/questions/11964552/finding-quartiles
		 auto const Q1 = predictions_copy.size() / 4;
		 auto const Q2 = predictions_copy.size() / 2;
		 auto const Q3 = Q1 + Q2;

		 std::nth_element(predictions_copy.begin(),          predictions_copy.begin() + Q1, predictions_copy.end());
		 std::nth_element(predictions_copy.begin() + Q1 + 1, predictions_copy.begin() + Q2, predictions_copy.end());
		 std::nth_element(predictions_copy.begin() + Q2 + 1, predictions_copy.begin() + Q3, predictions_copy.end());


		 double IQR=predictions_copy[Q3] - predictions_copy[Q1];
		 double h=2*(IQR/std::cbrt(frgs.size()));

		 float range=h;

		 double threshold=0;
		 for(double th=0.0; th < SelectNN::model_num+1 ; ++th){

		 threshold=(range*th)+minprobability;
		 thresholds.push_back(threshold);

		 }
		 SelectNN::removeuselessthresholds(frgs,thresholds);

}
 void SelectNN::minmaxprobability(std::vector<int> frgs,float &maxprobability,float &minprobability){
	  maxprobability=0;
	  minprobability=1;

	 				//find max and min
	 				for (int i = 0; i < frgs.size(); ++i){
	 					if(SelectNN::predictions[i] > maxprobability)
	 					maxprobability=SelectNN::predictions[i];
	 					if(SelectNN::predictions[i] < minprobability  )
	 					minprobability=SelectNN::predictions[i];

	 				}
}
void SelectNN::fixedthresholds(std::vector<int> frgs){
	std::vector<float> thresholds;
		float maxprobability=0;
		float minprobability=1;

		SelectNN::minmaxprobability(frgs,maxprobability,minprobability);
					float range=maxprobability-minprobability;
					double threshold=0;
					for(double th=0.0; th < SelectNN::model_num+1 ; ++th){
						threshold=((range/SelectNN::model_num)*th)+minprobability;

						thresholds.push_back(threshold);

					}
				SelectNN::removeuselessthresholds(frgs,thresholds);
}
 void SelectNN::removeuselessthresholds(std::vector<int> frgs,std::vector<float> thresholds){

	 // useless threshold that remove the same number of fragments as the previous threshold

	 int countusethresholds=0;
	 for(int i=0 ; i < thresholds.size();++i){


	 			int removedfrg=SelectNN::countremovedfrg(frgs,thresholds[i]);

	 			for(int nexth=i+1; nexth < thresholds.size();++nexth){
	 				int removedfrgnext=SelectNN::countremovedfrg(frgs,thresholds[nexth]);
	 				if(removedfrgnext!=removedfrg){
	 					moving_thresholds.push_back(thresholds[nexth-1]);
	 					i=nexth-1;

	 					countusethresholds++;
	 					break;
	 				}
	 			}


	 			}

}
int SelectNN::countremovedfrg(std::vector<int> frgs, double ths){
	int removedfrg=0;
	for (int p = 0; p < frgs.size(); ++p) {

	float pre = SelectNN::predictions[p];

	if (pre < ths ) {

	removedfrg++;
	}


	}

	return removedfrg;
}
SelectNN::~SelectNN() {
	// TODO Auto-generated destructor stub
}

