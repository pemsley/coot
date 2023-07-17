/*
 * buccaneer-NN-features.cpp
 *
 *  Created on: 18 Sep 2020
 *      Author: Emad Alharbi
 */

#include "buccaneer-NN-features.h"

#include "protein_db.h"
#include <limits>
#include <cmath>
Features::Features() {
	// TODO Auto-generated constructor stub

}

Features::~Features() {
	// TODO Auto-generated destructor stub
}
Score* Features::llk_score(clipper::MPolymer fragment, LLK_map_target llktgt,
		clipper::Xmap<float> xwrk) {
	double score = 0;
	double better = std::numeric_limits<double>::max(); // because llk in negative and lowest is better
	double worse = std::numeric_limits<double>::max() * -1;
	std::vector<float> num;
	for (int r = 0; r < fragment.size(); ++r) {

		Ca_group ca(fragment[r]);
		float ca_score = llktgt.llk(xwrk, ca.rtop_from_std_ori());

		score += ca_score;
		if (ca_score < better)
			better = ca_score;
		if (ca_score > worse)
			worse = ca_score;

		num.push_back(ca_score);
	}

	Score *st = new Score(better, (score / fragment.size()), score, worse,
			standard_deviation(num, score / fragment.size()));
	return st;
}
Score* Features::ramachandran_score(clipper::MPolymer fragment) {

	// determine if the fragment in allowed or favored
	double favored = 0;
	double allowed = 0;
	double outliers = 0;

	for (int r = 1; r < fragment.size(); ++r) {

		if (r + 1 < fragment.size()) {

			clipper::MMonomer mm(fragment[r]); //residue
			clipper::MMonomer nextmm(fragment[r + 1]); //residue
			clipper::MMonomer premm(fragment[r - 1]); //residue
			clipper::Ramachandran rama1 = clipper::Ramachandran(
					clipper::Ramachandran::All);
			int ca = mm.lookup(" CA ", clipper::MM::ANY);
			int n = mm.lookup(" N  ", clipper::MM::ANY);
			int c = mm.lookup(" C  ", clipper::MM::ANY);

			int canext = nextmm.lookup(" CA ", clipper::MM::ANY);
			int nnext = nextmm.lookup(" N  ", clipper::MM::ANY);
			int cnext = nextmm.lookup(" C  ", clipper::MM::ANY);

			int capre = premm.lookup(" CA ", clipper::MM::ANY);
			int npre = premm.lookup(" N  ", clipper::MM::ANY);
			int cpre = premm.lookup(" C  ", clipper::MM::ANY);

			double psi = clipper::Coord_orth::torsion(mm[n].coord_orth(),
					mm[ca].coord_orth(), mm[c].coord_orth(),
					nextmm[nnext].coord_orth());
			;
			double phi = clipper::Coord_orth::torsion(premm[c].coord_orth(),
					mm[n].coord_orth(), mm[ca].coord_orth(),
					mm[c].coord_orth());
			;

			if (rama1.favored(phi, psi)) {
				favored++;
			} else if (rama1.allowed(phi, psi)) {
				allowed++;
			} else {
				outliers++;
			}

		}
	}

	Score *st = new Score(favored, allowed, outliers, 0, 0);

	return st;

}

Score* Features::rmsd(clipper::MPolymer fragment, int nresf) {

	// calculate RMSD for fragment based on the top 500 fragments database.
	double score = 0.0;
	double highest = 0.0;
	double lowest = 100.0;
	const char *clibdptr = getenv("CLIBD");
	clipper::String prtdb = clipper::String(clibdptr) + "/prot500.db";
	const ProteinDB::ChainDB chaindb(prtdb);
	std::vector<float> num;

	for (int r = 0; r <= fragment.size() - nresf; r++) {
		ProteinDB::Chain frag;
		bool valid = true;
		for (int dr = 0; dr < nresf; dr++) {
			ProteinDB::Residue r1(fragment[r + dr]);
			frag.add_residue(r1);
			valid = valid && (r1.flag() == ProteinDB::Residue::NORMAL);
			if (dr > 0)
				valid =
						valid
								&& ((frag[dr].coord_n() - frag[dr - 1].coord_c()).lengthsq()
										< pow(2.0, 2));
		}
		double rmsd = -1;
		if (valid) {
			ProteinDB::ChainDB fragdb(frag);
			std::vector<ProteinDB::Chain> frags = chaindb.match_fragment(fragdb,
					50);
			rmsd = frag.rmsd(frags[0]);
		}
		score += rmsd;
		num.push_back(rmsd);
		if (rmsd > highest)
			highest = rmsd;
		if (rmsd < lowest)
			lowest = rmsd;
	}

	Score *st = new Score(highest, (score / fragment.size()), score, lowest,
			standard_deviation(num, score / fragment.size()));
	return st;
}
Score* Features::main_chain_densities_Score(clipper::MPolymer fragment,
		clipper::Xmap<float> xwrk) {
	double densitiy_score = 0.0;
	std::vector<float> main_chain_densities_Score =
			ProteinTools::main_chain_densities(fragment, xwrk);
	double max = 0.0;
	double min = std::numeric_limits<double>::max();
	;
	std::vector<float> num;
	for (int i = 0; i < main_chain_densities_Score.size(); ++i) {
		densitiy_score += main_chain_densities_Score[i];
		if (main_chain_densities_Score[i] > max)
			max = main_chain_densities_Score[i];
		if (main_chain_densities_Score[i] < min)
			min = main_chain_densities_Score[i];
		num.push_back(main_chain_densities_Score[i]);
	}

	Score *st = new Score(max, (densitiy_score / fragment.size()), densitiy_score,
			min, standard_deviation(num, densitiy_score / fragment.size()));

	return st;
}

double Features::standard_deviation(std::vector<float> num, double mean) {
// standard deviation
	double SD = 0;
	for (int n = 0; n < num.size(); ++n) {

		SD += pow(num[n] - mean, 2);
	}
	SD = SD / num.size();
	return sqrt(SD);
}

