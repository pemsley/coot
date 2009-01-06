
#include <math.h>
#include <iostream>
#include <fstream>

#include "coot-coord-utils.hh"
#include "bfkurt.hh"


coot_extras::b_factor_analysis::b_factor_analysis(const CMMDBManager *mol_in, bool is_mol_from_shelx_flag_in) { 

   // for each chain, we want a vector of residues, which contain the kertosis
   // 
   CMMDBManager *mol = (CMMDBManager *) mol_in; // ghastly.
   is_mol_from_shelx_flag = is_mol_from_shelx_flag_in;
   
   // recall kurtosis, $k$ of $N$ observations:
   // k = \frac{\Sigma(x_i - \mu)^4} {N \sigma^4} - 3    
   // (x_i - \mu)^4 = x_i^4 + 4x_i^3\mu + 6x_i^2\mu^2 + 4x_i\mu^3 + \mu^4


   if (mol) { 
      CModel *model_p = mol->GetModel(1);
      if (model_p) { 
	 CChain *chain_p;
	 int nchains = model_p->GetNumberOfChains();
	 if (nchains > 0) { 
	    for (int ichain=0; ichain<nchains; ichain++) {
	       chain_p = model_p->GetChain(ichain);
	       if (chain_p) {  
		  std::vector<my_stats_t> residue_prop;
		  std::string chain_id(chain_p->GetChainID());

		  int nres = chain_p->GetNumberOfResidues();
		  if (nres > 0) {
		     for (int ires=0; ires<nres; ires++) { 
			PCResidue residue_p = chain_p->GetResidue(ires);
			if (residue_p) { 
			   std::string res_name(residue_p->GetResName());
			   if (res_name != "HOH" && 
			       res_name != "WAT") { 
			      int seqnum = residue_p->GetSeqNum();
			      my_stats_t my_stats = stats(residue_p);
			      residue_prop.push_back(my_stats);
			   }
			}
		     }
		     // After all the residues in a chain:
		     // 
		     std::pair<std::string, std::vector<my_stats_t> > 
			named_residue_props(chain_id, residue_prop);
		     kurtoses.push_back(named_residue_props);
		  } else {
		     // we need a kurtoses for every chain to match
		     // display of graphs in graphics-info-graphs'
		     // b_factor_graphs()
		     std::pair<std::string, std::vector<my_stats_t> > 
			named_residue_props(chain_id, residue_prop);
		     kurtoses.push_back(named_residue_props);
		  }
	       }
	    }
	    set_questionable_flags(3.0); // uses kurtoses
	 }
      }
   }
}



void
coot_extras::b_factor_analysis::set_questionable_flags(float z) {

   double sum=0.0;
   double sum_sq=0.0;
   int n = 0;
   double statistic;
   
   for(unsigned int ichain=0; ichain<kurtoses.size(); ichain++) {
      for(unsigned int ires=0; ires<kurtoses[ichain].second.size(); ires++) {
	 if (kurtoses[ichain].second[ires].n > 1) { 
	    statistic = kurtoses[ichain].second[ires].std_dev;
	    sum += statistic;
	    sum_sq += statistic*statistic;
	    n++;
	 }
      }
   }

   if (n>1) {
      double mean = sum/double(n);
      double var  = sum_sq/double(n) - mean*mean;
      double std_dev = sqrt(var);

      // std::cout << "mean: " << mean << " from " << sum << "/" << n << std::endl;

      for(unsigned int ichain=0; ichain<kurtoses.size(); ichain++) {
	 for(unsigned int ires=0; ires<kurtoses[ichain].second.size(); ires++) {
	    if (kurtoses[ichain].second[ires].n > 1) { 
	       statistic = kurtoses[ichain].second[ires].std_dev;
	       if (statistic > mean + z * std_dev) {
		  kurtoses[ichain].second[ires].questionable_flag = 1;
	       }
	    } 
	 }
      }
   }
} 

coot_extras::my_stats_t 
coot_extras::b_factor_analysis::stats(CResidue *residue_p) const { 

   coot_extras::my_stats_t my_stats;
   my_stats.mean = 0;
   my_stats.std_dev = 0;
   my_stats.n = 0;

   PPCAtom residue_atoms;
   int nResidueAtoms;

   double running_sum; 
   double running_sum_2;
   double running_sum_3;
   double running_sum_4;
   
   double mean = 0, variance = 0, std_dev = 0, kt, kurtosis = 0;

   double bf;
   double occ;
   double bfo;
   double occ_sum = 0.0;
   
   running_sum   = 0.0;
   running_sum_2 = 0.0;
   running_sum_3 = 0.0;
   running_sum_4 = 0.0;
   
   residue_p->GetAtomTable(residue_atoms, nResidueAtoms);
   if (nResidueAtoms > 0) { 
      for (int i=0; i<nResidueAtoms; i++) {
	 std::string ele = residue_atoms[i]->element;
	 if ((ele != " H") && (ele != " D")) {
	    bf = residue_atoms[i]->tempFactor;
	    occ = residue_atoms[i]->occupancy;
	    // ignore atoms with silly (or shelx?) B factors and occs
	    if (((bf > 0.0) && (occ >= 0.0) && (occ <= 1.0)) ||
		(is_mol_from_shelx_flag && (occ < 11.001) && (occ>10.999))) {

	       if (is_mol_from_shelx_flag)
		  occ = 1.0; // only accepted atoms with occ 11.0
	       
	       occ_sum += occ;
	       bfo = bf * occ;
	       running_sum   += bfo;
	       running_sum_2 += bfo*bfo;
	       running_sum_3 += bfo*bfo*bfo;
	       running_sum_4 += bfo*bfo*bfo*bfo;
	    }
	 }
      }
      CAtom *intel_at = coot::util::intelligent_this_residue_mmdb_atom(residue_p);
      my_stats.atom_name = intel_at->name;
      double div = occ_sum;
      if (div > 0) { 
	 mean = running_sum / div;
	 variance = running_sum_2/div - mean*mean;
//  	 std::cout << "Variance: " << residue_p->GetSeqNum() << " "
//  		   << residue_p->GetChainID() << " " << variance
// 		   << " natoms: " << div << std::endl;
	 if (variance < 0)
	    variance = 0; // fixes numerical instability.
	 std_dev = sqrt(variance);

	 // Don't bother calculating kurtosis for less than 2 atoms
	 // 
	 if (nResidueAtoms > 1) { 

	    kt = running_sum_4 
	       - 4.0*running_sum_3*mean
	       + 6.0*running_sum_2*mean*mean
	       - 4.0*running_sum*mean*mean*mean
	       + mean*mean*mean*mean*div;
	 
	    kurtosis = kt/(div*variance*variance) 
	       - 3.0;

	 } else {
	    // some useless number
	    kurtosis = -999.9;
	 }
//       } else {
// 	 std::cout << "DEBUG:: ignoring this residue with no no-zero occ atoms"
// 		   << std::endl;
      } 
      my_stats.resname = residue_p->GetResName();
      my_stats.mean = mean;
      my_stats.std_dev = std_dev;
      my_stats.kurtosis = kurtosis;
      my_stats.n = nResidueAtoms; // should be div, perhaps.
      // and the data from the molecules residue:
      my_stats.resno = residue_p->GetSeqNum();
      my_stats.inscode = residue_p->GetInsCode();
   }
   return my_stats;
} 


short int
coot_extras::b_factor_analysis::write_table(const std::string &filename, 
					    const std::string &pdb_filename,
					    short int only_questionables_flag) const { 

   std::ofstream outfile(filename.c_str());

   if (!outfile) { 
      std::cout << "Cannot open output file" << std::endl;
   } else { 

      outfile << "<validation>\n";
      outfile << "   <date>20031029</date>\n";
      outfile << "   <validation-program>bfactan</validation-program>\n";
      outfile << "   <title>Validated by bfactan</title>\n";
      outfile << "   <bfactan-info version=\"0.0\" />\n";
      outfile << "   <chain-list>\n";

      for(unsigned int ichain=0; ichain<kurtoses.size(); ichain++) {
	 outfile << "      <chain>\n";
	 outfile << "          <chain-id>" 
		   << kurtoses[ichain].first << "</chain-id>\n";
	 outfile << "          <residue-list>\n";
	 for(unsigned int ires=0; ires<kurtoses[ichain].second.size(); ires++) {
	    if (kurtoses[ichain].second[ires].n > 0) {

	       if (kurtoses[ichain].second[ires].questionable_flag ||
		   only_questionables_flag == 0) { 
		  outfile << "             <residue>\n";
		  outfile << "                <sequence-number>"
			  << kurtoses[ichain].second[ires].resno
			  << "</sequence-number>\n";
		  // insertion code?
		  if (kurtoses[ichain].second[ires].inscode != "") {
		     outfile << "                <insertion-code>"
			     << kurtoses[ichain].second[ires].inscode
			     << "</insertion-code>\n";
		  } 
		  outfile << "                 <residue-temperature-factor-outlier>\n";
		  outfile << "                     <b-factor-mean>";
		  outfile << kurtoses[ichain].second[ires].mean;
		  outfile << "</b-factor-mean>\n";
		  
		  if (kurtoses[ichain].second[ires].n > 1) {
		     
		     outfile << "                     <b-factor-standard-deviation>";
		     outfile << kurtoses[ichain].second[ires].std_dev;
		     outfile << "</b-factor-standard-deviation>\n";
		     
		     outfile << "                     <b-factor-kurtosis>";
		     outfile << kurtoses[ichain].second[ires].kurtosis;
		     outfile << "</b-factor-kurtosis>\n";
		  }
		  outfile << "                 </residue-temperature-factor-outlier>\n";
		  outfile << "             </residue>\n";
	       }
	    }
	 }
	 outfile << "          </residue-list>\n";
	 outfile << "      </chain>\n";
      }
      outfile << "   </chain-list>\n";
      outfile << "</validation>\n";
   }
   return 1;
}


std::vector<coot_extras::my_chain_of_stats_t>
coot_extras::b_factor_analysis::chain_details() const {

   std::vector<coot_extras::my_chain_of_stats_t> v;
   coot_extras::my_chain_of_stats_t c;

   for(unsigned int ichain=0; ichain<kurtoses.size(); ichain++) {
      c.residue_properties = kurtoses[ichain].second;
      c.chain_id            =kurtoses[ichain].first;
      v.push_back(c);
   }
   return v;
}
