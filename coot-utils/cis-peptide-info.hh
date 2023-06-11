
#ifndef COOT_UTILS_CIS_PEPTIDE_INFO_HH
#define COOT_UTILS_CIS_PEPTIDE_INFO_HH

#include "mini-mol/atom-quads.hh"
#include "geometry/residue-and-atom-specs.hh"

namespace coot {
   namespace util {

      class cis_peptide_info_t {
      public:

	 int serial_number;
	 std::string chain_id_1;
	 std::string residue_name_1;
	 int resno_1;
	 std::string ins_code_1;
	 std::string chain_id_2;
	 std::string residue_name_2;
	 int resno_2;
	 std::string ins_code_2;
	 int model_number;
	 float omega_torsion_angle; // in degrees.

	 // normal constructor used by count_cis_peptides():
	 cis_peptide_info_t(const std::string &chain_id,
			    residue_spec_t res1,
			    residue_spec_t res2,
			    int model_number_in,
			    float tors_in) :
            chain_id_1(chain_id),
            ins_code_1(res1.ins_code),
            chain_id_2(chain_id),
            ins_code_2(res2.ins_code)
         {
	    model_number = model_number_in;
	    serial_number = -1; // unset
	    resno_1 = res1.res_no;
	    resno_2 = res2.res_no;
	    omega_torsion_angle = tors_in;
	 }

	 // Full constructor
	 cis_peptide_info_t(int serial_number_in,
			    const std::string &chain_id_1_in,
			    const std::string &residue_name_1_in,
			    int resno_1_in,
			    const std::string &ins_code_1_in,
			    const std::string &chain_id_2_in,
			    const std::string &residue_name_2_in,
			    int resno_2_in,
			    const std::string &ins_code_2_in,
			    int model_number_in,
			    float omega_torsion_angle_in) :
            serial_number(serial_number_in),
            chain_id_1(chain_id_1_in),
            residue_name_1(residue_name_1_in),
            ins_code_1(ins_code_1_in),
            chain_id_2(chain_id_2_in),
            residue_name_2(residue_name_2_in),
            ins_code_2(ins_code_2_in)
         {

	    serial_number = serial_number_in;
	    resno_1 = resno_1_in;
	    resno_2 = resno_2_in;
	    model_number = model_number_in;
	    omega_torsion_angle = omega_torsion_angle_in;
	 }

	 // Full from mmdb structure
	 explicit cis_peptide_info_t(mmdb::CisPep *cis) :
            serial_number(cis->serNum),
            chain_id_1(cis->chainID1),
            residue_name_1(cis->pep1),
            resno_1(cis->seqNum1),
            ins_code_1(cis->icode1),
            chain_id_2(cis->chainID2),
            residue_name_2(cis->pep2),
            resno_2(cis->seqNum2),
            ins_code_2(cis->icode2),
            model_number(cis->modNum),
            omega_torsion_angle(cis->measure) {}

	 std::string string() const;

	 bool operator==(const cis_peptide_info_t &a) {
	    bool r = 0;

	    // The model number in pdb files usually bogus because of badly formed CISPEP cards

//  	    std::cout << "comparing "  // << model_number
//  		      << " :" << chain_id_1 << ": " << resno_1 << " :" << ins_code_1 << ": :"
//  		      <<         chain_id_2 << ": " << resno_2 << " :" << ins_code_2 << ": "
//  		      << "\nto\n"
//  		      << "          " // << a.model_number
//  		      << " :" << a.chain_id_1 << ": " << a.resno_1 << " :" << a.ins_code_1 << ": :"
//  		      <<         a.chain_id_2 << ": " << a.resno_2 << " :" << a.ins_code_2 << ": "
//  		      << std::endl;
	    // if (a.model_number == model_number) {
	       if (a.chain_id_1 == chain_id_1) {
		  if (a.chain_id_2 == chain_id_2) {
		     if (a.resno_1 == resno_1) {
			if (a.resno_2 == resno_2) {
			   if (a.ins_code_1 == ins_code_1) {
			      if (a.ins_code_2 == ins_code_2) {
				 r = 1;
			      }
			   }
			}
		     }
		  }
	       }
            // }
	    return r;
	 }

      };

      // what type of cis-peptide is this?
      // twisted? pre-pro? non-pre-pro-cis
      // Of course "twisted" is not strictly cis.
      // perhaps this function should be called unorthodox_peptide_torsion_quad_info_t
      //
      class cis_peptide_quad_info_t {
      public:
	 enum type_t { UNSET_TYPE, CIS, PRE_PRO_CIS, TWISTED_TRANS };
	 atom_quad quad;
	 atom_index_quad index_quad;
	 type_t type;
	 cis_peptide_quad_info_t(const atom_quad &q, const atom_index_quad &iq, type_t t_in) :
	    quad(q), index_quad(iq), type(t_in) {}
      };



   }
}

#endif // COOT_UTILS_CIS_PEPTIDE_INFO_HH

