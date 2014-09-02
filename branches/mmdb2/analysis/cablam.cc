

#include <map>
#include <clipper/core/coords.h>

#include "coot-utils/coot-coord-utils.hh"
#include "cablam.hh"

clipper::Coord_orth
coot::cablam::get_closest_CA_CA_approach(const coot::torsion_atom_quad &quad) const {

   clipper::Coord_orth CA_p(quad.atom_1->x, quad.atom_1->y, quad.atom_1->z);
   clipper::Coord_orth CA_t(quad.atom_2->x, quad.atom_2->y, quad.atom_2->z);
   clipper::Coord_orth  O_t(quad.atom_4->x, quad.atom_4->y, quad.atom_4->z);

   clipper::Coord_orth PT = CA_t - CA_p;
   clipper::Coord_orth PT_unit(PT.unit());
   clipper::Coord_orth PO =  O_t - CA_p;

   double PO_length = sqrt(PO.lengthsq());
   
   double cos_alpha_p = (PT * PO)/(PO_length * sqrt(PT.lengthsq()));

   double PC_length = cos_alpha_p * PO_length;
   clipper::Coord_orth Closest_approach = CA_p + PC_length * PT_unit;
   return Closest_approach;
} 


coot::cablam::cablam(PCResidue *residues, int n_sel_residues) {

   std::map<CResidue *, torsion_atom_quad> residue_quads;

   for (unsigned int ires=1; ires<n_sel_residues; ires++) {

      if (ires < n_sel_residues) { 
	 CResidue *res_p = residues[ires-1];
	 CResidue *res_t = residues[ires];

	 // set the 4 atoms for a given residue n:  (n-1)Ca, nCA, nC, nO
	 // in that order
	 // 
	 CAtom *n_p_CA = res_p->GetAtom(" CA ");
	 CAtom *n_t_CA = res_t->GetAtom(" CA ");
	 CAtom *n_t_C  = res_t->GetAtom(" C  ");
	 CAtom *n_t_O  = res_t->GetAtom(" O  ");

	 if (n_p_CA && n_t_CA && n_t_C && n_t_O) {

	    // Is it plausibly a trans torsion? (measure CA-CA)
	    clipper::Coord_orth pt_1(n_p_CA->x, n_p_CA->y, n_p_CA->z);
	    clipper::Coord_orth pt_2(n_t_CA->x, n_t_CA->y, n_t_CA->z);
	    double d = clipper::Coord_orth::length(pt_1, pt_2);

	    // if so, add it then.
	    if (d<4 && d>3.5) { 
	       atom_quad a(n_p_CA,n_t_CA,n_t_C,n_t_O);
	       residue_quads[res_t] = torsion_atom_quad(n_p_CA,n_t_CA,n_t_C,n_t_O, 0, a.torsion(), 1);
	    }
	 }
      }
   }

   for (unsigned int ires=1; ires<n_sel_residues; ires++) {
      std::map<CResidue *, torsion_atom_quad>::const_iterator it_this;
      std::map<CResidue *, torsion_atom_quad>::const_iterator it_prev;

      it_prev = residue_quads.find(residues[ires-1]);
      it_this = residue_quads.find(residues[ires]);

      if (it_prev != residue_quads.end()) { 
	 if (it_this != residue_quads.end()) {
	    // make pseudo position Cl_ap
	    clipper::Coord_orth Cl_ap_p = get_closest_CA_CA_approach(it_prev->second);
	    clipper::Coord_orth Cl_ap_t = get_closest_CA_CA_approach(it_this->second);
	    clipper::Coord_orth O_prev(it_prev->second.atom_4->x,
				       it_prev->second.atom_4->y,
				       it_prev->second.atom_4->z);
	    clipper::Coord_orth O_this(it_this->second.atom_4->x,
				       it_this->second.atom_4->y,
				       it_this->second.atom_4->z);
	    double clipper_torsion = clipper::Coord_orth::torsion(O_prev, Cl_ap_p, Cl_ap_t, O_this);
	    int sse_prev = it_prev->second.atom_4->residue->SSE;
	    int sse_this = it_this->second.atom_4->residue->SSE;
	    std::string sp = coot::util::sse_to_string(sse_prev);
	    std::string st = coot::util::sse_to_string(sse_this);
	    std::cout << " torsion " << clipper::Util::rad2d(clipper_torsion) << " "
		      << sse_prev << " " << sse_this << " " << sp << " " << st << std::endl;
	 }
      }
   }
} 
