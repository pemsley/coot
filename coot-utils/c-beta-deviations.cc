
#include "compat/coot-sysdep.h"

#include "c-beta-deviations.hh"

#include "geometry/residue-and-atom-specs.hh"

#ifdef COOT_ENABLE_WINAPI_SUSPENSION
# undef GetAtomName
#endif // COOT_ENABLE_WINAPI_SUSPENSION

// return a map or a list at some stage.
// I am not sure that we want a map. Maybe we want a vector.
// We the value of the outside map is a map on the alt confs of the residue
// - or maybe we should use the residue spec?
//
std::map<mmdb::Residue *, std::map<std::string, coot::c_beta_deviation_t> >
coot::get_c_beta_deviations(mmdb::Manager *mol) {

   std::map<mmdb::Residue *, std::map<std::string, c_beta_deviation_t> > m; // return this
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
	 // std::cout << "ichain " << ichain << std::endl;
	 mmdb::Chain *chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 for (int ires=0; ires<nres; ires++) {
	    mmdb::Residue *residue_p = chain_p->GetResidue(ires);
	    int n_atoms = residue_p->GetNumberOfAtoms();
	    std::map<std::string, c_beta_deviation_t> cbdm = get_c_beta_deviations(residue_p);
	    if (cbdm.size()) {
	       m[residue_p] = cbdm;
	       if (false) { // debug
		  std::map<std::string, c_beta_deviation_t>::const_iterator it;
		  it = cbdm.begin();
		  std::cout << " " << residue_spec_t(residue_p) << " " << it->second.dist << std::endl;
	       }
	    }
	 }
      }
   }
   return m;
}

std::map<std::string, coot::c_beta_deviation_t>
coot::get_c_beta_deviations(mmdb::Residue *residue_p) {

   std::map<std::string, c_beta_deviation_t> m;
   std::string res_name(residue_p->GetResName());

   int n_atoms = residue_p->GetNumberOfAtoms();
   mmdb::Atom *at_0 = nullptr;
   mmdb::Atom *at_1 = nullptr;
   mmdb::Atom *at_2 = nullptr;
   mmdb::Atom *at_3 = nullptr;
   std::map<std::string, atom_quad> alt_conf_map;
   for (int iat=0; iat<n_atoms; iat++) {
      mmdb::Atom *at = residue_p->GetAtom(iat);
      std::string atom_name(at->GetAtomName());
      std::string alt_conf(at->altLoc);
      if (atom_name == " N  ") alt_conf_map[alt_conf].atom_1 = at;
      if (atom_name == " CA ") alt_conf_map[alt_conf].atom_2 = at;
      if (atom_name == " C  ") alt_conf_map[alt_conf].atom_3 = at;
      if (atom_name == " CB ") alt_conf_map[alt_conf].atom_4 = at;
   }

   // this doesn't work when CB has alt-conf but N (or CA) does not.
   // needs reworking... construct a quad for all CBs in residue, the
   // atoms of which fall back to '' if needed.

   std::map<std::string, atom_quad>::const_iterator it;
   for(it=alt_conf_map.begin(); it!=alt_conf_map.end(); it++) {
      const atom_quad &q = it->second;
      if (q.filled_p()) {
	 clipper::Coord_orth CB_real_pos(q.atom_4->x, q.atom_4->y, q.atom_4->z); // no coot-utils.h included
	 clipper::Coord_orth CB_ideal_pos = make_CB_ideal_pos(q, res_name);
	 double dsqrd = (CB_ideal_pos-CB_real_pos).lengthsq();
	 double d = sqrt(dsqrd);
	 c_beta_deviation_t cbd(q.atom_4, CB_ideal_pos, d);
	 const std::string &key = it->first;
	 m[key] = cbd;
      }
   }

   return m;
}

clipper::Coord_orth
coot::make_CB_ideal_pos(const coot::atom_quad &q, const std::string &res_name) {

// > median(a$V7[a$V5=='PRO'])
// [1] 115.23
// > median(a$V9[a$V5=='PRO'])
// [1] -119.849
// > median(a$V7[a$V5!='PRO'])
// [1] 122.902
// > median(a$V9[a$V5!='PRO'])
// [1] -122.602

   bool is_PRO = false;
   if (res_name == std::string("PRO"))
      is_PRO = true;

   clipper::Coord_orth pt_1(q.atom_1->x, q.atom_1->y, q.atom_1->z);
   clipper::Coord_orth pt_2(q.atom_2->x, q.atom_2->y, q.atom_2->z);
   clipper::Coord_orth pt_3(q.atom_3->x, q.atom_3->y, q.atom_3->z);

   double l = 1.53;
   double a1 = clipper::Util::d2rad(111.0);
   double a2 = clipper::Util::d2rad(111.7);
   double t1 = clipper::Util::d2rad(122.9);
   double t2 = clipper::Util::d2rad(-122.6);

   if (res_name == "ALA") l = 1.509;
   if (res_name == "ASP") l = 1.531;
   if (res_name == "ASN") l = 1.531;
   if (res_name == "CYS") l = 1.524;
   if (res_name == "GLU") l = 1.530;
   if (res_name == "PHE") l = 1.531;
   if (res_name == "HIS") l = 1.535;
   if (res_name == "ILE") l = 1.542;
   if (res_name == "LYS") l = 1.532;
   if (res_name == "LEU") l = 1.532;
   if (res_name == "MET") l = 1.532;
   if (res_name == "PRO") l = 1.534;
   if (res_name == "GLN") l = 1.530;
   if (res_name == "ARG") l = 1.532;
   if (res_name == "SER") l = 1.507;
   if (res_name == "THR") l = 1.534;
   if (res_name == "VAL") l = 1.541;
   if (res_name == "TRP") l = 1.534;
   if (res_name == "TYR") l = 1.531;

   if (res_name == "ALA") a1 = clipper::Util::d2rad(109.912);
   if (res_name == "ASP") a1 = clipper::Util::d2rad(111.338);
   if (res_name == "ASN") a1 = clipper::Util::d2rad(111.766);
   if (res_name == "CYS") a1 = clipper::Util::d2rad(110.827);
   if (res_name == "GLU") a1 = clipper::Util::d2rad(110.374);
   if (res_name == "PHE") a1 = clipper::Util::d2rad(110.494);
   if (res_name == "HIS") a1 = clipper::Util::d2rad(110.437);
   if (res_name == "ILE") a1 = clipper::Util::d2rad(110.820);
   if (res_name == "LYS") a1 = clipper::Util::d2rad(110.374);
   if (res_name == "LEU") a1 = clipper::Util::d2rad(108.955);
   if (res_name == "MET") a1 = clipper::Util::d2rad(110.906);
   if (res_name == "PRO") a1 = clipper::Util::d2rad(103.430); // interesting
   if (res_name == "GLN") a1 = clipper::Util::d2rad(110.374);
   if (res_name == "ARG") a1 = clipper::Util::d2rad(110.374);
   if (res_name == "SER") a1 = clipper::Util::d2rad(110.990);
   if (res_name == "THR") a1 = clipper::Util::d2rad(111.125);
   if (res_name == "VAL") a1 = clipper::Util::d2rad(111.441);
   if (res_name == "TRP") a1 = clipper::Util::d2rad(110.562);
   if (res_name == "TYR") a1 = clipper::Util::d2rad(110.494);

   if (res_name == "ALA") a2 = clipper::Util::d2rad(111.490);
   if (res_name == "ASP") a2 = clipper::Util::d2rad(111.804);
   if (res_name == "ASN") a2 = clipper::Util::d2rad(111.540);
   if (res_name == "CYS") a2 = clipper::Util::d2rad(109.612);
   if (res_name == "GLU") a2 = clipper::Util::d2rad(111.037);
   if (res_name == "PHE") a2 = clipper::Util::d2rad(111.331);
   if (res_name == "HIS") a2 = clipper::Util::d2rad(112.128);
   if (res_name == "ILE") a2 = clipper::Util::d2rad(111.764);
   if (res_name == "LYS") a2 = clipper::Util::d2rad(111.037);
   if (res_name == "LEU") a2 = clipper::Util::d2rad(111.075);
   if (res_name == "MET") a2 = clipper::Util::d2rad(109.344);
   if (res_name == "PRO") a2 = clipper::Util::d2rad(110.031);
   if (res_name == "GLN") a2 = clipper::Util::d2rad(111.037);
   if (res_name == "ARG") a2 = clipper::Util::d2rad(111.037);
   if (res_name == "SER") a2 = clipper::Util::d2rad(111.379);
   if (res_name == "THR") a2 = clipper::Util::d2rad(111.511);
   if (res_name == "VAL") a2 = clipper::Util::d2rad(111.388);
   if (res_name == "TRP") a2 = clipper::Util::d2rad(111.644);
   if (res_name == "TYR") a2 = clipper::Util::d2rad(111.331);
   if (is_PRO) {
     t1 = clipper::Util::d2rad( 115.20);
     t2 = clipper::Util::d2rad(-119.80);
   }
   clipper::Coord_orth pt_trial1 = clipper::Coord_orth(pt_1, pt_3, pt_2, l, a2, t1);
   clipper::Coord_orth pt_trial2 = clipper::Coord_orth(pt_3, pt_1, pt_2, l, a1, t2);

   double delta_trial_points = std::sqrt((pt_trial2-pt_trial1).lengthsq());

   if (false)
      std::cout << "delta-trial-pts " << delta_trial_points << " " << res_name << std::endl;

   clipper::Coord_orth pt_trial = 0.5 * (pt_trial1 + pt_trial2);

   return pt_trial;
}
