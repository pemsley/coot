
#include <map>
#include <iomanip>
#include <fstream>
#include "utils/coot-utils.hh"
#include "geometry/residue-and-atom-specs.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/atom-selection-container.hh"
#include "coot-utils/helix-like.hh"
#include "daca.hh"

coot::daca::box_index_t::box_index_t(const clipper::Coord_orth &pos) {
   double box_width = 1.0;
   idx_x = floor(pos.x()/box_width);
   idx_y = floor(pos.y()/box_width);
   idx_z = floor(pos.z()/box_width);

}

// the reverse of the above - make a point in the middle of the box
clipper::Coord_orth
coot::daca::box_index_t::coord_orth() const {
   double box_width = 1.0;
   double x = static_cast<double>(idx_x) * box_width + 0.5 * box_width;
   double y = static_cast<double>(idx_y) * box_width + 0.5 * box_width;
   double z = static_cast<double>(idx_z) * box_width + 0.5 * box_width;
   return clipper::Coord_orth(x,y,z);
}

bool
coot::daca::box_index_t::operator<(const coot::daca::box_index_t &other) const {
   if (other.idx_x < idx_x) return true;
   if (other.idx_x > idx_x) return false;
   if (other.idx_y < idx_y) return true;
   if (other.idx_y > idx_y) return false;
   if (other.idx_z < idx_z) return true;
   if (other.idx_z > idx_z) return false;
   return false;
}


// I can't get where this should go.
#if 0
std::ostream &
coot::daca::operator<<(std::ostream &s, const coot::daca::box_index_t &bi) {
   s << "[box " << bi.x << " " << bi.y << " " << bi.z << "]";
   return s;
}
#endif

void
coot::daca::fill_reference_fragments() {

   std::string pkg_data_dir = coot::package_data_dir();
   std::string fn = util::append_dir_file(pkg_data_dir, "standard-residues.pdb");
   if (file_exists(fn)) {
      atom_selection_container_t asc = get_atom_selection(fn, false, false);
      if (asc.read_success) {
         // make reference fragments for each of the residues
         // std::cout << "Now do things with residues in standard_residues\n";
         int imod = 1;
         mmdb::Model *model_p = asc.mol->GetModel(imod);
         if (model_p) {
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               int nres = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<nres; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     std::string res_name(residue_p->GetResName());
                     std::vector<std::vector<std::string> > atom_name_sets =
                        atom_names_for_fragments(res_name);

                     for(unsigned int iset=0; iset<atom_name_sets.size(); iset++) {
                        std::vector<clipper::Coord_orth> v; // this has special order
                        const std::vector<std::string> &us = atom_name_sets[iset];
                        std::vector<std::string>::const_iterator it;
                        for (it=us.begin(); it!=us.end(); it++) {
                           const std::string &frag_atom_name = *it;
                           int n_residue_atoms;
                           mmdb::PPAtom residue_atoms = 0;
                           residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
                           for (int iat=0; iat<n_residue_atoms; iat++) {
                              mmdb::Atom *at = residue_atoms[iat];
                              std::string this_atom_name(at->GetAtomName());
                              if (this_atom_name == frag_atom_name) {
                                 clipper::Coord_orth pos = co(at);
                                 v.push_back(pos);
                                 break;
                              }
                           }
                        }
                        if (v.size() == us.size()) {
                           // add another set of coordinates to reference_fragments[res_name]
                           clipper::Coord_orth sum(0,0,0); // for calculating the centre of the fragment
                           for (unsigned int ii=0; ii<v.size(); ii++)
                              sum += v[ii];
                           double m = 1.0/static_cast<double>(v.size());
                           clipper::Coord_orth fragment_centre(sum * m);
                           for (unsigned int ii=0; ii<v.size(); ii++)
                              v[ii] -= fragment_centre;
                           reference_fragments[res_name].push_back(v);
                           if (false)
                              std::cout << " filling " << residue_p << " " << res_name << " "
                                        << v.size() << std::endl;
                        }
                     }
                  }
               }
            }
         }
      }
   } else {
      std::cout << "File not found " << fn << std::endl;
   }

   if (false) { // debugging
      std::map<std::string, std::vector<std::vector<clipper::Coord_orth> > >::const_iterator it;
      for (it =reference_fragments.begin(); it!=reference_fragments.end(); it++){
         std::cout << "Reference Residue type " << it->first << "\n";
         std::vector<std::vector<clipper::Coord_orth> >::const_iterator itvv;
         for (itvv=it->second.begin(); itvv!=it->second.end(); itvv++) {
            const std::vector<clipper::Coord_orth> &v = *itvv;
            std::vector<clipper::Coord_orth>::const_iterator itv;
            for (itv=v.begin(); itv!=v.end(); itv++)
               std::cout << itv->format() << " " ;
            std::cout << std::endl;
         }
      }
   }

}

std::vector<std::pair<mmdb::Atom *, std::string> >
coot::daca::make_typed_atoms(mmdb::Model *model_p, const coot::protein_geometry &geom) const {

   std::vector<std::pair<mmdb::Atom *, std::string> > v;

   std::map<std::string, dictionary_residue_restraints_t> dictionary_map;

   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int nres = chain_p->GetNumberOfResidues();
         for (int ires=0; ires<nres; ires++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(ires);
            if (residue_p) {
               std::string res_type = residue_p->GetResName();
               std::map<std::string, dictionary_residue_restraints_t>::const_iterator it;
               it = dictionary_map.find(res_type);
               if (it == dictionary_map.end()) {
                  std::pair<bool, dictionary_residue_restraints_t> restraints =
                     geom.get_monomer_restraints(res_type, protein_geometry::IMOL_ENC_ANY);
                  if (restraints.first) {
                     dictionary_map[res_type] = restraints.second;
                  }
               }
            }
         }
      }
   }

   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int nres = chain_p->GetNumberOfResidues();
         for (int ires=0; ires<nres; ires++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(ires);
            if (residue_p) {
               int n_atoms = residue_p->GetNumberOfAtoms();
               if (n_atoms > 0) {
                  std::string res_type = residue_p->GetResName();
                  std::map<std::string, dictionary_residue_restraints_t>::const_iterator it;
                  it = dictionary_map.find(res_type);
                  if (it != dictionary_map.end()) {
                     for (int iat=0; iat<n_atoms; iat++) {
                        mmdb::Atom *at = residue_p->GetAtom(iat);
                        if (at) {
                           std::string atom_name(at->GetAtomName());
                           const std::string type = it->second.type_energy(atom_name);

                           // check  for atom name being "N" here //  C is Correct type
                           if (atom_name == " N  ") {
                              std::pair<mmdb::Atom *, std::string> p(at, "NH1");
                              v.push_back(p);
                           } else {
                              if (! type.empty()) {
                                 std::pair<mmdb::Atom *, std::string> p(at, type);
                                 v.push_back(p);
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   return v;
}

std::vector<std::vector<std::string> >
coot::daca::atom_names_for_fragments(const std::string &res_name) const {
   std::vector<std::vector<std::string> > v;
   if (res_name == "GLY") {
      std::vector<std::string> s{" N  ", " C  ", " CA "};
      v.push_back(s);
   }
   if (res_name == "ALA") {
      std::vector<std::string> s{" N  ", " C  ", " CA ", " CB "};
      v.push_back(s);
   }
   if (res_name == "CYS") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " SG "};
      v.push_back(s1);
      v.push_back(s2);
   }
   if (res_name == "ASP") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG "};
      std::vector<std::string> s3{" CB ", " CG ", " OD1", " OD2"};
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
   }
   if (res_name == "GLU") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG "};
      std::vector<std::string> s3{" CB ", " CG ", " CD "};
      std::vector<std::string> s4{" CG ", " CD ", " OE1", " OE2"};
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
      v.push_back(s4);
   }
   if (res_name == "PHE") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG "};
      std::vector<std::string> s3{" CB ", " CG ", " CD1", " CD2", " CE1", " CE2", " CZ "};
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
   }
   if (res_name == "HIS") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG "};
      std::vector<std::string> s3{" CB ", " CG ", " ND1", " CE1", " NE2", " CD2"};
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
   }
   if (res_name == "ILE") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG1", " CG2"};
      std::vector<std::string> s3{" CB ", " CG1", " CD1"};
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
   }
   if (res_name == "LYS") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG "};
      std::vector<std::string> s3{" CB ", " CG ", " CD "};
      std::vector<std::string> s4{" CG ", " CD ", " CE "};
      std::vector<std::string> s5{" CD ", " CE ", " NZ "};
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
      v.push_back(s4);
      v.push_back(s5);
   }
   if (res_name == "LEU") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG "};
      std::vector<std::string> s3{" CB ", " CG ", " CD1", " CD2"};
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
   }
   if (res_name == "MET") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG "};
      std::vector<std::string> s3{" CB ", " CG ", " SD "};
      std::vector<std::string> s4{" CG ", " SD ", " CE "};
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
      v.push_back(s4);
   }
   if (res_name == "ASN") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG "};
      std::vector<std::string> s3{" CB ", " CG ", " OD1", " ND2"};
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
   }
   if (res_name == "PRO") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB ", " CG ", " CD "};
      v.push_back(s1);
   }
   if (res_name == "GLN") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG "};
      std::vector<std::string> s3{" CB ", " CG ", " CD "};
      std::vector<std::string> s4{" CG ", " CD ", " OE1", " NE2"};
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
      v.push_back(s4);
   }
   if (res_name == "ARG") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG "};
      std::vector<std::string> s3{" CB ", " CG ", " CD "};
      std::vector<std::string> s4{" CG ", " CD ", " NE "};
      std::vector<std::string> s5{" NE ", " CZ ", " NH1", " NH2"}; // + CD?
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
      v.push_back(s4);
      v.push_back(s5);
   }
   if (res_name == "SER") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " OG "};
      v.push_back(s1);
      v.push_back(s2);
   }
   if (res_name == "THR") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " OG1", " CG2"};
      v.push_back(s1);
      v.push_back(s2);
   }
   if (res_name == "VAL") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG1", " CG2"};
      v.push_back(s1);
      v.push_back(s2);
   }
   if (res_name == "TRP") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG "};
      std::vector<std::string> s3{" CB ", " CG ", " CD1", " CD2", " NE1", " CE2", " CE3", " CZ2", " CZ3", " CH2"};
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
   }
   if (res_name == "TYR") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG "};
      std::vector<std::string> s3{" CB ", " CG ", " CD1", " CD2", " CE1", " CE2", " CZ ", " OH "};
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
   }
   return v;
}

std::vector<std::vector<mmdb::Atom *> >
coot::daca::get_daca_fragments(mmdb::Residue *reference_residue_p) const {

   std::vector<std::vector<mmdb::Atom *> > v;
   std::string res_name(reference_residue_p->GetResName());
   std::vector<std::vector<std::string> > atom_name_vec_vec = atom_names_for_fragments(res_name);
   std::vector<std::vector<std::string> >::const_iterator it_1;
   std::vector<std::string>::const_iterator it_2;
   for (it_1=atom_name_vec_vec.begin(); it_1!=atom_name_vec_vec.end(); it_1++){
      const std::vector<std::string> &atom_names = *it_1;
      std::vector<mmdb::Atom *> atom_vec;
      for (it_2=atom_names.begin(); it_2!=atom_names.end(); it_2++) {
         const std::string &atom_name = *it_2;
         mmdb::PPAtom residue_atoms = 0;
         int n_residue_atoms;
         reference_residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            mmdb::Atom *at = residue_atoms[iat];
            std::string this_atom_name(at->GetAtomName());
            if (atom_name == this_atom_name) {
               std::string alt_loc(at->altLoc);
               if (alt_loc.empty())
                  atom_vec.push_back(at);
            }
         }
      }
      if (atom_names.size() == atom_vec.size()) {
         v.push_back(atom_vec);
      } else {
         // debugging
         // we can get here when the residue has alt confs.
         if (false) // too noisy for now.
            std::cout << "INFO:: atom count mismatch in getting fragments for "
                      << residue_spec_t(reference_residue_p) << " "
                      << reference_residue_p->GetResName() << " "
                      << atom_vec.size() << std::endl;
      }
   }
   return v;
}

std::pair<bool, clipper::RTop_orth>
coot::daca::get_frag_to_reference_rtop(const std::string &res_name,
                                       const unsigned int &frag_idx,
                                       const std::vector<mmdb::Atom *> &fragment_atoms) const {
   clipper::RTop_orth rtop;
   bool status = true;
   std::map<std::string, std::vector<std::vector<clipper::Coord_orth> > >::const_iterator it;
   it = reference_fragments.find(res_name);
   if (it != reference_fragments.end()) {
      if (frag_idx < it->second.size()) {
         const std::vector<clipper::Coord_orth> &ref_atom_positions = it->second[frag_idx];
         // convert fragment_atoms to coordinates
         std::vector<clipper::Coord_orth> residue_fragment_atoms;
         for (unsigned int i=0; i<fragment_atoms.size(); i++) {
            clipper::Coord_orth pos = co(fragment_atoms[i]);
            residue_fragment_atoms.push_back(pos);
         }
         if (ref_atom_positions.size() == residue_fragment_atoms.size()) {
            clipper::RTop_orth rtop_1(residue_fragment_atoms, ref_atom_positions);
            rtop = rtop_1;
         } else {
            std::cout << "size error in get_frag_to_reference_rtop()" << std::endl;
            status = false;
         }
      } else {
         std::cout << "index vector error in get_frag_to_reference_rtop() " << frag_idx << " res_name"
                   << std::endl;
         status = false;
      }

   } else {
      std::cout << "index residue error in get_frag_to_reference_rtop() " << frag_idx << " res_name"
                << std::endl;
      status = false;
   }
   return std::pair<bool, clipper::RTop_orth> (status, rtop);
}

void
coot::daca::add_to_box(mode_t mode,
                       const std::string &residue_type,
                       bool is_helical_flag,
                       unsigned int frag_index,
                       const box_index_t &box_index,
                       const std::string &atom_type,
                       unsigned int counts) {

   std::string box_key = residue_type + "-non-helical";
   if (is_helical_flag) box_key = residue_type + "-helical";

   if (mode == REFERENCE) {
      boxes[box_key][frag_index][atom_type][box_index] += counts;
   }

   if (mode == ANALYSIS) {
      if (frag_index >= boxes_for_testing[residue_type].size())
         boxes_for_testing[box_key].resize(5);
      boxes_for_testing[box_key][frag_index][atom_type][box_index] += counts;
   }
}

void
coot::daca::debug_boxes(const std::string &prefix) const {

   std::map<std::string, std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > >::const_iterator it;
   for (it=boxes.begin(); it!=boxes.end(); it++) {
      const std::string residue_type = it->first;
      std::cout << "========== debug_boxes(): " << prefix << " Residue Type " << residue_type << std::endl;
      const std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > &frag_boxes = it->second;
      std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > >::const_iterator it_v;
      for (unsigned int ifrag=0; ifrag<frag_boxes.size(); ifrag++) {
         const std::map<std::string, std::map<box_index_t, unsigned int> > &typed_boxes = frag_boxes[ifrag];
         std::map<std::string, std::map<box_index_t, unsigned int> >::const_iterator it_typed_box;
         for (it_typed_box=typed_boxes.begin(); it_typed_box!=typed_boxes.end(); it_typed_box++) {
            std::string atom_type = it_typed_box->first;
            if (residue_type.substr(0,3) == "ARG") {
               if (ifrag == 0) {
                  std::cout << "========== debug_boxes(): " << prefix << " Residue Type " << residue_type << " frag index "
                            << ifrag << " atom_type " << atom_type << std::endl;
                  if (true) {
                     std::map<box_index_t, unsigned int>::const_iterator it_box;
                     for (it_box=it_typed_box->second.begin(); it_box!=it_typed_box->second.end(); it_box++) {
                        const box_index_t &bi = it_box->first;
                        unsigned int count = it_box->second;
                        std::cout << " "
                                  << std::setw(2) << bi.idx_x << " " << std::setw(2) << bi.idx_y << " " << std::setw(2) << bi.idx_z << " "
                                  << std::setw(3) << count << std::endl;
                     }
                  }
               }
            }
         }
      }
   }
}


void
coot::daca::write_tables(const std::string &dir) const {

   std::cout << "write_tables(): write " << boxes.size() << " boxes " << std::endl;
   coot::util::create_directory(dir);

   std::map<std::string, std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > >::const_iterator it;
   for (it=boxes.begin(); it!=boxes.end(); it++) {
      const std::string &residue_type = it->first;
      std::cout << "============= write_tables(): Residue Type " << residue_type << std::endl;
      const std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > &frag_boxes = it->second;
      for (unsigned int i=0; i<frag_boxes.size(); i++) {
         const std::map<std::string, std::map<box_index_t, unsigned int> > &typed_boxes = frag_boxes[i];
         std::map<std::string, std::map<box_index_t, unsigned int> >::const_iterator it_typed_box;
         for (it_typed_box=typed_boxes.begin(); it_typed_box!=typed_boxes.end(); it_typed_box++) {
            std::string atom_type = it_typed_box->first;
            if (false)
               std::cout << "----------------- write_tables(): Residue Type " << residue_type << " " << i << " atom type "
                         << atom_type << std::endl;
            std::string box_file_name = residue_type + "-" + util::int_to_string(i) + "-" + atom_type + ".table";
            std::string full_box_file_name = coot::util::append_dir_file(dir, box_file_name);
            std::ofstream f(full_box_file_name.c_str());
            if (f) {
               std::map<box_index_t, unsigned int>::const_iterator it_box;
               for (it_box=it_typed_box->second.begin(); it_box!=it_typed_box->second.end(); it_box++) {
                  const box_index_t &bi = it_box->first;
                  unsigned int count = it_box->second;
                  f << " "
                    << std::setw(2) << bi.idx_x << " " << std::setw(2) << bi.idx_y << " " << std::setw(2) << bi.idx_z << " "
                    << std::setw(3) << count << "\n";
               }
               f.close();
            }
         }
      }
   }
}

void
coot::daca::read_many_tables(const std::vector<std::string> &dirs) {
   presize_boxes();
   for (unsigned int i=0; i<dirs.size(); i++) {
      std::cout << "read tables directory " << dirs[i] << std::endl;
      read_tables(dirs[i]);
   }
}

void
coot::daca::read_tables(const std::string &dir) {

   if (! boxes_have_been_resized)
      presize_boxes();

   std::string glob_pattern = "*.table";
   std::vector<std::string> files = coot::util::glob_files(dir, glob_pattern);
   for (unsigned int i=0; i<files.size(); i++) {
      std::string file_name = files[i];
      // std::cout << "read table file " << file_name << std::endl;

      std::pair<std::string, std::string> z_parts = coot::util::split_string_on_last_slash(file_name);
      std::vector<std::string> fn_parts = coot::util::split_string(z_parts.second, "-");

      if (false) {
         std::cout << "fn_parts: " << std::endl;
         for (unsigned int i=0; i<fn_parts.size(); i++)
            std::cout << fn_parts[i] << " ";
         std::cout << std::endl;
      }

      if (fn_parts.size() == 4 || fn_parts.size() == 5) {
         try {
            std::string res_name = fn_parts[0];
            std::string ss_type = "helical";
            int ss_type_index = 0;
            unsigned int frag_string_index = 2;
            unsigned int atom_type_index = 3;
            bool is_helical_flag = true;
            if (fn_parts[1] == "non") {
               ss_type = "non-helical";
               ss_type_index = 1;
               frag_string_index = 3;
               atom_type_index = 4;
               is_helical_flag = false;
            }
            std::string frag_string = fn_parts[frag_string_index];
            int frag_index = coot::util::string_to_int(frag_string);
            const std::string &at_raw = fn_parts[atom_type_index];
            unsigned int l = at_raw.size();
            std::string atom_type = at_raw.substr(0,l-6);
            if (false)
               std::cout << " decoded: " << res_name << " " << ss_type << " " << frag_index
                         << " " << atom_type << std::endl;

            std::string line;
            std::vector<std::string> lines;
            std::ifstream f(files[i].c_str());
            while (std::getline(f, line)) {
               lines.push_back(line);
            }
            for (unsigned int j=0; j<lines.size(); j++) {
               const std::string &line = lines[j];
               std::vector<std::string> parts = coot::util::split_string_on_whitespace_no_blanks(line);
               if (parts.size() == 4) {
                  // .. x y z count
                  try {
                     int x = coot::util::string_to_int(parts[0]);
                     int y = coot::util::string_to_int(parts[1]);
                     int z = coot::util::string_to_int(parts[2]);
                     int c = coot::util::string_to_int(parts[3]);
                     box_index_t bi(x,y,z);
                     add_to_box(REFERENCE, res_name, is_helical_flag, frag_index, bi, atom_type, c);
                  }
                  catch (const std::runtime_error &rte) {
                     std::cout << "failed to parse " << line << " from " << files[i] << " " << rte.what() << std::endl;
                  }
               }
            }
         }
         catch (const std::runtime_error &rte) { }
      }
   }
}


void
coot::daca::fill_helix_flags(mmdb::Model *model_p, mmdb::Manager *mol) {

   std::vector<std::string> ch_ids;
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         ch_ids.push_back(chain_p->GetChainID());
      }
   }

   for (unsigned int ich=0; ich<ch_ids.size(); ich++) {
      int residue_selection_handle = mol->NewSelection();
      mol->Select (residue_selection_handle, mmdb::STYPE_RESIDUE, 0,
         ch_ids[ich].c_str(),
         mmdb::ANY_RES, "*",  // starting res
         mmdb::ANY_RES, "*",  // ending res
         "*",  // residue name
         "*",  // Residue must contain this atom name?
         "*",  // Residue must contain this Element?
         "*",  // altLocs
         mmdb::SKEY_NEW // selection key
      );

      std::vector<mmdb::Residue *> helical_residues_in_chain = like_a_helix(mol, residue_selection_handle);
      for (unsigned int i=0; i<helical_residues_in_chain.size(); i++)
      helical_residues.push_back(helical_residues_in_chain[i]);

      mol->DeleteSelection(residue_selection_handle);
   }
}

bool
coot::daca::atom_is_close_to_a_residue_atom(mmdb::Atom *at, mmdb::Residue *reference_residue_p) const {
   float d_close = 1.7 + 1.7 + 1.5; // or so
   float dd_close = d_close * d_close;
   bool status = false;
   int n_residue_atoms;
   mmdb::PPAtom residue_atoms = 0;
   reference_residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      mmdb::Atom *ref_at = residue_atoms[iat];
      float dd =
         (at->x - ref_at->x) * (at->x - ref_at->x) +
         (at->y - ref_at->y) * (at->y - ref_at->y) +
         (at->z - ref_at->z) * (at->z - ref_at->x);
      if (dd < dd_close) {
         status = true;
         break;
      }
   }
   return status;
}


bool
coot::daca::atom_is_neighbour_mainchain(mmdb::Atom *at, mmdb::Residue *reference_residue_p) const {
   bool status = false;
   int idx_res_1 = reference_residue_p->index;
   int idx_res_2 = at->residue->index;
   int idx_delta = abs(idx_res_2 - idx_res_1);
   if (idx_delta < 2) {
      std::string atom_name(at->GetAtomName());
      if (atom_name == " N  ") { status = true; }
      if (atom_name == " CA ") { status = true; }
      if (atom_name == " C  ") { status = true; }
      if (atom_name == " O  ") { status = true; }
   }
   return status;
}

void
coot::daca::presize_boxes(mode_t mode) {

   std::vector<std::string> types = { "GLY", "ALA", "CYS", "ASP", "GLU", "PHE", "HIS", "ILE", "LYS", "LEU",
                                      "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"};

   if (mode == REFERENCE) {
      boxes_have_been_resized = true;
      for (auto type : types) {
         const std::vector<std::string> h_types = {"-helical", "-non-helical"};
         for (auto h : h_types) {
            std::string key = type + h;
            boxes[key].resize(5);
         }
      }
   }
   if (mode == ANALYSIS) {
      for (auto type : types) {
         const std::vector<std::string> h_types = {"-helical", "-non-helical"};
         for (auto h : h_types) {
            std::string key = type + h;
            boxes_for_testing[key].resize(5);
         }
      }
   }
}

#include "geometry/main-chain.hh"

void
coot::daca::calculate_daca(mmdb::Residue *reference_residue_p,
                           const std::vector<std::pair<mmdb::Atom *, std::string> > &typed_atoms,
                           coot::daca::mode_t mode) {

   // brain-dead distance search (sad face)

   // has fill_helix_flags() been called before now?

   double d_crit = 8.0; // or something
   double dd_crit = d_crit * d_crit;

   presize_boxes();

   std::string res_name(reference_residue_p->GetResName());
   int reference_residue_seqnum = reference_residue_p->GetSeqNum();
   std::vector<std::vector<mmdb::Atom *> > fragments = get_daca_fragments(reference_residue_p);
   if (false)
      std::cout << "debug:: fragments.size() " << fragments.size() << " "
                << residue_spec_t(reference_residue_p)
                << " " << reference_residue_p->GetResName() << std::endl;
   for (unsigned int ifrag=0; ifrag<fragments.size(); ifrag++) {
      const std::vector<mmdb::Atom *> &atom_vec(fragments[ifrag]);
      std::vector<clipper::Coord_orth> reference_positions_vec;
      clipper::Coord_orth sum(0,0,0); // for calculating the centre of the fragment
      std::vector<mmdb::Atom *>::const_iterator it;
      for (it=atom_vec.begin(); it!=atom_vec.end(); it++) {
         clipper::Coord_orth pos = co(*it);
         reference_positions_vec.push_back(pos);
         sum += pos;
      }
      if (reference_positions_vec.size() > 2) {
         if (reference_positions_vec.size() == atom_vec.size()) {

            double m = 1.0/static_cast<double>(reference_positions_vec.size());
            clipper::Coord_orth frag_centre(sum * m);
            // Get the RTop that transforms the fragment to a reference
            // fragment at the origin.
            std::pair<bool, clipper::RTop_orth> frag_to_reference_rtop_pair = get_frag_to_reference_rtop(res_name, ifrag, atom_vec);
            if (! frag_to_reference_rtop_pair.first) continue;
            const clipper::RTop_orth &frag_to_reference_rtop = frag_to_reference_rtop_pair.second;
            for (unsigned int ita=0; ita<typed_atoms.size(); ita++) {
               mmdb::Atom *at = typed_atoms[ita].first;

               // don't consider atoms in this residue, of course
               if (at->residue == reference_residue_p)
                  continue;

               // don't consider peptide neighbour mainchain
               int res_no_delta = at->residue->GetSeqNum() - reference_residue_seqnum;
               if (fabs(res_no_delta) < 2)
                  if (at->residue->chain == reference_residue_p->chain)
                     if (is_main_chain_p(at))
                        continue;

               double dd =
                  (at->x - frag_centre.x()) * (at->x - frag_centre.x()) +
                  (at->y - frag_centre.y()) * (at->y - frag_centre.y()) +
                  (at->z - frag_centre.z()) * (at->z - frag_centre.z());
               if (dd < dd_crit) {
                  // Good, found something
                  if (atom_is_close_to_a_residue_atom(at, reference_residue_p)) {
                     if (! atom_is_neighbour_mainchain(at, reference_residue_p)) {
                        clipper::Coord_orth at_pos = co(at);
                        clipper::Coord_orth transformed_pos = frag_to_reference_rtop * at_pos;
                        box_index_t box_index(transformed_pos);
                        bool helical_flag = false;
                        if (std::find(helical_residues.begin(), helical_residues.end(), reference_residue_p) != helical_residues.end())
                        helical_flag = true;
                        add_to_box(mode, res_name, helical_flag, ifrag, box_index, typed_atoms[ita].second);
                     }
                  }
               }
            }
         } else {
            std::cout << "OOps in atom set vs reference set size test " << std::endl;
         }
      } else {
         std::cout << "ERROR:: in calculate_daca(): This can't happen. reference positions size "
                   << reference_positions_vec.size() << " " << residue_spec_t(reference_residue_p)
                   << std::endl;
      }
   }
}

void
coot::daca::write_tables_using_reference_structures_from_dir(const std::string &dir_name,
                                                             const std::string &output_tables_dir) {

   protein_geometry geom;
   geom.init_standard();
   std::vector<std::string> files = util::glob_files(dir_name, "*.pdb");

   for (unsigned int i=0; i<files.size(); i++) {
      std::string fn = files[i];
      atom_selection_container_t asc = get_atom_selection(fn, false, false);
      if (asc.read_success) {
         mmdb::Model *model_p = asc.mol->GetModel(1);
         if (model_p) {
            fill_helix_flags(model_p, asc.mol);
            std::vector<std::pair<mmdb::Atom *, std::string> > ta = make_typed_atoms(model_p, geom);
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               int nres = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<nres; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     std::string res_name(residue_p->GetResName());
                     if (res_name == "HOH") continue;
                     if (! util::is_standard_amino_acid_name(res_name)) continue;
                     calculate_daca(residue_p, ta, REFERENCE);
                  }
               }
            }
         }
      }
   }

   // debug_boxes("done-write-tables-using-reference-structures");

   write_tables(output_tables_dir);

}

void
coot::daca::score_molecule(const std::string &pdb_file_name) {

   std::cout << "score_molecule() " << pdb_file_name << std::endl;
   if (coot::file_exists(pdb_file_name)) {
      atom_selection_container_t asc = get_atom_selection(pdb_file_name, false, false);
      if (asc.read_success) {
         mmdb::Model *model_p = asc.mol->GetModel(1);
         if (model_p) {

            std::cout << "score " << pdb_file_name << std::endl;

            protein_geometry geom;
            geom.init_standard();
            presize_boxes(ANALYSIS);

            fill_helix_flags(model_p, asc.mol);
            std::vector<std::pair<mmdb::Atom *, std::string> > ta = make_typed_atoms(model_p, geom);
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               int nres = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<nres; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     std::string res_name(residue_p->GetResName());
                     if (res_name == "HOH") continue;
                     if (! util::is_standard_amino_acid_name(res_name)) continue;
                     calculate_daca(residue_p, ta, ANALYSIS);
                  }
               }
            }

            compare_boxes();
         }
      }
   } else {
      std::cout << "No such file " << pdb_file_name << std::endl;
   }

}

void
coot::daca::compare_boxes() const {

   unsigned int n_daca = 0;
   unsigned int n_hits = 0;
   unsigned int sum = 0;

   std::map<std::string, std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > >::const_iterator it;
   for (it =boxes_for_testing.begin(); it!=boxes_for_testing.end(); it++) {
      const std::string &res_name_with_ss(it->first);
      const std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > &v(it->second);
      for (unsigned int idx_frag=0; idx_frag<v.size(); idx_frag++) {
         const std::map<std::string, std::map<box_index_t, unsigned int> > &m1(v[idx_frag]);
         std::map<std::string, std::map<box_index_t, unsigned int> >::const_iterator it_1;
         for (it_1=m1.begin(); it_1!=m1.end(); it_1++) {
            const std::string &atom_type = it_1->first;
            const std::map<box_index_t, unsigned int> &m2(it_1->second);
            std::map<box_index_t, unsigned int>::const_iterator it_2;
            for (it_2=m2.begin(); it_2!=m2.end(); it_2++) {
               const box_index_t &bi = it_2->first;
               int count_analysis = it_2->second;
               n_daca++;

               // does there exist a reference count for that?
               std::map<std::string, std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > >::const_iterator it_ref;
               it_ref = boxes.find(res_name_with_ss);
               if (it_ref == boxes.end()) {
                  std::cout << "Failed to find reference for type " << res_name_with_ss << std::endl;
               } else {
                  const std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > &v_ref(it_ref->second);
                  if (! v_ref.empty()) {
                     const std::map<std::string, std::map<box_index_t, unsigned int> > &m1_ref(v_ref[idx_frag]); // vector sized correctly?
                     std::map<std::string, std::map<box_index_t, unsigned int> >::const_iterator it_1_ref;
                     it_1_ref = m1_ref.find(atom_type);
                     if (it_1_ref == m1_ref.end()) {
                        std::cout << "Failed to find reference for type " << res_name_with_ss
                                  << " frag-index " << idx_frag << " atom-type " << atom_type
                                  << " we have map size " << m1_ref.size() << std::endl;
                     } else {
                        const std::map<box_index_t, unsigned int> &m2_ref(it_1_ref->second);
                        std::map<box_index_t, unsigned int>::const_iterator it_2_ref;
                        it_2_ref = m2_ref.find(bi);
                        if (it_2_ref == m2_ref.end()) {
                           std::cout << "Failed to find reference for " << res_name_with_ss << " "
                                     << idx_frag << " " << atom_type << " box_index "
                                     << bi.idx_x << " " << bi.idx_y << " " << bi.idx_z << std::endl;

                        } else {
                           int count_ref = it_2_ref->second;
                           if (false)
                              std::cout << res_name_with_ss << " " << idx_frag << " " << atom_type << " "
                                        << bi.idx_x << " " << bi.idx_y << " " << bi.idx_z << " "
                                        << count_ref << "\n";
                           sum += count_ref;
                           // std::cout << "sum " << sum << "\n";
                           n_hits++;
                        }
                     }
                  } else {
                     std::cout << "v_ref is empty for " << it_ref->first << std::endl;
                  }
               }
            }
         }
      }
   }
   std::cout << "compare_boxes() n_daca " << n_daca << " n_hits " << n_hits
             << " sum " << sum << std::endl;
}

