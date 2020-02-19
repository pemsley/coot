
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

void
coot::daca::fill_reference_fragments() {

   std::string pkg_data_dir = coot::package_data_dir();
   std::string fn = util::append_dir_file(pkg_data_dir, "standard-residues.pdb");
   if (file_exists(fn)) {
      atom_selection_container_t asc = get_atom_selection(fn, false, false);
      if (asc.read_success) {
         // make reference fragments for each of the residues
         std::cout << "Now do things with residues in standard_residues\n";
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
                           if (true)
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
         std::cout << "INFO:: atom count mismatch in getting fragments for "
                   << residue_spec_t(reference_residue_p) << " "
                   << reference_residue_p->GetResName() << " "
                   << atom_vec.size() << std::endl;
      }
   }
   return v;
}

clipper::RTop_orth
coot::daca::get_frag_to_reference_rtop(const std::string &res_name,
                                       const unsigned int &frag_idx,
                                       const std::vector<mmdb::Atom *> &fragment_atoms) const {
   clipper::RTop_orth rtop;
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
         }
      } else {
         std::cout << "index vector error in get_frag_to_reference_rtop() " << frag_idx << " res_name"
                   << std::endl;
      }

   } else {
      std::cout << "index residue error in get_frag_to_reference_rtop() " << frag_idx << " res_name"
                << std::endl;
   }
   return rtop;
}

void
coot::daca::add_to_box(const std::string &residue_type,
                       bool is_helical_flag,
                       unsigned int frag_index,
                       const box_index_t &box_index,
                       const std::string &atom_type) {

   std::string box_key = residue_type + "-non-helical";
   if (is_helical_flag) box_key = residue_type + "-helical";

   if (frag_index < boxes[residue_type].size()) {
      boxes[box_key][frag_index][atom_type][box_index]++;
   } else {
      boxes[box_key].resize(frag_index+1);
      boxes[box_key][frag_index][atom_type][box_index]++;
   }
}

void
coot::daca::debug_boxes() const {

   std::map<std::string, std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > >::const_iterator it;
   for (it=boxes.begin(); it!=boxes.end(); it++) {
      const std::string residue_type = it->first;
      std::cout << "============= Residue Type " << residue_type << " helical ============ " << std::endl;
      const std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > &frag_boxes = it->second;
      std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > >::const_iterator it_v;
      for (it_v=frag_boxes.begin(); it_v!=frag_boxes.end(); it_v++) {
         const std::map<std::string, std::map<box_index_t, unsigned int> > &typed_boxes = *it_v;
         std::map<std::string, std::map<box_index_t, unsigned int> >::const_iterator it_typed_box;
         for (it_typed_box=typed_boxes.begin(); it_typed_box!=typed_boxes.end(); it_typed_box++) {
            std::string atom_type = it_typed_box->first;
            std::cout << "----------------- Residue Type " << residue_type << " atom_type " << atom_type << std::endl;
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


void
coot::daca::write_tables() const {

   std::cout << "write_tables(): write " << boxes.size() << " boxes " << std::endl;

   std::map<std::string, std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > >::const_iterator it;
   for (it=boxes.begin(); it!=boxes.end(); it++) {
      const std::string residue_type = it->first;
      std::cout << "============= write_tables(): Residue Type " << residue_type << " ============ " << std::endl;
      const std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > &frag_boxes = it->second;
      for (unsigned int i=0; i<frag_boxes.size(); i++) {
         const std::map<std::string, std::map<box_index_t, unsigned int> > &typed_boxes = frag_boxes[i];
         std::map<std::string, std::map<box_index_t, unsigned int> >::const_iterator it_typed_box;
         for (it_typed_box=typed_boxes.begin(); it_typed_box!=typed_boxes.end(); it_typed_box++) {
            std::string atom_type = it_typed_box->first;
            std::cout << "----------------- write_tables(): Residue Type " << residue_type << " atom type " << atom_type << std::endl;
            std::string box_file_name = residue_type + "-" + util::int_to_string(i) + "-" + atom_type + ".table";
            std::cout << "box_file_name: " << box_file_name << std::endl;
            std::ofstream f(box_file_name.c_str());
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
coot::daca::calculate_daca(mmdb::Residue *reference_residue_p,
                           const std::vector<std::pair<mmdb::Atom *, std::string> > &typed_atoms) {

   // brain-dead distance search (sad face)

   // has fill_helix_flags() been called before now?

   double d_crit = 8.0; // or something
   double dd_crit = d_crit * d_crit;

   std::string res_name(reference_residue_p->GetResName());
   std::vector<std::vector<mmdb::Atom *> > fragments = get_daca_fragments(reference_residue_p);
   if (false)
      std::cout << "debug:: fragments.size() " << fragments.size() << " "
                << residue_spec_t(reference_residue_p)
                << " " << reference_residue_p->GetResName() << std::endl;
   for (unsigned int i=0; i<fragments.size(); i++) {
      const std::vector<mmdb::Atom *> &atom_vec(fragments[i]);
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
            clipper::RTop_orth frag_to_reference_rtop = get_frag_to_reference_rtop(res_name, i, atom_vec);
            for (unsigned int ita=0; ita<typed_atoms.size(); ita++) {
               mmdb::Atom *at = typed_atoms[ita].first;
               // don't consider atoms in this residue, of course
               if (at->residue == reference_residue_p)
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
                        add_to_box(res_name, helical_flag, i, box_index, typed_atoms[ita].second);
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
coot::daca::write_tables_using_reference_structures_from_dir(const std::string &dir_name) {

   protein_geometry geom;
   geom.init_standard();
   std::vector<std::string> files = util::glob_files(dir_name, "*.pdb");

   for (unsigned int i=0; i<files.size(); i++) {
      std::string fn = files[i];
      std::cout << fn << std::endl;
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
                     calculate_daca(residue_p, ta);
                  }
               }
            }
         }
      }
   }

   debug_boxes();
   write_tables();

}
