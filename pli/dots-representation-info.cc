
#include "coot-utils/coot-coord-utils.hh"
#include "dots-representation-info.hh"
#include "api/coot-colour.hh"

// make a surface on mol
pli::dots_representation_info_t::dots_representation_info_t(mmdb::Manager *mol) {

   is_closed = 0;
   int SelHnd = mol->NewSelection();
   mol->SelectAtoms(SelHnd, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", "*", "*");
   mmdb::Manager *dum = NULL;
   coot::colour_t dummy_col;
   bool use_single_colour = false;
   // add_dots(SelHnd, mol, dum, 1.0, dummy_col, use_single_colour);
   mol->DeleteSelection(SelHnd);
}

pli::dots_representation_info_t::dots_representation_info_t(mmdb::Manager *mol,
                                                            mmdb::Manager *mol_exclude) {

   is_closed = 0;
   int SelHnd = mol->NewSelection();
   mol->SelectAtoms(SelHnd, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", "*", "*");
   coot::colour_t dummy_col;
   bool use_single_colour = false;
   // add_dots(SelHnd, mol, mol_exclude, 1.0, dummy_col, use_single_colour);
   mol->DeleteSelection(SelHnd);
}

double
pli::dots_representation_info_t::get_radius(const std::string &ele) const {

   double radius = 1.70;
   if (ele == " H")
      radius = 1.20;
   if (ele == " N")
      radius = 1.55;
   if (ele == " O")
      radius = 1.52;
   if (ele == " S")
      radius = 1.8;
   // PDBv3
   if (ele == "H")
      radius = 1.20;
   if (ele == "N")
      radius = 1.55;
   if (ele == "O")
      radius = 1.52;
   if (ele == "S")
      radius = 1.8;
   return radius;
}

// create (and later delete, of course) a new molecule by deep copying
// and assembling the passed residues.  Use that to make an atom
// selection which gets passed to
// dots_representation_info_t::solvent_exposure().  Notice that we
// pass back atom specs.
//
std::vector<std::pair<coot::atom_spec_t, float> >
pli::dots_representation_info_t::solvent_accessibilities(mmdb::Residue *res_ref,
                                                         const std::vector<mmdb::Residue *> &near_residues) const {


   std::vector<std::pair<coot::atom_spec_t, float> > v;

   // no solvent exposures if the ligand does not have *any*
   // neighbours (i.e. it's floating in space (probably the only
   // residue in the molecule).
   if (near_residues.size() == 0)
      return v;

   std::vector<mmdb::Residue *> residues = near_residues;
   residues.push_back(res_ref);

   std::pair<bool, mmdb::Manager *> mol =
      coot::util::create_mmdbmanager_from_residue_vector(residues, 0);

   if (mol.first) {

      int SelHnd = mol.second->NewSelection();
      mol.second->SelectAtoms(SelHnd, 0, res_ref->GetChainID(),
                              res_ref->GetSeqNum(), res_ref->GetInsCode(),
                              res_ref->GetSeqNum(), res_ref->GetInsCode(),
                              "*", "*", "*", "*");

      std::vector<std::pair<mmdb::Atom *, float> > se = solvent_exposure(SelHnd, mol.second);
      v.resize(se.size());
      for (unsigned int i=0; i<se.size(); i++) {
         v[i] = std::pair<coot::atom_spec_t, float> (coot::atom_spec_t(se[i].first), se[i].second);
      }

      mol.second->DeleteSelection(SelHnd);
      delete mol.second;
   }
   return v;
}


std::vector<pli::solvent_exposure_difference_helper_t>
pli::dots_representation_info_t::solvent_exposure_differences(mmdb::Residue *res_ref,
                                                              const std::vector<mmdb::Residue *> &near_residues) const {

   std::vector<pli::solvent_exposure_difference_helper_t> v;

   std::vector<mmdb::Residue *> residues = near_residues;
   residues.push_back(res_ref);

   std::pair<bool, mmdb::Manager *> mol_holo =
      coot::util::create_mmdbmanager_from_residue_vector(residues, 0);

   std::pair<bool, mmdb::Manager *> mol_apo =
      coot::util::create_mmdbmanager_from_residue_vector(near_residues, 0);

   if (mol_holo.first) {
      if (mol_apo.first) {

         for (unsigned int ir=0; ir<near_residues.size(); ir++) {
            std::string res_name = near_residues[ir]->GetResName();
            if (res_name != "HOH") {
               int SelHnd_holo = mol_holo.second->NewSelection();
               int SelHnd_apo  =  mol_apo.second->NewSelection();
               mol_holo.second->SelectAtoms(SelHnd_holo, 0, near_residues[ir]->GetChainID(),
                                            near_residues[ir]->GetSeqNum(), near_residues[ir]->GetInsCode(),
                                            near_residues[ir]->GetSeqNum(), near_residues[ir]->GetInsCode(),
                                            "*", "*", "!H", "*");
               mol_apo.second->SelectAtoms(SelHnd_apo, 0, near_residues[ir]->GetChainID(),
                                           near_residues[ir]->GetSeqNum(), near_residues[ir]->GetInsCode(),
                                           near_residues[ir]->GetSeqNum(), near_residues[ir]->GetInsCode(),
                                           "*", "*", "!H", "*");
               std::vector<std::pair<mmdb::Atom *, float> > se_holo = solvent_exposure(SelHnd_holo, mol_holo.second);
               std::vector<std::pair<mmdb::Atom *, float> > se_apo  = solvent_exposure(SelHnd_apo,   mol_apo.second);

               double se_frac_holo = 0.0;
               double se_frac_apo  = 0.0;

               for (unsigned int iah=0; iah<se_holo.size(); iah++) {
                  se_frac_holo += se_holo[iah].second;
               }
               for (unsigned int iaa=0; iaa<se_apo.size(); iaa++) {
                  se_frac_apo += se_apo[iaa].second;
               }

               if (0)
                  std::cout << "storing " << coot::residue_spec_t(near_residues[ir]) << " "
                            << near_residues[ir]->GetResName() << " "
                            << se_frac_holo << " " << se_frac_apo << std::endl;
               coot::residue_spec_t res_spec(near_residues[ir]);
               pli::solvent_exposure_difference_helper_t sed(res_spec, se_frac_holo, se_frac_apo);
               v.push_back(sed);

               mol_holo.second->DeleteSelection(SelHnd_holo);
               mol_apo.second->DeleteSelection(SelHnd_apo);
            }
         }

         delete mol_apo.second;
      }
      delete mol_holo.second;
   }
   return v;
}



// simply transfer the atoms of mol to the points vector
//
void
pli::dots_representation_info_t::pure_points(mmdb::Manager *mol) {

   is_closed = 0;
   int imod = 1;
   std::vector<clipper::Coord_orth> local_points;

   mmdb::Model *model_p = mol->GetModel(imod);
   mmdb::Chain *chain_p;
   // run over chains of the existing mol
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      mmdb::Residue *residue_p;
      mmdb::Atom *at;
      for (int ires=0; ires<nres; ires++) {
         residue_p = chain_p->GetResidue(ires);
         int n_atoms = residue_p->GetNumberOfAtoms();

         for (int iat=0; iat<n_atoms; iat++) {
            at = residue_p->GetAtom(iat);
            local_points.push_back(clipper::Coord_orth(at->x, at->y, at->z));
         }
      }
   }
   coot::colour_t col(0.3, 0.4, 0.5);
   std::pair<coot::colour_t, std::vector<clipper::Coord_orth> > p(col, local_points);
   points.push_back(p);
}


std::vector<std::pair<mmdb::Atom *, float> >
pli::dots_representation_info_t::solvent_exposure(int SelHnd_in, mmdb::Manager *mol) const {

   std::vector<std::pair<mmdb::Atom *, float> > v;
   if (mol) {

      double dot_density = 0.35;
      //
      double phi_step = 5.0 * (M_PI/180.0);
      double theta_step = 5.0 * (M_PI/180.0);
      if (dot_density > 0.0) {
         phi_step   /= dot_density;
         theta_step /= dot_density;
      }

      double water_radius = 1.4;
      double fudge = 1.0;
      mmdb::PPAtom atoms = 0;
      int n_atoms;
      mol->GetSelIndex(SelHnd_in, atoms, n_atoms);
      std::vector<double> radius(n_atoms);

      for (int iat=0; iat<n_atoms; iat++) {
         std::string ele(atoms[iat]->element);
         radius[iat] = get_radius(ele);
      }

      mmdb::PPAtom atoms_all = 0;
      int n_atoms_all;
      int SelHnd_all = mol->NewSelection();
      mol->SelectAtoms(SelHnd_all, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", "*", "*");
      mol->GetSelIndex(SelHnd_all, atoms_all, n_atoms_all);

      for (int iatom=0; iatom<n_atoms; iatom++) {
         if (! atoms[iatom]->isTer()) {
            clipper::Coord_orth centre(atoms[iatom]->x,
                                       atoms[iatom]->y,
                                       atoms[iatom]->z);
            bool even = 1;
            int n_points = 0;
            int n_sa = 0;
            for (double theta=0; theta<M_PI; theta+=theta_step) {
               double phi_step_inner = phi_step + 0.1 * pow(theta-0.5*M_PI, 2);
               for (double phi=0; phi<2*M_PI; phi+=phi_step_inner) {
                  if (even) {
                     double r = fudge * (radius[iatom] + water_radius);
                     clipper::Coord_orth pt(r*cos(phi)*sin(theta),
                                            r*sin(phi)*sin(theta),
                                            r*cos(theta));
                     pt += centre;
                     n_points++;

                     // now, is pt closer to (a water centre around)
                     // another atom?

                     bool is_solvent_accessible = 1;
                     for (int i_all=0; i_all<n_atoms_all; i_all++) {
                        // don't exclude from self
                        mmdb::Atom *other_at = atoms_all[i_all];
                        std::string other_res_name = other_at->GetResName();
                        if (other_res_name != "HOH") {
                           if (atoms[iatom] != other_at) {
                              std::string other_ele = other_at->element;
                              if (other_ele != " H") {
                                 double other_atom_r = fudge * (get_radius(other_ele) + water_radius);
                                 double other_atom_r_sq = other_atom_r * other_atom_r;
                                 clipper::Coord_orth pt_other(other_at->x, other_at->y, other_at->z);
                                 if ((pt-pt_other).lengthsq() < other_atom_r_sq) {
                                    is_solvent_accessible = 0;
                                    break;
                                 }
                              }
                           }
                        }
                     }
                     if (is_solvent_accessible)
                        n_sa++;
                  }
                  even = 1 - even;
               }
            }

            double exposure_frac = double(n_sa)/double(n_points);
            if (0)
               std::cout << "Atom " << atoms[iatom]->name << " has exposure " << n_sa << "/" << n_points
                         << " = " << exposure_frac << std::endl;
            std::pair<mmdb::Atom *, float> p(atoms[iatom], exposure_frac);
            v.push_back(p);
         }
      }
      mol->DeleteSelection(SelHnd_all); // presumably this was missing before... 20101230
   }
   return v;
}


#if 0
std::vector<coot::solvent_exposure_difference_helper_t>
pli::dots_representation_info_t::solvent_exposure_differences(mmdb::Residue *res_ref,
                                                               const std::vector<mmdb::Residue *> &near_residues) const {

   std::vector<coot::solvent_exposure_difference_helper_t> v;

   std::vector<mmdb::Residue *> residues = near_residues;
   residues.push_back(res_ref);

   std::pair<bool, mmdb::Manager *> mol_holo =
      coot::util::create_mmdbmanager_from_residue_vector(residues, 0);

   std::pair<bool, mmdb::Manager *> mol_apo =
      coot::util::create_mmdbmanager_from_residue_vector(near_residues, 0);

   if (mol_holo.first) {
      if (mol_apo.first) {

         for (unsigned int ir=0; ir<near_residues.size(); ir++) {
            std::string res_name = near_residues[ir]->GetResName();
            if (res_name != "HOH") {
               int SelHnd_holo = mol_holo.second->NewSelection();
               int SelHnd_apo  =  mol_apo.second->NewSelection();
               mol_holo.second->SelectAtoms(SelHnd_holo, 0, near_residues[ir]->GetChainID(),
                                            near_residues[ir]->GetSeqNum(), near_residues[ir]->GetInsCode(),
                                            near_residues[ir]->GetSeqNum(), near_residues[ir]->GetInsCode(),
                                            "*", "*", "!H", "*");
               mol_apo.second->SelectAtoms(SelHnd_apo, 0, near_residues[ir]->GetChainID(),
                                           near_residues[ir]->GetSeqNum(), near_residues[ir]->GetInsCode(),
                                           near_residues[ir]->GetSeqNum(), near_residues[ir]->GetInsCode(),
                                           "*", "*", "!H", "*");
               std::vector<std::pair<mmdb::Atom *, float> > se_holo = solvent_exposure(SelHnd_holo, mol_holo.second);
               std::vector<std::pair<mmdb::Atom *, float> > se_apo  = solvent_exposure(SelHnd_apo,   mol_apo.second);

               double se_frac_holo = 0.0;
               double se_frac_apo  = 0.0;

               for (unsigned int iah=0; iah<se_holo.size(); iah++) {
                  se_frac_holo += se_holo[iah].second;
               }
               for (unsigned int iaa=0; iaa<se_apo.size(); iaa++) {
                  se_frac_apo += se_apo[iaa].second;
               }

               if (0)
                  std::cout << "storing " << coot::residue_spec_t(near_residues[ir]) << " "
                            << near_residues[ir]->GetResName() << " "
                            << se_frac_holo << " " << se_frac_apo << std::endl;
               coot::residue_spec_t res_spec(near_residues[ir]);
               coot::solvent_exposure_difference_helper_t sed(res_spec, se_frac_holo, se_frac_apo);
               v.push_back(sed);

               mol_holo.second->DeleteSelection(SelHnd_holo);
               mol_apo.second->DeleteSelection(SelHnd_apo);
            }
         }

         delete mol_apo.second;
      }
      delete mol_holo.second;
   }
   return v;
}
#endif


coot::colour_t
pli::dots_representation_info_t::get_colour(const std::string &ele) const {

   coot::colour_t col(0.4, 0.4, 0.4);
   if (ele == " C")
      col = coot::colour_t(0.33, 0.4, 0.2);
   if (ele == " N")
      col = coot::colour_t(0.2, 0.2, 0.6);
   if (ele == " O")
      col = coot::colour_t(0.6, 0.2, 0.2);
   if (ele == " S")
      col = coot::colour_t(0.5, 0.5, 0.2);
   // PDB 3
   if (ele == "C")
      col = coot::colour_t(0.33, 0.4, 0.3);
   if (ele == "N")
      col = coot::colour_t(0.2, 0.2, 0.6);
   if (ele == "O")
      col = coot::colour_t(0.6, 0.2, 0.2);
   if (ele == "S")
      col = coot::colour_t(0.5, 0.5, 0.2);

   return col;
}



void
pli::dots_representation_info_t::add_dots(int SelHnd, mmdb::Manager *mol,
                                          mmdb::Manager *mol_exclude,
                                          double dot_density,
                                          const coot::colour_t &single_colour,
                                          bool use_single_colour) {

   mmdb::PPAtom atoms = NULL;
   int n_atoms;

   double phi_step = 5.0 * (M_PI/180.0);
   double theta_step = 5.0 * (M_PI/180.0);
   if (dot_density > 0.0) {
      phi_step   /= dot_density;
      theta_step /= dot_density;
   }
   mol->GetSelIndex(SelHnd, atoms, n_atoms);
   std::vector<double> radius(n_atoms);
   std::vector<double> radius_exclude;
   std::vector<coot::colour_t> colour(n_atoms);
   for (int iat=0; iat<n_atoms; iat++) {
      std::string ele(atoms[iat]->element);
      radius[iat] = get_radius(ele);
      if (use_single_colour)
         colour[iat] = single_colour;
      else
         colour[iat] = get_colour(ele);
   }

   int n_atoms_exclude = 0;
   mmdb::PPAtom atoms_exclude = NULL;
   int SelHnd_exclude = 0;
   if (mol_exclude) {
      SelHnd_exclude = mol_exclude->NewSelection();
      mol_exclude->SelectAtoms(SelHnd_exclude, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", "*", "*");
      mol_exclude->GetSelIndex(SelHnd_exclude, atoms_exclude, n_atoms_exclude);
      radius_exclude.resize(n_atoms_exclude);
      for (int iat=0; iat<n_atoms_exclude; iat++) {
         std::string ele(atoms_exclude[iat]->element);
         radius_exclude[iat] = get_radius(ele);
      }
   }

   for (int iatom=0; iatom<n_atoms; iatom++) {
      std::vector<clipper::Coord_orth> local_points;
      coot::colour_t col = colour[iatom];
      if (! atoms[iatom]->isTer()) {
         clipper::Coord_orth centre(atoms[iatom]->x,
                                    atoms[iatom]->y,
                                    atoms[iatom]->z);
         bool even = true;
         for (double theta=0; theta<M_PI; theta+=theta_step) {
            double phi_step_inner = phi_step + 0.1 * pow(theta-0.5*M_PI, 2);
            for (double phi=0; phi<2*M_PI; phi+=phi_step_inner) {
               if (even) {

                  // Is there another atom in the same residue as this
                  // atom, that is closer to pt than the centre atom?
                  // If so, don't draw this point.

                  double atom_radius = radius[iatom];

                  clipper::Coord_orth pt(atom_radius*cos(phi)*sin(theta),
                                         atom_radius*sin(phi)*sin(theta),
                                         atom_radius*cos(theta));
                  pt += centre;

                  bool draw_it = 1;

                  // it might be possible to speed this up by precalculating all dot
                  // points and then doing a findcontacts to the atoms of the atom
                  // selection. That is involved though.

                  for (int jatom=0; jatom<n_atoms; jatom++) {
                     if (jatom != iatom) {
                        if (! atoms[jatom]->isTer()) {
                           double radius_j = radius[jatom];
                           double radius_j_squared = radius_j * radius_j;
                           clipper::Coord_orth pt_j(atoms[jatom]->x, atoms[jatom]->y, atoms[jatom]->z);
                           if ((pt-pt_j).lengthsq() < radius_j_squared) {
                              draw_it = false;
                              break;
                           }
                        }
                     }
                  }

                  // and now, don't draw if far from exclude molecule
                  // atoms (if there are any).
                  //
                  if (n_atoms_exclude) {
                     if (draw_it) {
                        draw_it = false;
                        double dist_j = 4.0;
                        double dist_j_squared = dist_j * dist_j;
                        for (int jatom=0; jatom<n_atoms_exclude; jatom++) {
                           if (! atoms_exclude[jatom]->isTer()) {
                              clipper::Coord_orth pt_j(atoms_exclude[jatom]->x,
                                                       atoms_exclude[jatom]->y,
                                                       atoms_exclude[jatom]->z);
                              if ((pt-pt_j).lengthsq() < dist_j_squared) {
                                 draw_it = true;
                                 break;
                              }
                           }
                        }
                     }
                  }

                  if (draw_it) {
                     local_points.push_back(pt);
                  }
               }
               even = 1 - even;
            }
         }
      }
      std::pair<coot::colour_t, std::vector<clipper::Coord_orth> > p(col, local_points);
      points.push_back(p);
   }
   if (mol_exclude) {
      mol_exclude->DeleteSelection(SelHnd_exclude);
   }

   unsigned int n_balls = 0;
   for (unsigned int i=0; i<points.size(); i++)
      n_balls += points[i].second.size();

   // imm.setup_instancing_buffers(n_balls);

   // std::vector<> balls;
   // float scale = 0.05;
   // for (unsigned int i=0; i<points.size(); i++) {
   //    glm::vec4 col = points[i].first.to_glm();
   //    const auto &pt_vec = points[i].second;
   //    for (unsigned int j=0; j<pt_vec.size(); j++) {
   //       const auto &pt = pt_vec[j];
   //       glm::vec3 pos(pt.x(), pt.y(), pt.z());
   //       Instanced_Markup_Mesh_attrib_t ball(col, pos, scale);
   //       balls.push_back(ball);
   //    }
   // }

}
