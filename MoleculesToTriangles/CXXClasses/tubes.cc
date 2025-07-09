
#include <iostream>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>  // to_string()

#include <clipper/core/coords.h>
#include "MoleculesToTriangles/CXXClasses/NRStuff.h"
#include "coot-utils/cylinder.hh"

#include "MoleculesToTriangles/CXXSurface/CXXCoord.h"
#include "mmdb2/mmdb_selmngr.h"
#include "tubes.hh"

void secondary_structure_header_to_residue_sse(mmdb::Manager *mol);

class helix_residues_info_t {
public:
   std::string chain_id;
   std::vector<std::pair<mmdb::Residue *, float> > residues; // residue_p and bendiness
   helix_residues_info_t() {}
   helix_residues_info_t(const std::string &ch) : chain_id(ch) {}
   bool is_empty() const {
      return residues.empty();
   }
   void add_residue(mmdb::Residue *residue_p) {
      residues.push_back(std::make_pair(residue_p, -1.0f));
   }
   void clear() {
      chain_id = "";
      residues.clear();
   }
   unsigned int size() const {
      return residues.size();
   }
   void resize(unsigned int n) {
      residues.resize(n);
   }
   enum class end_cap_style { FLAT, HEMISPHERE };
};

class coil_residues_info_t : public helix_residues_info_t {
public:
   std::vector<clipper::Coord_orth> extra_points_for_front;
   std::vector<clipper::Coord_orth> extra_points_for_end;
   coil_residues_info_t() {}
   coil_residues_info_t(const std::string &ch) : helix_residues_info_t(ch) {}
};

std::vector<clipper::Coord_orth>
get_ref_coords_simple() {

   std::vector<clipper::Coord_orth> rc(12);
   rc[0] = clipper::Coord_orth (    -0.882,    -1.264,    -2.357);
   rc[1] = clipper::Coord_orth (    -0.278,    -2.269,     -1.47);
   rc[2] = clipper::Coord_orth (     0.556,    -1.575,    -0.386);
   rc[3] = clipper::Coord_orth (     0.453,      -1.9,     0.805);
   rc[4] = clipper::Coord_orth (      1.36,    -0.637,    -0.845);
   rc[5] = clipper::Coord_orth (      2.25,     0.152,     0.023);
   rc[6] = clipper::Coord_orth (     1.432,      0.88,     1.092);
   rc[7] = clipper::Coord_orth (     1.779,     0.862,     2.283);
   rc[8] = clipper::Coord_orth (     0.367,     1.501,     0.622);
   rc[9] = clipper::Coord_orth (    -0.557,     2.261,     1.476);
   rc[10] = clipper::Coord_orth (    -1.122,     1.354,     2.567);
   rc[11] = clipper::Coord_orth (    -1.156,     1.725,      3.75);
   return rc;
}

std::vector<clipper::Coord_orth>
get_coords_for_residue(mmdb::Residue *residue_p) {

   int n_atoms = residue_p->GetNumberOfAtoms();
   std::vector<clipper::Coord_orth> coords_for_residue;
   mmdb::Atom *at_N  = 0;
   mmdb::Atom *at_CA = 0;
   mmdb::Atom *at_C  = 0;
   mmdb::Atom *at_O  = 0;
   for (int iat=0; iat<n_atoms; iat++) {
      mmdb::Atom *at = residue_p->GetAtom(iat);
      if (! at->isTer()) {
         std::string atom_name = at->name;
         // PDBv3 FIXME
         if (! at_N) {
            if (atom_name == " N  ") {
               at_N = at;
            }
         }
         if (! at_CA) {
            if (atom_name == " CA ") {
               at_CA = at;
            }
         }
         if (! at_C) {
            if (atom_name == " C  ") {
               at_C = at;
            }
         }
         if (! at_O) {
            if (atom_name == " O  ") {
               at_O = at;
            }
         }
      }
   }
   if (at_N && at_CA && at_C && at_O) {
      clipper::Coord_orth con(at_N->x, at_N->y, at_N->z);
      coords_for_residue.push_back(con);
      clipper::Coord_orth coca(at_CA->x, at_CA->y, at_CA->z);
      coords_for_residue.push_back(coca);
      clipper::Coord_orth coc(at_C->x, at_C->y, at_C->z);
      coords_for_residue.push_back(coc);
      clipper::Coord_orth coo(at_O->x, at_O->y, at_O->z);
      coords_for_residue.push_back(coo);
   }
   return coords_for_residue;
 };

std::vector<helix_residues_info_t>
split_very_bent_helices(const std::vector<helix_residues_info_t> &hriv_in) {

   auto get_trace = [] (const clipper::Mat33<double> &mat) {
      double trace = mat(0,0) + mat(1,1) + mat(2,2);
      return trace;
   };

   auto strip_rot = [] (const clipper::RTop_orth &rtop_in) {
      clipper::RTop_orth rtop = rtop_in;
      rtop.trn() = clipper::Coord_orth(0,0,0);
      return rtop;
   };

   auto analyse_angle_between_residues = [strip_rot] (unsigned int helix_index,
                                                      helix_residues_info_t &helix,
                                                      unsigned int ihres,
                                                      const std::vector<clipper::Coord_orth> &ref_coords) {
      int status = 1;

      if ((ihres+3) < helix.size()) {
         std::vector<clipper::Coord_orth> coords_for_residue_0 = get_coords_for_residue(helix.residues[ihres].first);
         std::vector<clipper::Coord_orth> coords_for_residue_1 = get_coords_for_residue(helix.residues[ihres+1].first);
         std::vector<clipper::Coord_orth> coords_for_residue_2 = get_coords_for_residue(helix.residues[ihres+2].first);
         std::vector<clipper::Coord_orth> coords_for_residue_3 = get_coords_for_residue(helix.residues[ihres+3].first);
         std::vector<clipper::Coord_orth> coords_for_r_1;
         std::vector<clipper::Coord_orth> coords_for_r_2;
         if (coords_for_residue_0.size() == 4) {
            if (coords_for_residue_1.size() == 4) {
               if (coords_for_residue_2.size() == 4) {
                  if (coords_for_residue_3.size() == 4) {
                     for (unsigned int i=0; i<coords_for_residue_0.size(); i++) {
                        coords_for_r_1.push_back(coords_for_residue_0[i]);
                     }
                     for (unsigned int i=0; i<coords_for_residue_1.size(); i++) {
                        coords_for_r_1.push_back(coords_for_residue_1[i]);
                        coords_for_r_2.push_back(coords_for_residue_1[i]);
                     }
                     for (unsigned int i=0; i<coords_for_residue_2.size(); i++) {
                        coords_for_r_1.push_back(coords_for_residue_2[i]);
                        coords_for_r_2.push_back(coords_for_residue_2[i]);
                     }
                     for (unsigned int i=0; i<coords_for_residue_3.size(); i++) {
                        coords_for_r_2.push_back(coords_for_residue_3[i]);
                     }

                     clipper::RTop_orth rtop_1(ref_coords, coords_for_r_1);
                     clipper::RTop_orth rtop_2(ref_coords, coords_for_r_2);
                     clipper::RTop_orth rtop_1_rot = strip_rot(rtop_1);
                     clipper::RTop_orth rtop_2_rot = strip_rot(rtop_2);

                     clipper::Coord_orth z_pt(0,0,1);
                     clipper::Coord_orth pt_rot_1 = rtop_1_rot * z_pt;
                     clipper::Coord_orth pt_rot_2 = rtop_2_rot * z_pt;
                     double dp = clipper::Coord_orth::dot(pt_rot_1, pt_rot_2);
                     helix.residues[ihres+1].second = dp;
                     if (false)
                        std::cout << "set bendiness helix-index " << helix_index
                                  << " res-index " << ihres
                                  << " as " << dp << std::endl;
                     double theta = acos(dp);
                     if (dp < 0.6) {
                        status = 0;
                     }
                  }
               }
            }
         }
      }
      return status;
   };

   std::vector<clipper::Coord_orth> ref_coords =  get_ref_coords_simple();
   std::vector<helix_residues_info_t> hriv = hriv_in;
   std::vector<helix_residues_info_t> new_helices;
   for (unsigned int ih=0; ih<hriv.size(); ih++) {
      helix_residues_info_t &helix = hriv[ih];
      for (unsigned int ihres=0; ihres<helix.residues.size(); ihres++) {
         int angle_status = analyse_angle_between_residues(ih, helix, ihres, ref_coords);

         if (false) {
            if (angle_status == 0) { // need to split
               // So we have to split this helix here.
               // Chop residues out of this helix and make an new one.
               // Which we will add to new_helices.
               helix_residues_info_t new_helix(helix.chain_id);
               for (unsigned int jj=1; jj<helix.residues.size(); jj++) {
                  if (jj > ihres)
                     new_helix.add_residue(helix.residues[jj].first);
               }
               if (new_helix.size() > 2)
                  new_helices.push_back(new_helix);
               helix.resize(ih+1);
               break;
            }
         }
      }
   }
   if (! new_helices.empty())
      hriv.insert(hriv.end(), new_helices.begin(), new_helices.end());
   return hriv;
}

// getting the coil to match up with the helices is tricky.
// Let's try doing them both at the same time.
std::pair<std::vector<helix_residues_info_t>, std::vector<coil_residues_info_t> >
make_coil_splines_and_helicies(mmdb::Manager *mol, int atom_selection_handle, bool use_header,
                               helix_residues_info_t::end_cap_style end_cap_style) {

   std::vector<helix_residues_info_t> hriv;
   std::vector<coil_residues_info_t> coil_runs_of_residues;

   // --------------- using header mode ---------------------
   if (use_header) {

      // find the helices

      for (int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
         mmdb::Model *model_p = mol->GetModel(imod);
         if (model_p) {
            int n_helices = model_p->GetNumberOfHelices();
            for (int ih=1; ih<=n_helices; ih++) {
               mmdb:: Helix *helix_p = model_p->GetHelix(ih);
               if (helix_p) {
                  std::string helix_chain_id = helix_p->initChainID;
                  helix_residues_info_t running_helix(helix_chain_id);
                  int n_chains = model_p->GetNumberOfChains();
                  for (int ichain=0; ichain<n_chains; ichain++) {
                     mmdb::Chain *chain_p = model_p->GetChain(ichain);
                     if (chain_p) {
                        std::string this_chain_id = chain_p->GetChainID();
                        if (this_chain_id == helix_chain_id) {
                           int n_residues = chain_p->GetNumberOfResidues();
                           for (int ires=0; ires<n_residues; ires++) {
                              mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                              if (residue_p) {
                                 int res_no = residue_p->GetSeqNum();
                                 if (res_no >= helix_p->initSeqNum) {
                                    if (res_no <= helix_p->endSeqNum) {
                                       running_helix.add_residue(residue_p);
                                    }
                                 }
                                 if (res_no >= helix_p->endSeqNum) break;
                              }
                           }
                        }
                     }
                  }
                  if (! running_helix.is_empty()) hriv.push_back(running_helix);
               }
            }
         }
      }

      // find the coils

      for (int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
         mmdb::Model *model_p = mol->GetModel(imod);
         if (model_p) {
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               int n_res = chain_p->GetNumberOfResidues();
               mmdb::Residue *residue_prev_p = nullptr;
               std::string chain_id = chain_p->GetChainID();
               coil_residues_info_t running_coil(chain_id);

               for (int ires=0; ires<n_res; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  int res_no = residue_p->GetSeqNum();

                  // Is this residue in any of the helices?
                  // If not, then add it to the running coil
                  // (or start a new one using residue_prev_p)

                  bool residue_is_in_helix = false;
                  int n_helices = model_p->GetNumberOfHelices();
                  for (int ih=1; ih<=n_helices; ih++) {
                     mmdb:: Helix *helix_p = model_p->GetHelix(ih);
                     if (helix_p) {
                        std::string helix_chain_id = helix_p->initChainID;
                        if (helix_chain_id == chain_id) {
                           if (res_no >= helix_p->initSeqNum) {
                              if (res_no <= helix_p->endSeqNum) {
                                 // insertion codes make things more complex.
                                 residue_is_in_helix = true;
                                 break;
                              }
                           }
                        }
                     }
                  }
                  if (residue_is_in_helix) {
                     if (! running_coil.is_empty()) {
                        // end of a coil
                        running_coil.add_residue(residue_p);
                        coil_runs_of_residues.push_back(running_coil);
                        running_coil = coil_residues_info_t(chain_id);
                     }
                  } else {
                     if (running_coil.is_empty()) {
                        if (residue_prev_p)
                           running_coil.add_residue(residue_prev_p);
                        running_coil.add_residue(residue_p);
                     } else {
                        running_coil.add_residue(residue_p);
                     }
                  }

                  // for next round
                  residue_prev_p = residue_p;
               }
               if (! running_coil.is_empty()) coil_runs_of_residues.push_back(running_coil);
            }
         }
      }
   }

   // --------------- residue by residue mode ---------------------

   if (! use_header) {

      for (int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
         mmdb::Model *model_p = mol->GetModel(imod);
         if (model_p) {
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               int n_res = chain_p->GetNumberOfResidues();
               helix_residues_info_t running_helix;
               coil_residues_info_t running_coil;
               mmdb::Residue *residue_prev_p = nullptr;
               for (int ires=0; ires<n_res; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p->SSE == mmdb::SSE_Helix) {
                     if (running_helix.is_empty()) {
                        running_helix = helix_residues_info_t(chain_p->GetChainID());
                        if (residue_p->isInSelection(atom_selection_handle))
                           running_helix.add_residue(residue_p);
                        if (! running_coil.is_empty()) {
                           running_coil.add_residue(residue_p);
                           coil_runs_of_residues.push_back(running_coil);
                           running_coil.clear();
                        }
                     } else {
                        running_helix.add_residue(residue_p);
                     }
                     if (! running_coil.is_empty()) {
                        if (residue_p->isInSelection(atom_selection_handle)) {
                           running_coil.add_residue(residue_p);
                        }
                     }
                  } else {
                     if (! running_helix.is_empty()) {
                        hriv.push_back(running_helix);
                        running_helix.clear();
                     }
                     if (! running_coil.is_empty()) {
                        if (residue_p->isInSelection(atom_selection_handle)) {
                           running_coil.add_residue(residue_p);
                        }
                     } else {
                        running_coil = coil_residues_info_t(chain_p->GetChainID());
                        if (residue_prev_p)
                           running_coil.add_residue(residue_prev_p);
                        running_coil.add_residue(residue_p);
                     }
                  }
                  // for next round
                  residue_prev_p = residue_p;
               }
               if (! running_helix.is_empty()) hriv.push_back(running_helix);
               if (! running_coil.is_empty()) coil_runs_of_residues.push_back(running_coil);
            }
         }
      }
   }

   std::cout << ":::::::::::: coil_runs_of_residues() n coils " << coil_runs_of_residues.size() << std::endl;
   std::vector<helix_residues_info_t> clean_hriv = split_very_bent_helices(hriv);
   return std::make_pair(clean_hriv, coil_runs_of_residues);

}

coot::simple_mesh_t
make_mesh_for_helical_representation(const std::vector<helix_residues_info_t> &helices,
                                     mmdb::Manager *mol,
                                     float radius_for_helices,
                                     unsigned int n_slices_for_helices) {

   auto get_ref_coords = [] (mmdb::Manager *helix_mol) {

      std::vector<clipper::Coord_orth> coords;
      int imod = 1;
      mmdb::Model *model_p = helix_mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  int n_atoms = residue_p->GetNumberOfAtoms();
                  int res_no = residue_p->GetSeqNum();
                  if (res_no > 3) break;
                  for (int iat=0; iat<n_atoms; iat++) {
                     mmdb::Atom *at = residue_p->GetAtom(iat);
                     std::string atom_name = at->name;
                     // PDBv3 FIXME
                     if (atom_name == " N  " ||
                         atom_name == " CA " ||
                         atom_name == " C  " ||
                         atom_name == " O  ")  {
                        clipper::Coord_orth pt(at->x, at->y, at->z);
                        coords.push_back(pt);
                     }
                  }
               }
            }
         }
      }
      return coords;
   };

   // replace this at some stage
   auto hsv_to_rgb = [] (double hue, double sat, double val) {
      double c = val * sat;
      double x = c * (1 - fabs(fmod(hue / 60.0, 2) - 1));
      double m = val - c;
      double r, g, b;
      if (hue < 60) {
         r = c; g = x; b = 0;
      } else if (hue < 120) {
         r = x; g = c; b = 0;
      } else if (hue < 180) {
         r = 0; g = c; b = x;
      } else if (hue < 240) {
         r = 0; g = x; b = c;
      } else if (hue < 300) {
         r = x; g = 0; b = c;
      } else {
         r = c; g = 0; b = x;
      }
      return glm::vec4(r+m, g+m, b+m, 1.0f);
   };

   auto get_helix_colour = [hsv_to_rgb] (unsigned int ih) {
      double hue = static_cast<double>(ih) * 0.2695 *  360.0 *  0.1;
      double sat = 0.8;
      double val = 0.5;
      glm::vec4 col = hsv_to_rgb(hue, sat, val);
      return col;
   };

   auto get_helix_colour_by_bendiness = [] (int helix_index,
                                            const helix_residues_info_t &helix, unsigned int ih) {

      float bendiness = helix.residues[ih].second;
      if (bendiness == -1.0f) { // unset
         glm::vec4 unassigned_bendiness_colour(0.152, 0.12f, 1.0f, 1.0f);
         return unassigned_bendiness_colour;
      } else {
         // bendiness is typically between 0.98 and 0.999
         std::cout << "bendiness: helix-index " << helix_index << "  res-index " << ih
                   << " " << bendiness << std::endl;
         float v1 = 0.99 - bendiness;
         float v2 = v1 * 80.0f;
         float s = 0.5;
         float green = 0.5 - fabsf(v2 * 0.5f);
         float red  = 0.5 + s * v2;
         float blue = 0.5 - s * v2;
         if (red > 1.0) red = 1.0f;
         if (red < 0.0) red = 0.0f;
         if (blue > 1.0) blue = 1.0f;
         if (blue < 0.0) blue = 0.0f;
         return glm::vec4(red, green, blue, 1.0f);
      }
   };


      // the parameters are for the first residue in the triplet. If there are 2 more residues
   // then use the generate the coords for the rtop - and make the rtop
   //
   auto get_rtop = [] (unsigned int idx_res,
                       const helix_residues_info_t &helix,
                       std::vector<clipper::Coord_orth> &ref_coords) {

      bool status = false;
      clipper::RTop_orth rtop = clipper::RTop_orth::identity();

      if ((idx_res+2) < helix.residues.size()) {

         std::vector<clipper::Coord_orth> coords;
         std::vector<clipper::Coord_orth> coords_for_residue_0 = get_coords_for_residue(helix.residues[idx_res].first);
         std::vector<clipper::Coord_orth> coords_for_residue_1 = get_coords_for_residue(helix.residues[idx_res+1].first);
         std::vector<clipper::Coord_orth> coords_for_residue_2 = get_coords_for_residue(helix.residues[idx_res+2].first);
         if (coords_for_residue_0.size() == 4) {
            if (coords_for_residue_1.size() == 4) {
               if (coords_for_residue_2.size() == 4) {
                  for (unsigned int i=0; i<coords_for_residue_0.size(); i++)
                     coords.push_back(coords_for_residue_0[i]);
                  for (unsigned int i=0; i<coords_for_residue_1.size(); i++)
                     coords.push_back(coords_for_residue_1[i]);
                  for (unsigned int i=0; i<coords_for_residue_2.size(); i++)
                     coords.push_back(coords_for_residue_2[i]);
               }
            }
         }

         if (ref_coords.size() == coords.size()) {
            rtop = clipper::RTop_orth(ref_coords, coords);
            status = true;
         } else {
            std::cout << "debug:: get_rtop(): No size match for " << helix.chain_id << " "
                      << " : " << ref_coords.size() << " " << coords.size() << std::endl;
         }
      } else {
         std::cout << "ERROR:: bad residue index in get_rtop() " << idx_res << std::endl;
      }
      return std::make_pair(status, rtop);
   };

   auto strip_trn = [] (const clipper::RTop_orth &rtop_in) {
      clipper::RTop_orth rtop(rtop_in.rot(), clipper::Coord_orth(0,0,0));
      return rtop;
   };

   coot::simple_mesh_t m;

   if (! helices.empty()) {

      mmdb::Manager *helix_mol = new mmdb::Manager();
      helix_mol->ReadCoorFile("theor-helix-z-ori-v2.pdb");
      std::vector<clipper::Coord_orth> ref_coords = get_ref_coords(helix_mol);

      if (false) {
         for (unsigned int ii=0; ii<ref_coords.size(); ii++) {
            std::cout << "ref-coords: " << ii << " " << ref_coords[ii].format() << std::endl;
         }
      }

      // let's pre-calcuate the ring points
      std::vector<clipper::Coord_orth> ring_points;
      std::vector<clipper::Coord_orth> ring_point_normals;
      for (unsigned int i_pt=0; i_pt<n_slices_for_helices; i_pt++) {
         float theta = M_PI * 2.0 * static_cast<float>(i_pt) / static_cast<float>(n_slices_for_helices);
         float x = radius_for_helices * sinf(theta);
         float y = radius_for_helices * cosf(theta);
         float z = 0;
         float xn = sinf(theta);
         float yn = cosf(theta);
         float zn = 0;
         ring_points.push_back(clipper::Coord_orth(x,y,z));
         ring_point_normals.push_back(clipper::Coord_orth(xn,yn,zn));
      }

      // let's pre-calculate the end-cap points
      unsigned int n_end_cap_rings = 16;
      bool flat_end_caps = false; // so hemispheres
      std::vector<clipper::Coord_orth> end_cap_points;
      std::vector<clipper::Coord_orth> end_cap_point_normals;
      std::vector<g_triangle> end_cap_triangles;
      for (unsigned int ir=0; ir<n_end_cap_rings; ir++) {
         float scale_factor = static_cast<float>(ir+1)/static_cast<float>(n_end_cap_rings);
         // if (scale_factor == 0.0f) scale_factor = 0.0001f;
         for (unsigned int i_pt=0; i_pt<n_slices_for_helices; i_pt++) {
            float theta = M_PI * 2.0 * static_cast<float>(i_pt) / static_cast<float>(n_slices_for_helices);
            float rs = radius_for_helices * scale_factor;
            float x = rs * sinf(theta);
            float y = rs * cosf(theta);
            float zz = radius_for_helices * radius_for_helices - x * x - y * y;
            if (zz < 0.0) zz = 0.0;
            float z = - std::sqrt(zz);
            // z = 0.0;
            float xn = scale_factor * sinf(theta);
            float yn = scale_factor * cosf(theta);
            float zznn = 1.0f - xn * xn - yn * yn;
            if (zznn < 0.0) zznn = 0.0;
            float zn = std::sqrt(zznn);
            if (flat_end_caps) {
               z = 0.0f;
               xn = 0.0f; yn = 0.0f; zn = 1.0f;
            }
            end_cap_points.push_back(clipper::Coord_orth(x,y,z));
            end_cap_point_normals.push_back(clipper::Coord_orth(xn, yn, zn));
         }
      }
      // now fill end_cap_triangles
      for (unsigned int ir=0; ir<(n_end_cap_rings-1); ir++) {
         unsigned int ring_offset = ir * n_slices_for_helices;
         for (unsigned int idx=0; idx<(n_slices_for_helices-1); idx++) {
            int idx_0 = idx;
            int idx_1 = idx + 1;
            int idx_2 = idx + n_slices_for_helices;
            int idx_3 = idx + n_slices_for_helices + 1;
            g_triangle t1(idx_0, idx_1, idx_2);
            g_triangle t2(idx_1, idx_3, idx_2);
            t1.rebase(ring_offset);
            t2.rebase(ring_offset);
            end_cap_triangles.push_back(t1);
            end_cap_triangles.push_back(t2);
         }
         // and the triangles that connect the start the the end
         // (these seem to have "backward" winding compared to other start-to-end triangles - FWIW)
         g_triangle t1(0, n_slices_for_helices, n_slices_for_helices - 1);
         g_triangle t2(n_slices_for_helices, 2 * n_slices_for_helices - 1, n_slices_for_helices -1);
         t1.rebase(ring_offset);
         t2.rebase(ring_offset);
         end_cap_triangles.push_back(t1);
         end_cap_triangles.push_back(t2);
      }

      for (unsigned int ih=0; ih<helices.size(); ih++) {
         const auto &helix = helices[ih];
         if (helix.size() > 3) {
            glm::vec4 col = get_helix_colour(ih);
            std::vector<coot::api::vnc_vertex> vertices;
            std::vector<g_triangle> triangles;

            clipper::RTop_orth rtop_for_start_end_cap(clipper::RTop_orth::null());
            clipper::RTop_orth rtop_for_end_end_cap(clipper::RTop_orth::null());

            unsigned int n_rings_in_helix = 0;
            for (unsigned int ihres=0; ihres<helix.residues.size(); ihres++) {
               std::pair<bool, clipper::RTop_orth> rtop_pair = get_rtop(ihres, helix, ref_coords);
               if (rtop_pair.first) {
                  const clipper::RTop_orth &rtop = rtop_pair.second;
                  clipper::RTop_orth rtop_for_normals = strip_trn(rtop);
                  rtop_for_normals.trn() = clipper::Coord_orth(0,0,0);
                  if (ihres == 0) rtop_for_start_end_cap = rtop;
                  rtop_for_end_end_cap = rtop; // running replaced, until it isn't
                  std::vector<clipper::Coord_orth> transformed_points(ring_points.size());
                  std::vector<clipper::Coord_orth>    rotated_normals(ring_points.size());
                  col = get_helix_colour_by_bendiness(ih, helix, ihres);
                  for (unsigned int i_pt=0; i_pt<ring_points.size(); i_pt++) {
                     const clipper::Coord_orth &pt   = ring_points[i_pt];
                     const clipper::Coord_orth &n_pt = ring_point_normals[i_pt];
                     clipper::Coord_orth t_pt = rtop * pt;
                     clipper::Coord_orth t_n = rtop_for_normals * n_pt;
                     transformed_points[i_pt] = t_pt;
                     rotated_normals[i_pt] = t_n;
                     glm::vec3 v(t_pt.x(), t_pt.y(), t_pt.z());
                     glm::vec3 n(t_n.x(), t_n.y(), t_n.z());
                     coot::api::vnc_vertex vnc_v(v, n, col);
                     vertices.push_back(vnc_v);
                  }
                  n_rings_in_helix++;
               }
            }
            if (n_rings_in_helix > 1) {
               for (unsigned int i=0; i<(n_rings_in_helix-1); i++) {
                  unsigned int ring_offset = i * n_slices_for_helices;
                  const glm::vec3 &p0 = vertices[ring_offset].pos;
                  float dd_best = 9999999.9;
                  unsigned int j_best = n_slices_for_helices + 1;
                  for (unsigned int j=0; j<n_slices_for_helices; j++) {
                     const glm::vec3 &testing_pos = vertices[ring_offset + n_slices_for_helices + j].pos;
                     glm::vec3 d = testing_pos - p0;
                     float dd_test = d.x * d.x + d.y * d.y + d.z * d.z;
                     if (dd_test < dd_best) {
                        dd_best = dd_test;
                        j_best = j;
                     }
                  }
                  if (j_best != (n_slices_for_helices + 1)) {
                     for (unsigned int j=0; j<(n_slices_for_helices-1); j++) {
                        unsigned int jo_1 = j_best + j;
                        unsigned int jo_2 = j_best + j + 1;
                        if (jo_1 >= n_slices_for_helices) jo_1 -= n_slices_for_helices;
                        if (jo_2 >= n_slices_for_helices) jo_2 -= n_slices_for_helices;
                        // winding is important
                        g_triangle t1(j + 1, j, n_slices_for_helices + jo_1);
                        g_triangle t2(n_slices_for_helices + jo_1, n_slices_for_helices + jo_2, j + 1);
                        t1.rebase(ring_offset);
                        t2.rebase(ring_offset);
                        triangles.push_back(t1);
                        triangles.push_back(t2);
                     }
                     // now the join the end to the start: - I had to get a pen and paper out for this...
                     unsigned int jbo = j_best + n_slices_for_helices - 1;
                     if (j_best == 0) jbo = 2 * n_slices_for_helices - 1;
                     g_triangle t_end_1(0, n_slices_for_helices - 1, j_best + n_slices_for_helices);
                     g_triangle t_end_2(j_best + n_slices_for_helices, n_slices_for_helices -1, jbo);
                     t_end_1.rebase(ring_offset);
                     t_end_2.rebase(ring_offset);
                     triangles.push_back(t_end_1);
                     triangles.push_back(t_end_2);
                  }
               }
            }
            coot::simple_mesh_t helices_mesh(vertices, triangles);
            m.add_submesh(helices_mesh);

            // now add the end-caps:
            // start end-cap:
            glm::vec4 unassigned_bendiness_colour(0.152, 0.12f, 1.0f, 1.0f); // matches above in lambda
            col = unassigned_bendiness_colour;
            std::vector<g_triangle> end_cap_triangles_with_reversed_windings = end_cap_triangles;
            for (auto &tri : end_cap_triangles_with_reversed_windings)
               tri.reverse_winding();
            std::vector<coot::api::vnc_vertex> end_cap_vertices;
            // rtop_for_start_end_cap = clipper::RTop_orth::identity();
            std::cout << "debug:: rtop_for_start_end_cap:\n" << rtop_for_start_end_cap.format() << std::endl;
            clipper::RTop_orth rtop_for_start_end_cap_rot = strip_trn(rtop_for_start_end_cap);
            std::cout << "debug:: rtop_for_start_end_cap_rot:\n" << rtop_for_start_end_cap_rot.format() << std::endl;
            for (unsigned int i=0; i<end_cap_points.size(); i++) {
               clipper::Coord_orth ecp = end_cap_points[i];
               clipper::Coord_orth vc = rtop_for_start_end_cap * ecp;
               clipper::Coord_orth ecpn = end_cap_point_normals[i];
               ecpn = clipper::Coord_orth(ecpn.x(), ecpn.y(), -ecpn.z());
               clipper::Coord_orth nc = rtop_for_start_end_cap_rot * ecpn;
               glm::vec3 v(vc.x(), vc.y(), vc.z());
               glm::vec3 n(nc.x(), nc.y(), nc.z());
               float ss = n.x * n.x + n.y * n.y + n.z * n.z;
               coot::api::vnc_vertex vnc_v(v, n, col);
               end_cap_vertices.push_back(vnc_v);
            }
            coot::simple_mesh_t start_end_cap_sub_mesh(end_cap_vertices, end_cap_triangles_with_reversed_windings);
            m.add_submesh(start_end_cap_sub_mesh);

            // end end cap
            end_cap_vertices.clear();
            clipper::RTop_orth rtop_for_end_end_cap_rot = strip_trn(rtop_for_end_end_cap);
            for (unsigned int i=0; i<end_cap_points.size(); i++) {
               clipper::Coord_orth ecp = end_cap_points[i];
               ecp = clipper::Coord_orth(ecp.x(), ecp.y(), -ecp.z());
               clipper::Coord_orth vc = rtop_for_end_end_cap * ecp;
               clipper::Coord_orth nc = rtop_for_end_end_cap_rot * (-1.0f * end_cap_point_normals[i]);
               glm::vec3 v(vc.x(), vc.y(), vc.z());
               glm::vec3 n(-nc.x(), -nc.y(), -nc.z());
               coot::api::vnc_vertex vnc_v(v, n, col);
               end_cap_vertices.push_back(vnc_v);
            }
            coot::simple_mesh_t end_end_cap_sub_mesh(end_cap_vertices, end_cap_triangles);
            m.add_submesh(end_end_cap_sub_mesh);

         }
      }
   } // test for empty helices

   return m;
}


coot::simple_mesh_t
make_mesh_for_coil_representation(const std::vector<coil_residues_info_t> &coils,
                                  float radius,
                                  int Cn, int accuracy,
                                  unsigned int n_slices) {

   auto fcxx_to_glm = [] (const FCXXCoord &c) {
      return glm::vec3(c.x(), c.y(), c.z());
   };

   auto rotate_around_vector = [] (const glm::vec3 &direction,
                                   const glm::vec3 &position,
                                   const glm::vec3 &origin_shift,
                                   double angle) {

      glm::vec3 unit_vec = glm::normalize(direction);

      double l = unit_vec[0];
      double m = unit_vec[1];
      double n = unit_vec[2];

      double ll = l*l;
      double mm = m*m;
      double nn = n*n;
      double cosk = cos(angle);
      double sink = sin(angle);
      double I_cosk = 1.0 - cosk;

      glm::mat3 r(ll+(mm+nn)*cosk,    l*m*I_cosk-n*sink,  n*l*I_cosk+m*sink,
                  l*m*I_cosk+n*sink,  mm+(ll+nn)*cosk,    m*n*I_cosk-l*sink,
                  n*l*I_cosk-m*sink,  m*n*I_cosk+l*sink,  nn+(ll+mm)*cosk );

      glm::vec3 p1 = position - origin_shift;
      glm::vec3 p2 = r * p1;
      glm::vec3 p3 = p2 + origin_shift;
      return p3;
   };

   auto get_ring_around_mid_point = [rotate_around_vector] (const glm::vec3 &p1,
                                                            const glm::vec3 &p2,
                                                            const glm::vec3 &p3,
                                                            unsigned int n_slices,
                                                            float radius) {

      std::vector<std::pair<glm::vec3, glm::vec3> > v;
      glm::vec3 v1 = glm::normalize(p2-p1);
      glm::vec3 v2 = glm::normalize(p3-p2);
      glm::vec3 delta = (v2 - v1) * 0.5f;
      glm::vec3 delta_uv = glm::normalize(delta);
      glm::vec3 pos_start = p2 + radius *  delta_uv;
      glm::vec3 rotation_ori = v1 + v2;
      glm::vec3 origin(0,0,0);
      // this is the vector to spin around the point on the circle
      for (unsigned int i_slice=0; i_slice<n_slices; i_slice++) {
         double angle = 2 * M_PI * static_cast<double>(i_slice) / static_cast<double>(n_slices);
         glm::vec3 pt = rotate_around_vector(rotation_ori, pos_start, p2, angle);
         glm::vec3 n  = rotate_around_vector(rotation_ori, delta_uv,  origin, angle);
         v.push_back(std::make_pair(pt, n));
      }
      return v;
   };

   auto make_continuous_spline_mesh = [get_ring_around_mid_point, fcxx_to_glm]
      (const std::vector<FCXXCoord> &spline_points,
       float radius,
       unsigned int n_slices) {

      std::vector<std::vector<std::pair<glm::vec3, glm::vec3> > > rings;
      coot::simple_mesh_t m;
      if (spline_points.size() > 2) {
         std::size_t spline_points_end = spline_points.size() -1;
         for (std::size_t i=1; i<spline_points_end; i++) {
            glm::vec3 start_point = fcxx_to_glm(spline_points[i-1]);
            glm::vec3 mid_point   = fcxx_to_glm(spline_points[i]);
            glm::vec3 next_point  = fcxx_to_glm(spline_points[i+1]);
            std::vector<std::pair<glm::vec3, glm::vec3> > r =
               get_ring_around_mid_point(start_point, mid_point, next_point, n_slices, radius);
            rings.push_back(r);
         }

         glm::vec4 col(0.6, 0.6, 0.6, 1.0);
         for (const auto &ring : rings) {
            for (const auto &point : ring) {
               coot::api::vnc_vertex v(point.first, point.second, col);
               m.vertices.push_back(v);
            }
         }

         // simple minded first:
         bool simple_minded = false;
         if (simple_minded) {

            // so now we have the vertices - let's make some triangles!
            for (unsigned int i=0; i<(rings.size()-1); i++) {
               unsigned int ring_offset = i * n_slices;
               for (unsigned int j=0; j<n_slices-1; j++) {
                  // winding is important
                  g_triangle t1(ring_offset + j + 1, ring_offset + j, ring_offset + n_slices + j);
                  g_triangle t2(ring_offset + n_slices + j, ring_offset + n_slices + j + 1, ring_offset + j + 1);
                  m.triangles.push_back(t1);
                  m.triangles.push_back(t2);
               }
               // now the join the end to the start:
               g_triangle t_end_1(ring_offset, ring_offset + n_slices - 1, ring_offset + n_slices);
               g_triangle t_end_2(ring_offset + n_slices, ring_offset + n_slices -1, ring_offset + 2 * n_slices -1);
               m.triangles.push_back(t_end_1);
               m.triangles.push_back(t_end_2);
            }

         } else {

            for (unsigned int i=0; i<(rings.size()-1); i++) {
               unsigned int ring_offset = i * n_slices;
               const glm::vec3 &p0 = rings[i][0].first;
               float dd_best = 9999999.9;
               unsigned int j_best = 13;
               for (unsigned int j=0; j<n_slices; j++) {
                  const glm::vec3 &testing_pos = rings[i+1][j].first;
                  glm::vec3 d = testing_pos - p0;
                  float dd_test = d.x * d.x + d.y * d.y + d.z * d.z;
                  if (dd_test < dd_best) {
                     dd_best = dd_test;
                     j_best = j;
                  }
               }
               if (j_best != 13) {
                  for (unsigned int j=0; j<(n_slices-1); j++) {
                     unsigned jo_1 = j_best + j;
                     unsigned jo_2 = j_best + j + 1;
                     if (jo_1 >= n_slices) jo_1 -= n_slices;
                     if (jo_2 >= n_slices) jo_2 -= n_slices;
                     // winding is important
                     g_triangle t1(j + 1, j, n_slices + jo_1);
                     g_triangle t2(n_slices + jo_1, n_slices + jo_2, j + 1);
                     t1.rebase(ring_offset);
                     t2.rebase(ring_offset);
                     m.triangles.push_back(t1);
                     m.triangles.push_back(t2);
                  }
                  // now the join the end to the start: - I had to get a pen and paper out for this...
                  unsigned int jbo = j_best + n_slices - 1;
                  if (j_best == 0) jbo = 2 * n_slices - 1;
                  g_triangle t_end_1(0, n_slices - 1, j_best + n_slices);
                  g_triangle t_end_2(j_best + n_slices, n_slices -1, jbo);
                  t_end_1.rebase(ring_offset);
                  t_end_2.rebase(ring_offset);
                  m.triangles.push_back(t_end_1);
                  m.triangles.push_back(t_end_2);
               }
            }
         }
      }

      if (true) { // check mesh

         for (unsigned int ii=0; ii<m.triangles.size(); ii++) {
            const g_triangle &t1 = m.triangles[ii];
            if (t1.point_id[0]>m.vertices.size() ||
                t1.point_id[1]>m.vertices.size() ||
                t1.point_id[2]>m.vertices.size()) {
               std::cout << "ERROR:: make_continuous_spline_mesh() triangle t1 point id out of range "
                         << t1.point_id[0] << " " << t1.point_id[1] << " " << t1.point_id[2] << std::endl;
            }
         }
      }

      return m;
   };

   coot::simple_mesh_t m;
   for (unsigned int ic=0; ic<coils.size(); ic++) {
      // std::cout << "::::: coil " << ic << " of " << coils.size() << std::endl;
      const auto &coil = coils[ic];
      std::vector<FCXXCoord> ctlPts;
      for (unsigned int ir=0; ir<coil.size(); ir++) {
         mmdb::Residue *residue_p = coil.residues[ir].first;
         mmdb::Atom *ca_at = residue_p->GetAtom(" CA ");
         if (ca_at) {
            FCXXCoord fc(ca_at->x, ca_at->y, ca_at->z);
            ctlPts.push_back(fc);
         }
      }

      if (ctlPts.size() > 1) {
         CoordSpline cs;
         int nsteps =  accuracy * (ctlPts.size() - 1);
         int iinterp = 1;
         try {
            std::vector<FCXXCoord> v = cs.SplineCurve(ctlPts, nsteps, Cn, iinterp);
            coot::simple_mesh_t cs = make_continuous_spline_mesh(v, radius, n_slices);
            if (false)
               std::cout << ":::::::::::::::: continuous_spline: " << cs.vertices.size() << " " << cs.triangles.size()
                      << std::endl;
            m.add_submesh(cs);
         }
         catch (const std::runtime_error &e) {
            std::cout << "WARNING::" << e.what() << std::endl;
         }
      }
   }

   std::cout << "debug:: returning from make_mesh_for_coil_representation() " << std::endl;
   return m;
}

coot::simple_mesh_t
make_tubes_representation(mmdb::Manager *mol,
                          const std::string &atom_selection_str,
                          const std::string &colour_scheme,
                          float radius_for_coil,
                          int Cn_for_coil, int accuracy_for_coil,
                          unsigned int n_slices_for_coil,
                          int secondaryStructureUsageFlag) {

   std::cout << "---------------- start make_tubes_representation() " << std::endl;

   coot::simple_mesh_t m;
   float radius_for_helices = 2.5;
   unsigned int n_slices_for_helices = 16;
   helix_residues_info_t::end_cap_style end_cap_style = helix_residues_info_t::end_cap_style::FLAT;

   if (secondaryStructureUsageFlag == 2) { // CALC_SECONDARY_STRUCTURE in MyMolecule
      int nModels = mol->GetNumberOfModels();
      for (int iModel = 1; iModel <= nModels; iModel++){
         mmdb::Model *model = mol->GetModel(iModel);
         model->CalcSecStructure(true);
      }
   }

   if (secondaryStructureUsageFlag == 0)
      secondary_structure_header_to_residue_sse(mol);

   int sel_hnd = mol->NewSelection(); // d
   mol->Select(sel_hnd, mmdb::STYPE_RESIDUE, atom_selection_str.c_str(), mmdb::SKEY_NEW);

   bool use_header = true;
   std::pair<std::vector<helix_residues_info_t>, std::vector<coil_residues_info_t> >
      helices_and_coils_pair = make_coil_splines_and_helicies(mol, sel_hnd, use_header, end_cap_style);

   const auto &helices = helices_and_coils_pair.first;
   const auto &coils   = helices_and_coils_pair.second;

   coot::simple_mesh_t helices_mesh =
      make_mesh_for_helical_representation(helices, mol, radius_for_helices, n_slices_for_helices);

   coot::simple_mesh_t coils_mesh =
      make_mesh_for_coil_representation(coils, radius_for_coil, Cn_for_coil, accuracy_for_coil,
                                        n_slices_for_coil);

   m.add_submesh(helices_mesh);
   m.add_submesh(coils_mesh);

   std::cout << "---------------- done make_tubes_representation() " << std::endl;

   return m;
}
