/*
 * src/molecular-mesh-generator-mol-tris.cc
 *
 * Copyright 2020 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

// Mesh generation code for MolecularTriangles

#include <memory>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>  // to_string()

#include <MoleculesToTriangles/CXXClasses/RendererGL.h>
#include <MoleculesToTriangles/CXXClasses/Light.h>
#include <MoleculesToTriangles/CXXClasses/Camera.h>
#include <MoleculesToTriangles/CXXClasses/SceneSetup.h>
#include <MoleculesToTriangles/CXXClasses/ColorScheme.h>
#include <MoleculesToTriangles/CXXClasses/MyMolecule.h>
#include <MoleculesToTriangles/CXXClasses/RepresentationInstance.h>
#include <MoleculesToTriangles/CXXClasses/MolecularRepresentationInstance.h>
#include <MoleculesToTriangles/CXXClasses/VertexColorNormalPrimitive.h>
#include <MoleculesToTriangles/CXXClasses/BallsPrimitive.h>

#include "molecular-mesh-generator.hh"
#include "coot-utils/oct.hh"

std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
molecular_mesh_generator_t::molecular_representation_instance_to_mesh(std::shared_ptr<MolecularRepresentationInstance> molrepinst_1,
                                                                      const std::vector<std::pair<std::string, float> > &M2T_float_params,
                                                                      const std::vector<std::pair<std::string, int> >   &M2T_int_params) {

   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > vp;

   std::shared_ptr<Representation> r = molrepinst_1->getRepresentation();

   if (true) {
      if (! M2T_float_params.empty())
         for (const auto &par : M2T_float_params)
            std::cout << "            sending MT2 float " << par.first << " " << par.second << std::endl;
      if (! M2T_int_params.empty())
         for (const auto &par : M2T_int_params)
            std::cout << "            sending MT2 int " << par.first << " " << par.second << std::endl;
   }

   if (! M2T_float_params.empty())
      for (const auto &par : M2T_float_params)
         r->updateFloatParameter(par.first, par.second);
   if (! M2T_int_params.empty())
      for (const auto &par : M2T_int_params)
         r->updateIntParameter(par.first, par.second);

   r->redraw();
   std::vector<std::shared_ptr<DisplayPrimitive> > vdp = r->getDisplayPrimitives();
   auto displayPrimitiveIter = vdp.begin();
   std::vector<s_generic_vertex> vertices;
   std::vector<g_triangle> triangles;

#if 0
   for (displayPrimitiveIter=vdp.begin(); displayPrimitiveIter != vdp.end(); displayPrimitiveIter++) {
      DisplayPrimitive &displayPrimitive = **displayPrimitiveIter;
      if (displayPrimitive.type() == DisplayPrimitive::PrimitiveType::BallsPrimitive){
         BallsPrimitive &ballsCollection = dynamic_cast<BallsPrimitive &>(displayPrimitive);
         const std::vector<Ball> &balls = ballsCollection.getBalls();
         std::vector<Ball>::const_iterator ball_iter;
         for (ball_iter=balls.begin(); ball_iter!=balls.end(); ++ball_iter) {
            float radius = ball_iter->radius;
            std::cout << " radii " << radius << " " << ball_iter->radiusAlongNormal << std::endl;
            glm::vec3 p(ball_iter->centre.x(), ball_iter->centre.y(), ball_iter->centre.z());
            glm::vec3 n(ball_iter->normal.x(), ball_iter->normal.y(), ball_iter->normal.z());
            glm::vec4 c(ball_iter->color.x(),  ball_iter->color.y(),  ball_iter->color.z(), 1.0f);
            std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > ball_info =
               make_octasphere(3, p, radius, ball_iter->radiusAlongNormal, n, c);
            add_to_mesh(&vp, ball_info);
         }
      }
   }
#endif

   for (displayPrimitiveIter=vdp.begin(); displayPrimitiveIter != vdp.end(); displayPrimitiveIter++) {
      DisplayPrimitive &displayPrimitive = **displayPrimitiveIter;
      if (displayPrimitive.type() == DisplayPrimitive::PrimitiveType::SurfacePrimitive    ||
          displayPrimitive.type() == DisplayPrimitive::PrimitiveType::BoxSectionPrimitive ||
          displayPrimitive.type() == DisplayPrimitive::PrimitiveType::BallsPrimitive      ||
          displayPrimitive.type() == DisplayPrimitive::PrimitiveType::CylinderPrimitive ){
         displayPrimitive.generateArrays();

         VertexColorNormalPrimitive &surface = dynamic_cast<VertexColorNormalPrimitive &>(displayPrimitive);
         vertices.resize(surface.nVertices());

         auto vcnArray = surface.getVertexColorNormalArray();
         for (unsigned int iVertex=0; iVertex < surface.nVertices(); iVertex++){
            s_generic_vertex &gv = vertices[iVertex];
            VertexColorNormalPrimitive::VertexColorNormal &vcn = vcnArray[iVertex];
            for (int i=0; i<3; i++) {
               gv.pos[i]    = vcn.vertex[i];
               gv.normal[i] = vcn.normal[i];
               gv.color[i]  = 0.0037f * vcn.color[i];
            }
            gv.color[3] = 1.0;
            if ((gv.color[0] + gv.color[1] + gv.color[2]) > 10.0) {
               gv.color *= 0.00392; // /255
            }
         }

         auto indexArray = surface.getIndexArray();
         triangles.resize(surface.nTriangles());
         for (unsigned int iTriangle=0; iTriangle<surface.nTriangles(); iTriangle++){
            g_triangle &gt = triangles[iTriangle];
            for (int i=0; i<3; i++)
               gt[i] = indexArray[3*iTriangle+i];
         }
         if (false)
            std::cout << "in molecular_representation_instance_to_mesh(): Here B adding "
                      << vertices.size() << " vertices and " << triangles.size() << " triangles" << std::endl;
         add_to_mesh(&vp, vertices, triangles);
      }
   }
   return vp;
}

void
molecular_mesh_generator_t::add_selection_and_colour(const std::string &sel, const std::string &col) {

   //    selection_colours.push_back(std::pair<std::string, std::string>(sel, col));
}


int
molecular_mesh_generator_t::get_max_resno_for_polymer(mmdb::Chain *chain_p) const {

   int res_no_max = -1;
   int nres = chain_p->GetNumberOfResidues();
   for (int ires=nres-1; ires>=0; ires--) {
      mmdb::Residue *residue_p = chain_p->GetResidue(ires);
      if (residue_p) {
         int seq_num = residue_p->GetSeqNum();
         if (seq_num > res_no_max) {
            if (residue_p->isAminoacid() || residue_p->isNucleotide()) {
               res_no_max = seq_num;
            }
         }
      }
   }
   return res_no_max;
}

std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
molecular_mesh_generator_t::get_molecular_triangles_mesh(mmdb::Manager *mol,
                                                         mmdb::Chain *chain_p,
                                                         const std::string &colour_scheme,
                                                         const std::string &style,
                                                         int secondary_structure_usage_flag,
                                                         const std::vector<std::pair<std::string, float> > &M2T_float_params,
                                                         const std::vector<std::pair<std::string, int> >   &M2T_int_params) {

   if (true)
      std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ get_molecular_triangles_mesh() "
                << " chain-id: " << chain_p->GetChainID() << " colour_scheme: "
                << colour_scheme << " style " << style << std::endl;

   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > vp;

   if (! mol) {
      std::cout << "ERROR:: null mol " << __FUNCTION__ << "()" << std::endl;
      return vp;
   }

   std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ get_molecular_triangles_mesh() Here A " << std::endl;
   auto my_mol = std::make_shared<MyMolecule>(mol, secondary_structure_usage_flag);
   auto ss_cs = ColorScheme::colorBySecondaryScheme();
   auto ribbon_ramp_cs = ColorScheme::colorRampChainsScheme();
   std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ get_molecular_triangles_mesh() Here B " << std::endl;
   auto chain_cs = ColorScheme::colorChainsScheme();
   auto this_cs = chain_cs;

   std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ get_molecular_triangles_mesh() Here C " << std::endl;

   if (colour_scheme == "colorRampChainsScheme" || colour_scheme == "Ramp") {
      this_cs = ribbon_ramp_cs;
      std::cout << "mmgt::get_molecular_triangles_mesh(): with colorRampChainsScheme" << std::endl;
      int nres = chain_p->GetNumberOfResidues();
      if (nres > 0) {
         std::string atom_selection_str = "//" + std::string(chain_p->GetChainID());
         int min_resno = chain_p->GetResidue(0)->GetSeqNum();
         int max_resno = get_max_resno_for_polymer(chain_p);
         std::cout << "mmgt::get_molecular_triangles_mesh(): with min_resno " << min_resno << std::endl;
         std::cout << "mmgt::get_molecular_triangles_mesh(): with max_resno " << max_resno << std::endl;
         if (max_resno > 0) {
            std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ get_molecular_triangles_mesh() Here D " << std::endl;
            AtomPropertyRampColorRule apcrr;
            apcrr.setNumberOfRampPoints(max_resno-min_resno);
            apcrr.setStartValue(min_resno);
            apcrr.setEndValue(max_resno);
            auto apcrr_p = std::make_shared<AtomPropertyRampColorRule> (apcrr);
            std::cout << "DEBUG:: get_molecular_triangles_mesh() A calls addRule" << std::endl;
            ribbon_ramp_cs->addRule(apcrr_p);
            std::cout << "this_cs " << this_cs << std::endl;
            std::shared_ptr<MolecularRepresentationInstance> molrepinst =
               MolecularRepresentationInstance::create(my_mol, this_cs, atom_selection_str, style);
            vp = molecular_representation_instance_to_mesh(molrepinst, M2T_float_params, M2T_int_params);
            std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ get_molecular_triangles_mesh() Here E "
                      << atom_selection_str << " " << style << std::endl;
         }
      }
   }
   if (true)
      std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ get_molecular_triangles_mesh() --- done ---"
                << std::endl;
   return vp;
}


std::vector<molecular_triangles_mesh_t>
molecular_mesh_generator_t::get_molecular_triangles_mesh(mmdb::Manager *mol,
                                                         const std::string &selection_string, // mmdb-format
                                                         const std::string &colour_scheme,
                                                         const std::string &style,
                                                         int secondary_structure_usage_flag,
                                                         const std::vector<std::pair<std::string, float> > &M2T_float_params,
                                                         const std::vector<std::pair<std::string, int> > &M2T_int_params
                                                         ) {

   std::vector<molecular_triangles_mesh_t> mtm;

   // std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > vp;

   if (! mol) {
      std::cout << "ERROR:: null mol " << __FUNCTION__ << "()" << std::endl;
      return mtm;
   }

   int imodel = 1;
   mmdb::Model *model_p = mol->GetModel(imodel);
   if (model_p) {
      auto my_mol = std::make_shared<MyMolecule>(mol, secondary_structure_usage_flag);
      auto ss_cs = ColorScheme::colorBySecondaryScheme();
      auto ribbon_ramp_cs = ColorScheme::colorRampChainsScheme();
      auto chain_cs = ColorScheme::colorChainsScheme();
      auto ele_cs   = ColorScheme::colorByElementScheme();
      auto this_cs = chain_cs;

      if (colour_scheme == "colorChainsScheme" || colour_scheme == "Chain")
         this_cs = chain_cs;
      if (colour_scheme == "colorRampChainsScheme" || colour_scheme == "Ramp") {
         this_cs = ribbon_ramp_cs;
      }
      if (colour_scheme == "colorByElementScheme" || colour_scheme == "Element") {
         this_cs = ele_cs;
         for (unsigned int i=0; i<selection_colours.size(); i++) {
            const std::string ss = selection_colours[i].first;
            const std::string cs = selection_colours[i].second;
            std::cout << "DEBUG:: get_molecular_triangles_mesh() B calls addRule" << std::endl;
            this_cs->addRule(SolidColorRule::colorRuleForSelectionAndName(std::shared_ptr<CompoundSelection>(new CompoundSelection(ss,ss)), cs));
         }
      }
      if (colour_scheme == "colorBySecondaryScheme" || colour_scheme == "Secondary") {
         this_cs = ss_cs;
      }

      if (colour_scheme == "colorRampChainsScheme" || colour_scheme == "Ramp") {

         std::cout << "debug:: in get_molecular_triangles_mesh() with colorRampChainsScheme" << std::endl;
         for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
            mmdb::Model *model_inner_p = mol->GetModel(imod);
            if (model_p) {
               int n_chains = model_inner_p->GetNumberOfChains();
               for (int ichain=0; ichain<n_chains; ichain++) {
                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                  int nres = chain_p->GetNumberOfResidues();
                  if (nres > 0) {
                     int min_resno = chain_p->GetResidue(0)->GetSeqNum();
                     int max_resno = chain_p->GetResidue(nres-1)->GetSeqNum();
                     AtomPropertyRampColorRule apcrr;
                     apcrr.setStartValue(min_resno);
                     apcrr.setEndValue(max_resno);
                     auto apcrr_p = std::make_shared<AtomPropertyRampColorRule> (apcrr);
                     std::cout << "DEBUG:: get_molecular_triangles_mesh() C calls addRule" << std::endl;
                     ribbon_ramp_cs->addRule(apcrr_p);
                  }
               }
            }
         }
      }

      std::string atom_selection_str = selection_string;

      std::shared_ptr<MolecularRepresentationInstance> molrepinst_1 =
         MolecularRepresentationInstance::create(my_mol, this_cs, atom_selection_str, style);

      // now fill vp.first and vp.second
      // use fcxxcoord_to_glm() if needed.

      std::shared_ptr<Representation> r = molrepinst_1->getRepresentation();
      for (const auto &par : M2T_float_params)
         r->updateFloatParameter(par.first, par.second);
      for (const auto &par : M2T_int_params)
         r->updateIntParameter(par.first, par.second);
      r->redraw();
      std::vector<std::shared_ptr<DisplayPrimitive> > vdp = r->getDisplayPrimitives();
      auto displayPrimitiveIter = vdp.begin();
      // int i=0;
      std::vector<s_generic_vertex> vertices;
      std::vector<g_triangle> triangles;

#if 0
      for (displayPrimitiveIter=vdp.begin(); displayPrimitiveIter != vdp.end(); displayPrimitiveIter++) {
         DisplayPrimitive &displayPrimitive = **displayPrimitiveIter;
         if (displayPrimitive.type() == DisplayPrimitive::PrimitiveType::BallsPrimitive){
            BallsPrimitive &ballsCollection = dynamic_cast<BallsPrimitive &>(displayPrimitive);
            const std::vector<Ball> &balls = ballsCollection.getBalls();
            std::vector<Ball>::const_iterator ball_iter;
            for (ball_iter=balls.begin(); ball_iter!=balls.end(); ++ball_iter) {
               float radius = ball_iter->radius;
               std::cout << " radii " << radius << " " << ball_iter->radiusAlongNormal << std::endl;
               glm::vec3 p(ball_iter->centre.x(), ball_iter->centre.y(), ball_iter->centre.z());
               glm::vec3 n(ball_iter->normal.x(), ball_iter->normal.y(), ball_iter->normal.z());
               glm::vec4 c(ball_iter->color.x(),  ball_iter->color.y(),  ball_iter->color.z(), 1.0f);
               std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > ball_info =
                  make_octasphere(3, p, radius, ball_iter->radiusAlongNormal, n, c);
               add_to_mesh(&vp, ball_info);
            }
         }
      }
#endif

      unsigned int idx = 0;
      for (displayPrimitiveIter=vdp.begin(); displayPrimitiveIter != vdp.end(); displayPrimitiveIter++, ++idx) {

         molecular_triangles_mesh_t current_primitives;

         DisplayPrimitive &displayPrimitive = **displayPrimitiveIter;
         if (displayPrimitive.type() == DisplayPrimitive::PrimitiveType::SurfacePrimitive    ||
             displayPrimitive.type() == DisplayPrimitive::PrimitiveType::BoxSectionPrimitive ||
             displayPrimitive.type() == DisplayPrimitive::PrimitiveType::BallsPrimitive      ||
             displayPrimitive.type() == DisplayPrimitive::PrimitiveType::CylinderPrimitive
             ) {

            current_primitives.type_index = displayPrimitive.type();
            displayPrimitive.generateArrays();

            VertexColorNormalPrimitive &surface = dynamic_cast<VertexColorNormalPrimitive &>(displayPrimitive);
            vertices.resize(surface.nVertices());

            auto vcnArray = surface.getVertexColorNormalArray();
            for (unsigned int iVertex=0; iVertex < surface.nVertices(); iVertex++) {
               s_generic_vertex &gv = vertices[iVertex];
               VertexColorNormalPrimitive::VertexColorNormal &vcn = vcnArray[iVertex];
               for (int i=0; i<3; i++) {
                  gv.pos[i]    = vcn.vertex[i];
                  gv.normal[i] = vcn.normal[i];
                  gv.color[i]  = 0.0037f * vcn.color[i];
               }
               gv.color[3] = 1.0;

               // 20220401-PE This is what the code used to say
               //     // BoxSectionPrimitive have not had their colours scaled for some reason.
               //     if (displayPrimitive.type() == DisplayPrimitive::PrimitiveType::BoxSectionPrimitive) {
               //        // std::cout << "vertex iVertex " << iVertex << " " << glm::to_string(gv.color) << std::endl;
               //        for (int ii=0; ii<3; ii++) gv.color[ii] /= 255.0;
               //     }

               // But on Mac, we saw black strands. So maybe they *have* been scaled now, but not on my PC
               // (where I am writng this).
               // So add in a hack, looking for colours more than 2 (say). If they don't exist, then
               // we don't need to scale down the colours.

               bool do_scale_colours = false;
               if (displayPrimitive.type() == DisplayPrimitive::PrimitiveType::BoxSectionPrimitive)
                  for (int ii=0; ii<3; ii++)
                     if (gv.color[ii]  > 2.0) {
                        do_scale_colours = true;
                        break;
                     }

               if (do_scale_colours) {
                  // BoxSectionPrimitive have not had their colours scaled for some reason.
                  if (displayPrimitive.type() == DisplayPrimitive::PrimitiveType::BoxSectionPrimitive) {
                     // std::cout << "vertex iVertex " << iVertex << " " << glm::to_string(gv.color) << std::endl;
                     for (int ii=0; ii<3; ii++) gv.color[ii] /= 255.0;
                  }
               }
            }

            auto indexArray = surface.getIndexArray();
            triangles.clear();
            triangles.resize(surface.nTriangles());
            for (unsigned int iTriangle=0; iTriangle<surface.nTriangles(); iTriangle++){
               g_triangle &gt = triangles[iTriangle];
               for (int ii=0; ii<3; ii++)
                  gt[ii] = indexArray[3*iTriangle+ii];
            }
            current_primitives.add_to_mesh(vertices, triangles);
         }

         if (! current_primitives.vertices.empty())
            mtm.push_back(current_primitives);

      } // end of displayPrimitive loop

   } // valid model_p test

   // std::cout << "INFO:: " << __FUNCTION__  << "() n_primitives: " << mtm.size() << std::endl;
   return mtm;

}

#include "ud-colour-rule.hh"

std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
molecular_mesh_generator_t::get_molecular_triangles_mesh_for_ribbon_with_user_defined_residue_colours(mmdb::Manager *mol, mmdb::Chain *chain_p,
                                                                                                      const std::vector<std::pair<unsigned int, coot::colour_holder> > &user_defined_colours,
                                                                                                      int secondary_structure_usage_flag,
                                                                                                      const std::vector<std::pair<std::string, float> > &M2T_float_params,
                                                                                                      const std::vector<std::pair<std::string, int> >   &M2T_int_params) {

   auto debug_the_colours = [] (mmdb::Manager *mol, ud_colour_rule &cr) {

      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  int n_atoms = residue_p->GetNumberOfAtoms();
                  for (int iat=0; iat<n_atoms; iat++) {
                     mmdb::Atom *at = residue_p->GetAtom(iat);
                     if (! at->isTer()) {
                        FCXXCoord fc_col = cr.colorForAtom(at);
                        std::cout << "debug:: get_molecular_triangles_mesh_for_ribbon_with_user_defined_residue_colours() atom colour "
                                  << coot::atom_spec_t(at) << " " << fc_col[0] << " " << fc_col[1] << " " << fc_col[2] << "\n";
                     }
                  }
               }
            }
         }
      }
   };

   auto debug_the_colour_table = [] (const std::vector<std::pair<unsigned int, coot::colour_holder> > &user_defined_colours) {

      std::cout << "-------------- mmg get_molecular_triangles_mesh_for_ribbon_with_user_defined_residue_colours debug_the_colour_table start --- " << std::endl;
      for (unsigned int i=0; i<user_defined_colours.size(); i++) {
         unsigned int idx = user_defined_colours[i].first;
         const auto &col  = user_defined_colours[i].second;
         std::cout << "debug colour " << idx << " " << col << std::endl;
      }
      std::cout << "-------------- mmg get_molecular_triangles_mesh_for_ribbon_with_user_defined_residue_colours debug_the_colour_table done --- " << std::endl;
   };

   // -------------------------------------------------------------------------

   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > vp;

   auto my_mol = std::make_shared<MyMolecule>(mol, secondary_structure_usage_flag);
   auto ribbon_ramp_cs = ColorScheme::colorRampChainsScheme();
   auto this_cs = ribbon_ramp_cs;

   int udd_handle = mol->GetUDDHandle(mmdb::UDR_ATOM, "user-defined-atom-colour-index");

   if (false) {
      debug_the_colour_table(user_defined_colours);
   }

   ud_colour_rule cr(udd_handle, mol, user_defined_colours);
   std::string atom_selection_str = "//" + std::string(chain_p->GetChainID());
   std::shared_ptr<CompoundSelection> comp_sel = std::make_shared<CompoundSelection>(atom_selection_str);
   cr.setCompoundSelection(comp_sel);
   auto udcr_p = std::make_shared<ud_colour_rule>(cr);
   std::cout << "DEBUG:: get_molecular_triangles_mesh_for_ribbon_with_user_defined_residue_colours() calls addRule" << std::endl;
   ribbon_ramp_cs->addRule(udcr_p);

   if (false)
      debug_the_colours(mol, cr);

   int nres = chain_p->GetNumberOfResidues();
   if (nres > 0) {
      std::string style = "Ribbon";
      std::shared_ptr<MolecularRepresentationInstance> molrepinst =
         MolecularRepresentationInstance::create(my_mol, this_cs, atom_selection_str, style);
      vp = molecular_representation_instance_to_mesh(molrepinst, M2T_float_params, M2T_int_params);
   }

   return vp;
}

