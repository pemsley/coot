
// Mesh generation code for MolecularTriangles

#ifdef USE_MOLECULES_TO_TRIANGLES

#include <memory>
#include <CXXClasses/RendererGL.h>
#include <CXXClasses/Light.h>
#include <CXXClasses/Camera.h>
#include <CXXClasses/SceneSetup.h>
#include <CXXClasses/ColorScheme.h>
#include <CXXClasses/MyMolecule.h>
#include <CXXClasses/RepresentationInstance.h>
#include <CXXClasses/MolecularRepresentationInstance.h>
#include <CXXClasses/VertexColorNormalPrimitive.h>
#include <CXXClasses/BallsPrimitive.h>

#include "molecular-mesh-generator.hh"
#include "oct.hh"

std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
molecular_mesh_generator_t::molecular_representation_instance_to_mesh(std::shared_ptr<MolecularRepresentationInstance>  molrepinst_1) {

   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > vp;

   std::shared_ptr<Representation> r = molrepinst_1->getRepresentation();
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
      if (displayPrimitive.type() == DisplayPrimitive::PrimitiveType::SurfacePrimitive ||
          displayPrimitive.type() == DisplayPrimitive::PrimitiveType::BoxSectionPrimitive ||
          displayPrimitive.type() == DisplayPrimitive::PrimitiveType::BallsPrimitive ||
          displayPrimitive.type() == DisplayPrimitive::PrimitiveType::CylinderPrimitive ){
         displayPrimitive.generateArrays();

         VertexColorNormalPrimitive &surface = dynamic_cast<VertexColorNormalPrimitive &>(displayPrimitive);
         vertices.resize(surface.nVertices());

         auto vcnArray = surface.getVertexColorNormalArray();
         for (unsigned int iVertex=0; iVertex < surface.nVertices(); iVertex++){
            s_generic_vertex &gv = vertices[iVertex];
            VertexColorNormalPrimitive::VertexColorNormal &vcn = vcnArray[iVertex];
            for (int i=0; i<3; i++) {
               gv.pos[i] = vcn.vertex[i];
               gv.normal[i] = vcn.normal[i];
               gv.color[i]  = 0.0037f * vcn.color[i];
            }
            gv.color[3] = 1.0;
         }

         auto indexArray = surface.getIndexArray();
         unsigned long nIndices = 3 * surface.nTriangles();
         triangles.resize(surface.nTriangles());
         for (unsigned int iTriangle=0; iTriangle<surface.nTriangles(); iTriangle++){
            g_triangle &gt = triangles[iTriangle];
            for (int i=0; i<3; i++)
               gt[i] = indexArray[3*iTriangle+i];
         }
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


std::vector<molecular_triangles_mesh_t>
molecular_mesh_generator_t::get_molecular_triangles_mesh(mmdb::Manager *mol,
                                                         const std::string &selection_string, // mmdb-format
                                                         const std::string &colour_scheme,
                                                         const std::string &style) {

   std::vector<molecular_triangles_mesh_t> mtm;

   if (! mol) {
      std::cout << "ERROR:: null mol " << __FUNCTION__ << "()" << std::endl;
      return mtm;
   }

   int imodel = 1;
   mmdb::Model *model_p = mol->GetModel(imodel);
   if (model_p) {
      auto my_mol = std::make_shared<MyMolecule>(mol);
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
            this_cs->addRule(SolidColorRule::colorRuleForSelectionAndName(std::shared_ptr<CompoundSelection>(new CompoundSelection(ss,ss)), cs));
         }
      }
      if (colour_scheme == "colorBySecondaryScheme" || colour_scheme == "Secondary") {
         this_cs = ss_cs;
      }

      if (colour_scheme == "colorRampChainsScheme" || colour_scheme == "Ramp") {

         std::cout << "here with colorRampChainsScheme" << std::endl;
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
      r->redraw();
      std::vector<std::shared_ptr<DisplayPrimitive> > vdp = r->getDisplayPrimitives();
      auto displayPrimitiveIter = vdp.begin();
      int i=0;
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
         if (displayPrimitive.type() == DisplayPrimitive::PrimitiveType::SurfacePrimitive ||
             displayPrimitive.type() == DisplayPrimitive::PrimitiveType::BoxSectionPrimitive ||
             displayPrimitive.type() == DisplayPrimitive::PrimitiveType::BallsPrimitive ||
             displayPrimitive.type() == DisplayPrimitive::PrimitiveType::CylinderPrimitive ){
            displayPrimitive.generateArrays();

            VertexColorNormalPrimitive &surface = dynamic_cast<VertexColorNormalPrimitive &>(displayPrimitive);
            vertices.resize(surface.nVertices());

            auto vcnArray = surface.getVertexColorNormalArray();
            for (unsigned int iVertex=0; iVertex < surface.nVertices(); iVertex++){
               s_generic_vertex &gv = vertices[iVertex];
               VertexColorNormalPrimitive::VertexColorNormal &vcn = vcnArray[iVertex];
               for (int ii=0; ii<3; ii++) {
                  gv.pos[ii]    = vcn.vertex[ii];
                  gv.normal[ii] = vcn.normal[ii];
                  gv.color[ii]  = 0.0037f * vcn.color[ii];
               }
               gv.color[3] = 1.0;
            }

            auto indexArray = surface.getIndexArray();
            unsigned long nIndices = 3 * surface.nTriangles();
            triangles.resize(surface.nTriangles());
            for (unsigned int iTriangle=0; iTriangle<surface.nTriangles(); iTriangle++){
               g_triangle &gt = triangles[iTriangle];
               for (int ii=0; ii<3; ii++)
                  gt[ii] = indexArray[3*iTriangle+ii];
            }
            // add_to_mesh(&vp, vertices, triangles);
            // mtm.add_to_mesh(vertices, triangles);
            std::string prim_name = "displayPrimitiveNameHere";
            molecular_triangles_mesh_t prim(vertices, triangles, prim_name);
            mtm.push_back(prim);
         }
      }
   }

   std::cout << "INFO:: " << __FUNCTION__  << "() n_primitives: " << mtm.size() << std::endl;
   return mtm;

}

#endif

