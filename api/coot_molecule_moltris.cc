
#include "coot_molecule.hh"

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


//! Add a colour rule: eg. ("//A", "red")
void
coot::molecule_t::add_colour_rule(const std::string &selection, const std::string &colour_name) {


   colour_rules.push_back(std::make_pair(selection, colour_name));

}


//! delete all the colour rules
void
coot::molecule_t::delete_colour_rules() {
   colour_rules.clear();
}

coot::simple_mesh_t
coot::molecule_t::get_molecular_representation_mesh(const std::string &atom_selection_str,
                                                    const std::string &colour_scheme,
                                                    const std::string &style) const {

   auto get_max_resno_for_polymer = [] (mmdb::Chain *chain_p) {
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
   };

   auto molecular_representation_instance_to_mesh = [] (std::shared_ptr<MolecularRepresentationInstance> molrepinst) {
      coot::simple_mesh_t mesh;

      std::shared_ptr<Representation> r = molrepinst->getRepresentation();
      r->redraw(); // 20221207-PE this was missing!
      std::vector<std::shared_ptr<DisplayPrimitive> > vdp = r->getDisplayPrimitives();

      auto displayPrimitiveIter = vdp.begin();
      for (displayPrimitiveIter=vdp.begin(); displayPrimitiveIter != vdp.end(); displayPrimitiveIter++) {
         DisplayPrimitive &displayPrimitive = **displayPrimitiveIter;
         if (displayPrimitive.type() == DisplayPrimitive::PrimitiveType::SurfacePrimitive    ||
             displayPrimitive.type() == DisplayPrimitive::PrimitiveType::BoxSectionPrimitive ||
             displayPrimitive.type() == DisplayPrimitive::PrimitiveType::BallsPrimitive      ||
             displayPrimitive.type() == DisplayPrimitive::PrimitiveType::CylinderPrimitive ){
            displayPrimitive.generateArrays();

            coot::simple_mesh_t submesh;
            VertexColorNormalPrimitive &surface = dynamic_cast<VertexColorNormalPrimitive &>(displayPrimitive);
            submesh.vertices.resize(surface.nVertices());

            auto vcnArray = surface.getVertexColorNormalArray();
            for (unsigned int iVertex=0; iVertex < surface.nVertices(); iVertex++){
               coot::api::vnc_vertex &gv = submesh.vertices[iVertex];
               VertexColorNormalPrimitive::VertexColorNormal &vcn = vcnArray[iVertex];
               for (int i=0; i<3; i++) {
                  gv.pos[i]    = vcn.vertex[i];
                  gv.normal[i] = vcn.normal[i];
                  gv.color[i]  = 0.0037f * vcn.color[i];
               }
               if ((gv.color[0] + gv.color[1] + gv.color[2]) > 10.0) {
                  gv.color *= 0.00392; // /255
               }
               gv.color[3] = 1.0;
            }

            auto indexArray = surface.getIndexArray();
            submesh.triangles.resize(surface.nTriangles());
            for (unsigned int iTriangle=0; iTriangle<surface.nTriangles(); iTriangle++){
               g_triangle &gt = submesh.triangles[iTriangle];
               for (int i=0; i<3; i++)
                  gt[i] = indexArray[3*iTriangle+i];
            }
            if (false)
               std::cout << "in molecular_representation_instance_to_mesh(): Here B adding "
                         << submesh.vertices.size() << " vertices and "
                         << submesh.triangles.size() << " triangles" << std::endl;
            mesh.add_submesh(submesh);
         }
      }
      
      return mesh;
   };

   class chain_info_t {
   public:
      mmdb::Chain *chain_p;
      int resno_min;
      int resno_max;
      chain_info_t(mmdb::Chain *c, int min, int max) : chain_p(c), resno_min(min), resno_max(max) {}
   };

   auto get_chains_in_selection = [] (std::shared_ptr<MyMolecule> my_mol,
                                const std::string &atom_selection_str) {
      std::vector<chain_info_t> ci;
      MyMolecule *mm = my_mol.get();
      mmdb::Manager *mol = my_mol.get()->mmdb;

      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            int resno_max = -999999;
            int resno_min = 999999;
            if (n_res > 0) {
               for (int ires=0; ires<n_res; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     int res_no = residue_p->GetSeqNum();
                     std::string res_name(residue_p->GetResName());
                     if (res_name != "HOH") {
                        if (res_no > resno_max) resno_max = res_no;
                        if (res_no < resno_min) resno_min = res_no;
                     }
                  }
               }
            }
            ci.push_back(chain_info_t(chain_p, resno_min, resno_max));
         }
      }

      return ci;
   };

   auto ramp_chains = [get_chains_in_selection,
                       molecular_representation_instance_to_mesh]
      (std::shared_ptr<MyMolecule> my_mol,
       const std::string &atom_selection_str,
       const std::string &style) {

      coot::simple_mesh_t mesh;
      auto ramp_cs  = ColorScheme::colorRampChainsScheme();

      std::vector<chain_info_t> ci = get_chains_in_selection(my_mol, atom_selection_str);

      for (const auto &ch : ci) {

         std::string atom_selection_str = "//" + std::string(ch.chain_p->GetChainID());
         AtomPropertyRampColorRule apcrr;
         apcrr.setStartValue(ch.resno_min);
         apcrr.setEndValue(ch.resno_max);
         auto apcrr_p = std::make_shared<AtomPropertyRampColorRule> (apcrr);
         ramp_cs->addRule(apcrr_p);
         std::shared_ptr<MolecularRepresentationInstance> molrepinst =
            MolecularRepresentationInstance::create(my_mol, ramp_cs, atom_selection_str, style);
         coot::simple_mesh_t submesh = molecular_representation_instance_to_mesh(molrepinst);
         mesh.add_submesh(submesh);
      }
      return mesh;
   };

   coot::simple_mesh_t mesh;

   auto my_mol = std::make_shared<MyMolecule>(atom_sel.mol);
   // auto chain_cs = ColorScheme::colorChainsScheme();
   auto chain_cs = ColorScheme::colorChainsSchemeWithColourRules(colour_rules);
   auto ele_cs   = ColorScheme::colorByElementScheme();
   auto ss_cs    = ColorScheme::colorBySecondaryScheme();
   auto bf_cs    = ColorScheme::colorBFactorScheme();
   auto this_cs  = chain_cs; // default
   if (colour_scheme == "Chains")    this_cs = chain_cs;
   if (colour_scheme == "Element")   this_cs = ele_cs;
   if (colour_scheme == "BFactor")   this_cs = bf_cs;
   if (colour_scheme == "Secondary") this_cs = ss_cs;
   if (colour_scheme == "RampChains") {
      mesh = ramp_chains(my_mol, atom_selection_str, style);
   } else {
      std::shared_ptr<MolecularRepresentationInstance> molrepinst =
         MolecularRepresentationInstance::create(my_mol, this_cs, atom_selection_str, style);
      mesh = molecular_representation_instance_to_mesh(molrepinst);

      if (false) {
         for (unsigned int i=0; i<mesh.vertices.size(); i++) {
            const auto &vertex = mesh.vertices[i];
            std::cout << i << " " << glm::to_string(vertex.pos) << " " << glm::to_string(vertex.color) << std::endl;
         }
      }
   }

   mesh.fill_colour_map(); // for blendering

   return mesh;
}
