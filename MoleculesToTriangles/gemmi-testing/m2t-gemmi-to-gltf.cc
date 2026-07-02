/*
 * MoleculesToTriangles/gemmi-testing/m2t-gemmi-to-gltf.cc
 *
 * Stand-alone test driver for the gemmi-native (coot::m2t) MoleculesToTriangles
 * twins. Reads coordinates with gemmi (no mmdb), builds a MolecularRepresentation,
 * converts the DisplayPrimitives to a coot::simple_mesh_t and writes a glTF file.
 *
 * The conversion mirrors coot::molecule_t::get_molecular_representation_mesh_gemmi()
 * in api/coot-molecule-moltris.cc, but stays free of mmdb and molecule_t so that it
 * can be linked against cootmoleculestotriangles alone.
 *
 * Usage:
 *    m2t-gemmi-to-gltf <coords-file> [selection-cid] [colour-scheme] [style] [out.glb]
 *
 *    colour-scheme: Element | Chains | RampChains | Secondary | BFactor   (default Chains)
 *    style:         Ribbon | MolecularSurface | VdWSurface | Cylinders |
 *                   Calpha | Sticks | Spheres | ...                       (default Ribbon)
 *
 * This program is part of Coot. Licensed under LGPL v3+.
 */

#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "MoleculesToTriangles/CXXClasses/MyMolecule-gemmi.hh"
#include "MoleculesToTriangles/CXXClasses/ColorScheme-gemmi.hh"
#include "MoleculesToTriangles/CXXClasses/MolecularRepresentationInstance-gemmi.hh"
#include "MoleculesToTriangles/CXXClasses/MolecularRepresentation-gemmi.hh"
#include "MoleculesToTriangles/CXXClasses/DisplayPrimitive.h"
#include "MoleculesToTriangles/CXXClasses/VertexColorNormalPrimitive.h"

#include "coot-utils/simple-mesh.hh"

static std::shared_ptr<coot::m2t::ColorScheme>
make_colour_scheme(const std::string &colour_scheme) {

   if (colour_scheme == "Chains")     return coot::m2t::ColorScheme::colorChainsScheme();
   if (colour_scheme == "RampChains") return coot::m2t::ColorScheme::colorRampChainsScheme();
   if (colour_scheme == "Secondary")  return coot::m2t::ColorScheme::colorBySecondaryScheme();
   if (colour_scheme == "BFactor")    return coot::m2t::ColorScheme::colorBFactorScheme();
   return coot::m2t::ColorScheme::colorByElementScheme();
}

// Append a primitive's geometry into mesh, offsetting the triangle indices.
static void
add_primitive_to_mesh(coot::simple_mesh_t &mesh, VertexColorNormalPrimitive &prim) {

   unsigned int idx_base = mesh.vertices.size();
   unsigned int n_vertices = prim.nVertices();
   unsigned int n_triangles = prim.nTriangles();
   auto vcnArray = prim.getVertexColorNormalArray();
   auto indexArray = prim.getIndexArray();

   mesh.vertices.reserve(idx_base + n_vertices);
   for (unsigned int iv = 0; iv < n_vertices; iv++) {
      VertexColorNormalPrimitive::VertexColorNormal &vcn = vcnArray[iv];
      coot::api::vnc_vertex gv;
      for (int i = 0; i < 3; i++) {
         gv.pos[i]    = vcn.vertex[i];
         gv.normal[i] = vcn.normal[i];
         gv.color[i]  = 0.0037f * vcn.color[i];  // /255 (with a touch of headroom, as in moltris)
      }
      gv.color[3] = 0.00392f * vcn.color[3];
      if ((gv.color[0] + gv.color[1] + gv.color[2]) > 10.0) {
         gv.color *= 0.00392f; // /255
         if (gv.color[3] > 0.99) gv.color[3] = 1.0;
      }
      mesh.vertices.push_back(gv);
   }

   mesh.triangles.reserve(mesh.triangles.size() + n_triangles);
   for (unsigned int it = 0; it < n_triangles; it++) {
      g_triangle gt(indexArray[3*it+0] + idx_base,
                    indexArray[3*it+1] + idx_base,
                    indexArray[3*it+2] + idx_base);
      mesh.triangles.push_back(gt);
   }
}

int main(int argc, char **argv) {

   if (argc < 2) {
      std::cout << "Usage: " << argv[0]
                << " <coords-file> [selection-cid] [colour-scheme] [style] [out.glb]\n"
                << "   colour-scheme: Element | Chains | RampChains | Secondary | BFactor   (default Chains)\n"
                << "   style:         Ribbon | MolecularSurface | VdWSurface | Cylinders | Calpha ...\n";
      return 1;
   }

   std::string coords_file    = argv[1];
   // gemmi::Selection select-all is "/*/*/*/*" (NOT the mmdb "//").
   std::string selection_cid  = (argc > 2) ? argv[2] : "/*/*/*/*";
   std::string colour_scheme  = (argc > 3) ? argv[3] : "Chains";
   std::string style          = (argc > 4) ? argv[4] : "Ribbon";
   std::string out_file       = (argc > 5) ? argv[5] : "m2t-gemmi.glb";

   int secondaryStructureUsageFlag = coot::m2t::CALC_SECONDARY_STRUCTURE;

   std::cout << "Reading " << coords_file << " ..." << std::endl;
   auto my_mol = std::make_shared<coot::m2t::MyMolecule>(coords_file, secondaryStructureUsageFlag);

   std::shared_ptr<coot::m2t::ColorScheme> cs = make_colour_scheme(colour_scheme);

   auto molrepinst = coot::m2t::MolecularRepresentationInstance::create(my_mol, cs, selection_cid, style);
   std::shared_ptr<coot::m2t::MolecularRepresentation> r = molrepinst->getRepresentation();
   if (! r) {
      std::cout << "ERROR:: null MolecularRepresentation" << std::endl;
      return 1;
   }
   r->redraw();

   coot::simple_mesh_t mesh;
   std::vector<std::shared_ptr<DisplayPrimitive> > vdp = r->getDisplayPrimitives();
   std::cout << "Got " << vdp.size() << " display primitives" << std::endl;
   for (auto it = vdp.begin(); it != vdp.end(); ++it) {
      DisplayPrimitive &dp = **it;
      DisplayPrimitive::PrimitiveType t = dp.type();
      if (t == DisplayPrimitive::PrimitiveType::SurfacePrimitive    ||
          t == DisplayPrimitive::PrimitiveType::BoxSectionPrimitive ||
          t == DisplayPrimitive::PrimitiveType::BallsPrimitive      ||
          t == DisplayPrimitive::PrimitiveType::CylinderPrimitive) {
         dp.generateArrays();
         VertexColorNormalPrimitive &prim = dynamic_cast<VertexColorNormalPrimitive &>(dp);
         add_primitive_to_mesh(mesh, prim);
      }
   }

   std::cout << "Mesh: " << mesh.vertices.size() << " vertices, "
             << mesh.triangles.size() << " triangles" << std::endl;

   if (mesh.vertices.empty()) {
      std::cout << "WARNING:: empty mesh - nothing to write" << std::endl;
      return 1;
   }

   // export_to_gltf(file, roughness, metalicity, binary): the 2nd float is the
   // PBR metallicFactor. 0.0 = dielectric (clay/plastic look), 1.0 = full metal.
   float roughness  = 0.1f;
   float metalicity = 0.9f;
   mesh.export_to_gltf(out_file, roughness, metalicity, true);
   std::cout << "Wrote " << out_file << std::endl;
   return 0;
}
