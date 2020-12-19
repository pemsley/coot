#ifndef GRAPHICS_LIGAND_MESH_MOLECULE_HH
#define GRAPHICS_LIGAND_MESH_MOLECULE_HH

#include "lidia-core/lig-build.hh"
#include "lidia-core/lbg-molfile.hh"
#include "coot-colour.hh"
#include "Shader.hh"
#include "LigandViewMesh.hh"
#include "HUDTextureMesh.hh"

class graphics_ligand_mesh_atom : public lig_build::atom_t {
public:
   graphics_ligand_mesh_atom(const lig_build::pos_t &pos_in, const std::string &ele_in, int formal_charge_in) :
      lig_build::atom_t(pos_in, ele_in, formal_charge_in) { }
   coot::colour_t colour;
   coot::colour_t get_colour(bool against_a_dark_background) const;
};

class graphics_ligand_mesh_bond : public lig_build::bond_t {
public:
   graphics_ligand_mesh_bond(int first, int second, lig_build::bond_t::bond_type_t type) :
      lig_build::bond_t(first, second, type) {}
   void gl_bond(const lig_build::pos_t &pos_1, const lig_build::pos_t &pos_2,
		bool shorten_first, bool shorten_second, lig_build::bond_t::bond_type_t bt);
   void gl_bond_double_aromatic_bond(const lig_build::pos_t &pos_1, const lig_build::pos_t &pos_2,
				     bool shorten_first, bool shorten_second);
   void gl_bond_double_bond(const lig_build::pos_t &pos_1, const lig_build::pos_t &pos_2,
                            bool shorten_first, bool shorten_second);
};

class graphics_ligand_mesh_molecule_t : public lig_build::molecule_t<graphics_ligand_mesh_atom,
                                                                     graphics_ligand_mesh_bond>  {


   class position_triple_t {
      glm::vec2 pos_t_to_glm(const lig_build::pos_t &p) {
         return glm::vec2(p.x, p.y);
      }
   public:
      glm::vec2 positions[3];
      position_triple_t(const lig_build::pos_t &p1,
                        const lig_build::pos_t &p2,
                        const lig_build::pos_t &p3) {
         positions[0] = pos_t_to_glm(p1);
         positions[1] = pos_t_to_glm(p2);
         positions[2] = pos_t_to_glm(p3);
      }
   };

   void init_from_molfile_molecule(const lig_build::molfile_molecule_t &mol);
   void fill_mesh();
   std::pair<std::vector<glm::vec2>, std::vector<position_triple_t> > fill_mesh_bonds();
   void fill_mesh_atoms();
   
public:
   // this is a static in graphics_info_t
   graphics_ligand_mesh_molecule_t() { imol = -1;
      scale_correction.first = true;
      scale_correction.second = 1.0;
   }
   ~graphics_ligand_mesh_molecule_t();
   LigandViewMesh mesh; // contains vectors for both lines and triangles. Don't use indexing to draw.
   HUDTextureMesh hud_texture_tmesh; // for (non-Carbon) atoms (text)
   std::pair<bool, double> scale_correction;
   bool setup_from(int imol, mmdb::Residue *residue_p, const std::string &alt_conf, coot::protein_geometry *geom_p);
   int imol;
   void draw(Shader *lines_shader_p, Shader *hud_text_shader_p,
             float aspect_ratio, const std::map<GLchar, FT_character> &ft_characters);

};


#endif // GRAPHICS_LIGAND_MESH_MOLECULE_HH
