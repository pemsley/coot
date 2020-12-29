
#include <iostream>


// Do we need this block?
//
#include <iomanip>
#ifdef USE_PYTHON
#include "Python.h"
#endif
#include "graphics-info.h" // because gtk_gl_area_attach_buffers()
///////////////////////////

#include "graphics-ligand-mesh-molecule.hh"



#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#include "lidia-core/rdkit-interface.hh"
#endif

graphics_ligand_mesh_molecule_t::~graphics_ligand_mesh_molecule_t() {}

// I used to understand polymorphism. Now this looks to me like a mystical incantation.
template<class graphics_ligand_mesh_atom, class graphics_ligand_mesh_bond> lig_build::molecule_t<graphics_ligand_mesh_atom, graphics_ligand_mesh_bond>::~molecule_t() {}

bool
graphics_ligand_mesh_molecule_t::setup_from(int imol_in, mmdb::Residue *residue_p,
                                            const std::string &alt_conf,
                                            coot::protein_geometry *geom_p) {

   std::cout << "graphics_ligand_mesh_molecule_t::setup_from() !!!!!!!!! " << std::endl;

   bool status = false;

   imol = imol_in;

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

   imol = imol_in;
   if (residue_p) {
      try {
	 std::string res_name = residue_p->GetResName();
	 std::pair<bool, coot::dictionary_residue_restraints_t> p =
	    geom_p->get_monomer_restraints_at_least_minimal(res_name, imol);
	 if (! p.first) {
	    std::cout << "DEBUG:: graphics_ligand_molecule: No restraints for \""
		      << res_name << "\"" << std::endl;
	 } else {
	    const coot::dictionary_residue_restraints_t &restraints = p.second;
	    RDKit::RWMol rdkm = coot::rdkit_mol(residue_p, restraints, alt_conf);

	    // std::cout << "--------------------- graphics-ligand-view setup_from() 1 " << std::endl;
	    // coot::debug_rdkit_molecule(&rdkm);

	    unsigned int n_atoms = rdkm.getNumAtoms();
	    if (n_atoms > 1) {
	       // return a kekulize mol
	       RDKit::RWMol rdk_mol_with_no_Hs = coot::remove_Hs_and_clean(rdkm);

	       double weight_for_3d_distances = 0.005;
	       int mol_2d_depict_conformer = coot::add_2d_conformer(&rdk_mol_with_no_Hs,
                                                                    weight_for_3d_distances);

	       // why is there no connection between a lig_build molecule_t
	       // and a rdkit molecule conformer?

	       // For now hack around using a molfile molecule...
	       //
	       // I think I should have a rdkit_mol->lig_build::molecule_t converter
	       // (for later).

	       lig_build::molfile_molecule_t m =
		  coot::make_molfile_molecule(rdk_mol_with_no_Hs, mol_2d_depict_conformer);

	       if (false) {
		  std::cout << "make_molfile_molecule() makes molecule: " << std::endl;
		  m.debug();
	       }

	       init_from_molfile_molecule(m);

	       status = true; // OK, if we got to here...
               
            }
         }
      }
      catch (const std::runtime_error &coot_error) {
	 std::cout << coot_error.what() << std::endl;
      }
      catch (const std::exception &rdkit_error) {
	 std::cout << rdkit_error.what() << std::endl;
      }
   }

#endif

   return status;

}

void
graphics_ligand_mesh_molecule_t::init_from_molfile_molecule(const lig_build::molfile_molecule_t &mol_in) {

   atoms.clear();
   bonds.clear();

   for (unsigned int iat=0; iat<mol_in.atoms.size(); iat++) {
      const lig_build::molfile_atom_t &at_in = mol_in.atoms[iat];
      graphics_ligand_mesh_atom at(lig_build::pos_t(at_in.atom_position.x(), at_in.atom_position.y()),
                                   at_in.element, at_in.formal_charge);
      at.atom_name = at_in.name;
      at.aromatic  = at_in.aromatic;
      // what about chiral here (a lig_build::atom_t does not have chiral information).
      atoms.push_back(at);
   }

   for (unsigned int ib=0; ib<mol_in.bonds.size(); ib++) {
      const lig_build::molfile_bond_t &bond_in = mol_in.bonds[ib];
      graphics_ligand_mesh_bond b(bond_in.index_1, bond_in.index_2, bond_in.bond_type);
      bonds.push_back(b);
   }
   assign_ring_centres();
   scale_correction = mol_in.get_scale_correction(); // so that the median bond length is 1.0

   fill_mesh();

   hud_texture_tmesh.setup_quad();
}

void
graphics_ligand_mesh_molecule_t::fill_mesh() {

   gtk_gl_area_attach_buffers(GTK_GL_AREA(graphics_info_t::glareas[0])); // needed?

   std::pair<std::vector<glm::vec2>, std::vector<position_triple_t> > p = fill_mesh_bonds();

   std::vector<glm::vec2> triangle_vertices(3*p.second.size());
   for (unsigned int i=0; i<p.second.size(); i++) {
      triangle_vertices[3*i  ] = p.second[i].positions[0];
      triangle_vertices[3*i+1] = p.second[i].positions[1];
      triangle_vertices[3*i+2] = p.second[i].positions[2];
   }

   mesh.import(p.first, triangle_vertices);

   fill_mesh_atoms();
}

std::pair<std::vector<glm::vec2>, std::vector<graphics_ligand_mesh_molecule_t::position_triple_t> >
graphics_ligand_mesh_molecule_t::fill_mesh_bonds() {

   auto pos_t_to_glm = [] (const lig_build::pos_t &p) {
                          return glm::vec2(p.x, p.y);
                       };

   // add to line_vertices and polygon_vertices
   auto gl_bond = [pos_t_to_glm] (const graphics_ligand_mesh_bond &bond,
                                  const lig_build::pos_t &pos_1_raw,
                                  const lig_build::pos_t &pos_2_raw,
                                  bool shorten_first, bool shorten_second,
                                  std::vector<glm::vec2> &line_vertices,
                                  std::vector<position_triple_t> &polygon_vertices) {

                     lig_build::bond_t::bond_type_t bt = bond.get_bond_type();
                     double shorten_fraction = 0.8;

                     lig_build::pos_t pos_1 = pos_1_raw;
                     lig_build::pos_t pos_2 = pos_2_raw;

                     // fraction_point() returns a point that is (say) 0.8 of the way
                     // from p1 (first arg) to p2 (second arg).
                     //
                     if (shorten_first)
                        pos_1 = lig_build::pos_t::fraction_point(pos_2_raw, pos_1_raw, shorten_fraction);
                     if (shorten_second)
                        pos_2 = lig_build::pos_t::fraction_point(pos_1_raw, pos_2_raw, shorten_fraction);

                     switch (bt) {
                     case lig_build::bond_t::SINGLE_BOND:
                        // Add new cases, a bit of a hack of course.
                     case lig_build::bond_t::SINGLE_OR_DOUBLE:
                     case lig_build::bond_t::SINGLE_OR_AROMATIC:
                     case lig_build::bond_t::AROMATIC_BOND:
                     case lig_build::bond_t::DELOC_ONE_AND_A_HALF:
                     case lig_build::bond_t::BOND_ANY:
                        {
                           line_vertices.push_back(pos_t_to_glm(pos_1));
                           line_vertices.push_back(pos_t_to_glm(pos_2));
                        }
                        break;

                     case lig_build::bond_t::DOUBLE_BOND:
                     case lig_build::bond_t::DOUBLE_OR_AROMATIC:
                        {
                           if (bond.have_centre_pos()) {
                              lig_build::pos_t pos_1_local = pos_1;
                              lig_build::pos_t pos_2_local = pos_2;
                              double shorten_fraction_for_double = 0.82;
                              if (shorten_first)
                                 pos_1 = lig_build::pos_t::fraction_point(pos_2_local, pos_1_local, shorten_fraction_for_double);
                              if (shorten_second)
                                 pos_2 = lig_build::pos_t::fraction_point(pos_1_local, pos_2_local, shorten_fraction_for_double);

                              std::pair<lig_build::pos_t, lig_build::pos_t> p =
                                 bond.make_double_aromatic_short_stick(pos_1_local, pos_2_local, shorten_first, shorten_second);

                              line_vertices.push_back(pos_t_to_glm(pos_1));
                              line_vertices.push_back(pos_t_to_glm(pos_2));
                              line_vertices.push_back(pos_t_to_glm(p.first));
                              line_vertices.push_back(pos_t_to_glm(p.second));

                           } else {
                              std::pair<std::pair<lig_build::pos_t, lig_build::pos_t>, std::pair<lig_build::pos_t, lig_build::pos_t> > p =
                                 bond.make_double_bond(pos_1, pos_2, shorten_first, shorten_second);
                              line_vertices.push_back(pos_t_to_glm(p.first.first));
                              line_vertices.push_back(pos_t_to_glm(p.first.second));
                              line_vertices.push_back(pos_t_to_glm(p.second.first));
                              line_vertices.push_back(pos_t_to_glm(p.second.second));
                           }
                        }
                        break;

                     case lig_build::bond_t::TRIPLE_BOND:
                        {
                           lig_build::pos_t buv = (pos_2-pos_1).unit_vector();
                           lig_build::pos_t buv_90 = buv.rotate(90);
                           double small = 0.25;
                           lig_build::pos_t p1 = pos_1 + buv_90 * small;
                           lig_build::pos_t p2 = pos_2 + buv_90 * small;
                           lig_build::pos_t p3 = pos_1;
                           lig_build::pos_t p4 = pos_2;
                           lig_build::pos_t p5 = pos_1 - buv_90 * small;
                           lig_build::pos_t p6 = pos_2 - buv_90 * small;

                           line_vertices.push_back(pos_t_to_glm(p1));
                           line_vertices.push_back(pos_t_to_glm(p2));
                           line_vertices.push_back(pos_t_to_glm(p3));
                           line_vertices.push_back(pos_t_to_glm(p4));
                           line_vertices.push_back(pos_t_to_glm(p5));
                           line_vertices.push_back(pos_t_to_glm(p6));
                        }
                        break;

                     case lig_build::bond_t::IN_BOND:
                        {
                           // set of lines
                           std::vector<std::pair<lig_build::pos_t, lig_build::pos_t> > vp =
                              lig_build::pos_t::make_wedge_in_bond(pos_1, pos_2);
                           if (vp.size()) {
                              for (unsigned int i=0; i<vp.size(); i++) {
                                 line_vertices.push_back(pos_t_to_glm(vp[i].first));
                                 line_vertices.push_back(pos_t_to_glm(vp[i].second));
                              }
                           }
                        }
                        break;

                     case lig_build::bond_t::OUT_BOND:
                        {
                           // filled shape
                           std::vector<lig_build::pos_t> v =
                              lig_build::pos_t::make_wedge_out_bond(pos_1, pos_2);
                           if (v.size() == 4) {
                              // make 2 triangle from this
                              polygon_vertices.push_back(position_triple_t(v[0], v[1], v[2]));
                              polygon_vertices.push_back(position_triple_t(v[2], v[3], v[0]));
                           }
                        }
                        break;
                     case lig_build::bond_t::BOND_UNDEFINED:
                        break;
                     }
                  };

   std::vector<glm::vec2> line_vertices;
   std::vector<position_triple_t> polygon_vertices; // for wedge bonds;

   for (unsigned int ib=0; ib<bonds.size(); ib++) {
      int idx_1 = bonds[ib].get_atom_1_index();
      int idx_2 = bonds[ib].get_atom_2_index();
      if ((idx_1 != UNASSIGNED_INDEX) && (idx_2 != UNASSIGNED_INDEX)) {
         const graphics_ligand_mesh_bond &bond = bonds[ib];

	 // c.f. widgeted_bond_t::construct_internal()
	 bool shorten_first = false;
	 bool shorten_second = false;
	 if (atoms[idx_1].atom_id != "C") {
	    shorten_first = true;
	 }
	 if (atoms[idx_2].atom_id != "C") {
	    shorten_second = true;
	 }
	 lig_build::pos_t pos_1 =  atoms[idx_1].atom_position;
	 lig_build::pos_t pos_2 =  atoms[idx_2].atom_position;

         // now the modern OpenGL equivalent of gl_bond()

         gl_bond(bond, pos_1, pos_2, shorten_first, shorten_second, line_vertices, polygon_vertices);

      }
   }

   return std::pair<std::vector<glm::vec2>, std::vector<graphics_ligand_mesh_molecule_t::position_triple_t> > (line_vertices, polygon_vertices);

}

void
graphics_ligand_mesh_molecule_t::fill_mesh_atoms() {

   // this function should be called assign_atom_colours() we don't do anything with the mesh.
   // All complexities happen at draw time.

   bool background_black_flag = true;
   for (unsigned int iat=0; iat<atoms.size(); iat++) {
      const std::string &ele = atoms[iat].element;
      if (ele != "C") {
	 std::vector<unsigned int> local_bonds = bonds_having_atom_with_atom_index(iat);
	 bool gl_flag = true;
	 lig_build::atom_id_info_t atom_id_info = make_atom_id_by_using_bonds(iat, ele, local_bonds, gl_flag);
	 coot::colour_t col = atoms[iat].get_colour(background_black_flag); // using ele
         atoms[iat].colour = col;
	 if (true)
	    std::cout << "in fill_mesh_atoms() atom_index " << iat << " with charge " << atoms[iat].charge
		      << " made atom_id_info:\n" << atom_id_info << " assign colour " << col << std::endl;
      }
   }


}

coot::colour_t 
graphics_ligand_mesh_atom::get_colour(bool dark_bg) const {

   coot::colour_t col;

   if (element == "Br") {
      col.col[0] = 0.66;
      col.col[1] = 0.2;
      col.col[2] = 0.2;
   }
   if (element == "I") {
      col.col[0] = 0.42;
      col.col[1] = 0.1;
      col.col[2] = 0.8;
   }
   if ((element == "F") || (element == "Cl")) {
      col.col[0] = 0.3;
      col.col[1] = 0.7;
      col.col[2] = 0.3;
   }
   if (element == "O") {
      col.col[0] = 0.9;
      col.col[1] = 0.3;
      col.col[2] = 0.3;
   }
   if (element == "P")  {
      col.col[0] = 0.7;
      col.col[1] = 0.3;
      col.col[2] = 0.9;
   }
   if ((element == "S") || (element == "Se")) {
      col.col[0] = 0.76;
      col.col[1] = 0.76;
      col.col[2] = 0.2;
   }
   if (element == "N")  {
      col.col[0] = 0.5;
      col.col[1] = 0.5;
      col.col[2] = 1.0;
   }
   return col;
}


void
graphics_ligand_mesh_molecule_t::draw(Shader *shader_p, Shader *hud_text_shader_p,
                                      float widget_height, float widget_width, const std::map<GLchar, FT_character> &ft_characters) {

   auto pos_t_to_glm = [] (const lig_build::pos_t &p) {
                          return glm::vec2(p.x, p.y);
                       };

   // draw line bonds and wedge bonds
   //
   mesh.draw(shader_p, widget_height, widget_width);

   // draw atoms (text)
   //
   float gl_flag = true;
   for (unsigned int iat=0; iat<atoms.size(); iat++) {
      graphics_ligand_mesh_atom &atom = atoms[iat];
      const std::string &ele = atom.element;
      if (ele != "C") {
	 std::vector<unsigned int> local_bonds = bonds_having_atom_with_atom_index(iat);
         lig_build::atom_id_info_t atom_id_info = make_atom_id_by_using_bonds(iat, ele, local_bonds, gl_flag);
         for (unsigned int i=0; i<atom_id_info.n_offsets(); i++) {
            const lig_build::offset_text_t &offset = atom_id_info.offsets[i];
            std::string label = offset.text;
            glm::vec2 pos(-0.61, -0.61);
            pos += 0.05 * pos_t_to_glm(atom.atom_position);
            if (offset.text_pos_offset == lig_build::offset_text_t::UP)
               pos.y += 0.03;
            if (offset.text_pos_offset == lig_build::offset_text_t::DOWN)
               pos.y -= 0.03;
            if (offset.subscript)
               pos.y -= 0.012;
            if (offset.superscript)
               pos.y -= 0.012;
            pos.x += 0.03 * 0.08 * offset.tweak.x;
            pos.y += 0.03 * 0.08 * offset.tweak.y;
            float sc = 0.000184;
            if (offset.subscript)   sc *= 0.8;
            if (offset.superscript) sc *= 0.8;
            glm::vec2 scales(sc, sc);
            hud_texture_tmesh.set_position_and_scales(pos, scales);
            std::cout << "draw() calling draw_label(): iat " << iat << " ioff " << i
                      << " \"" << offset << " " << offset.text << "\" "
                      << offset.text.length() << " colour " << atom.colour << std::endl;
            glm::vec4 colour = atom.colour.to_glm();
            hud_texture_tmesh.draw_label(label, colour, hud_text_shader_p, ft_characters);
         }
      }
   }
}
