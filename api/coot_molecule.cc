
#include <clipper/ccp4/ccp4_map_io.h>

#include "coot-utils/coot-coord-utils.hh"
#include "coot_molecule.hh"
#include "ideal/pepflip.hh"
#include "rama-plot-phi-psi.hh"

bool
coot::molecule_t::is_valid_model_molecule() const {

   bool status = false;
   if (atom_sel.mol)
      status = true;
   return status;

}

bool
coot::molecule_t::is_valid_map_molecule() const {

   bool status = false;
   if (! xmap.is_null()) {
      status = true;
   }
   return status;
}

std::pair<bool, coot::residue_spec_t>
coot::molecule_t::cid_to_residue_spec(const std::string &cid) {

   bool status = false;
   coot::residue_spec_t rs;
   if (atom_sel.mol) {
      int selHnd = atom_sel.mol->NewSelection(); // d
      mmdb::Residue **SelResidues;
      int nSelResidues = 0;
      atom_sel.mol->Select(selHnd, mmdb::STYPE_CHAIN, cid.c_str(), mmdb::SKEY_NEW);
      atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);
      atom_sel.mol->DeleteSelection(selHnd);
      if (nSelResidues > 0) {
         mmdb::Residue *residue_p = SelResidues[0];
         coot::residue_spec_t rs_inner(residue_p);
         rs = rs_inner;
         status = true;
      }
   }
   return std::make_pair(status, rs);
}

int coot::molecule_t::flipPeptide(const coot::residue_spec_t &rs, const std::string &alt_conf) {

   int result = coot::pepflip(atom_sel.mol, rs.chain_id, rs.res_no, rs.ins_code, alt_conf);
   return result;

}


void
coot::molecule_t::make_bonds_type_checked(coot::protein_geometry *geom_p, const char *caller) {

   bool debug = false;

   // Note caller can be 0 (e.g. with clang) - so be aware of that when debugging.

   std::string caller_s("NULL");
   if (caller) caller_s = std::string(caller);

   bool is_intermediate_atoms_molecule = false; // 20221005-PE IMPORT-HACK
   if (debug)
      std::cout << "debug:: plain make_bonds_type_checked() --------start--------- called by "
                << caller_s << "() with is_intermediate_atoms_molecule: " << is_intermediate_atoms_molecule
                << std::endl;
   if (debug)
      std::cout << "--------- make_bonds_type_checked() called with bonds_box_type "
                << bonds_box_type << " vs "
                << "NORMAL_BONDS " << coot::NORMAL_BONDS << " "
                << "BONDS_NO_HYDROGENS " << coot::BONDS_NO_HYDROGENS << " "
                << "COLOUR_BY_CHAIN_BONDS " << coot::COLOUR_BY_CHAIN_BONDS << " "
                << "COLOUR_BY_MOLECULE_BONDS " << coot::COLOUR_BY_MOLECULE_BONDS << " "
                << "CA_BONDS " << coot::CA_BONDS << " "
                << "CA_BONDS_PLUS_LIGANDS " << coot::CA_BONDS_PLUS_LIGANDS << " "
                << "COLOUR_BY_USER_DEFINED_COLOURS___BONDS " << coot::COLOUR_BY_USER_DEFINED_COLOURS____BONDS << " "
                << std::endl;

   // Delete this in due course
   // graphics_info_t g; // urgh!  (But the best solution?)

   bool force_rebonding = true; // if we get here, this must be true (?)

   // coot::protein_geometry *geom_p = g.Geom_p();

   std::set<int> dummy;

   if (bonds_box_type == coot::NORMAL_BONDS)
      makebonds(geom_p, dummy);

#if 0
   if (bonds_box_type == coot::BONDS_NO_HYDROGENS)
      makebonds(geom_p, dummy);
   if (bonds_box_type == coot::CA_BONDS)
      make_ca_bonds();
   if (bonds_box_type == coot::COLOUR_BY_CHAIN_BONDS || bonds_box_type == coot::COLOUR_BY_CHAIN_GOODSELL) {
      // Baah, we have to use the static in graphics_info_t here as it
      // is not a per-molecule property.
      std::set<int> s;
      bool goodsell_mode = false;
      if (bonds_box_type == coot::COLOUR_BY_CHAIN_GOODSELL)
         goodsell_mode = true;
      make_colour_by_chain_bonds(s, g.rotate_colour_map_on_read_pdb_c_only_flag, goodsell_mode, force_rebonding);
   }
   if (bonds_box_type == coot::COLOUR_BY_MOLECULE_BONDS)
      make_colour_by_molecule_bonds(force_rebonding);
   if (bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS)
      make_ca_plus_ligands_bonds(g.Geom_p());
   if (bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS_AND_SIDECHAINS)
      make_ca_plus_ligands_and_sidechains_bonds(g.Geom_p());
   if (bonds_box_type == coot::BONDS_NO_WATERS)
      bonds_no_waters_representation();
   if (bonds_box_type == coot::BONDS_SEC_STRUCT_COLOUR)
      bonds_sec_struct_representation();
   if (bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR)
      ca_plus_ligands_sec_struct_representation(g.Geom_p());
   if (bonds_box_type == coot::COLOUR_BY_RAINBOW_BONDS)
      ca_plus_ligands_rainbow_representation(g.Geom_p());
   if (bonds_box_type == coot::COLOUR_BY_OCCUPANCY_BONDS)
      occupancy_representation();
   if (bonds_box_type == coot::COLOUR_BY_B_FACTOR_BONDS)
      b_factor_representation();
   if (bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS_B_FACTOR_COLOUR)
      b_factor_representation_as_cas();
   if (bonds_box_type == coot::COLOUR_BY_USER_DEFINED_COLOURS____BONDS)
      user_defined_colours_representation(g.Geom_p(), true, g.draw_missing_loops_flag); // hack,
                                                             // because we need to remeber somehow
                                                             // if this was called with all-atom or CA-only.
                                                             // See c-interface.cc
                                                             // graphics_to_user_defined_atom_colours_representation()
                                                             // Perhaps we need two functions
                                                             // user_defined_colours_representation_all()
                                                             // user_defined_colours_representation_Calpha() [+ ligands]

   if (bonds_box_type == coot::COLOUR_BY_USER_DEFINED_COLOURS_CA_BONDS)
      user_defined_colours_representation(g.Geom_p(), false, g.draw_missing_loops_flag); // hack,

#endif

   // bleugh. But if we don't do this here, where *do* we do it?
   // Should the glci be passed to make_bonds_type_checked()?  Urgh.
   // That is called from many places....
   //
#ifndef EMSCRIPTEN

#if 0 // 20221005-PE not sure what these are
   // all these will need to be changed or removed
   update_additional_representations(glci, g.Geom_p());
   update_fixed_atom_positions();
   update_ghosts();
   update_extra_restraints_representation();
#endif

#endif

   if (debug) {
      std::cout << "debug:: -------------- make_bonds_type_checked() done " << std::endl;
   }
}

std::vector<std::pair<coot::Cartesian, coot::util::phi_psi_t> >
coot::molecule_t::ramachandran_validation() const {

   std::vector<std::pair<coot::Cartesian, coot::util::phi_psi_t> > v;

   rama_plot::phi_psis_for_model_t ppm(atom_sel.mol);
   // This: std::map<coot::residue_spec_t, phi_psi_t> ppm.phi_psi  is now filled

   std::map<residue_spec_t, rama_plot::phi_psi_t>::const_iterator it;
   for (it=ppm.phi_psi.begin(); it!=ppm.phi_psi.end(); ++it) {
      const auto &phi_psi(it->second);
      mmdb::Residue *rp = phi_psi.residue_prev;
      mmdb::Residue *rt = phi_psi.residue_this;
      mmdb::Residue *rn = phi_psi.residue_next;
      if (rp && rt && rn) {
         mmdb::Atom *at = rt->GetAtom(" CA "); // 20221006-PE alt-confs another day
         if (at) {
            coot::Cartesian pos(at->x, at->y, at->z);
            coot::util::phi_psi_t cupp(rp, rt, rn);
            auto p = std::make_pair(pos, cupp);
            v.push_back(p);
         }
      }
   }

   return v;
}


// returns either the specified atom or null if not found
mmdb::Atom *
coot::molecule_t::get_atom(const coot::atom_spec_t &atom_spec) const {

   mmdb::Atom *at = coot::util::get_atom(atom_spec, atom_sel.mol);
   return at;
}


// returns either the specified residue or null if not found
mmdb::Residue *
coot::molecule_t::get_residue(const coot::residue_spec_t &residue_spec) const {

   mmdb::Residue *r = coot::util::get_residue(residue_spec, atom_sel.mol);
   return r;

}


int
coot::molecule_t::writeMap(const std::string &file_name) const {

   int status = 0;

   if (! xmap.is_null()) {
      clipper::CCP4MAPfile mapout;
      mapout.open_write(file_name);
      mapout.export_xmap(xmap);
      mapout.close_write();
      status = 1;
   }

   return status;

}

coot::simple_mesh_t
coot::molecule_t::get_map_contours_mesh(clipper::Coord_orth position, float radius, float contour_level) const {

   coot::simple_mesh_t m;

   // should the contouring be threaded?

   return m;
}
