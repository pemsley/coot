#ifndef COOT_PLI_RESIDUE_CIRCLE_HH
#define COOT_PLI_RESIDUE_CIRCLE_HH

#include <string>
#include <vector>
#include <clipper/core/coords.h>

#include "lidia-core/lig-build.hh"
#include "geometry/residue-and-atom-specs.hh"
#include "bond-to-ligand.hh"

class residue_circle_t {
   bool   se_diff_set_;
   double se_holo;
   double se_apo;
   int stacking_type;
   std::vector<std::string> ligand_ring_atom_names;
   std::string ligand_cation_atom_name;
   bool is_a_primary_residue_; // primary residues are stacking or
   // bond and are placed first looking
   // for non-overlaping
   // interaction-bond lines.
public:
   // CATION_PI_STACKING sets ligand_cationic_atom_name, not the
   // ligand_ring_atom_names vector.
   //
   enum { PI_PI_STACKING,
      PI_CATION_STACKING, // for cations on the protein residues (ligand pi)
      CATION_PI_STACKING, // for cations on the ligand (protein TRY, PRO, TRP)
   };

   clipper::Coord_orth trans_rel_pos_3d;
   clipper::Coord_orth residue_centre_real;

   lig_build::pos_t pos; // coordinate system of the ligand atoms
   coot::residue_spec_t spec;
   std::string residue_type;
   std::string residue_label;
   std::vector<bond_to_ligand_t> bonds_to_ligand;
   double water_dist_to_protein;
   residue_circle_t(const clipper::Coord_orth &pos_in,
                    const clipper::Coord_orth &click_pos_in,
                    coot::residue_spec_t spec_in,
                    const std::string &type_in,
                    const std::string &label_in) :
      trans_rel_pos_3d(pos_in),
      residue_centre_real(click_pos_in),
      spec(spec_in),
      residue_type(type_in),
      residue_label(label_in) {
      se_holo = 0.0;
      se_apo  = 0.0;
      se_diff_set_ = 0;
      stacking_type = -1; // should have enumerated value
      is_a_primary_residue_ = 0;
      water_dist_to_protein = 100;
   }
   void set_canvas_pos(const lig_build::pos_t &pos_in) {
      pos = pos_in;
   }
   void add_bond_to_ligand(const bond_to_ligand_t &bl) {
      bonds_to_ligand.push_back(bl);
      is_a_primary_residue_ = 1;
   }
   void set_solvent_exposure_diff(double se_holo_in, double se_apo_in) {
      se_holo = se_holo_in;
      se_apo = se_apo_in;
      se_diff_set_ = 1;
   }
   void set_stacking(const std::string &type_string,
                     const std::vector<std::string> &ligand_ring_atom_names_in,
                     const std::string &ligand_cation_atom_name_in) {
      is_a_primary_residue_ = 1;
      if (type_string == "pi-pi")
         stacking_type = PI_PI_STACKING;
      if (type_string == "pi-cation")
         stacking_type = PI_CATION_STACKING;
      if (type_string == "cation-pi")
         stacking_type = CATION_PI_STACKING;
      if (stacking_type == PI_PI_STACKING)
         ligand_ring_atom_names = ligand_ring_atom_names_in;
      if (stacking_type == PI_CATION_STACKING)
         ligand_ring_atom_names = ligand_ring_atom_names_in;
      if (stacking_type == CATION_PI_STACKING)
         ligand_cation_atom_name = ligand_cation_atom_name_in;
   }
   int get_stacking_type() const {
      return stacking_type;
   }
   std::string get_ligand_cation_atom_name() const {
      return ligand_cation_atom_name;
   }
   bool se_diff_set() const {
      return se_diff_set_;
   }
   std::pair<double, double> solvent_exposures() const {
      return std::pair<double, double> (se_holo, se_apo);
   }
   void set_water_dist_to_protein(double d) {
      water_dist_to_protein = d;
   }
   std::vector<std::string> get_ligand_ring_atom_names() const {
      return ligand_ring_atom_names;
   }
   bool has_ring_stacking_interaction() const {
      return (ligand_ring_atom_names.size() > 0);
   }
   bool is_a_primary_residue() const {
      return is_a_primary_residue_;
   }

   // std::vector<std::pair<lig_build::pos_t, double> >
   // get_attachment_points(const lig_build::molecule_t &mol) const;

   // friend std::ostream& operator<<(std::ostream &s, residue_circle_t r);
};
// std::ostream& operator<<(std::ostream &s, residue_circle_t r);


#endif // COOT_PLI_RESIDUE_CIRCLE_HH
