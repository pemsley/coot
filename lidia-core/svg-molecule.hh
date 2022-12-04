#ifndef SVG_MOLECULE_HH
#define SVG_MOLECULE_HH

#include "lig-build.hh"
#include "use-rdkit.hh"

class svg_atom_t : public lig_build::atom_t {
   // use self element to set the colour
   void set_colour();
   std::string colour;
public:
   svg_atom_t(const lig_build::pos_t &pos_in, const std::string &ele_in, int formal_charge_in) :
      lig_build::atom_t(pos_in, ele_in, formal_charge_in) { set_colour(); }

   std::string make_text_item(const lig_build::atom_id_info_t &atom_id_info_in,
                              const lig_build::pos_t &centre,
                              double scale, double median_bond_length) const;
};

class svg_bond_t : public lig_build::bond_t {
   std::string draw_sheared_or_darted_wedge_bond(const lig_build::pos_t &pos_1,
                                                 const lig_build::pos_t &pos_2,
                                                 const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
                                                 const lig_build::pos_t &centre,
                                                 double scale) const;
   std::string make_bond_line_string(const lig_build::pos_t &p1, const lig_build::pos_t &p2) const;
public:
   svg_bond_t(int first, int second, lig_build::bond_t::bond_type_t type) :
      lig_build::bond_t(first, second, type) {}

   std::string draw_bond(const svg_atom_t &at_1, const svg_atom_t &at_2,
                  bool at_1_in_ring_flag, bool at_2_in_ring_flag,
                  lig_build::bond_t::bond_type_t bt,
                  bool shorten_first, bool shorten_second,
                  const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_first_atom,
                  const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
                  const lig_build::pos_t &centre, double scale);
   std::string draw_double_in_ring_bond(const lig_build::pos_t &pos_1,
                                        const lig_build::pos_t &pos_2,
                                        bool shorten_first,
                                        bool shorten_second,
                                        const lig_build::pos_t &centre,
                                        double scale, bool dashed_inner=false);
   std::string draw_double_bond(const lig_build::atom_t &at_1,
                                const lig_build::atom_t &at_2,
                                bool shorten_first, bool shorten_second,
                                const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_first_atom,
                                const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
                                const lig_build::pos_t &centre, double scale);
   
};

class svg_molecule_t : public lig_build::molecule_t<svg_atom_t, svg_bond_t> {
public:
   svg_molecule_t() { median_bond_length_ = 1.0; }
   void import_rdkit_mol(RDKit::ROMol *mol, int iconf);
   std::string render_to_svg_string();
   double median_bond_length_;
   double get_scale() const;

   static lig_build::pos_t mol_coords_to_svg_coords(const lig_build::pos_t &pos_1,
                                                    const lig_build::pos_t &centre,
                                                    double scale);
};


#endif // SVG_MOLECULE_HH
