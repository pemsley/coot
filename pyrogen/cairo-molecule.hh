
#include <iostream> // for lig-build.hh - hmm.
#include <cairo/cairo.h>
#include "lidia-core/lig-build.hh"
#include "utils/coot-utils.hh"

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#include "lidia-core/rdkit-interface.hh"
#endif // MAKE_ENHANCED_LIGAND_TOOLS

namespace coot {

   class cairo_atom_t : public lig_build::atom_t {
   
   public:
      cairo_atom_t(lig_build::pos_t pos_in, std::string ele_in, int formal_charge_in) :
	 lig_build::atom_t(pos_in, ele_in, formal_charge_in) {
	 font_colour = "";
      }
      std::string font_colour;

      // needs to pass a colour attribute too
      void make_text_item(cairo_t *cr, const lig_build::atom_id_info_t &atom_id_info_in,
			  const lig_build::pos_t &centre, double scale) const;
      void set_colour(cairo_t *cr) const; // depending on ele

   };

   class cairo_bond_t : public lig_build::bond_t {
   public:
      cairo_bond_t(int first, int second, lig_build::bond_t::bond_type_t type) :
	 lig_build::bond_t(first, second, type) {}
      //
      // maybe we want yet more sophisticated bond_t constructors? (like wmolecule.hh)

      void draw_bond(cairo_t *cr, const lig_build::pos_t &pos_1, const lig_build::pos_t &pos_2,
		     bool shorten_first, bool shorten_second, lig_build::bond_t::bond_type_t bt,
		     const lig_build::pos_t &centre, double scale);
      void draw_double_in_ring_bond(cairo_t *cr,
				    const lig_build::pos_t &pos_1,
				    const lig_build::pos_t &pos_2,
				    bool shorten_first,
				    bool shorten_second,
				    const lig_build::pos_t &centre,
				    double scale, bool dashed_inner=false);
      void draw_double_bond(cairo_t *cr,
			    const lig_build::pos_t &pos_1, const lig_build::pos_t &pos_2,
			    const lig_build::pos_t &centre, double scale);
   };

   class cairo_molecule_t : public lig_build::molecule_t<cairo_atom_t,
							 cairo_bond_t> {

      std::pair<bool, double> scale_correction;
      void draw_bonds();

   public:
      cairo_molecule_t() {}
      ~cairo_molecule_t() {}
#ifdef MAKE_ENHANCED_LIGAND_TOOLS
      // the passed mol can't be const because we use getConformer() which is not const.
      cairo_molecule_t(RDKit::ROMol *mol, int iconf) { import_rdkit_mol(mol, iconf); }
      void import_rdkit_mol(RDKit::ROMol *mol, int iconf);
#endif // MAKE_ENHANCED_LIGAND_TOOLS

      void render(const std::string &png_file_name, unsigned int npx=300);

      // helper function
      static
      lig_build::pos_t mol_coords_to_cairo_coords(const lig_build::pos_t &pos_1,
						  const lig_build::pos_t &centre,
						  double scale);

   };

   void cairo_png_depict(const std::string &mmcif_file_name,
			 const std::string &comp_id,
			 const std::string png_file_name,
			 unsigned int npx=300);

}
