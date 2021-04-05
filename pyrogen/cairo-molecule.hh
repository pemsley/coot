#ifndef CAIRO_MOLECULE_HH
#define CAIRO_MOLECULE_HH

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#include "Python.h"
#endif

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
			  const lig_build::pos_t &centre,
			  double scale, double median_bond_length) const;
      void set_colour(cairo_t *cr) const; // depending on ele

   };

   class cairo_bond_t : public lig_build::bond_t {
      void draw_sheared_or_darted_wedge_bond(cairo_t *cr,
					     const lig_build::pos_t &pos_1,
					     const lig_build::pos_t &pos_2,
					     const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
					     const lig_build::pos_t &centre,
					     double scale) const;
   public:
      cairo_bond_t(int first, int second, lig_build::bond_t::bond_type_t type) :
	 lig_build::bond_t(first, second, type) {}
      //
      // maybe we want yet more sophisticated bond_t constructors? (like wmolecule.hh)

      void draw_bond(cairo_t *cr, const cairo_atom_t &at_1, const cairo_atom_t &at_2,
		     bool at_1_in_ring_flag, bool at_2_in_ring_flag,
		     lig_build::bond_t::bond_type_t bt,
		     bool shorten_first, bool shorten_second,
		     const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_first_atom,
		     const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
		     const lig_build::pos_t &centre, double scale);
      void draw_double_in_ring_bond(cairo_t *cr,
				    const lig_build::pos_t &pos_1,
				    const lig_build::pos_t &pos_2,
				    bool shorten_first,
				    bool shorten_second,
				    const lig_build::pos_t &centre,
				    double scale, bool dashed_inner=false);
      void draw_double_bond(cairo_t *cr,
			    const lig_build::atom_t &at_1,
			    const lig_build::atom_t &at_2,
			    bool shorten_first, bool shorten_second,
			    const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_first_atom,
			    const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
			    const lig_build::pos_t &centre, double scale);
   };

   class cairo_molecule_t : public lig_build::molecule_t<cairo_atom_t,
							 cairo_bond_t> {

      std::pair<bool, double> scale_correction;
      void draw_bonds();
      void render(cairo_t *cr);
      double median_bond_length_;
      void set_highlight_colour(cairo_t *cr, unsigned int idx);
      double get_scale() const;
      std::vector<unsigned int> find_bonds_for_atoms(const std::vector<unsigned int> &highlight_atom_indices) const;
      void debug_box(cairo_t *cr);

   public:
      cairo_molecule_t() { median_bond_length_ = 1.5; }
      ~cairo_molecule_t() {}
#ifdef MAKE_ENHANCED_LIGAND_TOOLS
      // the passed mol can't be const because we use getConformer() which is not const.
      cairo_molecule_t(RDKit::ROMol *mol, int iconf) { import_rdkit_mol(mol, iconf); }
      void import_rdkit_mol(RDKit::ROMol *mol, int iconf);
#endif // MAKE_ENHANCED_LIGAND_TOOLS

      void draw_atom_highlights(cairo_t *cr, const lig_build::pos_t &centre,
				double scale,
				const std::vector<unsigned int> &highlight_atom_indices,
				const std::vector<unsigned int> &highlight_bond_indices,
				bool use_highlight_bond_indices_flag);
      // render to file
      void render_to_file(const std::string &png_file_name, unsigned int npx=300,
			  const std::pair<bool, colour_holder> &bg_col=std::pair<bool, colour_holder>(false, colour_holder()));
      // render to string
      std::string render_to_png_string(const std::vector<unsigned int> &atom_highlight_list,
				       const std::vector<unsigned int> &bond_highlight_list,
				       bool use_highlight_bond_indices_flag,
				       unsigned int npx=300);

      // render to string
      std::string render_to_svg_string(const std::vector<unsigned int> &atom_highlight_list,
				       const std::vector<unsigned int> &bond_highlight_list,
				       bool use_highlight_bond_indices_flag,
				       unsigned int npx=300);

      // helper function
      static
      lig_build::pos_t mol_coords_to_cairo_coords(const lig_build::pos_t &pos_1,
						  const lig_build::pos_t &centre,
						  double scale);

      static cairo_status_t png_stream_writer(void *closure_in,
					      const unsigned char *data,
					      unsigned int length);
   };

   void cairo_png_depict_from_mmcif(const std::string &mmcif_file_name,
				    const std::string &comp_id,
				    const std::string png_file_name,
				    unsigned int npx=300,
				    PyObject *background_colour=0);

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
   std::string cairo_png_string_from_mol(RDKit::ROMol *m, int iconf = -1,
					 PyObject *highlight_atom_list=0,
					 PyObject *highlight_bond_list=0,
					 PyObject *highlight_atom_colours_dict=0,
					 PyObject *highlight_bond_colours_dict=0,
					 unsigned int npx=300);
   std::string cairo_svg_string_from_mol(RDKit::ROMol *m, int iconf = -1,
					 PyObject *highlight_atom_list=0,
					 PyObject *highlight_bond_list=0,
					 PyObject *highlight_atom_colours_dict=0,
					 PyObject *highlight_bond_colours_dict=0,
					 unsigned int npx=300);
   // both of the above are trivial wrappers for
   std::string cairo_image_string_from_mol(RDKit::ROMol *m, int iconf = -1,
					   PyObject *highlight_atom_list=0,
					   PyObject *highlight_bond_list=0,
					   PyObject *highlight_atom_colours_dict=0,
					   PyObject *highlight_bond_colours_dict=0,
					   bool png_vs_svg_mode=true,
					   unsigned int npx=300);
#endif // MAKE_ENHANCED_LIGAND_TOOLS

}




#endif // CAIRO_MOLECULE_HH
