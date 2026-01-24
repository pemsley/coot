
#include "graphics-info.h"

//   static
coot::protein_geometry *
graphics_info_t::Geom_p() { return geom_p; }

// static
bool
graphics_info_t::is_valid_model_molecule(int imol) {

   bool v = 0;
   if (imol >= 0) {
      if (imol < n_molecules()) {
         if (molecules[imol].has_model()) {
            v = 1;
         }
      }
   }
   return v;
}

// static
bool
graphics_info_t::is_valid_map_molecule(int imol) {

   bool v = 0;
   if (imol >= 0) {
      if (imol < n_molecules()) {
	 if (molecules[imol].has_xmap()) {
            v = 1;
	 }
	 // NXMAP-FIXME // do I want to check for nxmap here too?
      }
   }
   return v;
}

// static
bool
graphics_info_t::is_difference_map(int imol) {
   bool status = false;
   if (imol >= 0) {
      if (imol < n_molecules()) {
         status = molecules[imol].is_difference_map_p();
      }
   }
   return status;
}

// static
void
graphics_info_t::add_to_rotation_centre(const glm::vec3 &offset) {
   rotation_centre_x += offset.x;
   rotation_centre_y += offset.y;
   rotation_centre_z += offset.z;
}

// static
clipper::Coord_orth
graphics_info_t::get_rotation_centre_co() {
   return clipper::Coord_orth(rotation_centre_x, rotation_centre_y, rotation_centre_z);
}

// static
glm::vec3
graphics_info_t::get_rotation_centre() {
   return glm::vec3(rotation_centre_x, rotation_centre_y, rotation_centre_z);
}


// static
bool
graphics_info_t::display_mode_use_secondary_p() {

   bool r = false;
   if ((display_mode == coot::SIDE_BY_SIDE_STEREO) ||
       (display_mode == coot::SIDE_BY_SIDE_STEREO_WALL_EYE) ||
       (display_mode == coot::DTI_SIDE_BY_SIDE_STEREO)) {
      r = true;
   }
   return r;
}

// static
void
graphics_info_t::erase_last_molecule() {
     // std::vector<molecule_class_info_t>::iterator it = molecules.end();
     // std << "DEBUG:: Erasing molecule number " << it->MoleculeNumber() << std::endl;
/*      std::cout << "DEBUG:: Erasing the back molecule " << molecules.size() - 1  */
/* 	       << " which says that it has molecule number "  */
/* 	       << molecules[molecules.size() -1].MoleculeNumber() << std::endl; */
     molecules.pop_back();
}

std::string
graphics_info_t::get_directory_for_fileselection() const {
   return directory_for_fileselection;
}

std::string
graphics_info_t::get_directory_for_filechooser() const {
   return directory_for_filechooser;
}

float
graphics_info_t::get_clipping_plane_front() const {
   if (perspective_projection_flag)
      return screen_z_near_perspective;
   else
      return clipping_front;
}

float
graphics_info_t::get_clipping_plane_back() const {
   if (perspective_projection_flag)
      return screen_z_far_perspective;
   else
      return clipping_back;
}

// static
int
graphics_info_t::n_map_molecules() {
   int n = 0;
   for (unsigned int i=0; i<molecules.size(); i++) {
      if (is_valid_map_molecule(i))
         n++;
   }
   return n;
}

short int
graphics_info_t::GetActiveMapDrag() const { return active_map_drag_flag; };

void
graphics_info_t::set_rotation_centre_cross_hairs_colour(const glm::vec4 &c) { rotation_centre_cross_hairs_colour = c; }

// static
void
graphics_info_t::Increment_Frames() {
   Frames++;
}

void
graphics_info_t::save_display_control_widget_in_graphics(GtkWidget *widget) {
   display_control_window_ = widget;
}

GtkWidget *
graphics_info_t::get_display_control_window() {
   return display_control_window_;
}

GtkWidget *
graphics_info_t::display_control_window() {
   return get_display_control_window();
}

void
graphics_info_t::set_use_harmonic_approximations_for_nbcs(bool flag) {
   use_harmonic_approximation_for_NBCs = flag;
}

// static
glm::vec4
graphics_info_t::get_background_colour() {
   return glm::vec4(background_colour, 1.0f);
}

void
graphics_info_t::quanta_buttons() {
   button_1_mask_ = GDK_BUTTON2_MASK;
   button_2_mask_ = GDK_BUTTON1_MASK;
}

// mouse buttons
GdkModifierType
graphics_info_t::gdk_button1_mask() { return button_1_mask_;}
GdkModifierType
graphics_info_t::gdk_button2_mask() { return button_2_mask_;}
GdkModifierType
graphics_info_t::gdk_button3_mask() { return button_3_mask_;}

// static
coot::Cartesian
graphics_info_t::to_cartesian(const clipper::Coord_orth &co) {
   return coot::Cartesian(co.x(), co.y(), co.z());
}

// static
clipper::Coord_orth
graphics_info_t::to_coord_orth(const coot::Cartesian &c) {
   return clipper::Coord_orth(c.x(), c.y(), c.z());
}

void
graphics_info_t::set_find_ligands_mols(int map, int protein,
                                       const std::vector<std::pair<int, bool> > &ligand_wiggly_info) {
   find_ligand_map_mol_ = map;
   find_ligand_protein_mol_ = protein;
   // *find_ligand_ligand_mols_ = ligand_wiggly_info; // lets add some protection..
   find_ligand_ligand_mols_->clear();
   for (unsigned int ilig=0; ilig<ligand_wiggly_info.size(); ilig++) {
      int il=ligand_wiggly_info[ilig].first;
      if (il < n_molecules()) {
         if (molecules[il].atom_sel.n_selected_atoms > 0) {
	    find_ligand_ligand_mols_->push_back(ligand_wiggly_info[ilig]);
         }
      }
   }
}

int
graphics_info_t::find_ligand_map_mol() const {
   return find_ligand_map_mol_;
}

int
graphics_info_t::find_ligand_protein_mol() const {
   return find_ligand_protein_mol_;
}

// static
void
graphics_info_t::set_ligand_protein_mol(int imol) {
   if (imol >=0)
      if (imol < n_molecules())
	 find_ligand_protein_mol_ = imol;
}
// static
void
graphics_info_t::set_ligand_map_mol(int imol) {
   if (imol >=0)
      if (imol < n_molecules())
         find_ligand_map_mol_ = imol;
}

// static
void
graphics_info_t::find_ligand_add_rigid_ligand(int imol) {
   if (imol >=0) {
      if (imol < n_molecules()) {
	 if (molecules[imol].has_model()) {
            find_ligand_ligand_mols_->push_back(std::pair<int, bool>(imol, 0));
	 }
      }
   }
}
// static
void
graphics_info_t::find_ligand_add_flexible_ligand(int imol) {
   if (imol >=0) {
      if (imol < n_molecules()) {
	 if (molecules[imol].has_model()) {
            find_ligand_ligand_mols_->push_back(std::pair<int, bool>(imol, 1));
	 }
      }
   }
}

void
graphics_info_t::set_find_ligand_do_real_space_refine_state(bool state) {
   find_ligand_do_real_space_refine_ = state;
}

bool
graphics_info_t::find_ligand_do_real_space_refine_state() {
   return find_ligand_do_real_space_refine_;
}

std::vector<std::pair<int, bool> >
graphics_info_t::find_ligand_ligand_mols() const {
   return *find_ligand_ligand_mols_;
}


void
graphics_info_t::find_ligand_clear_ligand_mols() {
   find_ligand_ligand_mols_->clear();
}

int
graphics_info_t::last_restraints_size() const {
   // It's OK to call this when there are no restraints - e.g. we move by rotate/translate
   // rather than during a refinement.
   if (! last_restraints) {
      return 0;
   } else {
      return last_restraints->size();
   }
}

// static
void
graphics_info_t::fill_rotamer_probability_tables() {

   if (! rot_prob_tables.tried_and_failed()) {
      rot_prob_tables.fill_tables();
   }
}

// static
void
graphics_info_t::set_model_display_radius(bool on_off, float radius_in) {
   model_display_radius.first  = on_off;
   model_display_radius.second = radius_in;
}


void
graphics_info_t::set_convert_dictionary_planes_to_improper_dihedrals(bool state) {
   convert_dictionary_planes_to_improper_dihedrals_flag = state;
}

// static
bool
graphics_info_t::have_user_defined_colours() {
   return ! user_defined_colours.empty();
}

void
graphics_info_t::set_python_draw_function(const std::string &f) {
   python_draw_function_string = f;
}

// static
void graphics_info_t::read_inchikeys() {
   inchikey_store.init();
}

int graphics_info_t::intelligent_get_scroll_wheel_map() const {

   int swm = -1;
   if (is_valid_map_molecule(scroll_wheel_map))
      if (molecules[scroll_wheel_map].get_map_is_displayed())
         swm = scroll_wheel_map;
   if (swm == -1) {
      // find the first visible map then
      int molecules_size = molecules.size();
      for(int imol=0; imol<molecules_size; imol++) {
         if (is_valid_map_molecule(imol)) {
            if (molecules[imol].get_map_is_displayed()) {
               swm = imol;
               break;
            }
         }
      }
   }
   return swm;
}
