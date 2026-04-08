#ifndef COOT_SRC_TRANSLATION_GIZMO_HH
#define COOT_SRC_TRANSLATION_GIZMO_HH

#include "coot-utils/simple-mesh.hh"
#include "coords/Cartesian.hh"

class translation_gizmo_t {
   void init(); // just make the mesh - don't do OpenGL things
public:
   translation_gizmo_t() { init(); }
   coot::simple_mesh_t mesh;
   float scale_factor; // used for picking
   coot::Cartesian gizmo_origin;
   enum pick_info_t {NONE, X_AXIS, Y_AXIS, Z_AXIS };
   pick_info_t pick(const coot::Cartesian &pt_front, const coot::Cartesian &pt_back);
   void translate(const coot::Cartesian &t); // move the mesh (relative)
   void set_position(const coot::Cartesian &p);
   void scale(float scale_factor); // scale the mesh. Zero or negative scale_factor is ignored.
   void set_scale_absolute(float scale_factor); // scale the mesh. Zero or negative scale_factor is ignored.
   void reset_scale();
   // either one or the other. The gizmo should not be displayed if the object
   // (or molecule) is not displayed. -1 means "not attached"
   enum { UNATTACHED = -1 };
   int attached_to_generic_display_object_number;
   int attached_to_molecule_number;
};


#endif // COOT_SRC_TRANSLATION_GIZMO_HH

