
#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

// #include <GL/glu.h> OpenGLv1

#include "old-generic-display-object.hh"
#include "graphics-info.h"
#include "c-interface-widgets.hh"

// ---------------------- generic objects -----------------------------

void
coot::old_generic_display_object_t::add_line(const coot::colour_holder &colour_in,
                                             const std::string &colour_name,
                                             const int &width_in, 
                                             const std::pair<clipper::Coord_orth, clipper::Coord_orth> &coords_in) {

   int lines_set_index = -1; // magic unset value
   
   for (unsigned int ils=0; ils<lines_set.size(); ils++) {
      if (lines_set[ils].colour_name == colour_name) {
	 if (lines_set[ils].width == width_in) {
	    lines_set_index = ils;
	    break;
	 }
      }
   }

   if (lines_set_index == -1) {
      old_generic_display_line_set_t t(colour_in, colour_name, width_in);
      lines_set.push_back(t);
      lines_set_index = lines_set.size() -1;
   }

   old_generic_display_line_t line(coords_in);
   lines_set[lines_set_index].add_line(line);

}


void coot::old_generic_display_object_t::add_point(const coot::colour_holder &colour_in,
                                                   const std::string &colour_name,
                                                   const int &size_in, 
                                                   const clipper::Coord_orth &coords_in) {

   int points_set_index = -1; // magic unset number
   for (unsigned int ips=0; ips<points_set.size(); ips++) {
      if (points_set[ips].colour_name == colour_name) {
	 if (points_set[ips].size == size_in) {
	    points_set_index = ips;
	    break;
	 }
      }
   }
   if (points_set_index == -1) {
      coot::old_generic_display_point_set_t point_set(colour_in, colour_name, size_in);
      point_set.add_point(coords_in);
      points_set.push_back(point_set);
   } else {
      // normal case
      points_set[points_set_index].add_point(coords_in);
   }
   
}

void
coot::old_generic_display_object_t::add_dodecahedron(const colour_holder &colour_in,
						 const std::string &colour_name,
						 double radius,
						 const clipper::Coord_orth &pos) {

   dodec d;
   dodec_t dod(d, radius, pos);
   dod.col = colour_in;
   dodecs.push_back(dod);
}

void
coot::old_generic_display_object_t::add_pentakis_dodecahedron(const colour_holder &colour_in,
							  const std::string &colour_name,
							  double stellation_factor,
							  double radius,
							  const clipper::Coord_orth &pos) {

   pentakis_dodec d(stellation_factor);
   pentakis_dodec_t pdod(d, radius, pos);
   pdod.col = colour_in;

   pentakis_dodecs.push_back(pdod);
}


// static
void
graphics_info_t::draw_generic_objects() {
   graphics_info_t g;
   g.draw_generic_objects_solid();
}



// static
void
graphics_info_t::draw_generic_objects_simple() {

   // std::cout << "debug:: drawing " << generic_objects_p->size()
   // << " generic objects" << std::endl;

   unsigned int n_points = 0;
   for (unsigned int i=0; i<generic_objects_p->size(); i++) {

      if ((*generic_objects_p)[i].is_displayed_flag) {

	 // if this is attached to a molecule that is not displayed, skip it.
	 if (generic_objects_p->at(i).is_valid_imol()) { // i.e. is not UNDEFINED
	    int imol = generic_objects_p->at(i).get_imol();
	    if (is_valid_model_molecule(imol))
	       if (! graphics_info_t::molecules[imol].is_displayed_p()) {
		  continue;
	       }
	 } else {
	    if (generic_objects_p->at(i).is_intermediate_atoms_object()) {
	       if (! moving_atoms_asc)
		  continue;
	       if (! moving_atoms_asc->mol)
		  continue;
	    }
	 }

	 // Lines
	 for (unsigned int ils=0; ils< (*generic_objects_p)[i].lines_set.size(); ils++) {
	    glLineWidth((*generic_objects_p)[i].lines_set[ils].width);
	    glColor3f((*generic_objects_p)[i].lines_set[ils].colour.red,
		      (*generic_objects_p)[i].lines_set[ils].colour.green,
		      (*generic_objects_p)[i].lines_set[ils].colour.blue);
	    glBegin(GL_LINES);
	    unsigned int s = (*generic_objects_p)[i].lines_set[ils].lines.size();
	    for (unsigned int iline=0; iline<s; iline++) {
	       glVertex3f((*generic_objects_p)[i].lines_set[ils].lines[iline].coords.first.x(),
			  (*generic_objects_p)[i].lines_set[ils].lines[iline].coords.first.y(),
			  (*generic_objects_p)[i].lines_set[ils].lines[iline].coords.first.z());
	       glVertex3f((*generic_objects_p)[i].lines_set[ils].lines[iline].coords.second.x(),
			  (*generic_objects_p)[i].lines_set[ils].lines[iline].coords.second.y(),
			  (*generic_objects_p)[i].lines_set[ils].lines[iline].coords.second.z());
	    }
	    glEnd();
	 }

	 // Points
	 for (unsigned int ips=0; ips<(*generic_objects_p)[i].points_set.size(); ips++) {
	    glPointSize((*generic_objects_p)[i].points_set[ips].size);
	    glColor3f((*generic_objects_p)[i].points_set[ips].colour.red,
		      (*generic_objects_p)[i].points_set[ips].colour.green,
		      (*generic_objects_p)[i].points_set[ips].colour.blue);
	    glBegin(GL_POINTS);
	    unsigned int npoints = (*generic_objects_p)[i].points_set[ips].points.size();
	    n_points += npoints;
	    for (unsigned int ipoint=0; ipoint<npoints; ipoint++) { 
	       glVertex3f((*generic_objects_p)[i].points_set[ips].points[ipoint].x(),
			  (*generic_objects_p)[i].points_set[ips].points[ipoint].y(),
			  (*generic_objects_p)[i].points_set[ips].points[ipoint].z());
	    }
	    glEnd();
	 }

         // Display lists
	 for (unsigned int idl=0; idl<(*generic_objects_p)[i].GL_display_list_handles.size(); idl++) {
             glCallList((*generic_objects_p)[i].GL_display_list_handles[idl]);
         }
      }
   }
   // VRCoot info
   // std::cout << "drew " << n_points << " points" << std::endl; -> ~2000
}

void
graphics_info_t::draw_generic_objects_solid() {

   if (generic_objects_p->size()) {
      for (unsigned int i=0; i<generic_objects_p->size(); i++) {
         const coot::old_generic_display_object_t &obj = generic_objects_p->at(i);
	 if (obj.is_displayed_flag) {
            // std::cout << "draw_generic_objects_solid() " << i << std::endl;
         }
      }
   }
}

#if 0
// static
void
graphics_info_t::draw_generic_objects_solid_old() {

   graphics_info_t g;
   double radius = 0.02;

   // Don't mess with the lighting if we aren't drawing anything
   // 
   if (generic_objects_p->size()) {

      glEnable(GL_LIGHTING);
      glEnable(GL_LIGHT1);
      glEnable(GL_LIGHT0);
      glDisable(GL_LIGHT2); // only for cut-glass mode
      glDisable(GL_COLOR_MATERIAL);
      glEnable(GL_NORMALIZE); // slows things, but makes the shiny nice

      for (unsigned int i=0; i<generic_objects_p->size(); i++) {

	 if ((*generic_objects_p)[i].is_displayed_flag) {


	    // if this is attached to a molecule that is not displayed, skip it.
	    if ((*generic_objects_p)[i].is_valid_imol()) { // i.e. is not UNDEFINED
	       int imol = (*generic_objects_p)[i].get_imol();
	       if (is_valid_model_molecule(imol))
		  if (! graphics_info_t::molecules[imol].is_displayed_p()) {
		     continue;
		  } 
	    } 
	    

	    // Previously (r4209) I had noted that
	    // glEnable(GL_COLOR_MATERIAL) needed for correct tube
	    // colours.
	    // 
	    // 20120903 but this creates problems when displaying
	    // chemical features (they become solid yellow (with nvidia
	    // drivers?)).  I currently don't understand what problems
	    // result if we do not enable GL_COLOR_MATERIAL here. (Mogul
	    // output markup seems OK?).  So comment out.
	    // 
	    // glEnable(GL_COLOR_MATERIAL);
	 
	    GLfloat  mat_specular[]  = {0.9, 0.9, 0.9, 1};
	    GLfloat  mat_shininess[] = {80};
	 
	    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  mat_specular);
	    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
	    
	    // Lines
	    for (unsigned int ils=0; ils< (*generic_objects_p)[i].lines_set.size(); ils++) {
	    
	       glLineWidth((*generic_objects_p)[i].lines_set[ils].width);
	       // 	    glColor3f((*generic_objects_p)[i].lines_set[ils].colour.red,
	       // 		      (*generic_objects_p)[i].lines_set[ils].colour.green,
	       // 		      (*generic_objects_p)[i].lines_set[ils].colour.blue);

	       GLfloat  mat_diffuse[]  = {float((*generic_objects_p)[i].lines_set[ils].colour.red   * 0.8),
					  float((*generic_objects_p)[i].lines_set[ils].colour.green * 0.8),
					  float((*generic_objects_p)[i].lines_set[ils].colour.blue  * 0.8),
					  1.0};
	       glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,   mat_diffuse);
	       glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   mat_diffuse);
	    
	       unsigned int s = (*generic_objects_p)[i].lines_set[ils].lines.size();
	       for (unsigned int iline=0; iline<s; iline++) {

		  g.graphics_object_internal_single_tube((*generic_objects_p)[i].lines_set[ils].lines[iline].coords.first,
							 (*generic_objects_p)[i].lines_set[ils].lines[iline].coords.second,
							 (*generic_objects_p)[i].lines_set[ils].width * radius,
							 coot::ROUND_ENDS);
	       }
	    }
	 
	    // Points
	    for (unsigned int ips=0; ips<(*generic_objects_p)[i].points_set.size(); ips++) {
	       // 	    glColor3f((*generic_objects_p)[i].points_set[ips].colour.red,
	       // 		      (*generic_objects_p)[i].points_set[ips].colour.green,
	       // 		      (*generic_objects_p)[i].points_set[ips].colour.blue);
	    
	       int sphere_slices = 5;
	       int sphere_stacks = 5;
	       float feature_opacity = 1.0;
	       unsigned int npoints = (*generic_objects_p)[i].points_set[ips].points.size();
	       for (unsigned int ipoint=0; ipoint<npoints; ipoint++) {

		  const coot::old_generic_display_object_t &obj = (*generic_objects_p)[i];
		  GLfloat  mat_specular[]  = {obj.points_set[ips].colour.red,
					      obj.points_set[ips].colour.green,
					      obj.points_set[ips].colour.blue,
					      feature_opacity};
		  GLfloat  mat_shininess[] = {15};
		  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  mat_specular);
		  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
		  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,   mat_specular);
		  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   mat_specular);

	       
		  GLUquadric* sphere_quad = gluNewQuadric();
		  glPushMatrix();
		  glTranslatef((*generic_objects_p)[i].points_set[ips].points[ipoint].x(),
			       (*generic_objects_p)[i].points_set[ips].points[ipoint].y(),
			       (*generic_objects_p)[i].points_set[ips].points[ipoint].z());	 
		  gluSphere(sphere_quad,
			    (*generic_objects_p)[i].points_set[ips].size * radius,
			    sphere_slices, sphere_stacks);
		  gluDeleteQuadric(sphere_quad);
		  glPopMatrix();	 
	       }
	    }

	    // Other stuff:
	    float feature_opacity = 0.6;
	    // feature_opacity = 1.0; // hack for to fix shininess

	    // spheres

	    if ((*generic_objects_p)[i].spheres.size()) {
	       glEnable (GL_BLEND); // these 2 lines are needed to make the transparency work.
	       glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	       for (unsigned int isphere=0; isphere<(*generic_objects_p)[i].spheres.size(); isphere++) { 

		  const coot::old_generic_display_object_t &obj = (*generic_objects_p)[i];
		  GLfloat  mat_specular[]  = {obj.spheres[isphere].col.col[0],
					      obj.spheres[isphere].col.col[1],
					      obj.spheres[isphere].col.col[2], 
					      feature_opacity};
		  GLfloat  mat_shininess[] = {15};
		  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  mat_specular);
		  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
		  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,   mat_specular);
		  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   mat_specular);
	       
		  int sphere_slices = 10;
		  int sphere_stacks = 10;
		  GLUquadric* sphere_quad = gluNewQuadric();
		  glPushMatrix();
		  glTranslatef((*generic_objects_p)[i].spheres[isphere].centre.x(),
			       (*generic_objects_p)[i].spheres[isphere].centre.y(),
			       (*generic_objects_p)[i].spheres[isphere].centre.z());
		  gluSphere(sphere_quad,
			    (*generic_objects_p)[i].spheres[isphere].radius,
			    sphere_slices, sphere_stacks);
		  gluDeleteQuadric(sphere_quad);
		  glPopMatrix();
	       }
	    }

	    // arrows
	    if ((*generic_objects_p)[i].arrows.size()) {
	       glEnable (GL_BLEND);
	       glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	       for (unsigned int iarrow=0; iarrow<(*generic_objects_p)[i].arrows.size(); iarrow++) {
		  const coot::old_generic_display_object_t &obj = (*generic_objects_p)[i];
		  GLfloat  mat_specular[]  = {obj.arrows[iarrow].col.col[0],
					      obj.arrows[iarrow].col.col[1],
					      obj.arrows[iarrow].col.col[2], 
					      feature_opacity};
		  GLfloat  mat_shininess[] = {15};
		  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  mat_specular);
		  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
		  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,   mat_specular);
		  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   mat_specular);
		  g.graphics_object_internal_arrow((*generic_objects_p)[i].arrows[iarrow].start_point,
						   (*generic_objects_p)[i].arrows[iarrow].end_point,
						   0.3, 0.1); 
	       }
	    }

	    // tori
	    if ((*generic_objects_p)[i].tori.size()) {
	       glEnable (GL_BLEND);
	       glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	       for (unsigned int itor=0; itor<(*generic_objects_p)[i].tori.size(); itor++) {
		  const coot::old_generic_display_object_t &obj = (*generic_objects_p)[i];
		  GLfloat  mat_diffuse[]  = {obj.tori[itor].col.col[0],
					     obj.tori[itor].col.col[1],
					     obj.tori[itor].col.col[2], 
					     feature_opacity};
		  GLfloat  mat_specular[]  = {0.1, 0.1, 0.1, 1.0};
		  GLfloat  mat_shininess[] = {1};
		  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  mat_specular);
		  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
		  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,   mat_diffuse);
		  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   mat_diffuse);

		  g.graphics_object_internal_torus(obj.tori[itor].start_point,
						   obj.tori[itor].end_point,
						   obj.tori[itor].radius_1,
						   obj.tori[itor].radius_2,
						   obj.tori[itor].n_ring_atoms);
	       }
	    }

	    // arcs
	    if ((*generic_objects_p)[i].arcs.size()) {
	    
	       for (unsigned int iarc=0; iarc<(*generic_objects_p)[i].arcs.size(); iarc++) {
		  const coot::old_generic_display_object_t &obj = (*generic_objects_p)[i];

		  // glEnable(GL_COLOR_MATERIAL);
		  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	       
		  GLfloat  mat_diffuse[]  = {float(obj.arcs[iarc].col.col[0] * 0.8),
					     float(obj.arcs[iarc].col.col[1] * 0.8),
					     float(obj.arcs[iarc].col.col[2] * 0.8),
					     1.0};
		  // 	       GLfloat  mat_specular[]  = {obj.arcs[iarc].col.col[0],
		  // 					   obj.arcs[iarc].col.col[1],
		  // 					   obj.arcs[iarc].col.col[2], 
		  // 					   1.0};
		  GLfloat  mat_specular[]  = {0.6, 0.6, 0.6, 1};
		  GLfloat  mat_shininess[] = {35};
		  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  mat_specular);
		  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
		  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,   mat_diffuse);
		  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   mat_diffuse);

		  g.graphics_object_internal_arc(obj.arcs[iarc].start_angle,
						 obj.arcs[iarc].end_angle,
						 obj.arcs[iarc].start_point,
						 obj.arcs[iarc].start_dir,
						 obj.arcs[iarc].normal,
						 obj.arcs[iarc].radius,
						 obj.arcs[iarc].radius_inner);
	       }
	    }


	    // dodecahdrons
	    //
	    if ((*generic_objects_p)[i].dodecs.size()) {
	       float feature_opacity = 0.7;
	       glEnable (GL_BLEND);
	       glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	       for (unsigned int idodec=0; idodec<(*generic_objects_p)[i].dodecs.size(); idodec++) {
		  const coot::old_generic_display_object_t &obj = (*generic_objects_p)[i];
		  GLfloat  mat_diffuse[]  = {obj.dodecs[idodec].col.red,
					     obj.dodecs[idodec].col.green,
					     obj.dodecs[idodec].col.blue, 
					     feature_opacity};
		  GLfloat  mat_specular[]  = {0.3, 0.3, 0.3, 1.0};
		  GLfloat  mat_shininess[] = {1};
		  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  mat_specular);
		  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
		  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,   mat_diffuse);
		  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   mat_diffuse);

		  g.graphics_object_internal_dodec(obj.dodecs[idodec]);
	       }
	    }

	    // pentakis dodecahdrons
	    //
	    if ((*generic_objects_p)[i].pentakis_dodecs.size()) {
	       float feature_opacity = 0.93;
	       glEnable (GL_BLEND);
	       glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	       for (unsigned int idodec=0; idodec<(*generic_objects_p)[i].pentakis_dodecs.size(); idodec++) {
		  const coot::old_generic_display_object_t &obj = (*generic_objects_p)[i];
		  GLfloat  mat_diffuse[]  = {obj.pentakis_dodecs[idodec].col.red,
					     obj.pentakis_dodecs[idodec].col.green,
					     obj.pentakis_dodecs[idodec].col.blue, 
					     feature_opacity};
		  GLfloat  mat_specular[]  = {0.5, 0.5, 0.5, 1.0};
		  GLfloat  mat_shininess[] = {40};
		  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  mat_specular);
		  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
		  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,   mat_diffuse);
		  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   mat_diffuse);

		  g.graphics_object_internal_pentakis_dodec(obj.pentakis_dodecs[idodec]);
	       }
	    }
	 }
      }
      glDisable(GL_LIGHTING);
   }
}
#endif // draw_generic_objects_solid_old



/*! \brief return the number of generic display objects */
int number_of_generic_objects() {

   graphics_info_t g;
   return g.generic_objects_p->size();
}

std::pair<short int, std::string>
is_interesting_dots_object_next_p(const std::vector<std::string> &vs) {

   std::pair<short int, std::string> r(0, "");

   if (vs.size() == 3) {
//       std::cout << "Looking at bits:  \n  "; 
//       for (unsigned int i=0; i<3; i++) { 
// 	 std::cout << ":" << vs[i] << ": ";
//       }
//       std::cout << "\n"; 
      if ((vs[1] == "wide") && (vs[2] == "contact)")) {
	 r.first = 1;
	 r.second = "wide contact";
      }
      if ((vs[1] == "close") && (vs[2] == "contact)")) {
	 r.first = 1;
	 r.second = "close contact";
      }
      if ((vs[1] == "small") && (vs[2] == "overlap)")) {
	 r.first = 1;
	 r.second = "small overlap";
      }
      if ((vs[1] == "bad") && (vs[2] == "overlap)")) {
	 r.first = 1;
	 r.second = "bad overlap";
      }
      if (vs[1] == "H-bonds)") { 
	 r.first = 1;
	 r.second = "H-bonds";
      }
   }
   return r;
}

std::string probe_dots_short_contact_name_to_expanded_name(const std::string &short_name) {

   std::vector<std::pair<std::string, std::string> > names;
   names.push_back(std::pair<std::string, std::string>("wc", "wide contact"));
   names.push_back(std::pair<std::string, std::string>("cc", "close contact"));
   names.push_back(std::pair<std::string, std::string>("so", "small overlap"));
   names.push_back(std::pair<std::string, std::string>("bo", "bad overlap"));
   names.push_back(std::pair<std::string, std::string>("hb", "H-bonds"));

   std::string r = "unknown";
   for (int i=0; i<5; i++) {
      if (names[i].first == short_name) {
	 r = names[i].second;
	 break;
      }
   } 
   return r;
} 
