/* src/graphics-info-render.cc
 * 
 * Copyright 2004, 2006, 2007 by The University of York.
 * Author Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1335, USA.
 */

#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"

#if defined _MSC_VER
#include <windows.h>
#endif

#include <iostream>
#include <fstream>

#include "clipper/core/coords.h"
#include "clipper/core/rotation.h"

#include "gl-matrix.h"
#include "graphics-info.h"
#include "coot-render.hh"

#include "ccp4mg-utils/ppmutil.h"

#include "old-generic-display-object.hh"


void
coot::raytrace_info_t::add_geometry_objects(const std::vector<coot::simple_distance_object_t> &sdov) {

   int ndist = sdov.size();
   double dist;
   graphics_info_t g;

   if (ndist > 0) {
      meshed_generic_display_object mgdo;
      mgdo.mesh.draw_this_mesh = true;
      mgdo.mesh.set_name("distance-geometry");
      for (int i=0; i<ndist; i++) {
	 if (g.is_valid_model_molecule(sdov[i].imol_start)) {
	    if (g.is_valid_model_molecule(sdov[i].imol_end)) {
	       if (g.molecules[sdov[i].imol_start].is_displayed_p()) {
		  if (g.molecules[sdov[i].imol_end].is_displayed_p()) {

                     clipper::Coord_orth p_start(sdov[i].start_pos.x(),
                                                 sdov[i].start_pos.y(),
                                                 sdov[i].start_pos.z());
                     clipper::Coord_orth p_end(sdov[i].end_pos.x(),
                                               sdov[i].end_pos.y(),
                                               sdov[i].end_pos.z());

                     float dash_density = 5.0;
                     float ll = clipper::Coord_orth::length(p_start, p_end);
                     int n_dashes = int(dash_density * ll);
                     bool visible = true;

                     for (int idash=0; idash<(n_dashes-1); idash++) {
                        if (visible) {
                           float fracs = float(idash)/float(n_dashes);
                           float fracn = float(idash+1)/float(n_dashes);
                           clipper::Coord_orth p1 = p_start + fracs * (p_end - p_start);
                           clipper::Coord_orth p2 = p_start + fracn * (p_end - p_start);
                           std::pair<clipper::Coord_orth, clipper::Coord_orth> from_to(p1, p2);
                           colour_holder col(0.4, 0.65, 0.4);
                           std::string colour_name = "lightblue";
                           int width = 4;
                           mgdo.add_line(col, colour_name, width, from_to);
                           display_objects.push_back(mgdo);
                        }
                        visible = !visible;
                     }
		  }
	       }
	    }
	 }
      }

      for (int i=0; i<ndist; i++) {
	 if (g.is_valid_model_molecule(sdov[i].imol_start)) {
	    if (g.is_valid_model_molecule(sdov[i].imol_end)) {
	       if (g.molecules[sdov[i].imol_start].is_displayed_p()) {
		  if (g.molecules[sdov[i].imol_end].is_displayed_p()) {
                     clipper::Coord_orth text_pos = sdov[i].start_pos +
			0.5 * (sdov[i].end_pos - sdov[i].start_pos +
                               clipper::Coord_orth(0.0, 0.1, 0.1));
		     dist = clipper::Coord_orth::length(sdov[i].start_pos, sdov[i].end_pos);
                     std::string s = coot::util::float_to_string(dist);
                     std::pair<std::string, clipper::Coord_orth> p(s, text_pos);
                     add_label(p); // This is drawn in atom colour - which is not ideal.
                                   // Previous to the addition of this function,
                                   // the only thing with labels were atoms.
                                   // Fix another time.
		  }
	       }
	    }
	 }
      }
   }
}


// raster3d
short int
graphics_info_t::raster3d(std::string filename) {

   short int istate = 0;
   coot::colour_t background;
   background.col.resize(3);
   background.col[0] = background_colour[0];
   background.col[1] = background_colour[1];
   background.col[2] = background_colour[2];
   int width  = 600;
   int height = 600;
   if (glareas.size() > 0) {
      GtkAllocation allocation = get_glarea_allocation();
      width  = allocation.width;
      height = allocation.height;
   } 
   coot::raytrace_info_t rt(RotationCentre(), zoom, background,
                            width, height,
                            clipping_front,
                            raster3d_bond_thickness,
                            raster3d_bone_thickness,
                            raster3d_atom_radius,
                            raster3d_density_thickness);
   GL_matrix m;
   m.from_quaternion(quat);
   rt.set_view_matrix(m);

   rt.add_generic_display_objects(generic_display_objects);
   rt.set_raster3d_enable_shadows(raster3d_enable_shadows);
   bool is_bb = background_is_black_p();
   coot::colour_t atom_label_colour(font_colour.red, font_colour.green, font_colour.blue);
   rt.set_atom_label_colour(atom_label_colour);
   rt.set_font_size(raster3d_font_size);

   std::cout << "Generating raytrace molecule objects..." << std::endl;
   for (int imol=0; imol<n_molecules(); imol++) {
      std::cout << " molecule " << imol << " in  raytrace" << std::endl;

      if (molecules[imol].is_displayed_p()) {
	 if (molecules[imol].has_model()) {
	    rt.rt_mol_info.push_back(molecules[imol].fill_raster_model_info(is_bb));
	 }
	 coot::ray_trace_molecule_info rai = molecules[imol].fill_raster_additional_info();
	 rt.rt_mol_info.push_back(rai);

	 if (molecules[imol].has_xmap()) {  // NXMAP-FIXME
	    // map and skeleton
	    rt.rt_mol_info.push_back(molecules[imol].fill_raster_map_info(1));
	    if (molecules[imol].is_difference_map_p()) {
	       rt.rt_mol_info.push_back(molecules[imol].fill_raster_map_info(-1));
	    }
	 }
	 const std::vector<int> &l = molecules[imol].labelled_atom_index_list;
	 unsigned int n_atoms_to_label = l.size();
	 for (std::size_t i=0; i<l.size(); i++) {
	    mmdb::Atom *at = molecules[imol].atom_sel.atom_selection[l[i]];
	    // Look these up if anyone ever complains about "incorrect" labels
	    bool al_flag = false; // brief atom labels
	    bool sid_flag = false; // seg_ids flag
	    std::pair<std::string, clipper::Coord_orth> s = molecules[imol].make_atom_label_string(i, al_flag, sid_flag);
	    rt.add_label(s);
	 }
      }
   }
   std::cout << "Rendering raytrace..." << std::endl;
   rt.render_ray_trace(filename);
   std::cout << "done raytrace." << std::endl;
   return istate;
}



// povray
short int
graphics_info_t::povray(std::string filename) {

   short int istate=0;
   coot::colour_t background;
   background.col.resize(3);
   background.col[0] = background_colour[0];
   background.col[1] = background_colour[1];
   background.col[2] = background_colour[2];
   GtkAllocation allocation = get_glarea_allocation();
   coot::raytrace_info_t rt(RotationCentre(), zoom, background,
			    allocation.width,
			    allocation.height, 
			    clipping_front,
			    raster3d_bond_thickness,
			    raster3d_bone_thickness,
			    raster3d_atom_radius,
			    raster3d_density_thickness);
   GL_matrix m;
   m.from_quaternion(quat);
   rt.set_view_matrix(m);

   rt.add_generic_display_objects(generic_display_objects);

   // So where is the "eye"? We have to do an unproject:
   int x0 = allocation.width/2;
   int y0 = allocation.height/2;
   glm::vec3 glm_eye = get_world_space_eye_position();
   coot::Cartesian eye(glm_eye.x, glm_eye.y, glm_eye.z);

   // It seems that for raster3d, this eye position is good, but for
   // povray, we are very close in.  So, let's try moving the eye back
   // a scaled amount

   coot::Cartesian eye_centre = eye - rt.view_centre;
   eye_centre *= 7.5;
   coot::Cartesian far_eye = rt.view_centre - eye_centre;
   //std::cout <<"BL DEBUG:: eye and far eye all   " << eye << " " << far_eye<<std::endl;
   //std::cout <<"BL DEBUG:: view centre all       " << rt.view_centre <<std::endl;
   eye *= 1.015;
   rt.set_front_clipping_plane_point(eye);
   rt.set_camera_location(far_eye);
   bool is_bb = background_is_black_p();
   
   for (int imol=0; imol<n_molecules(); imol++) {
      std::cout << " molecule " << imol << " in  raytrace" << std::endl;
      
      if (molecules[imol].has_model()) {
         if (molecules[imol].is_displayed_p()) { 
            rt.rt_mol_info.push_back(molecules[imol].fill_raster_model_info(is_bb));
            //coot::Cartesian eye_centre = eye - rt.view_centre;
            //eye_centre *= 7.5;
            //coot::Cartesian far_eye = rt.view_centre - eye_centre;
            //std::cout <<"BL DEBUG:: eye and far eye model " << eye << " " << far_eye<<std::endl;
            //std::cout <<"BL DEBUG:: view centre model     " << rt.view_centre <<std::endl;
            //eye *= 1.015;
            //rt.set_front_clipping_plane_point(eye);
            //rt.set_camera_location(far_eye);
         }
      }
      if (molecules[imol].has_xmap()) {  // NXMAP-FIXME.
         rt.rt_mol_info.push_back(molecules[imol].fill_raster_map_info(1));
         if (molecules[imol].is_difference_map_p()) { 
            rt.rt_mol_info.push_back(molecules[imol].fill_raster_map_info(-1));
         }
      }
   }
   std::cout << "INFO Raytracing..." << std::endl;
   rt.povray_ray_trace(filename); // write the povray input file
   std::cout << "INFO:: Wrote raytrace input file." << std::endl;

   return istate;
}


// Renderman format:  (Fun 20090215)
// 
short int
graphics_info_t::renderman(std::string filename) {

   short int istate = 0;

   coot::colour_t background;
   background.col.resize(3);
   background.col[0] = background_colour[0];
   background.col[1] = background_colour[1];
   background.col[2] = background_colour[2];
   int width = 600;
   int height = 600;
   GtkAllocation allocation = get_glarea_allocation();
   width  = allocation.width;
   height = allocation.height;
   bool is_bb = background_is_black_p(); 

   coot::raytrace_info_t rt(RotationCentre(), zoom, background,
                            width, height,
                            clipping_front,
                            raster3d_bond_thickness,
                            raster3d_bone_thickness,
                            raster3d_atom_radius,
                            raster3d_density_thickness);
   GL_matrix m;
   m.from_quaternion(quat);
   rt.set_view_matrix(m);
   rt.set_quaternion(quat);

   float aspect_ratio = 1.3;
   if (use_graphics_interface_flag) { 
      aspect_ratio = float (allocation.width)/float (allocation.height);
   }
   
   rt.set_ortho_params(-0.3*zoom*aspect_ratio,
                       0.3*zoom*aspect_ratio,
                       -0.3*zoom, 0.3*zoom);

   rt.add_generic_display_objects(generic_display_objects);

   std::cout << "Generating raytrace molecule objects..." << std::endl;
   for (int imol=0; imol<n_molecules(); imol++) {
      std::cout << " molecule " << imol << " in  raytrace" << std::endl;
      
      if (molecules[imol].has_model()) {
         if (molecules[imol].is_displayed_p()) { 
            rt.rt_mol_info.push_back(molecules[imol].fill_raster_model_info(is_bb));
         }
      }
      if (molecules[imol].has_xmap()) { // NXMAP-FIXME
         rt.rt_mol_info.push_back(molecules[imol].fill_raster_map_info(1));
         if (molecules[imol].is_difference_map_p()) { 
            rt.rt_mol_info.push_back(molecules[imol].fill_raster_map_info(-1));
         }
      }
   }
   std::cout << "Rendering with renderman output..." << std::endl;
   rt.renderman_render(filename);
   std::cout << "done raytrace." << std::endl;
   return istate;

}

int 
coot::raytrace_info_t::render_ray_trace(std::string filename) {

   return render_ray_trace(filename, 1);
}
   
int
coot::raytrace_info_t::render_ray_trace(std::string filename, int reso_multiplier) {

   int istat = 1;
   std::ofstream render_stream(filename.c_str());

   if (!render_stream) {
      std::cout << "WARNING:: can't open file " << filename << std::endl;
      istat = 1;
      
   } else {

      // std::cout << "clipping: " << clipping << std::endl;

      // how do we preserve the aspect ratio?
      // adjust the tiles ratio according the the window height & width.

      int nxtiles = 80;
      int nytiles = 64;

      nxtiles = int(float(window_width)/8.0);
      nytiles = int(float(window_height)/8.0);

      if (reso_multiplier != 1) { 
         nxtiles *= reso_multiplier;
         nytiles *= reso_multiplier;
      }

      std::cout << "using tiles: " << nxtiles << " " << nytiles << std::endl;
   
      render_stream << "coot " << VERSION << ": Raster3D output\n";
      render_stream << nxtiles << " " << nytiles << "   tiles in x,y\n";
      render_stream << quality << " " << quality << "  pixels (x,y) per tile\n";
      render_stream << "4        anti-aliasing 3x3\n";
      //       render_stream << "0 0 0    background\n";
      render_stream << background.col[0] << " ";
      render_stream << background.col[1] << " ";
      render_stream << background.col[2] << "    background\n";
      if (raster3d_enable_shadows) 
         render_stream << "T";
      else 
         render_stream << "F";
      render_stream << "        shadows\n";
      render_stream << "25       Phong power\n";
      render_stream << "0.15     secondary light contribution\n";
      render_stream << "0.05     ambient light contribution\n";
      render_stream << "0.25     specular light contribution\n";
      render_stream << "0        eye position\n";
      render_stream << "1 1 1    main light source position\n";
//        render_stream << "1 0 0 0  input co-ordinates + radius transformation\n";
//        render_stream << "0 1 0 0  \n";
//        render_stream << "0 0 1 0  \n";

      // ------

      const float *mat = view_matrix.get();

      clipper::RTop_orth rtop;
      clipper::Mat33<double> clipper_mat(mat[0], mat[1], mat[ 2],
                                         mat[4], mat[5], mat[ 6],
                                         mat[8], mat[9], mat[10]);
      clipper::Coord_orth cco(-view_centre.x(), -view_centre.y(), -view_centre.z());
      rtop = clipper::RTop_orth(clipper_mat.inverse(), clipper::Coord_orth(0,0,0));

      clipper::Coord_orth new_centre = cco.transform(rtop);

      std::cout << "view centre:     " << cco.format()        << std::endl;
      std::cout << "render recentre: " << new_centre.format() << std::endl;

      // ------ 

      render_stream << mat[0] << " " << mat[1] << " " << mat[ 2] << " 0\n";
      render_stream << mat[4] << " " << mat[5] << " " << mat[ 6] << " 0\n";
      render_stream << mat[8] << " " << mat[9] << " " << mat[10] << " 0\n";
      render_stream << " " << new_centre.x() << " "
                    << new_centre.y() << " "
                    << new_centre.z() << " " << zoom*0.62 << "\n";
      //                    << new_centre.z() << " " << zoom*0.6 << "\n";
      render_stream << "3         mixed object types\n";
      render_stream << "*\n*\n*\n";

      // more head: fog & clipping planes
      // frontclip
      float frontclip =  0.1;
      float backclip  = -0.3;

      // fiddle
      // clipping 0: 
      // clipping 6: 
      // clipping -6: 
//       frontclip *= zoom/(1.0 + 0.1*clipping);
//       backclip  *= zoom/(1.0 + 0.1*clipping);

      // 1.0 is the "correct" addative factor, but lets skin it down a triffle.
      //      frontclip *= zoom*(0.85 - 0.1*clipping);
      //      backclip  *= zoom*(0.85 - 0.1*clipping);
      frontclip *= zoom*(1.0 - 0.1*clipping);
      backclip  *= zoom*(1.0 - 0.1*clipping);

      std::cout << "FRONTCLIP " << frontclip << std::endl;
      std::cout << "BACKCLIP "  << backclip << std::endl;

      render_stream << "16\n";
      render_stream << "FRONTCLIP " << frontclip << "\n";
      render_stream << "16\n";
      render_stream << "BACKCLIP " << backclip << "\n";
      render_stream << "16\n";
      render_stream << "FOG 0 1.0 1.0 1.0\n";
         
      render_molecules(render_stream);

      render_generic_objects(render_stream);

      render_labels(render_stream);
         
      render_stream.close();
      istat = 0;
   }
   return istat;
}

int coot::raytrace_info_t::renderman_render(std::string filename) {

   int istate = 0;
   int istat = 1;
   std::ofstream render_stream(filename.c_str());

   if (!render_stream) {
      std::cout << "WARNING:: can't open file " << filename << std::endl;
      istat = 1;
      
   } else {

      render_stream << "##RenderMan RIB-Structure 1.0\n";
      render_stream << "\n";
      render_stream << "FrameBegin 1\n";
      render_stream << "\n";
      render_stream << "Display \"" << filename << ".tif\" \"file\" \"rgba\"\n";
      render_stream << "Format 640 480 -1\n";
      render_stream << "ShadingRate 1\n";
      // render_stream << "Projection \"perspective\" \"fov\" [30]\n";
      render_stream << "Projection \"orthographic\"\n";

      // Get the camera orientation from the quaterions.

      // You can get the glOrtho from the GL context: Glget(GL_PROJECTION_MATRIX, &some_variable);
      
      render_stream << "ScreenWindow " << ortho_left << " " << ortho_right << " "
                    << ortho_bottom << " " << ortho_top << "\n";

      render_stream << "Exposure 1.0 1.3\n";
      // Atmosphere "mgfog" "background" [0.0 0.0 0.0] "mindistance" 196.67 "maxdistance" 206.66
      render_stream << "Translate 0 0 200\n";

      // render_stream << "FrameAspectRatio 1.33\n";
      render_stream << "Identity\n";
      render_stream << "\n";
      render_stream << "# Default distant headlight\n";
      render_stream << "LightSource \"distantlight\" 1\n";
      render_stream << "# Camera transformation\n"; 
      render_stream << "Translate 0 0 20\n"; 
      render_stream << "WorldBegin\n";
      
      render_stream << "Attribute \"visibility\"  # make objects visible to eye\n";
      render_stream << "Attribute \"trace\" \"bias\" 0.1\n";


      // mg code:

      //         def getCameraRotations(self,quat):
      //                 import pygl_coord
      //                 import math
      // 
      //                 dp = quat.Getdval()
      //                 tmp = pygl_coord.doublea(1)
      //                 quatp = tmp.frompointer(dp)
      // 
      //                 angle = 2*math.acos(quatp[0])
      //                 sina = math.sin(angle/2.)
      //                 x = quatp[1]/sina
      //                 y = quatp[2]/sina
      //                 z = quatp[3]/sina
      //                 return angle, x,y,z
      // 
      // ax,rx,ry,rz = self.getCameraRotations(quat)
      // outstream.Print("                Rotate "+str(-ax/math.pi*180.)+" "+str(rx)+" "+str(ry)+" "+str(-rz)+"\n")
      // outstream.Print("                Translate "+str(rpos[0])+" "+str(rpos[1])+" "+str(-rpos[2])+"\n");

      double angle = 2.0 * acos(view_quat.q0);
      double sina = sin(angle/2.0);
      double qx = view_quat.q1/sina;
      double qy = view_quat.q2/sina;
      double qz = view_quat.q3/sina;
      
      render_stream << "Rotate " // -317.207100403 0.242010379709 -0.95971828126 -0.142729803226
                    << angle*M_PI/180.0 << " " << qx << " " << qy << " " << qz << "\n";
      render_stream << "Translate " << view_centre.x() << " "
                    << view_centre.y() << " " << view_centre.z() 
                    << "\n";
      
                   // render_stream << "Identity\n";

      // molecule attributes
      
      for (unsigned int i=0; i<rt_mol_info.size(); i++) {
         std::cout << "rendman output for molecule : " << i << std::endl;
         rt_mol_info[i].renderman_molecule(render_stream, bond_thickness,
                                           atom_radius, density_thickness,
                                           bone_thickness);
      }
      
      render_stream << "WorldEnd\n";
      render_stream << "FrameEnd\n";

   } 
   return istate;
}



void
coot::raytrace_info_t::render_molecules(std::ofstream &render_stream) {

   for (unsigned int i=0; i<rt_mol_info.size(); i++) {
      std::cout << "rendering ray trace number: " << i << std::endl;
      render_stream << "# render for molecule no and name:"
                    << rt_mol_info[i].molecule_number << " "
                    << rt_mol_info[i].molecule_name   << "\n";
      rt_mol_info[i].render_molecule(render_stream, bond_thickness,
                                     atom_radius, density_thickness,
                                     bone_thickness);
   }
}

void
coot::ray_trace_molecule_info::render_molecule(std::ofstream &render_stream,
                                               float bond_thickness,
                                               float atom_radius,
                                               float density_thickness,
                                               float bone_thickness) {

   for(unsigned int id=0; id<density_lines.size(); id++) {
      render_stream << "5" << "\n";
      // coord1 radius coord2 dummy colour
      render_stream << "  "
                    << density_lines[id].first.x() << " "
                    << density_lines[id].first.y() << " "
                    << density_lines[id].first.z() << " "
                    << density_thickness << " "
                    << density_lines[id].second.x() << " "
                    << density_lines[id].second.y() << " "
                    << density_lines[id].second.z() << " "
                    << density_thickness << " "
                    << density_colour.col[0] << " "
                    << density_colour.col[1] << " "
                    << density_colour.col[2] << "\n";
   }

   for (unsigned int iset=0; iset<bond_lines.size(); iset++) {
      std::cout << "set " << iset << " colour " << bond_lines[iset].colour << std::endl;
      for (unsigned int ib=0; ib<bond_lines[iset].bonds.size(); ib++) {
         render_stream << "3" << "\n";
         // coord1 radius coord2 dummy colour
         render_stream << "  " 
                       << bond_lines[iset].bonds[ib].begin_pos.x() << " "
                       << bond_lines[iset].bonds[ib].begin_pos.y() << " "
                       << bond_lines[iset].bonds[ib].begin_pos.z() << " "
                       << bond_lines[iset].bonds[ib].bond_thickness << " "
                       << bond_lines[iset].bonds[ib].end_pos.x() << " "
                       << bond_lines[iset].bonds[ib].end_pos.y() << " "
                       << bond_lines[iset].bonds[ib].end_pos.z() << " "
                          << bond_lines[iset].bonds[ib].bond_thickness << " "
                       << bond_lines[iset].colour[0] << " "
                       << bond_lines[iset].colour[1] << " "
                       << bond_lines[iset].colour[2] << "\n";
      }
   }

   if (graphics_info_t::renderer_show_atoms_flag) { 
      for (unsigned int iat=0; iat<balls.size(); iat++) {
         double r = balls[iat].radius;
         render_stream << "2" << "\n";
         render_stream << balls[iat].pos.x() << " "
                       << balls[iat].pos.y() << " "
                       << balls[iat].pos.z() << " "
                       << r
                       << " " << balls[iat].colour[0]
                       << " " << balls[iat].colour[1]
                       << " " << balls[iat].colour[2]
                       << "\n";
      }
   }

   for (unsigned int ib=0; ib<bone_lines.size(); ib++) {
      render_stream << "5\n";
      // coord1 radius coord2 dummy colour
      render_stream << "  " 
                    << bone_lines[ib].first.x() << " "
                    << bone_lines[ib].first.y() << " "
                    << bone_lines[ib].first.z() << " "
                    << bone_thickness << " "
                    << bone_lines[ib].second.x() << " "
                    << bone_lines[ib].second.y() << " "
                    << bone_lines[ib].second.z() << " "
                    << bone_thickness << " "
                    << bones_colour.col[0] << " "
                    << bones_colour.col[1] << " "
                    << bones_colour.col[2] << "\n";
   }

   // Extra Restraints (Bond lines)
   // 
   for (unsigned int i=0; i<velr.size(); i++) {
      const extra_line_representation &l = velr[i];
      render_stream << "5\n";
      render_stream << l.p1.x() << " "
                    << l.p1.y() << " "
                    << l.p1.z() << " "
                    << l.thickness << " "
                    << l.p2.x() << " "
                    << l.p2.y() << " "
                    << l.p2.z() << " "
                    << l.thickness << " "
                    << l.c.col[0] << " "
                    << l.c.col[1] << " "
                    << l.c.col[2] << "\n";
   }

   for (unsigned int i=0; i<balls.size(); i++) { 
      const ball_t &ball = balls[i];
      render_stream << "2" << "\n";
      render_stream << ball.pos.x() << " "
                    << ball.pos.y() << " "
                    << ball.pos.z() << " "
                    << ball.radius 
                    << " " << ball.colour.col[0]
                    << " " << ball.colour.col[1]
                    << " " << ball.colour.col[2]
                    << "\n";
   }
}

void
coot::raytrace_info_t::render_generic_objects(std::ofstream &render_stream) const {

   for (unsigned int i=0; i<display_objects.size(); i++) {
      display_objects[i].raster3d(render_stream);
   }
}

void
coot::raytrace_info_t::render_labels(std::ofstream &s) const {

   if (labels.size() > 0) {
      s << "10\n";
      s << "\"Sans\" ";
      s << font_size_string;
      s << " \"Left-align\"\n";
      for (std::size_t i=0; i<labels.size(); i++) {
         s << "11\n  ";
         s << labels[i].second.x() << " "
           << labels[i].second.y() << " "
           << labels[i].second.z() << " "
           << atom_label_colour.col[0] << " "
           << atom_label_colour.col[1] << " "
           << atom_label_colour.col[2] << "\n"
           << labels[i].first << "\n";
      }
   }
}

void
coot::ray_trace_molecule_info::renderman_molecule(std::ofstream &render_stream,
                                                  float bond_thickness,
                                                  float atom_radius,
                                                  float density_thickness,
                                                  float bone_thickness) {

   for (unsigned int iset=0; iset<bond_lines.size(); iset++) {
      for (unsigned int ib=0; ib<bond_lines[iset].bonds.size(); ib++) {
         render_stream << "# render a bond\n";
         render_stream << "AttributeBegin\n";
         render_stream << "   Color ["
                       << bond_lines[iset].colour[0] << " "
                       << bond_lines[iset].colour[1] << " "
                       << bond_lines[iset].colour[2] << "]\n";
         render_stream << "   Surface \"plastic\" \"Ka\" [1] \"Kd\" [0.5] \"Ks\" 1 \"roughness\" 0.1\n";
         //                     << bond_lines[ib].first.x() << " "
         //                     << bond_lines[ib].first.y() << " "
         //                     << bond_lines[ib].first.z() << " "
         //                     << bond_lines[ib].second.x() << " "
         //                     << bond_lines[ib].second.y() << " "
         //                     << bond_lines[ib].second.z() << "\n";

         //       render_stream << "   TransformBegin\n"; // no need
         render_stream << "   Translate "
                       << bond_lines[iset].bonds[ib].begin_pos.x() << " "
                       << bond_lines[iset].bonds[ib].begin_pos.y() << " "
                       << bond_lines[iset].bonds[ib].begin_pos.z() << "\n";
         double l = (bond_lines[iset].bonds[ib].begin_pos - bond_lines[iset].bonds[ib].end_pos).amplitude();

         coot::Cartesian v = (bond_lines[iset].bonds[ib].begin_pos - bond_lines[iset].bonds[ib].end_pos);
         v.unit_vector_yourself();
         coot::Cartesian axis = coot::Cartesian::CrossProduct(v,Cartesian(0,0,1));
         double dp = coot::dot_product(v,coot::Cartesian(0,0,-1));
         // std::cout << " dot product: " << dp << std::endl;
         if (dp > 1.0)  dp =  1.0;
         if (dp < -1.0) dp = -1.0;
         double angle = -180.0/M_PI*acos(dp);
         if(fabs(axis.length())<1e-7) axis = coot::Cartesian(0,1,0);
         render_stream << "   Rotate "
                       << angle << " " << axis.x() << " " << axis.y() << " " << axis.z() << "\n";
         // Think about scaling the cylinder so that far away bonds are not tiny thin.
      
      
         render_stream << "   Cylinder 0.15 0 " << l << "  360\n";
         //       render_stream << "   TransformEnd\n";  // no need
         render_stream << "AttributeEnd\n";
      }
   }

}

void
coot::raytrace_info_t::set_ortho_params(float left, float right, float bottom, float top) {

   ortho_left = left;
   ortho_right = right;
   ortho_bottom = bottom;
   ortho_top = top;
} 


void 
meshed_generic_display_object::raster3d(std::ofstream &render_stream) const {

   std::cout << "soemthing here: raster3d()" << std::endl;
}


void 
coot::old_generic_display_object_t::raster3d(std::ofstream &render_stream) const {

   // Make sure that you have set is_displayed_flag for your generic_display_object_t

   if (! is_closed_flag) {
      if (is_displayed_flag) {
         for (unsigned int ils=0; ils<lines_set.size(); ils++) {
            float w = float(lines_set[ils].width)/80.0;
            for (unsigned int il=0; il<lines_set[ils].lines.size(); il++) {
               render_stream << "3" << "\n";
               render_stream << "   "
                             << lines_set[ils].lines[il].coords.first.x() << " " 
                             << lines_set[ils].lines[il].coords.first.y() << " " 
                             << lines_set[ils].lines[il].coords.first.z() << " "
                             << w << " "
                             << lines_set[ils].lines[il].coords.second.x() << " " 
                             << lines_set[ils].lines[il].coords.second.y() << " " 
                             << lines_set[ils].lines[il].coords.second.z() << " " 
                             << w << " "
                             << lines_set[ils].colour.red << " " 
                             << lines_set[ils].colour.green << " " 
                             << lines_set[ils].colour.blue << "\n";
            }
         } 
         for (unsigned int ips=0; ips<points_set.size(); ips++) {
            for (unsigned int ip=0; ip<points_set[ips].points.size(); ip++) {
               render_stream << "2" << "\n"
                             << "   "
                             << points_set[ips].points[ip].x() << " " 
                             << points_set[ips].points[ip].y() << " " 
                             << points_set[ips].points[ip].z() << " "
                             << 0.05 << " " 
                             << points_set[ips].colour.red << " " 
                             << points_set[ips].colour.green << " " 
                             << points_set[ips].colour.blue << "\n";
            }
         } 
      }
   } 
} 



// ---------------------------------------------------------------------
//                   povray
// ---------------------------------------------------------------------



int 
coot::raytrace_info_t::povray_ray_trace(std::string filename) {

   int ipov = 0;

   // setup the view and header:
   std::ofstream render_stream(filename.c_str());

   if (!render_stream) {
      std::cout << "WARNING:: can't open file " << filename << std::endl;
      ipov = 1;
      
   } else {

      // We have the eye position and the look at point, but we need
      // to rotate round the eye->centre_point vector.  In more
      // conventional terms this is a kappa (of polar phi, psi, kappa)
      // rotation.  Use the Rotation class of clipper to do the
      // conversion from a 3x3 matrix to polar coords.
      //
      clipper::Mat33<double> view_matrix_cl(view_matrix.matrix_element(0,0),
                                            view_matrix.matrix_element(0,1),
                                            view_matrix.matrix_element(0,2),
                                            view_matrix.matrix_element(1,0),
                                            view_matrix.matrix_element(1,1),
                                            view_matrix.matrix_element(1,2),
                                            -view_matrix.matrix_element(2,0),
                                            -view_matrix.matrix_element(2,1),
                                            -view_matrix.matrix_element(2,2));
                                            
      clipper::Mat33<double> view_matrix_inv_cl(view_matrix.matrix_element(0,0),
                                                view_matrix.matrix_element(1,0),
                                                view_matrix.matrix_element(2,0),
                                                view_matrix.matrix_element(0,1),
                                                view_matrix.matrix_element(1,1),
                                                view_matrix.matrix_element(2,1),
                                                -view_matrix.matrix_element(0,2),
                                                -view_matrix.matrix_element(1,2),
                                                -view_matrix.matrix_element(2,2));

      clipper::Polar_ccp4 polar = clipper::Rotation(view_matrix_cl).polar_ccp4();
      std::cout << "kappa: " << polar.kappa() << std::endl;

      GtkAllocation allocation = graphics_info_t::get_glarea_allocation();
      int x0 = allocation.width/2;
      int y0 = allocation.height/2;
      graphics_info_t g;
      glm::vec4 glm_screen_edge1 = g.unproject(0, 0, 0);
      glm::vec4 glm_screen_edge2 = g.unproject(x0*2, 0, 0);
      coot::Cartesian screen_edge1(glm_screen_edge1.x, glm_screen_edge1.y, glm_screen_edge1.z);
      coot::Cartesian screen_edge2(glm_screen_edge2.x, glm_screen_edge2.y, glm_screen_edge2.z);

      coot::Cartesian v1_2 = screen_edge2 - screen_edge1;

      clipper::Vec3<double> camera_location_cl(camera_location.x(),
					       camera_location.y(),
					       camera_location.z());

      clipper::Vec3<double> view_centre_cl(view_centre.x(),
                                           view_centre.y(),
                                           view_centre.z());

      float dir_len = (camera_location - view_centre).amplitude();

      clipper::Vec3<double> direction_cl(view_matrix.matrix_element(2,0),
                                         view_matrix.matrix_element(2,1),
                                         view_matrix.matrix_element(2,2));

      float tmp_len = view_centre_cl * direction_cl;
      
      float angle_factor = abs((v1_2.amplitude()/2)/(dir_len+tmp_len));
      //std::cout <<"BL DEBUG:: view centre 2 : " << view_centre_cl[0] << " " <<view_centre_cl[1] << " " << view_centre_cl[2] <<std::endl;
      //std::cout <<"BL DEBUG:: anglefactor   : " << angle_factor <<std::endl;
      //std::cout <<"BL DEBUG::   v1_2 ampl   : " << v1_2.amplitude() <<std::endl;
      //std::cout <<"BL DEBUG::   dir_len     : " << dir_len <<std::endl;
      //std::cout <<"BL DEBUG::   tmp_len     : " << tmp_len <<std::endl;
      //std::cout <<"BL DEBUG::   dir+tmp_len : " << dir_len+tmp_len <<std::endl;
      if (angle_factor > 1.99) {
        // simple protection, so that povray doesnt fail if angle get's too
        // large
        angle_factor = 1.99;
      }

      clipper::Vec3<double> tt_cl;
      for (int i=0; i<3; i++) {
        tt_cl[i] = dir_len * direction_cl[i];
      }
      clipper::Vec3<double> camera_translation_cl = view_centre_cl + tt_cl;

      // BL insert for spec
      //int x0 = graphics_info_t::glarea->allocation.width;
      //int y0 = graphics_info_t::glarea->allocation.height;
      float ratio = (float)x0/(float)y0;
      //

      double mv[16];
      glGetDoublev(GL_MODELVIEW_MATRIX, mv);
      
      // Still dont get the proper field of view angle properly done.
      // Just set it to 47 which is almost correct for all views I tested
      // But here an attempt to get FOV based on the projection matrix
      // assuming p[0] is 1/width of screen in orthographic projection
      /*
      double p[16];
      glGetDoublev(GL_PROJECTION_MATRIX, p);
      // FOV ?
      double width = 1 / p[0];
      double fov = 360 * atan(width / (dir_len * ratio)) / M_PI; // or dir_len * ratio?
      std::cout<< "BL DEBUG:: FOV? " << fov <<std::endl;
      */

      //

      render_stream << "#include \"colors.inc\"\n";
      render_stream << "background { color rgb <"
                    << graphics_info_t::background_colour[0] << ", "
                    << graphics_info_t::background_colour[1] << ", "
                    << graphics_info_t::background_colour[2] << "> }\n";
      /*      render_stream << "camera { orthographic\n"
                    << "        location  <"
                    << camera_translation_cl[0] << ", "
                    << camera_translation_cl[1] << ", "
                    << camera_translation_cl[2] << ">\n";
      */
      render_stream << "camera { orthographic\n"
                    << "         transform  {\n"
                    << "         matrix  <\n"
                    << "           "
                    << mv[0]  << ", " << mv[1]  << ", " << mv[2] << ",\n"
                    << "           "
                    << mv[4]  << ", " << mv[5]  << ", " << mv[6] << ",\n"
                    << "           "
                    << mv[8]  << ", " << mv[9]  << ", " << mv[10] << ",\n"
                    << "           "
                    << mv[12] << ", " << mv[13] << ", " << mv[14] << "\n"
                    << "         >\n"
                    << "         inverse }\n";


      render_stream << "         direction <0,0,-1>  \n";
      render_stream << "         location  <0,0," << dir_len * ratio << ">  \n";
      render_stream << "         angle     47  \n";
      render_stream << "         right     <"<< ratio <<",0,0> \n"
                    << "         up        <0,1,0>\n";
      /*
      render_stream << "        right     < " 
                    << view_matrix.matrix_element(0,0) * ratio << ", "
                    << view_matrix.matrix_element(0,1) * ratio << ", "
                    << view_matrix.matrix_element(0,2) * ratio << "> \n"
                    << "        up        < "
                    << view_matrix.matrix_element(1,0) << ", "
                    << view_matrix.matrix_element(1,1) << ", "
                    << view_matrix.matrix_element(1,2) << "> \n"
                    << "        direction < "
                    << -view_matrix.matrix_element(2,0) << ", "
                    << -view_matrix.matrix_element(2,1) << ", "
                    << -view_matrix.matrix_element(2,2) << "> \n";
      render_stream << "        angle  90* "<< angle_factor << " \n";
      */
      //render_stream << "         look_at  <"
      //<< view_centre.x() << ", "
      //<< view_centre.y() << ", "
      //<< view_centre.z() << ">}\n";
      render_stream << "}\n";
      render_stream << "light_source{<1500,  2500, -2500> colour White}\n";
      render_stream << "light_source{<1500, -2500,  1500> colour White}\n";
      render_stream << "light_source{<-1500, 2500,  1500> colour White}\n";

      povray_molecules(render_stream);
   }
   return ipov;
}


void
coot::raytrace_info_t::povray_molecules(std::ofstream &render_stream) {

   for (unsigned int i=0; i<rt_mol_info.size(); i++) {
      std::cout << "rendering povray ray trace number: " << i << std::endl;
      rt_mol_info[i].povray_molecule(render_stream,
                                     bond_thickness,
                                     density_thickness,
                                     atom_radius,
                                     zoom, view_centre, front_clipping_plane_point);
   }
}

void
coot::ray_trace_molecule_info::povray_molecule(std::ofstream &render_stream,
                                               float bond_thickness,
                                               float density_thickness,
                                               float atom_radius,
                                               float zoom,
                                               const coot::Cartesian &view_centre,
                                               const coot::Cartesian &front_clipping_plane_point) {

   // Clipping:
   //
   // We need to apply clipping planes.  How do we do that?
   //
   // We have the view centre to front clipping plane vector.
   // Let's look 90 degrees to that:
   //
   //           view            front GL
   //          centre           clipping
   //                            plane
   //             * -------------- * 
   //                             / \       .
   //                            /   \      .
   //                           /     \     .
   //        not clipped -->   o       \                         .
   //                                   o  <-- clipped out
   //
   // So the dot product of v1 (view centre to front clipping plane) and 
   // v2 (front clipping plane to point) should not be negative
   //
   // A similar argument applies to the back clipping plane.

   coot::Cartesian front_clip_to_centre_vec = view_centre - front_clipping_plane_point;
   coot::Cartesian back_clipping_plane_point = view_centre + front_clip_to_centre_vec;
   coot::Cartesian  back_clip_to_centre_vec = front_clipping_plane_point - view_centre;

   coot::Cartesian front_clipping_plane_point_new;
   coot::Cartesian tmp;
   // swap front and back
   //front_clipping_plane_point_new = back_clipping_plane_point;
   //back_clipping_plane_point = front_clipping_plane_point;
   //tmp = front_clip_to_centre_vec;
   //front_clip_to_centre_vec = back_clip_to_centre_vec;
   //back_clip_to_centre_vec = tmp;
   // now invert z
   //front_clipping_plane_point_new.invert_z();
   //back_clipping_plane_point.invert_z();
   //front_clip_to_centre_vec.invert_z();
   //back_clip_to_centre_vec.invert_z();

   // BL says:: for povray we have to invert z, i.e. z = -z; dont ask me why, but that's
   // how it is. We convert it here for now, but I guess there may be a better place and
   // way how to do it
   // we just do it before everything else get's done

   // elements of the form:
   // cylinder{<0,0,0>, <0,1,0>, 0.1 pigment {colour <0.1,0.2,0.3>} }
   for(unsigned int id=0; id<density_lines.size(); id++) {

      Cartesian v1 = density_lines[id].first  - front_clipping_plane_point;
      Cartesian v2 = density_lines[id].second - front_clipping_plane_point;
      Cartesian v3 = density_lines[id].first  -  back_clipping_plane_point;
      Cartesian v4 = density_lines[id].second -  back_clipping_plane_point;
      float dp1 = coot::dot_product(v1, front_clip_to_centre_vec);
      float dp2 = coot::dot_product(v2, front_clip_to_centre_vec);
      float dp3 = coot::dot_product(v3,  back_clip_to_centre_vec);
      float dp4 = coot::dot_product(v4,  back_clip_to_centre_vec);
      if ((dp1 > 0.0) && (dp2 > 0.0) && (dp3 > 0.0) && (dp4 > 0.0)) { 
        Cartesian line = density_lines[id].first - density_lines[id].second;
        if (line.length() > 0.001) {
         render_stream << "cylinder{<"
                       << density_lines[id].first.x() << ", "
                       << density_lines[id].first.y() << ", "
                       << density_lines[id].first.z() << ">\n "
                       << "         <"
                       << density_lines[id].second.x() << ", "
                       << density_lines[id].second.y() << ", "
                       << density_lines[id].second.z() << ">\n"
                       << "         " << density_thickness
                       << "   pigment { color <"
                       << density_colour.col[0] <<", "
                       << density_colour.col[1] <<", "
                       << density_colour.col[2] <<"> " << "} "
                       << "scale " << 1.0 
                       << "}\n";
        } else {
          render_stream<< "sphere{ <"
                       << density_lines[id].first.x() << ", "
                       << density_lines[id].first.y() << ", "
                       << density_lines[id].first.z() << "> "
                       << density_thickness
                       << "   pigment { color <"
                       << density_colour.col[0] <<", "
                       << density_colour.col[1] <<", "
                       << density_colour.col[2] <<">} "
                       << ""
                       << "}\n";
        }
      }
   }

   for (unsigned int iset=0; iset<bond_lines.size(); iset++) {
      for (unsigned int ib=0; ib<bond_lines[iset].bonds.size(); iset++) {
      
         Cartesian v1 = bond_lines[iset].bonds[ib].begin_pos  - front_clipping_plane_point;
         Cartesian v2 = bond_lines[iset].bonds[ib].end_pos    - front_clipping_plane_point;
         Cartesian v3 = bond_lines[iset].bonds[ib].begin_pos  -  back_clipping_plane_point;
         Cartesian v4 = bond_lines[iset].bonds[ib].end_pos    -  back_clipping_plane_point;
         float dp1 = coot::dot_product(v1, front_clip_to_centre_vec);
         float dp2 = coot::dot_product(v2, front_clip_to_centre_vec);
         float dp3 = coot::dot_product(v3,  back_clip_to_centre_vec);
         float dp4 = coot::dot_product(v4,  back_clip_to_centre_vec);
         if ((dp1 > 0.0) && (dp2 > 0.0) && (dp3 > 0.0) && (dp4 > 0.0)) { 
            render_stream << "cylinder{ <"
                          << bond_lines[iset].bonds[ib].begin_pos.x() << ", "
                          << bond_lines[iset].bonds[ib].begin_pos.y() << ", "
                          << bond_lines[iset].bonds[ib].begin_pos.z() << ">\n "
                          << "         <"
                          << bond_lines[iset].bonds[ib].end_pos.x() << ", "
                          << bond_lines[iset].bonds[ib].end_pos.y() << ", "
                          << bond_lines[iset].bonds[ib].end_pos.z() << ">\n"
                          << "         " << bond_thickness
                          << "   pigment { color <"
                          << bond_lines[iset].colour[0] <<", "
                          << bond_lines[iset].colour[1] <<", "
                          << bond_lines[iset].colour[2] <<"> " << "} "
               //                     << "scale " << 1.0*(100.0/zoom)
                          << "scale " << 1.0
                          << "}\n";
         }
      }
   }

   //                        << " " << atom[iat].second.col[0]
   // cylinder{<0,0,0>, <0,1,0>, 0.1 pigment {colour <0.1,0.2,0.3>} }
   if (graphics_info_t::renderer_show_atoms_flag) { 
      for (unsigned int iat=0; iat<balls.size(); iat++) {
      
         Cartesian v1 = balls[iat].pos  - front_clipping_plane_point;
         Cartesian v2 = balls[iat].pos  -  back_clipping_plane_point;
         float dp1 = coot::dot_product(v1, front_clip_to_centre_vec);
         float dp2 = coot::dot_product(v2,  back_clip_to_centre_vec);
         if ((dp1 > 0) && (dp2 > 0)) { 
            render_stream << "sphere{ <"
                          << balls[iat].pos.x() << ","
                          << balls[iat].pos.y() << ","
                          << balls[iat].pos.z() << ">"
                          << "0.3   pigment { color <"
                          << balls[iat].colour[0] <<", "
                          << balls[iat].colour[1] <<", "
                          << balls[iat].colour[2] <<">} "
                          << ""
                          << "} "
                          << "\n";
         }
      }
   }

}


// ---------------------------------------------------------------------
//                                 movies
// ---------------------------------------------------------------------
// static
void
graphics_info_t::dump_a_movie_image() {

   std::string number_str =
      coot::util::int_to_string(graphics_info_t::movie_frame_number);

   if (graphics_info_t::movie_frame_number < 10000)
      number_str = "0" + number_str;
   if (graphics_info_t::movie_frame_number < 1000)
      number_str = "0" + number_str;
   if (graphics_info_t::movie_frame_number < 100)
      number_str = "0" + number_str;
   if (graphics_info_t::movie_frame_number < 10)
      number_str = "0" + number_str;

   std::string file_name = graphics_info_t::movie_file_prefix;
   file_name += number_str;
   file_name += ".ppm";
   graphics_info_t::screendump_image(file_name);
   graphics_info_t::movie_frame_number++;

}


// static
int
graphics_info_t::screendump_image(const std::string &file_name) {

   GLint viewport[4];
   glGetIntegerv(GL_VIEWPORT, viewport);
   glPixelTransferi(GL_MAP_COLOR, GL_FALSE);
   glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
   glPixelStorei(GL_PACK_ALIGNMENT, 1);
   
   unsigned char* pixels = new unsigned char[viewport[2]*viewport[3]*IMAGEINFO_RGBA_SIZE];
   glReadPixels(0, 0, viewport[2], viewport[3], GL_RGBA, GL_UNSIGNED_BYTE, pixels);
   image_info iinfo(viewport[2], viewport[3], pixels, IMAGEINFO_RGBA);

   // file should be a ppm file for now, png when we add it to configure
   iinfo.invert();
   int istatus = 0;
   try { 
      istatus = iinfo.write(file_name.c_str());
   }
   catch (...) {
      std::string s("Can't write that image format at the moment.\n");
      s += "ppm is suggested instead.";
      wrapped_nothing_bad_dialog(s);
   }
   delete [] pixels; // does iinfo copy the data or the pointer? possible crash.
   return istatus;
}
