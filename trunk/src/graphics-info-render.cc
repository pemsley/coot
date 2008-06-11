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
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

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

#include "ppmutil.h"


// raster3d
short int
graphics_info_t::raster3d(std::string filename) {

   short int istate = 0;
   coot::colour_t background;
   background.col.resize(3);
   background.col[0] = background_colour[0];
   background.col[1] = background_colour[1];
   background.col[2] = background_colour[2];
   int width = 600;
   int height = 600;
   if (glarea) {
      width = glarea->allocation.width;
      height = glarea->allocation.height;
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

   std::cout << "Generating raytrace molecule objects..." << std::endl;
   for (int imol=0; imol<n_molecules(); imol++) {
      std::cout << " molecule " << imol << " in  raytrace" << std::endl;
      
      if (molecules[imol].has_model()) {
	 if (molecules[imol].is_displayed_p()) { 
	    rt.rt_mol_info.push_back(molecules[imol].fill_raster_model_info());
	 }
      }
      if (molecules[imol].has_map()) {
	 rt.rt_mol_info.push_back(molecules[imol].fill_raster_map_info(1));
	 if (molecules[imol].is_difference_map_p()) { 
	    rt.rt_mol_info.push_back(molecules[imol].fill_raster_map_info(-1));
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
   coot::raytrace_info_t rt(RotationCentre(), zoom, background,
			    glarea->allocation.width,
			    glarea->allocation.height, 
			    clipping_front,
			    raster3d_bond_thickness,
			    raster3d_bone_thickness,
			    raster3d_atom_radius,
			    raster3d_density_thickness);
   GL_matrix m;
   m.from_quaternion(quat);
   rt.set_view_matrix(m);
   // So where is the "eye"? We have to do an unproject:
   int x0 = glarea->allocation.width/2;
   int y0 = glarea->allocation.height/2;
   coot::Cartesian eye = unproject_xyz(x0, y0, 0);

   // It seems that for raster3d, this eye position is good, but for
   // povray, we are very close in.  So, let's try moving the eye back
   // a scaled amount

   
   for (int imol=0; imol<n_molecules(); imol++) {
      std::cout << " molecule " << imol << " in  raytrace" << std::endl;
      
      if (molecules[imol].has_model()) {
	 if (molecules[imol].is_displayed_p()) { 
	    rt.rt_mol_info.push_back(molecules[imol].fill_raster_model_info());
	    coot::Cartesian eye_centre = eye - rt.view_centre;
	    eye_centre *= 7.5;
	    coot::Cartesian far_eye = rt.view_centre - eye_centre;
	    std::cout <<"BL DEBUG:: eye is before " << eye << std::endl;
	    eye *= 1.02;
	    rt.set_front_clipping_plane_point(eye);
	    rt.set_camera_location(far_eye);
	 }
      }
      if (molecules[imol].has_map()) {
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


int 
coot::raytrace_info_t::render_ray_trace(std::string filename) {

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


      std::cout << "using tiles: " << nxtiles << " " << nytiles << std::endl;
   
      render_stream << "coot " << VERSION << ": Raster3D output\n";
      render_stream << nxtiles << " " << nytiles << "   tiles in x,y\n";
      render_stream << quality << " " << quality << "  pixels (x,y) per tile\n";
      render_stream << "4        anti-aliasing 3x3\n";
      //       render_stream << "0 0 0    background\n";
      render_stream << background.col[0] << " ";
      render_stream << background.col[1] << " ";
      render_stream << background.col[2] << "    background\n";
      render_stream << "T        shadows\n";
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
		    << new_centre.z() << " " << zoom*0.42 << "\n";
      //		    << new_centre.z() << " " << zoom*0.6 << "\n";
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
	 
      render_stream.close();
      istat = 0;
   }
   return istat;
} 


void
coot::raytrace_info_t::render_molecules(std::ofstream &render_stream) {

   for (unsigned int i=0; i<rt_mol_info.size(); i++) {
      std::cout << "rendering ray trace number: " << i << std::endl;
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

   for (unsigned int ib=0; ib<bond_lines.size(); ib++) {
      render_stream << "3" << "\n";
      // coord1 radius coord2 dummy colour
      render_stream << "  " 
		    << bond_lines[ib].first.x() << " "
		    << bond_lines[ib].first.y() << " "
		    << bond_lines[ib].first.z() << " "
		    << bond_thickness << " "
		    << bond_lines[ib].second.x() << " "
		    << bond_lines[ib].second.y() << " "
		    << bond_lines[ib].second.z() << " "
		    << bond_thickness << " "
		    << bond_colour[ib].col[0] << " "
		    << bond_colour[ib].col[1] << " "
		    << bond_colour[ib].col[2] << "\n";
   }

   if (graphics_info_t::renderer_show_atoms_flag) { 
      for (unsigned int iat=0; iat<atom.size(); iat++) {
	 render_stream << "2" << "\n";
	 render_stream << atom[iat].first.x() << " "
		       << atom[iat].first.y() << " "
		       << atom[iat].first.z() << " "
		       << atom_radius
		       << " " << atom[iat].second.col[0]
		       << " " << atom[iat].second.col[1]
		       << " " << atom[iat].second.col[2]
		       << "\n";
      }
   }

   for (unsigned int ib=0; ib<bone_lines.size(); ib++) {
      render_stream << "5" << "\n";
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

      int x0 = graphics_info_t::glarea->allocation.width/2;
      int y0 = graphics_info_t::glarea->allocation.height/2;
      coot::Cartesian screen_edge1 = unproject_xyz(0, 0, 0);
      coot::Cartesian screen_edge2 = unproject_xyz(x0*2, 0, 0);

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

      render_stream << "#include \"colors.inc\"\n";
      render_stream << "camera { orthographic\n            location <"
		    << camera_translation_cl[0] << ", "
		    << camera_translation_cl[1] << ", "
		    << camera_translation_cl[2] << ">\n";
      // BL insert for spec
      //int x0 = graphics_info_t::glarea->allocation.width;
      //int y0 = graphics_info_t::glarea->allocation.height;
      float ratio = (float)x0/(float)y0;
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
      render_stream << "         angle  90* "<< angle_factor << " \n";
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
	if (density_lines[id].first.x() != density_lines[id].second.x() &&
	    density_lines[id].first.y() != density_lines[id].second.y() &&
	    density_lines[id].first.z() != density_lines[id].second.z()) {	 
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
   for (unsigned int ib=0; ib<bond_lines.size(); ib++) {
      
      Cartesian v1 = bond_lines[ib].first  - front_clipping_plane_point;
      Cartesian v2 = bond_lines[ib].second - front_clipping_plane_point;
      Cartesian v3 = bond_lines[ib].first  -  back_clipping_plane_point;
      Cartesian v4 = bond_lines[ib].second -  back_clipping_plane_point;
      float dp1 = coot::dot_product(v1, front_clip_to_centre_vec);
      float dp2 = coot::dot_product(v2, front_clip_to_centre_vec);
      float dp3 = coot::dot_product(v3,  back_clip_to_centre_vec);
      float dp4 = coot::dot_product(v4,  back_clip_to_centre_vec);
      if ((dp1 > 0.0) && (dp2 > 0.0) && (dp3 > 0.0) && (dp4 > 0.0)) { 
	 render_stream << "cylinder{ <"
		       << bond_lines[ib].first.x() << ", "
		       << bond_lines[ib].first.y() << ", "
		       << bond_lines[ib].first.z() << ">\n "
		       << "         <"
		       << bond_lines[ib].second.x() << ", "
		       << bond_lines[ib].second.y() << ", "
		       << bond_lines[ib].second.z() << ">\n"
		       << "         " << bond_thickness
		       << "   pigment { color <"
		       << bond_colour[ib].col[0] <<", "
		       << bond_colour[ib].col[1] <<", "
		       << bond_colour[ib].col[2] <<"> " << "} "
	    // 		    << "scale " << 1.0*(100.0/zoom)
		       << "scale " << 1.0
		       << "}\n";
      }
   }

   // 		       << " " << atom[iat].second.col[0]
   // cylinder{<0,0,0>, <0,1,0>, 0.1 pigment {colour <0.1,0.2,0.3>} }
   if (graphics_info_t::renderer_show_atoms_flag) { 
      for (unsigned int iat=0; iat<atom.size(); iat++) {
      
	 Cartesian v1 = atom[iat].first  - front_clipping_plane_point;
	 Cartesian v2 = atom[iat].first  -  back_clipping_plane_point;
	 float dp1 = coot::dot_product(v1, front_clip_to_centre_vec);
	 float dp2 = coot::dot_product(v2,  back_clip_to_centre_vec);
	 if ((dp1 > 0) && (dp2 > 0)) { 
	    render_stream << "sphere{ <"
			  << atom[iat].first.x() << ","
			  << atom[iat].first.y() << ","
			  << atom[iat].first.z() << ">"
			  << "0.3   pigment { color <"
			  << atom[iat].second.col[0] <<", "
			  << atom[iat].second.col[1] <<", "
			  << atom[iat].second.col[2] <<">} "
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

   if (graphics_info_t::movie_frame_number < 1000)
      number_str = "0" + number_str;
   if (graphics_info_t::movie_frame_number < 100)
      number_str = "0" + number_str;
   if (graphics_info_t::movie_frame_number < 10)
      number_str = "0" + number_str;

   std::string file_name = graphics_info_t::movie_file_prefix;
   file_name += file_name + number_str;
   file_name += ".ppm";
   graphics_info_t::screendump_image(file_name);

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
