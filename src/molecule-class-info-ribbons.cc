/* src/molecule-class-info-ribbons.cc
 * 
 * Copyright 2008 by The University of York
 * Author: Bernhard Lohkamp
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */
#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#ifdef _MSC_VER
#include <windows.h>
#endif

#include <iostream>
#include <vector>

#include "mmdb-extras.h"
#include "mmdb.h"
#include "mmdb-crystal.h"

// guessing for now as I lost the orig file....
#include "mman_manager.h"
#include "splineinfo.h"
#include "atom_util.h"
 
#include "cbuild.h"
#include "cdisplayobject.h"

#include "CParamsManager.h"
#include "mg_colour.h"
#include "mmut_connectivity.h"

#include "graphics-info.h"
#include "globjects.h"

#include "molecule-class-info.h"
#include "sphere.h"


void
molecule_class_info_t::make_ribbons() {

  InitMatType();
  graphics_info_t g;
  std::string mon_lib_dir = g.master_mon_lib_dir;
  std::cout <<"BL DEBUG:: mon_lib_dir is " << mon_lib_dir<<std::endl;
  // here we dont know if it contains "data/monomers/", so if it doesnt
  // we add it -> move to coot-util at some point
  // and add windows "\" although maybe done already!?
  std::string mon_lib_stub;
  std::pair<std::string, std::string> tmp_pair;
  mon_lib_stub = mon_lib_dir;
  if (mon_lib_dir.substr(mon_lib_dir.size()-1) == "/") {
     mon_lib_stub = mon_lib_dir.substr(0, mon_lib_dir.size()-1);
  }
  tmp_pair = coot::util::split_string_on_last_slash(mon_lib_stub);
  if (tmp_pair.second != "monomers") {
    // we dont have monomers, but do we have at least data?
    if (tmp_pair.second != "data") {
      // we dont have data and monomers
      mon_lib_dir = mon_lib_stub;
      mon_lib_dir += "/data/monomers/";
    } else {
      // no monomers but data
      mon_lib_dir = mon_lib_stub;
      mon_lib_dir += "/monomers/";
    }
  }
  // now mon_lib_dir should be wit data and monomers

  std::string ener_lib = coot::util::append_dir_file(mon_lib_dir, "ener_lib.cif");
  if (!coot::file_exists(ener_lib)) {
    ener_lib = "";
    std::cout << "BL WARNING:: couldnt find ener_lib.cif, ribbons may not work" <<std::endl;
  }

  std::string mon_lib_list = coot::util::append_dir_file(mon_lib_dir, "mon_lib_list");
  if (!coot::file_exists(mon_lib_list)) {
    mon_lib_list = "";
    std::cout << "BL WARNING:: couldnt find mon_lib.list, ribbons may not work" <<std::endl;
  }

  std::string elements = coot::util::append_dir_file(mon_lib_dir, "elements.cif");
  if (!coot::file_exists(elements)) {
    elements = "";
    std::cout << "BL INFO:: couldnt find elements.cif, ribbons may not work but should [SMc]" <<std::endl;
  }

  CMMANManager* molHnd=0;
  //CMMDBManager* molHnd=0;
  CMGSBase *s_base = new CMGSBase((char *)mon_lib_dir.c_str(),
                                  (char *)mon_lib_dir.c_str(),
                                  "",
                                  (char *)ener_lib.c_str(),
                                  (char *)mon_lib_list.c_str(),
                                  (char *)elements.c_str());

  CMolBondParams *bond_params = new CMolBondParams(s_base);

  //molHnd = new CMMANManager(s_base, bond_params);
  molHnd = new CMMANManager();
  //molHnd = new CMMDBManager();
  g_print("BL DEBUG:: reading file %s\n", name_.c_str());
  //molHnd->ReadCoorFile(name_.c_str());
  //molHnd->(CMMANManager *)atom_sel.mol;
  //molHnd = atom_sel.mol;
  //CMMDBManager* newmol = get_residue_range_as_mol("A", 10, 90);
  molHnd = (CMMANManager*)atom_sel.mol;
  //molHnd = (CMMANManager*)newmol;
  //molHnd->SetSBaseAndBondParams(s_base, bond_params);
  int selHnd = molHnd->NewSelection();

  molHnd->SelectAtoms(selHnd, 0,"*",ANY_RES,"*",ANY_RES,"*","*","*","*","*",SKEY_OR );
  std::cout<<"BL DEBUG:: selHnd " << selHnd <<std::endl;
  //std::cout<<"BL DEBUG:: no of atoms in molHnd, i.e. selHnd" << molHnd->NumberOfAtoms(selHnd) <<std::endl;
  //molHnd->GetMolBonds();
  //molHnd->GetModel(1)->CalcSecStructure(0);
  molHnd->GetModel(1)->CalcSecStructure(0);
  std::cout<<"BL DEBUG:: gto to here " << selHnd <<std::endl;
  int nSelAtoms;
  PPCAtom selAtoms=0;

  molHnd->GetSelIndex(selHnd,selAtoms,nSelAtoms);
  std::cout<<"BL DEBUG:: gto to here 2 with no sel " << nSelAtoms <<std::endl;

  AtomColourVector cv;// = new AtomColourVector();

  int *colours = new int[nSelAtoms];
  int is, j;
  for(int i=0;i<nSelAtoms;i++) {
     colours[i] = 3;
  }
  PCColourSchemes schemes = new CColourSchemes();
  CMolColour cmc = CMolColour((CMMANManager*)molHnd, -1, schemes);
  cmc.SetMode (1, SECSTR, SECSTR, -1, 0);
  cmc.ReColour();
  //cmc.SetMode (1, BYATOMTYPE, BYATOMTYPE, -1, 0);
  //cmc.SetOneColour (5);
  //cmc.SetMode (1, ONECOLOUR, ONECOLOUR, -1, 0);
  cv = cmc.GetAtomColourVector();
  std::cout<<"BL DEBUG:: got to here 3 " <<std::endl;
  for (int i=0;i<10;i++)
     std::cout << "BL DEBUG:: col i before " << cv.GetRGB(i)[0] << " " <<  cv.GetRGB(i)[1] << " "<< cv.GetRGB(i)[2]<< std::endl;
  
  /*
  for (int i=0;i<cv.size();i++) {
     cv.GetRGB(i)[0] = cv.GetRGB(i)[0]/255.;
     cv.GetRGB(i)[1] = cv.GetRGB(i)[1]/255.;
     cv.GetRGB(i)[2] = cv.GetRGB(i)[2]/255.;
  }
  for (int i=0;i<10;i++) 
  std::cout << "BL DEBUG:: col i after " << cv.GetRGB(i)[0] << " " <<  cv.GetRGB(i)[1] << " "<< cv.GetRGB(i)[2]<< std::endl;
  */
  
  int spline_accu = 4 + ccp4mg_global_params.GetInt("solid_quality") * 4;
  //SplineInfo sinfo = GetSplineInfo(molHnd, selHnd, cv,
  //CMMANManager *manMolHnd = new CMMANManager(molHnd);
  SplineInfo sinfo = GetSplineInfo(molHnd, selHnd, cv,
                                   spline_accu, -1, -1,
                                   ribbon_params.GetInt("flatten_beta"),
                                   ribbon_params.GetInt("flatten_loop"),
                                   ribbon_params.GetInt("smooth_helix"));

  int mode=SPLINE;
  Ribbons.clear_prims();
  build_spline(sinfo, Ribbons, mode, ribbon_params, ccp4mg_global_params, "", "");

}

void
molecule_class_info_t::draw_ribbons() {

 if (cootribbons) {
//   Ribbons.set_transparent(1);
//   Ribbons.SetAlpha(0.5);
   // Maybe should be somewhere else?!
   glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
   glEnable(GL_NORMALIZE);
   glEnable(GL_LIGHTING);
   glShadeModel(GL_SMOOTH);
   glEnable(GL_DEPTH_TEST);
   glPolygonMode(GL_FRONT, GL_FILL);
   glPolygonMode(GL_BACK,GL_FILL);
   glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
   glEnable(GL_COLOR_MATERIAL);
   glEnable(GL_LINE_SMOOTH);
   glEnable(GL_BLEND);
   //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   glBlendFunc(GL_SRC_ALPHA, GL_ZERO);
   glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
   glHint(GL_FOG_HINT, GL_NICEST);
   glEnable(GL_LIGHT0);
   Ribbons.draw_solids();
   // disable again (for surface?!)
   glDisable(GL_LINE_SMOOTH);
   glDisable(GL_BLEND);
   glDisable(GL_LIGHTING);
 }

}

void
molecule_class_info_t::set_ribbon_param(const std::string &name,
					int value) {
   
  ribbon_params.SetInt(name.c_str(), value);
  make_ribbons(); //refresh
  graphics_info_t::graphics_draw();

}

void
molecule_class_info_t::set_ribbon_param(const std::string &name,
					float value) {

   ribbon_params.SetFloat(name.c_str(), value);
   make_ribbons(); //refresh
   graphics_info_t::graphics_draw();

}

void
molecule_class_info_t::set_global_ccp4mg_param(const std::string &name,
					       int value) {

  ccp4mg_global_params.SetInt(name.c_str(), value);
  make_ribbons(); // refresh?!
  graphics_info_t::graphics_draw();

}
 
void
molecule_class_info_t::set_global_ccp4mg_param(const std::string &name,
					       float value) {

  ccp4mg_global_params.SetFloat(name.c_str(), value);
  make_ribbons(); // refresh?!
  graphics_info_t::graphics_draw();

}


// for Aniso spheroids
// probably shoudlnt be here.... will move at some point
void
molecule_class_info_t::make_aniso_spheroids() {

  InitMatType();
  graphics_info_t g;
  std::string mon_lib_dir = g.master_mon_lib_dir;
  std::cout <<"BL DEBUG:: mon_lib_dir is " << mon_lib_dir<<std::endl;
  // here we dont know if it contains "data/monomers/", so if it doesnt
  // we add it -> move to coot-util at some point
  // and add windows "\" although maybe done already!?
  std::string mon_lib_stub;
  std::pair<std::string, std::string> tmp_pair;
  mon_lib_stub = mon_lib_dir;
  if (mon_lib_dir.substr(mon_lib_dir.size()-1) == "/") {
     mon_lib_stub = mon_lib_dir.substr(0, mon_lib_dir.size()-1);
  }
  tmp_pair = coot::util::split_string_on_last_slash(mon_lib_stub);
  if (tmp_pair.second != "monomers") {
    // we dont have monomers, but do we have at least data?
    if (tmp_pair.second != "data") {
      // we dont have data and monomers
      mon_lib_dir = mon_lib_stub;
      mon_lib_dir += "/data/monomers/";
    } else {
      // no monomers but data
      mon_lib_dir = mon_lib_stub;
      mon_lib_dir += "/monomers/";
    }
  }
  // now mon_lib_dir should be wit data and monomers

  std::string ener_lib = coot::util::append_dir_file(mon_lib_dir, "ener_lib.cif");
  if (!coot::file_exists(ener_lib)) {
    ener_lib = "";
    std::cout << "BL WARNING:: couldnt find ener_lib.cif, ribbons may not work" <<std::endl;
  }

  std::string mon_lib_list = coot::util::append_dir_file(mon_lib_dir, "mon_lib.list");
  if (!coot::file_exists(mon_lib_list)) {
    mon_lib_list = "";
    std::cout << "BL WARNING:: couldnt find mon_lib.list, ribbons may not work" <<std::endl;
  }

  std::string elements = coot::util::append_dir_file(mon_lib_dir, "elements.cif");
  if (!coot::file_exists(elements)) {
    elements = "";
    std::cout << "BL INFO:: couldnt find elements.cif, ribbons may not work but should [SMc]" <<std::endl;
  }

  AnisoSpheroids.clear_prims();
  CMMANManager* molHnd=0;
  CMGSBase *s_base = new CMGSBase((char *)mon_lib_dir.c_str(),
                                  (char *)mon_lib_dir.c_str(),
                                  "",
                                  (char *)ener_lib.c_str(),
                                  (char *)mon_lib_list.c_str(),
                                  (char *)elements.c_str());

  CMolBondParams *bond_params = new CMolBondParams(s_base);

  //molHnd = new CMMANManager(s_base, bond_params);
  molHnd = new CMMANManager();
  g_print("BL DEBUG:: reading file %s\n", name_.c_str());
  //molHnd->ReadCoorFile(name_.c_str());
  molHnd = (CMMANManager*)atom_sel.mol;
  int selHnd = molHnd->NewSelection();
  // PCColourSchemes schemes;

  molHnd->SelectAtoms(selHnd, 0,"*",ANY_RES,"*",ANY_RES,"*","*","*","*","*",SKEY_OR );
  int nSelAtoms;
  PPCAtom selAtoms=0;

  molHnd->GetSelIndex(selHnd,selAtoms,nSelAtoms);
  g_print ("BL DEBUG:: no of atoms: %i\n", nSelAtoms);

  AtomColourVector cv;

  PCColourSchemes schemes = new CColourSchemes();
  CColourScheme *AtomType = schemes->GetScheme("atomtype");
  char *atmtyps[7] = { "*"," C"," O", " N", " S"," H"," P" };
  char *atmcols[7] = { "grey", "yellow", "red", "blue" , "green", "grey","magenta" };
  //char *atmcols_coot[7] = {};
  std::vector<float> c_col = get_atom_colour_from_element(" C");
  std::cout << "BL DEBUG:: have atom col C " << c_col[0]<< " " << c_col[1] << " " << c_col[2]<< std::endl;
  std::vector<float> o_col = get_atom_colour_from_element(" O");
  std::cout << "BL DEBUG:: have atom col O " << o_col[0]<< " " << o_col[1] << " " << o_col[2]<< std::endl;

  int RC;
  RC = AtomType->SetSchemeString(7, atmtyps, atmcols);
  schemes->AtomType->GetSchemeCodes();
  CMolColour cmc = CMolColour(molHnd, -1, schemes);
  //cmc.SetMode (1, BYATOMTYPE, BYATOMTYPE,-1,0);
  //cmc.SetMode (1, BYRESTYPE, BYRESTYPE,-1,0);
  //cmc.SetMode (1, BYATOMTYPE, BYATOMTYPE,-1,0);
  cmc.SetMode (1, BYATOMTYPE, BYATOMTYPE);
  cmc.ReColour();

  std::vector<double> atom_radii;
  atom_radii = std::vector<double> (nSelAtoms);
//  for(int i=0;i<nSelAtoms;i++)
//     colours[i] = 0;

  for(int i=0;i<nSelAtoms;i++)
    atom_radii[i] = 1;
  //atom_radii = molHnd->GetAtomRadii(selHnd, VDWRADIUS, 1.);
  //cmc.Print();
  cv = cmc.GetAtomColourVector();
  for (int i=0;i<10;i++) 
    std::cout << "BL DEBUG:: col i " << cv.GetRGB(i)[0] << " " <<  cv.GetRGB(i)[1] << " "<< cv.GetRGB(i)[2]<< std::endl;

  // 1.2 is arbitrary scale, appears to be needed to get to same level... FIXME
  float spheroid_scale = graphics_info_t::show_aniso_atoms_probability / 50. * 1.2;
  
  g_print ("BL DEBUG:: no of atoms: %i\n", nSelAtoms);
  int aniso_style;
  aniso_style = SPHEROID_SOLID;
  //aniso_style = SPHEROID_AXES;
  //aniso_style = SPHEROID_SOLIDAXES;


  /*   std::vector<Cartesian>  cavertices;
   double *colour_array=0;
   double red[] = {1.0,0.0,0.0,1.0};
   double blue[] = {0.0,0.0,1.0,1.0};
   double green[] = {0.0,1.0,0.0,1.0};

   PolyCollection *polys = new PolyCollection();
   colour_array = green;
   cavertices.push_back(Cartesian(8.,22.,8.));
   SphereElement *sphere = new SphereElement(cavertices.back(),colour_array,cavertices.back(),0.8,1.0,1);


   polys->add_primitive(sphere);
   AnisoSpheroids.add_primitive(polys);
  */

  ccp4mg_global_params.SetInt("solid_quality", 2);

  DrawAnisoU(AnisoSpheroids,
	     selHnd, selAtoms, nSelAtoms,
	     cv, atom_radii,
	     aniso_style, spheroid_scale,
	     ccp4mg_global_params);

  /*colour_array = blue;
   cavertices.push_back(Cartesian(8.5,19.,8.5));
   sphere = new SphereElement(cavertices.back(),colour_array,cavertices.back(),0.8,1.0,1);
   PolyCollection *poly2 = new PolyCollection();
   poly2->add_primitive(sphere);
   AnisoSpheroids.add_primitive(poly2);
  */

}
/*
  ConnectivityDraw connd = ConnectivityDraw();
  Connectivity Conn = Connectivity();
  Conn.AddBonds(molHnd, selHnd, selAtoms, nSelAtoms, 0);
  std::cout<<"BL DEBUG:: no atom new " <<Conn.GetNumberOfAtoms()<<std::endl;;
  connd.SetParametersAndCalculate(Conn,
				  molHnd, AnisoSpheroids, 
				  SPHERES,   // mode,
				  //				  ANISOU,   // mode,
				  ribbon_params, // parameters.CParamsManager,
				  ccp4mg_global_params,  // global_parameters.CParamsManager,
				  nSelAtoms,
				  cv,   // atomColourVector,
				  atom_radii,  // atomRad,
				  "",
				  "",   // bumpmap,
				  -1, // stick_colour,
				  0,  // side_to_ribbon,
				  0,  // side_to_worm,
				  DRAW_ALL_BONDS //external_bonds_mode
				  );
}
*/

void
molecule_class_info_t::draw_aniso_spheroids() {

 if (cootanisospheroids) {
   // Maybe should be somewhere else?!
   // 21/3/12
   glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
   glEnable(GL_NORMALIZE);
   glEnable(GL_LIGHTING);
   glShadeModel(GL_SMOOTH);
   glEnable(GL_DEPTH_TEST);
   glEnable(GL_COLOR_MATERIAL);
   glEnable(GL_LINE_SMOOTH);
   glEnable(GL_BLEND);
   glPolygonMode(GL_FRONT, GL_FILL);
   glPolygonMode(GL_BACK,GL_FILL);
   glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
   glHint(GL_FOG_HINT, GL_NICEST);
   glEnable(GL_LIGHT0);
   AnisoSpheroids.draw_solids();
   glDisable(GL_LIGHTING);
   // disable again (for surface?!)
   glDisable(GL_LINE_SMOOTH);
   glDisable(GL_BLEND);

   /*
   std::vector<Cartesian>  cavertices;
   double *colour_array=0;
   double red[] = {1.0,0.0,0.0,1.0};
   double blue[] = {0.0,0.0,1.0,1.0};

   colour_array = red;
   cavertices.push_back(Cartesian(9.,20.,9.));
   SphereElement *sphere = new SphereElement(cavertices.back(),colour_array,cavertices.back(),0.8,1.0,1);
   matrix B(4,4);
   double Bfac = 1.;
   B(0,0) = Bfac;
   B(1,1) = B(2,2) = Bfac * 2.;
   B(3,3) = 1.0;
   cavertices.push_back(Cartesian(8.,20.,8.));
   Spheroid *spheroid = new Spheroid(cavertices.back(),colour_array,cavertices.back(),
                                     B,
                                     1.0, 2 , 0, 1);
   Displayobject test;
   test = Displayobject();
   test.add_primitive(sphere);
   test.add_primitive(spheroid);
   */

   /* 21/3/12
   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);

   glEnable(GL_LIGHTING);

   GLfloat bgcolor[4]={1.0, 1.0, 0.3, 1.0};
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   glMaterialfv(GL_FRONT, GL_SPECULAR, bgcolor);
   glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 128); 
   glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
   glEnable(GL_COLOR_MATERIAL);

   GLfloat  mat_specular[]  = {1.0, 0.3, 0.2, 1.0};
   GLfloat  mat_ambient[]   = {0.8, 0.1, 0.1, 1.0};
   GLfloat  mat_diffuse[]   = {0.2, 0.2, 0.2, 0.5};
   GLfloat  mat_shininess[] = {50.0};

   glClearColor(0.0, 0.0, 0.0, 0.0);
   glShadeModel(GL_SMOOTH);
   glMaterialfv(GL_FRONT, GL_SPECULAR,  mat_specular);
   glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
   glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
   glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

   glEnable(GL_DEPTH_TEST);
   glEnable(GL_NORMALIZE);

   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   glMaterialfv(GL_FRONT, GL_SPECULAR, bgcolor);

   glPushMatrix();
   */

   /*

   test.draw_solids();
   //
   cavertices.push_back(Cartesian(9.623,22.655,7.804));
   colour_array = blue;
   sphere = new SphereElement(cavertices.back(),colour_array,cavertices.back(),0.8,1.0,2);
   
   PolyCollection *polys = new PolyCollection();
   polys->add_primitive(sphere);

   AnisoSpheroids.add_primitive(polys);

   */
   /* 21/3/12
   AnisoSpheroids.draw_solids();

   glPopMatrix();

   glDisable(GL_LIGHTING);

   glDisable(GL_LIGHT0);
   glDisable(GL_LIGHTING);
   */

 }

}

std::vector<float>
molecule_class_info_t::get_atom_colour_from_element(const char *element) {

   std::vector<float> rgb(3);
   float rotation_size = bonds_colour_map_rotation/360.0;
   while (rotation_size > 1.0) { // no more black bonds?
      rotation_size -= 1.0;
   }
   int i_element;
   i_element = atom_colour(element);
   switch (i_element) {
   case YELLOW_BOND: 
      rgb[0] = 0.8; rgb[1] =  0.8; rgb[2] =  0.3;
      break;
   case BLUE_BOND: 
      rgb[0] = 0.5; rgb[1] =  0.5; rgb[2] =  1.0;
      break;
   case RED_BOND: 
      rgb[0] = 1.0; rgb[1] =  0.3; rgb[2] =  0.3;
      break;
   case GREEN_BOND:
      rgb[0] = 0.1; rgb[1] =  0.99; rgb[2] =  0.1;
      break;
   case GREY_BOND: 
      rgb[0] = 0.7; rgb[1] =  0.7; rgb[2] =  0.7;
      break;
   case MAGENTA_BOND:
      rgb[0] = 0.99; rgb[1] =  0.2; rgb[2] = 0.99;
      break;
   case ORANGE_BOND:
      rgb[0] = 0.89; rgb[1] =  0.89; rgb[2] = 0.1;
      break;
   case CYAN_BOND:
      rgb[0] = 0.1; rgb[1] =  0.89; rgb[2] = 0.89;
      break;
      
   default:
      rgb[0] = 0.8; rgb[1] =  0.2; rgb[2] =  0.2;
      rgb = rotate_rgb(rgb, float(imol_no*26.0/360.0));

   }
   // "correct" for the +1 added in the calculation of the rotation
   // size.
   // 21. is the default colour map rotation
   rgb = rotate_rgb(rgb, float(1.0 - 21.0/360.0));
   if (graphics_info_t::rotate_colour_map_on_read_pdb_c_only_flag) {
     if (i_element == YELLOW_BOND) { 
       rgb = rotate_rgb(rgb, rotation_size);
     } 
   } else {
     //       std::cout << "DEBUG: rotating coordinates colour map by "
     //                 << rotation_size * 360.0 << " degrees " << std::endl;
     rgb = rotate_rgb(rgb, rotation_size);
   }

   return rgb;

}
