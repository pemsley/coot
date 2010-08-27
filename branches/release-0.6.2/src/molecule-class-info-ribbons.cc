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

#ifdef _MSC_VER
#include <windows.h>
#endif

#include <iostream>
#include <vector>

#include "mmdb_manager.h"
#include "mmdb-extras.h"
#include "mmdb.h"
#include "mmdb-crystal.h"

// guessing for now as I lost the orig file....
#include "mattype_.h"
#include "mman_manager.h"
#include "splineinfo.h"
#include "atom_util.h"
 
#include "cbuild.h"
#include "cdisplayobject.h"

#include "CParamsManager.h"
#include "mmut_colour.h"


#include "graphics-info.h"

#include "molecule-class-info.h"


void
molecule_class_info_t::make_ribbons() {

  InitMatType();
  graphics_info_t g;
  coot::protein_geometry* gp = g.Geom_p();
  std::string mon_lib_dir = gp->get_mon_lib_dir();
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

  CMMANManager* molHnd=0;
  CMGSBase *s_base = new CMGSBase((char *)mon_lib_dir.c_str(),
                                  (char *)mon_lib_dir.c_str(),
                                  (char *)"",
                                  (char *)ener_lib.c_str(),
                                  (char *)mon_lib_list.c_str(),
                                  (char *)elements.c_str());

  CMolBondParams *bond_params = new CMolBondParams(s_base);

  molHnd = new CMMANManager(s_base, bond_params);
  g_print("BL DEBUG:: reading file %s\n", name_.c_str());
  molHnd->ReadCoorFile(name_.c_str());
  int selHnd = molHnd->NewSelection();

  molHnd->SelectAtoms(selHnd, 0,"*",ANY_RES,"*",ANY_RES,"*","*","*","*","*",SKEY_OR );
  molHnd->GetMolBonds();
  molHnd->GetModel(1)->CalcSecStructure(0);
  int nSelAtoms;
  PPCAtom selAtoms=0;

  molHnd->GetSelIndex(selHnd,selAtoms,nSelAtoms);


  AtomColourVector* cv = new AtomColourVector();
  int *colours = new int[nSelAtoms];
  for(int i=0;i<nSelAtoms;i++)
//     colours[i] = 0;
     colours[i] = 4;

  cv->SetAtomColours(nSelAtoms,colours);

  //  cmc->SetMode(SECSTR);
  //std:: cout << "BL DEBUG:: nSelAtoms " << nSelAtoms <<std::endl;
  //cmc->GetAtomColourVector(nSelAtoms, (int*)cv);

  int spline_accu = 4 + ribbon_global_params.GetInt("solid_quality") * 4;
  SplineInfo sinfo = GetSplineInfo(molHnd, selHnd, cv,
				   spline_accu, -1, -1, 1, 0, 0);

//  Displayobject Ribbons;
//  CParamsManager params;
//  CParamsManager global_params;

  printf("BL DEBUG:: params for ribbons:\n");
  printf("BL DEBUG:: solid_quality %i\n", ribbon_global_params.GetInt("solid_quality"));
  printf("BL DEBUG:: ribbons_style %i\n", ribbon_params.GetInt("ribbon_style"));
  printf("BL DEBUG:: flatten_loop %i\n", ribbon_params.GetInt("flatten_loop"));

//  global_params.SetInt("solid_quality", ribbon_settings->solid_quality);
//  global_params.SetInt("solid_quality", 4);

//  params.SetFloat("cylinder_width", ribbon_settings->cylinder_width);
//  params.SetFloat("ribbon_width", ribbon_settings->ribbon_width);
//  params.SetFloat("alpha_helix_width", ribbon_settings->alpha_helix_width);
//  params.SetFloat("arrow_width",ribbon_settings->arrow_width);
//  params.SetFloat("arrow_length",ribbon_settings->arrow_length);
//  params.SetFloat("worm_width", ribbon_settings->worm_width);
//  params.SetInt("ribbon_style", ribbon_settings->ribbon_style);
//  params.SetInt("two_colour_ribbon", ribbon_settings->two_colour_ribbon);
//  params.SetInt("helix_style", ribbon_settings->helix_style);
//  params.SetInt("flatten_loop", ribbon_settings->flatten_loop);
//  params.SetInt("flatten_beta", ribbon_settings->flatten_beta);
//  params.SetInt("spline_beta_flat", ribbon_settings->spline_beta_flat);

  int mode=SPLINE;
  build_spline(sinfo, Ribbons, mode, ribbon_params, ribbon_global_params, "", "");

}

void
molecule_class_info_t::draw_ribbons() {

 if (cootribbons) {
   Ribbons.draw_solids();
   Ribbons.set_transparent(0);
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
molecule_class_info_t::set_global_ribbon_param(const std::string &name,
					       int value) {

  ribbon_global_params.SetInt(name.c_str(), value);
  make_ribbons(); // refresh?!
  graphics_info_t::graphics_draw();

}
