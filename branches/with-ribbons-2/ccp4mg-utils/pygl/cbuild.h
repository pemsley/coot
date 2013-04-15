/*
     pygl/cbuild.h: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009-2010 University of York
     Copyright (C) 2012 STFC

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
*/


#ifndef _CCP4MG_BUILD_
#define _CCP4MG_BUILD_

#include "mgtree.h"
#include "cartesian.h"
#include "cdisplayobject.h"
#include "CParamsManager.h"
#include <vector>
#include <string>
#include "connect.h"
#include "splineinfo.h"
#include "rgbreps.h"
#include "atom_util.h"
#include <mmut_connectivity.h>
#include <mmut_basepairs.h>
#include <mmut_lipids.h>

enum { BONDS, CYLINDERS, SPHERES, BALLSTICK, SPLINE, WORM, FATWORM, FATBONDS, THINBONDS, NUCLEICBASEPAIRS, PYRAMIDS, LIPIDS, ANISOU, CIRCLES, IMPOSTER_SPHERES, VARIABLEWORM };
enum { NOLINE, LINE, DASHLINE, ARROW, DASHARROW, CYLINDER, CYLINDERARROW, CONE, DASHCYLINDER, DASHCYLINDERARROW };
enum { NOTLABELLED, LABELLEDCENTRE,  LABELLEDSTART, LABELLEDEND};

enum {SPHEROID_AXES,SPHEROID_SOLID,SPHEROID_SOLIDAXES};

enum {DRAW_ALL_BONDS, DRAW_INTERNAL_BONDS, DRAW_EXTERNAL_BONDS};

class ClickedLine {
  public:
   int first;
   int second;
   double dist;
   int symm;
};

void DrawLipids(Displayobject &obj, const std::vector<MMUTLipid> &lipids, CMMANManager *molHnd, const AtomColourVector &atom_colour_vector, const CParamsManager &params, const CParamsManager &global_params);

void DrawBasePairs(Displayobject &obj, const CNABasePairs &bp, const SplineInfo &splineinfo, const CParamsManager &params);
void DrawBaseBlocks(Displayobject &obj, CMMANManager *molHnd, int selHnd, PPCAtom selAtoms, int nSelAtoms, const AtomColourVector &atom_colour_vector, const CParamsManager &params);

void DrawAnisoU(Displayobject &obj, int selHndin, PPCAtom selAtoms, int nSelAtoms, const AtomColourVector &atom_colour_vector, const std::vector<double> &atomRadii,int style, double scale,const CParamsManager &global_params);

void DrawCircles(Displayobject &obj, const Connectivity &connectivity, const AtomColourVector &atom_colour_vector,const CParamsManager &params,const CParamsManager &global_params, int stick_colour);

void DrawImposterSpheres(Displayobject &obj, const Connectivity &connectivity, const AtomColourVector &atom_colour_vector, const CParamsManager &params, const CParamsManager &global_params, const std::vector<double> &atomRadii);

void DrawSimpleConnection(Displayobject &obj, const std::vector<SimpleConnection> &conn, const std::vector<double> &colour, int style=DASHLINE, int width=1, int labelstyle=NOTLABELLED, const std::string &labelcolour="", const std::string &family="",  const std::string &weight="", const std::string &slant="", const std::string &size="");

void DrawSimpleConnection(Displayobject &obj, const std::vector<SimpleConnection> &conn, const std::vector<double> &colour, int style=DASHLINE, int width=1, int labelstyle=NOTLABELLED, const std::vector<double> &labelcolf=std::vector<double>(), const std::string &family="",  const std::string &weight="", const std::string &slant="", const std::string &size="");

void DrawSimpleConnection(Displayobject &obj,
  const std::vector<SimpleConnection> &conn,
  const std::vector<double> &colour, int style, int width,
  int labelstyle, const std::string &labelcolour,
  const std::string &family,  const std::string &weight, 
  const std::string &slant, const std::string &size,
  const std::vector<int> &tags, const std::vector<int> &selTags);

void DrawSimpleConnection(Displayobject &obj,
  const std::vector<SimpleConnection> &conn,
  const std::vector<double> &colour, int style, int width,
  int labelstyle, const std::vector<double> &labelcolf,
  const std::string &family,  const std::string &weight, 
  const std::string &slant, const std::string &size,
  const std::vector<int> &tags, const std::vector<int> &selTags);

void build_spline(const SplineInfo &splineinfo, Displayobject &obj, int mode, const CParamsManager &params,  const CParamsManager &global_params, const std::string &texture, const std::string &bumpmap);
std::vector<int> GetPointsInVolume(Displayobject &obj, const std::vector<Cartesian> &atoms, const Volume &volume);
ClickedLine FindPoint(Displayobject &obj, const std::vector<Cartesian> &primorigin, const std::vector<Cartesian> &xyzbox, int symmetry=1);
ClickedLine FindLine(Displayobject &obj, const std::vector<Cartesian> &primorigin, const std::vector<Cartesian> &xyzbox, const std::vector<std::vector<int> > &conn_lists, int symmetry);
ClickedLine FindLine(Displayobject &obj, const std::vector<Cartesian> &primorigin, const std::vector<Cartesian> &xyzbox, const Connectivity &conn, int symmetry=1);
ClickedLine FindLine(Displayobject &obj, const std::vector<SimpleConnection> &conn, const std::vector<Cartesian> &xyzbox, int symmetry=1);

int AddTextLabel(Displayobject &obj, double x, double y, double z, const char *label);
int AddTextLabel(Displayobject &obj, double x, double y, double z, const std::string &label);
void DeleteTextLabel(Displayobject &obj, int text_id);
void AddBillBoardImage(Displayobject &obj, const char *filename, double *pcarts);
void AddBillBoardTextLabel(Displayobject &obj, double x, double y, const char *label);

int *GetTextIDS(Displayobject &obj);
int GetNumberOfTextIDS(Displayobject &obj);
void SetTextString(Displayobject &obj,int text_id, const char* new_string);
void SetTextString(Displayobject &obj,int text_id, const std::string &new_string);
const char* GetTextString(Displayobject &obj,int text_id);

class ConnectivityDraw{

  public:
    CMMANManager *molhnd;
    //Connectivity connectivity;
    std::vector<Cartesian> all_cart;
    ConnectivityDraw();
    void SetParametersAndCalculate(const Connectivity &connectivity_in, CMMANManager *molhnd, Displayobject &obj, int mode, const CParamsManager &params,const CParamsManager &global_params,  int nSelAtoms, const AtomColourVector &at_col_vec, const std::vector<double> &atomRadii, const std::string &texture, const std::string &bumpmap , int stick_colour=-1,int side_to_ribbon=0,int side_to_worm=0,int bonds_mode=DRAW_ALL_BONDS);
    ConnectivityDraw(const Connectivity &connectivity_in, CMMANManager *molhnd, Displayobject &obj, int mode, const CParamsManager &params,const CParamsManager &global_params,  int nSelAtoms, const AtomColourVector &at_col_vec, const std::vector<double> &atomRadii, const std::string &texture, const std::string &bumpmap , int stick_colour=-1,int side_to_ribbon=0,int side_to_worm=0,int bonds_mode=DRAW_ALL_BONDS);
    void RedrawPrimitives(Displayobject &obj, const Connectivity &connectivity, int mode, const CParamsManager &params, const CParamsManager &global_params, int nSelAtoms,  const AtomColourVector &at_col_vec, const std::vector<double> &atomRadii, const std::string &texture, const std::string &bumpmap, int stick_colour=-1,int side_to_ribbon=0,int side_to_worm=0,int bonds_mode=DRAW_ALL_BONDS);
	
};

void build_beta_surface(CMMANManager *molH, int atom_selHnd_in, Displayobject &obj, const CParamsManager &params, const AtomColourVector &atom_colour_vector);

#endif
