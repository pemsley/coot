/*
     pygl/build_tree_primitives.cc: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009 University of York

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

#ifdef _WIN32
#include <windows.h>
#endif

#include <iostream>
#include <stdlib.h>
#include "mgtree.h"
#include "cdisplayobject.h"
#include "cprimitive.h"
#include "CParamsManager.h"
#include "cartesian.h"
#include "cbuild.h"
#include "rgbreps.h"
#include "help.h"
#include "catmull.h"
#include <vector>
#include <utility>
#include <algorithm>
#include <math.h>
#include <string>
#include <string.h>
#include "cartesian.h"
#include "texture.h"
#include "matrix.h"
#include <stdlib.h>
#include "connect.h"
#include "splineinfo.h"
#include <mmut_connectivity.h>
#include <mmut_basepairs.h>
#include <mmut_lipids.h>
#include <mmut_util.h>
#include <mmdb_manager.h>
#include <atom_util.h>
#include <sstream>
#include <iostream>
#include <iomanip>

#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif

enum enum_SecStr { NOSECSTR, BETA, BULGE, TURN3, TURN4, TURN5, ALPHA }; // We want a #include enum_secstr from mmut_sec ....

ClickedLine FindLine(Displayobject &obj, const std::vector<Cartesian> &primorigin, const std::vector<Cartesian> &xyzbox, const std::vector<std::vector<int> > &conn_lists, int symmetry){
  /* Find lines nearest to, eg., a mouse click. */

  int clicked_symm = -1;

  std::vector<int> lines;
  lines.push_back(-1);
  lines.push_back(-1);

  matrix objrotmat = obj.quat.getMatrix();

  Cartesian front = Cartesian::MidPoint(objrotmat*xyzbox[2],objrotmat*xyzbox[6]);
  Cartesian back  = Cartesian::MidPoint(objrotmat*xyzbox[3],objrotmat*xyzbox[7]);

  std::vector<Plane> planes;
  std::vector<Cartesian> points;

  //planes.push_back(Plane(objrotmat*xyzbox[0],objrotmat*xyzbox[2],objrotmat*xyzbox[6])); // Front clipping plane
  //planes.push_back(Plane(objrotmat*xyzbox[1],objrotmat*xyzbox[7],objrotmat*xyzbox[3])); // Back clipping plane

  Volume v = GetClippingPlanes();
  planes.push_back(v.GetPlane(5));
  planes.push_back(v.GetPlane(4));
  points.push_back(planes[0].find_points_on_plane()[0]);
  points.push_back(planes[1].find_points_on_plane()[0]);

  //points.push_back(objrotmat*xyzbox[0]);
  //points.push_back(objrotmat*xyzbox[0]);

  double mindist = 1.0e+8;
  Cartesian prim0;
  Cartesian prim1;

  int i = 0;
  std::vector<std::vector<int> >::const_iterator conn_iter = conn_lists.begin();

  std::vector<Cartesian> unconnected;
  std::vector<int> unconnected_map;

  int nsym = obj.GetNumSymmetryMatrices();
  if(!symmetry) nsym = 0;
  int isym = 0;

  while(conn_iter!=conn_lists.end()){
    std::vector<int>::const_iterator k=conn_iter->begin();
    /* If k has no connections then we have a problem */
    if(k==conn_iter->end()){
      unconnected.push_back(primorigin[i]);
      unconnected_map.push_back(i);
    }
    while(k!=conn_iter->end()){
      isym = -1;
      do{
      prim0 = primorigin[i];
      prim1 = primorigin[*k];
      if(nsym>0&&isym>-1){
	matrix T = obj.GetSymmetryMatrix(isym);
	prim0 = T*prim0;
	prim1 = T*prim1;
      }
      std::vector<Plane>::iterator plane = planes.begin();
      std::vector<Cartesian>::const_iterator point = points.begin();
      int in_clip_planes0 = 1;
      int in_clip_planes1 = 1;
      while(plane!=planes.end()){
        Cartesian n = plane->get_normal();
        n.normalize();
        Cartesian p2prim = *point-prim0;
        p2prim.normalize();
        if(n.DotProduct(n,p2prim)>1e-3)
	  in_clip_planes0 = 0;
        Cartesian p2prim2 = *point-prim1;
        p2prim2.normalize();
        if(n.DotProduct(n,p2prim2)>1e-3)
	  in_clip_planes1 = 0;
        point++;
        plane++;
      }
      if(in_clip_planes0||in_clip_planes1){
	std::vector<double> linedist = DistanceBetweenTwoLines(front,back,prim0,prim1);
	double dist = fabs(linedist[0]);
	double u = linedist[2];
        if(dist<1.5){
        }
	if(u<0.25&&u>-0.25&&dist<0.5&&dist<mindist&&in_clip_planes0){
          mindist = dist;
	  lines[0] = i;
	  lines[1] = -1;
	  clicked_symm = isym;
	}
	if(u>0.75&&u<1.25&&dist<0.5&&dist<mindist&&in_clip_planes1){
          mindist = dist;
	  lines[0] = -1;
	  lines[1] = *k;
	  clicked_symm = isym;
	}
        if(u>=0.25&&u<=0.75&&dist<0.5&&dist<mindist&&in_clip_planes0&&in_clip_planes1) {
          mindist = dist;
	  lines[0] = i;
	  lines[1] = *k;
	  clicked_symm = isym;
        }
      }
      isym++;
      }while(isym<nsym);
      k++;
    }
    conn_iter++; i++;
  }

  ClickedLine nearprim = FindPoint(obj,unconnected,xyzbox,symmetry);
  if(nearprim.first>-1){
    Cartesian prim = primorigin[unconnected_map[nearprim.first]];
    if(nearprim.symm>-1){
      matrix T = obj.GetSymmetryMatrix(nearprim.symm);
      prim = T*prim;
    }
    std::vector<double> linedisttmp = DistanceBetweenTwoLines(front,back,prim,prim);
    double dist = fabs(linedisttmp[0]);
    if(dist<mindist){
      lines[0] = unconnected_map[nearprim.first];
      lines[1] = -1;
      mindist  = dist;
      clicked_symm = nearprim.symm;
    }
  }

  //cout << "Minimum distance: " << mindist << ", between: " << lines[0] << ", " << lines[1] << "\n";
  ClickedLine line;
  line.first  = lines[0];
  line.second = lines[1];
  line.dist   = mindist;
  line.symm = clicked_symm;
  return line;
}

ClickedLine FindLine(Displayobject &obj, const std::vector<SimpleConnection> &conn, const std::vector<Cartesian> &xyzbox, int symmetry){

  std::vector<int> lines;
  int clicked_symm = -1;

  lines.push_back(-1);
  lines.push_back(-1);

  matrix objrotmat = obj.quat.getMatrix();
  Cartesian dum;
  Cartesian front = dum.MidPoint(objrotmat*xyzbox[2],objrotmat*xyzbox[6]);
  Cartesian back  = dum.MidPoint(objrotmat*xyzbox[3],objrotmat*xyzbox[7]);

  std::vector<Plane> planes;
  std::vector<Cartesian> points;

  //planes.push_back(Plane(objrotmat*xyzbox[0],objrotmat*xyzbox[2],objrotmat*xyzbox[6])); // Front clipping plane
  //planes.push_back(Plane(objrotmat*xyzbox[1],objrotmat*xyzbox[7],objrotmat*xyzbox[3])); // Back clipping plane
  points.push_back(objrotmat*xyzbox[0]);
  points.push_back(objrotmat*xyzbox[1]);

  Volume v = GetClippingPlanes();
  planes.push_back(v.GetPlane(5));
  planes.push_back(v.GetPlane(4));

  double mindist = 1.0e+8;
  Cartesian prim0;
  Cartesian prim1;

  int i = 0;
  std::vector<SimpleConnection>::const_iterator conn_iter = conn.begin();

  int nsym = obj.GetNumSymmetryMatrices();
  if(!symmetry) nsym = 0;
  int isym = 0;

  while(conn_iter!=conn.end()){
      isym = -1;
      do{
      prim0 = conn_iter->first;
      prim1 = conn_iter->second;
      if(nsym>0&&isym>-1){
	matrix T = obj.GetSymmetryMatrix(isym);
	prim0 = T*prim0;
	prim1 = T*prim1;
      }
      std::vector<Plane>::iterator plane = planes.begin();
      std::vector<Cartesian>::const_iterator point = points.begin();
      int in_clip_planes = 1;
      while(plane!=planes.end()){
        Cartesian n = plane->get_normal();
        n.normalize();
        Cartesian p2prim = *point-prim0;
        p2prim.normalize();
        if(n.DotProduct(n,p2prim)>0.0)
	  in_clip_planes = 0;
        p2prim = *point-prim1;
        p2prim.normalize();
        if(n.DotProduct(n,p2prim)>0.0)
	  in_clip_planes = 0;
        point++;
        plane++;
      }
      if(in_clip_planes){
	std::vector<double> linedist = DistanceBetweenTwoLines(front,back,prim0,prim1);
	double dist = linedist[0];
	double u = linedist[2];
	if(u<0.25&&u>-0.25&&dist<0.5&&dist<mindist){
          mindist = dist;
	  lines[0] = i;
	  lines[1] = -1;
	  clicked_symm = isym;
	}
	if(u>0.75&&u<1.25&&dist<0.5&&dist<mindist){
          mindist = dist;
	  lines[0] = -1;
	  lines[1] = i;
	  clicked_symm = isym;
	}
        if(u>=0.25&&u<=0.75&&dist<0.5&&dist<mindist) {
          mindist = dist;
	  lines[0] = i;
	  lines[1] = i;
	  clicked_symm = isym;
        }
      }
      isym++;
      }while(isym<nsym);
    conn_iter++; i++;
  }

  ClickedLine line;
  line.first  = lines[0];
  line.second = lines[1];
  line.dist   = mindist;
  line.symm = clicked_symm;
  return line;
}



void DrawSimpleConnection(Displayobject &obj, const std::vector<SimpleConnection> &conn, const std::vector<double> &colour, int style, int width, int labelstyle, const std::string &labelcolour,
  const std::string &family,  const std::string &weight, 
  const std::string &slant, const std::string &size){
  double col[3] ={double(colour[0]),double(colour[1]),double(colour[2])};
  double alpha = double(colour[3]);
  int cylinders_accu = 8;
  double cylinder_radius_scale = 0.02;

  std::string label;
  std::vector<Cartesian> carts(2);
  Text *text;

  LineCollection *lines = new LineCollection();
  PolyCollection *polys = new PolyCollection();

  bool have_lines = false;
  bool have_polys = false;
  std::vector<SimpleConnection>::const_iterator i = conn.begin();
  while(i!=conn.end()){
    carts[0] = i->first;
    carts[1] = i->second;
    // style == NOLINE => do nothing
    if(style==DASHLINE){
      DashLineElement *line;
      line = new DashLineElement(carts,col,carts[0],double(width),alpha);
      lines->add_primitive(line);
      have_lines=true;
    }
    if(style==LINE){
      LineElement *line;
      line = new LineElement(carts,col,carts[0],double(width),alpha);
      lines->add_primitive(line);
      have_lines=true;
    }
    if(style==ARROW){
      Arrow *line;
      line = new Arrow(carts,col,carts[0],double(width),alpha);
      lines->add_primitive(line);
      have_lines=true;
    }
    if(style==DASHARROW){
      DashArrow *line;
      line = new DashArrow(carts,col,carts[0],double(width),alpha);
      lines->add_primitive(line);
      have_lines=true;
    }
    if(style==DASHCYLINDER){
      DashCylinderElement *line;
      line = new DashCylinderElement(carts,col,carts[0],double(width)*cylinder_radius_scale,
                                  alpha,cylinders_accu);
      line->SetDashLength(0.2);
      line->SetDashEnd(1);
      polys->add_primitive(line);
      have_polys=true;
    }
    if(style==CYLINDER){
      Cylinder *line;
      line = new Cylinder(carts,col,carts[0],double(width)*cylinder_radius_scale,cylinders_accu);
      polys->add_primitive(line);
      have_polys=true;
    }
    if(style==CYLINDERARROW){
      Cartesian p = 0.7 * carts[1] + 0.3 * carts[0];
      Cartesian p0 = carts[0];
      carts[0] = p;
      Cone *cone;
      cone = new Cone(carts,col,carts[0],double(width)*cylinder_radius_scale*2.0,alpha,cylinders_accu);
      polys->add_primitive(cone);
      carts[0] = p0;
      carts[1] = p;
      Cylinder *line;
      line = new Cylinder(carts,col,carts[0],double(width)*cylinder_radius_scale,alpha,cylinders_accu);
      polys->add_primitive(line);
      have_polys=true;
    }
    i++;
  }


  lines->SetSize((double)width);
  if(have_lines)obj.add_primitive(lines);
  if(have_polys)obj.add_primitive(polys);
  // Add text label
  if (labelstyle == NOTLABELLED) return;

  i = conn.begin();

  while(i!=conn.end()){
    if (labelcolour != "") 
    /*
     label = "<colour=\""+labelcolour+"\">"+ i->label+"</colour>";
    else
    */
      label =  i->label;
    if ( labelstyle == LABELLEDCENTRE )
      text = new Text ( i->first.MidPoint(i->first,i->second) , label,
              i->first.MidPoint(i->first,i->second));
    else if ( labelstyle == LABELLEDSTART)
      text = new Text ( i->first , label,i->first);
    else if ( labelstyle == LABELLEDEND)
      text = new Text ( i->second ,label,i->second);
    else
      text = new Text ( i->first , label,i->first);

    text->SetFontSize(18);
    if(family!="") text->SetFontFamily(family);
    if(weight!="") text->SetFontWeight(weight);
    if(slant!="") text->SetFontSlant(slant);
    if(size!="") text->SetFontSize(atoi(size.c_str()));
    if(labelcolour!=""&&labelcolour!="default"&&labelcolour!="complement"){
       std::vector<double> defcol = RGBReps::GetColour(labelcolour);
       text->SetColour(defcol[0],defcol[1],defcol[2],1.0);
    } else {
       text->SetDefaultColour();
    }
    //text->initialize();
    obj.add_text_primitive(text);
    i++;
  }
}


void DrawSimpleConnection(Displayobject &obj,
  const std::vector<SimpleConnection> &conn,
  const std::vector<double> &colour, int style, int width,
  int labelstyle, const std::string &labelcolour,
  const std::string &family,  const std::string &weight, 
  const std::string &slant, const std::string &size,
  const std::vector<int> &tags,const std::vector<int> &selTags){

  double col[3] ={double(colour[0]),double(colour[1]),double(colour[2])};
  double alpha = double(colour[3]);
  int cylinders_accu = 8;
  double cylinder_radius_scale = 0.02;
  bool apply_selection = false;
  std::string label;
  std::vector<Cartesian> carts(2);
  Text *text;

  if (selTags.size()>=1)apply_selection = true; 
  std::vector<SimpleConnection>::const_iterator i = conn.begin();
  std::vector<int>::const_iterator j = tags.begin();
  while(i!=conn.end()){
    if (!apply_selection || 
         std::find(selTags.begin(),selTags.end(),*j)!=selTags.end()) {

      carts[0] = i->first;
      carts[1] = i->second;
      if(style==DASHLINE){
        DashLine *line;
        line = new DashLine(carts,col,carts[0],double(width),alpha);
        obj.add_primitive(line);
      }
      if(style==LINE){
        Line *line;
        line = new Line(carts,col,carts[0],double(width),alpha);
        obj.add_primitive(line);
      }
      if(style==ARROW){
        Arrow *line;
        line = new Arrow(carts,col,carts[0],double(width),alpha);
        obj.add_primitive(line);
      }
      if(style==DASHARROW){
        DashArrow *line;
        line = new DashArrow(carts,col,carts[0],double(width),alpha);
        obj.add_primitive(line);
      }

      if(style==DASHCYLINDER){
        DashCylinderElement *line;
        line = new DashCylinderElement(carts,col,carts[0],double(width)*cylinder_radius_scale,
                                  alpha,cylinders_accu);
        line->SetDashLength(0.2);
        line->SetDashEnd(1);
        obj.add_primitive(line);
      }
 

      if(style==CYLINDER){
        Cylinder *line;
        line = new Cylinder(carts,col,carts[0],double(width)*cylinder_radius_scale,alpha,cylinders_accu);
        obj.add_primitive(line);
      }
      if(style==CYLINDERARROW){
        Cartesian p = 0.7 * carts[1] + 0.3 * carts[0];
        Cartesian p0 = carts[0];
        carts[0] = p;
        Cone *cone;
        cone = new Cone(carts,col,carts[0],double(width)*cylinder_radius_scale*2.0,alpha,cylinders_accu);
        obj.add_primitive(cone);
        carts[0] = p0;
        carts[1] = p;
        Cylinder *line;
        line = new Cylinder(carts,col,carts[0],double(width)*cylinder_radius_scale,alpha,cylinders_accu);
        obj.add_primitive(line);
      }
    }
    i++;
    j++;
  }

  // Add text label
  if (labelstyle == NOTLABELLED) return;

  i = conn.begin();
  j = tags.begin();
  while(i!=conn.end()){
    if (!apply_selection || 
           std::find(selTags.begin(),selTags.end(),*j)!=selTags.end()) {
      if (labelcolour != "") 
        label = "<colour=\""+labelcolour+"\">"+ i->label+"</colour>";
      else
        label =  i->label;
      if ( labelstyle == LABELLEDCENTRE )
        text = new Text ( i->first.MidPoint(i->first,i->second) , label, i->first.MidPoint(i->first,i->second));
      else if ( labelstyle == LABELLEDSTART)
        text = new Text ( i->first , label,i->first);
      else if ( labelstyle == LABELLEDEND)
        text = new Text ( i->second ,label,i->second);
      else 
        text = new Text ( i->first , label,i->first);
  
      obj.add_text_primitive(text);
      text->SetFontFamily(family);
      text->SetFontWeight(weight);
      text->SetFontSlant(slant);
      text->SetFontSize(atoi(size.c_str()));
    }
    i++;
    j++;
  }

}

std::vector<int> GetPointsInVolume(Displayobject &obj, const std::vector<Cartesian> &atoms, const Volume &volume){

  Cartesian atom_origin;
  int clicked;
  Plane plane;
  std::vector <Cartesian> points;
  Cartesian n;
  Cartesian p2atom;
  std::vector<int> clicked_atoms;

  for(unsigned int j=0;j<atoms.size();j++){
     clicked = 1;
     atom_origin = atoms[j];
     matrix mat = obj.quat.getInvMatrix();
     atom_origin = mat*atom_origin;
     atom_origin += obj.origin;
     for(int ii=0;ii<volume.GetNumberOfPlanes();ii++){
       plane = volume.GetPlane(ii);
       n = plane.get_normal();
       points = plane.find_points_on_plane();
       p2atom = points[0] - atom_origin;
       n.normalize();
       p2atom.normalize();
       if(n.DotProduct(n,p2atom)<0.0)
         clicked = 0;
     }
     if(clicked) clicked_atoms.push_back(j);
  }

  return clicked_atoms;

}

ClickedLine FindPoint(Displayobject &obj, const std::vector<Cartesian> &primorigin, const std::vector<Cartesian> &xyzbox, int symmetry){

  int clicked_symm = -1;
  int nearprim = findprimc(xyzbox,primorigin,obj.origin,obj.quat.getInvMatrix());
  std::vector<Cartesian> primorigin2 = primorigin;

  unsigned int nsym = obj.GetNumSymmetryMatrices();
  if(!symmetry) nsym = 0;

  if(nearprim==-1){
    for(unsigned int j=0;j<nsym;j++){
      matrix T = obj.GetSymmetryMatrix(j);
      primorigin2.clear();
      for(unsigned int i=0;i<primorigin.size();i++){
        Cartesian cart = primorigin[i];
        cart = T*cart;
        primorigin2.push_back(cart);
      }
      nearprim = findprimc(xyzbox,primorigin2,obj.origin,obj.quat.getInvMatrix());
      //std::cout << "For symmetry " << j << " found " << nearprim << std::endl;
      if(nearprim>-1){
	 clicked_symm = j;
         break;
      }
    }
  }

  ClickedLine line;
  line.first  = nearprim;
  line.second = -1;
  line.dist   = -1.0; // Calculate this later...
  line.symm = clicked_symm;

  return line;

}          

int *GetTextIDS(Displayobject &obj){
  return obj.GetTextIDS();
}

int GetNumberOfTextIDS(Displayobject &obj){
  return obj.GetNumberOfTextIDS();
}

void SetTextString(Displayobject &obj,int text_id, const char* new_string){
  obj.SetTextString(text_id,new_string);
}

void SetTextString(Displayobject &obj,int text_id, const std::string &new_string){
  obj.SetTextString(text_id,new_string);
}

const char* GetTextString(Displayobject &obj,int text_id){
  return obj.GetTextString(text_id);
}

void DeleteTextLabel(Displayobject &obj, int text_id){
  obj.DeleteTextPrimitive(text_id);
}

int AddTextLabel(Displayobject &obj, double x, double y, double z, const std::string &label){
  Text *text;
  Cartesian primorigin = Cartesian(x,y,z,1.0);
  text = new Text(primorigin,label,primorigin);
  obj.add_text_primitive(text);
  return text->GetID();
}

int AddTextLabel(Displayobject &obj, double x, double y, double z, const char *label){
  int newtextid = AddTextLabel(obj,x,y,z,std::string(label));
  return newtextid;
}

void AddBillBoardTextLabel(Displayobject &obj, double x, double y, const std::string &label){
  BillBoardText *text;
  Cartesian primorigin = Cartesian(x,y,0);
  text = new BillBoardText(primorigin,label,primorigin);
  obj.add_text_primitive(text);
}

void AddBillBoardTextLabel(Displayobject &obj, double x, double y, const char *label){
  AddBillBoardTextLabel(obj,x,y,std::string(label));
}

void FitToPolynomial(std::vector<Cartesian> &carts, int pass){
  std::vector <Cartesian> spline = SplineCurve(carts,(carts.size()-1)*4,2,pass);
  for(unsigned ii=1;ii<carts.size()-2;ii++){
    carts[ii] = spline[ii*(spline.size()+1)/(carts.size()-1)];
   }
}


void build_beta_surface(CMMANManager *molH, int atom_selHnd_in, Displayobject &obj, const CParamsManager &params, AtomColourVector *atom_colour_vector){


  std::cout << "build_beta_surface\n";
  int sec_str_mask[] = {1,1,1,0,0,0,0,0};
  std::string sec_strucs = molH->ListSecStructure(sec_str_mask);

  std::cout << sec_strucs << "\n";
  std::cout << "Done build_beta_surface\n";
  return;

  int CAselHnd;
  PPCAtom atomTable;
  int nAtoms;

  int atom_selHnd = molH->NewSelection();
  molH->Select(atom_selHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","*","*","*",SKEY_NEW);
  molH->Select(atom_selHnd,STYPE_ATOM,atom_selHnd_in,SKEY_AND);
  molH->ExcludeOverlappedAtoms(atom_selHnd,0.8);

  // Find all CA - to use as quick check if atom is in this set
  CAselHnd = molH->NewSelection();
  molH->Select(CAselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","CA","C","*",SKEY_NEW);
  molH->Select(CAselHnd,STYPE_ATOM,atom_selHnd_in,SKEY_AND);
  molH->ExcludeOverlappedAtoms(CAselHnd,0.8);
  molH->GetSelIndex ( atom_selHnd, atomTable, nAtoms );
  std::vector<Cartesian>  cavertices;
  double *colour_array=0;
  double red[] = {1.0,0.0,0.0,1.0};

  /* We won't try anything fancy with colours just yet. */
  //if(!atom_colour_vector)
     colour_array = red;

  double min_x = 1e+8;
  double min_y = 1e+8;
  double max_x = 1e-8;
  double max_y = 1e-8;

  std::cout << "\n";
  //PolyCollection *polys = new PolyCollection();
  for(int j=0;j<nAtoms;j++){
    if(atomTable[j]->isInSelection(CAselHnd) && molH->isAminoacid(atomTable[j]->residue)){
      PCAtom pCA = atomTable[j];
      PCResidue pRes = pCA->residue;
      // Save the CA pointer
      if(int(pRes->SSE)== SSE_Strand || int(pRes->SSE)== SSE_Bulge){
        cavertices.push_back(Cartesian(pCA->x,pCA->y,pCA->z));
        //SphereElement *sphere = new SphereElement(cavertices.back(),colour_array,cavertices.back(),0.4,1.0,2);
        //polys->add_primitive(sphere);
        if(pCA->x<min_x) min_x = pCA->x;
        if(pCA->y<min_y) min_y = pCA->y;
        if(pCA->x>max_x) max_x = pCA->x;
        if(pCA->y>max_y) max_y = pCA->y;
        std::cout << cavertices.back() << "\n";
        //if( atom_colour_vector ){
          //colour_array = atom_colour_vector->GetRGB(j);
        //}
      }
    }
  }
  //obj.add_primitive(polys);

  std::cout << "\n";
  std::cout << min_x << " " << max_x << "\n";
  std::cout << min_y << " " << max_y << "\n";
  std::cout << "\n";

  //std::cout << cavertices.size() << "\n";
  if(cavertices.size()>5){//Need at least 6 sets of coords to satisfy our 6 unknown coeffs.
    std::vector<Cartesian> carts(2);
    double width = 2.0;
    LineCollection *lines = new LineCollection();
    min_x -= 3;
    min_y -= 3;
    max_x += 3;
    max_y += 3;
    std::vector<double> poly_params = LeastSquaresQuadraticFit3D(cavertices);
    //std::cout << "Draw function " << poly_params[0] << "x^2 + " << poly_params[1] << "y^2 + " << poly_params[2] << "xy + " << poly_params[3] << "x + " << poly_params[4] << "y + " << poly_params[5] << "\n";
    //std::cout << "In range: " << min_x << " -> " << max_x << ", " << min_y << " -> " << max_y << "\n";
    double x = min_x;
    double delta_x = (max_x-min_x)/30.0;
    double delta_y = (max_y-min_y)/30.0;
    while(x<max_x){
      double y = min_y;
      while(y<max_y){
        double z = poly_params[0]*x*x + poly_params[1]*y*y + poly_params[2]*x*y + poly_params[3]*x + poly_params[4]*y + poly_params[5];
        carts[0] = Cartesian(x,y,z);
        z = poly_params[0]*x*x + poly_params[1]*(y+delta_y)*(y+delta_y) + poly_params[2]*x*(y+delta_y) + poly_params[3]*x + poly_params[4]*(y+delta_y) + poly_params[5];
        carts[1] = Cartesian(x,y+delta_y,z);
        //std::cout << carts[1].get_x()-carts[0].get_x() << " " << carts[1].get_y()-carts[0].get_y() << "\n";
        LineElement *line = new LineElement(carts,colour_array,carts[0],width,1.0);
        lines->add_primitive(line);
        z = poly_params[0]*(x+delta_x)*(x+delta_x) + poly_params[1]*y*y + poly_params[2]*(x+delta_x)*y + poly_params[3]*(x+delta_x) + poly_params[4]*y + poly_params[5];
        carts[1] = Cartesian(x+delta_x,y,z);
        line = new LineElement(carts,colour_array,carts[0],width,1.0);
        lines->add_primitive(line);
        y += delta_y;
      }
      x += delta_x;
    }
    lines->SetSize(width);
    obj.add_primitive(lines);
  }

}

void build_spline(const SplineInfo &splineinfo, Displayobject &obj, int mode, const CParamsManager &params,  const CParamsManager &global_params, const std::string &texture, const std::string &bumpmap){

  unsigned int i;
  int ribbon_accus[] = {4, 9, 12, 18, 30, 36};

  //std::cout << "into build_spline" << std::endl;

  int multicolour = 1;
  int spline_accu = 4+4*global_params.GetInt("solid_quality");
  

  if(bumpmap!=""&&mode!=BONDS&&mode!=FATBONDS&&mode!=THINBONDS){
    image_info iinfo = image_info(bumpmap.c_str());
    load_texture(iinfo,MIPMAP);
    glEnable(GL_TEXTURE_2D);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glEnable(GL_TEXTURE_GEN_S);
    glEnable(GL_TEXTURE_GEN_T);
    glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
    glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
  }


  std::vector<std::vector<Cartesian> >exp_atomColourVector;
  if(multicolour&&(int)splineinfo.colours.size()>0) {
    for(i=0;i<splineinfo.colours.size();i++){
      exp_atomColourVector.push_back(std::vector<Cartesian>(0));
      for(int j=0;j<int(splineinfo.colours[i].size())-1;j++){
	for(int k=0;k<spline_accu;k++){
          /* Need to be more clever since multicolour is too much of a catchall */
          double frac = 0.0;// double(k)/spline_accu;
          Cartesian frac_col = (1.0-frac) * splineinfo.colours[i][j] + frac * splineinfo.colours[i][j+1];
	  exp_atomColourVector[i].push_back(frac_col);
        }
      }
      for(int k=0;k<spline_accu;k++){
        double frac = 0.0;// double(k)/spline_accu;
        Cartesian frac_col = frac * splineinfo.colours[i][splineinfo.colours[i].size()-1] + (1.0-frac) * splineinfo.colours[i][splineinfo.colours[i].size()-2];
        exp_atomColourVector[i].push_back(frac_col);
      } 
    }
  }

  std::vector<Cartesian> cartesians;
  std::vector<Cartesian> n1_cartesians;
  std::vector<Cartesian> n2_cartesians;


  std::vector<std::vector<Cartesian> >::const_iterator splines_iter=splineinfo.splines.begin();
  std::vector<Cartesian>::const_iterator spline_iter;
  std::vector<std::vector<Cartesian> >::const_iterator n1_splines_iter=splineinfo.n1_splines.begin();
  std::vector<Cartesian>::const_iterator n1_spline_iter;
  std::vector<std::vector<Cartesian> >::const_iterator n2_splines_iter=splineinfo.n2_splines.begin();
  std::vector<Cartesian>::const_iterator n2_spline_iter;
  std::vector<std::vector<Cartesian> >::const_iterator colour_vecs=exp_atomColourVector.begin();
  std::vector<Cartesian>::const_iterator colour_vec;
  std::vector<Cartesian> colours_vec;

  std::vector<std::vector<Cartesian> >::const_iterator nasplines_iter=splineinfo.nasplines.begin();
  std::vector<std::vector<Cartesian> >::const_iterator n1_nasplines_iter=splineinfo.n1_nasplines.begin();
  std::vector<std::vector<Cartesian> >::const_iterator n2_nasplines_iter=splineinfo.n2_nasplines.begin();
  std::vector<std::vector<Cartesian> > naexp_atomColourVector;
  if(multicolour&&(int)splineinfo.nacolours.size()>0) {
    for(i=0;i<splineinfo.nacolours.size();i++){
      naexp_atomColourVector.push_back(std::vector<Cartesian>(0));
      for(unsigned int j=0;j<splineinfo.nacolours[i].size();j++)
	for(int k=0;k<spline_accu;k++)
	   naexp_atomColourVector[i].push_back(splineinfo.nacolours[i][j]);
    }
  }
  std::vector<std::vector<Cartesian> >::const_iterator nacolours_iter=naexp_atomColourVector.begin();
  while(nasplines_iter!=splineinfo.nasplines.end()){
    if(nasplines_iter->size()>2){
    float worm_width = params.GetFloat("worm_width");
    float ribbon_width = params.GetFloat("ribbon_width");
    int ribbon_style = params.GetInt("ribbon_style");
    int ribbon_accu = ribbon_accus[global_params.GetInt("solid_quality")];
    int spline_accu = 4+4*global_params.GetInt("solid_quality");
    if (mode == SPLINE|| mode == FATWORM){
      double *col_tmp = RGBReps::GetColourP(1);
      Ribbon *ribbon = new Ribbon(*nasplines_iter,*n1_nasplines_iter,*n2_nasplines_iter,col_tmp,(*nasplines_iter)[0],*nacolours_iter,worm_width,ribbon_width,worm_width,1.0,-((ribbon_style<<16)|ribbon_accu*2),spline_accu);
      obj.add_primitive(ribbon);
      delete [] col_tmp;
    }else{
      double *col_tmp = RGBReps::GetColourP(1);
      Ribbon *ribbon = new Ribbon(*nasplines_iter,*n1_nasplines_iter,*n2_nasplines_iter,col_tmp,(*nasplines_iter)[0],*nacolours_iter,worm_width,worm_width,worm_width,1.0,ribbon_accu*2,spline_accu);
      obj.add_primitive(ribbon);
      delete [] col_tmp;
    }
    }
    nasplines_iter++;
    n1_nasplines_iter++;
    n2_nasplines_iter++;
    if(multicolour) nacolours_iter++;
  }

  int totalpoints = 0;
  int nchain = 0;
  while(splines_iter!=splineinfo.splines.end()){
    std::vector<std::vector<int> > secstr_indices = splineinfo.secstr_indices[nchain];
    std::vector<std::vector<int> >::const_iterator secstr_iter=secstr_indices.begin();
    nchain++;
    if(splines_iter->size()<2){
      splines_iter++;
      n1_splines_iter++;
      n2_splines_iter++;
      if(multicolour) colour_vecs++;
      secstr_indices = splineinfo.secstr_indices[nchain];
      secstr_iter=secstr_indices.begin();
      nchain++;
    }
    spline_iter=(*splines_iter).begin();
    n1_spline_iter=(*n1_splines_iter).begin();
    n2_spline_iter=(*n2_splines_iter).begin();

    totalpoints += (*splines_iter).size();
    if(multicolour) colour_vec=(*colour_vecs).begin();
    //std::cout << "New chain, size: " <<  (*splines_iter).size() << "\n";
    while(secstr_iter!=secstr_indices.end()){
      if(secstr_iter<secstr_indices.end()-1) {
        //std::cout << (*secstr_iter)[0] << " to " << (*(secstr_iter+1))[0] << "\n";
        int begin = (*secstr_iter)[0]*spline_accu;
        int end = (*(secstr_iter+1))[0]*spline_accu;
        for(int i=begin;i<=end&&spline_iter!=(*splines_iter).end();i++){ 
	  cartesians.push_back(*spline_iter);
	  n1_cartesians.push_back(*n1_spline_iter);
	  n2_cartesians.push_back(*n2_spline_iter);
	  if(multicolour) colours_vec.push_back(*colour_vec);
	  n1_spline_iter++;
	  n2_spline_iter++;
	  spline_iter++;
	  if(multicolour) colour_vec++;
        }
        bool doOneMore = (spline_iter+1==(*splines_iter).end()||(spline_iter+2!=(*splines_iter).end()&&(*(spline_iter+2)-*(spline_iter+2)).length()>15./spline_accu));
	if(doOneMore&&spline_iter!=(*splines_iter).end()&&n1_spline_iter!=(*n1_splines_iter).end()&&n2_spline_iter!=(*n2_splines_iter).end()){
	  cartesians.push_back(*spline_iter);
	  n1_cartesians.push_back(*n1_spline_iter);
	  n2_cartesians.push_back(*n2_spline_iter);
	  if(multicolour&&(colour_vec!=(*colour_vecs).end())) colours_vec.push_back(*colour_vec);
        }
        //std::cout << "cartesians.size() " << cartesians.size() << "\n"; std::cout.flush();
        totalpoints++;
        spline_iter--;
        n1_spline_iter--;
        n2_spline_iter--;
	if(multicolour) colour_vec--;
      }else{
	while(spline_iter!=(*splines_iter).end()){
	  cartesians.push_back(*spline_iter);
	  n1_cartesians.push_back(*n1_spline_iter);
	  n2_cartesians.push_back(*n2_spline_iter);
	  if(multicolour) colours_vec.push_back(*colour_vec);
	  n1_spline_iter++;
	  n2_spline_iter++;
	  spline_iter++;
	  if(multicolour) colour_vec++;
	}
        spline_iter--;
        n1_spline_iter--;
        n2_spline_iter--;
	if(multicolour) colour_vec--;
      }
      if(cartesians.size()>2) {
        if((cartesians[cartesians.size()-1]-cartesians[cartesians.size()-2]).length()<1e-4){
           cartesians.pop_back();
           n1_cartesians.pop_back();
           n2_cartesians.pop_back();
        }
        //std::cout << cartesians.size() << " " << n1_cartesians.size() << " " << n2_cartesians.size() << "\n";
        //std::cout << cartesians[0] << " " << cartesians[cartesians.size()-1] << "\n";
        //if(cartesians.size()>1) cartesians.pop_back();
        float worm_width = params.GetFloat("worm_width");
        float arrow_length = params.GetFloat("arrow_width");
        float arrow_width = params.GetFloat("arrow_width");
        int ribbon_accu = ribbon_accus[global_params.GetInt("solid_quality")];
        int ribbon_style = params.GetInt("ribbon_style");
        int helix_style = params.GetInt("helix_style");
        int spline_accu = 4+4*global_params.GetInt("solid_quality");
        int two_colour_ribbon = params.GetInt("two_colour_ribbon");
        float alpha_helix_width;
        float beta_sheet_width;
	float helix_tube_diameter;
        if (mode == SPLINE || mode == FATWORM) {
          alpha_helix_width = params.GetFloat("alpha_helix_width");
          helix_tube_diameter =  params.GetFloat("helix_tube_diameter");
          beta_sheet_width = params.GetFloat("alpha_helix_width");
          //beta_sheet_width = params.GetFloat("beta_sheet_width");
        } else {
          alpha_helix_width = params.GetFloat("worm_width");
          beta_sheet_width = params.GetFloat("worm_width");
        }
        if(fabs(arrow_length)<1e-2) arrow_length = 1.0;
        if(arrow_width<beta_sheet_width) arrow_width = beta_sheet_width;
        if(colours_vec.size()>1){
          colours_vec.pop_back();
          colours_vec.push_back(colours_vec.back());
        }
        if((*secstr_iter)[1]==ALPHA&&cartesians.size()>4){
          if (mode == SPLINE||mode==WORM) {
	    Ribbon *ribbon;
            double *col_tmp = RGBReps::GetColourP((*secstr_iter)[1]+1);
            if(mode==WORM)
	      ribbon = new Ribbon(cartesians,n1_cartesians,n2_cartesians,col_tmp,cartesians[0],colours_vec,worm_width,alpha_helix_width,worm_width,1.0,ribbon_accu*2,spline_accu);
            else
	      ribbon = new Ribbon(cartesians,n1_cartesians,n2_cartesians,col_tmp,cartesians[0],colours_vec,worm_width,alpha_helix_width,worm_width,1.0,-((two_colour_ribbon<<20)|(helix_style<<16)|ribbon_accu*2),spline_accu);
	    obj.add_primitive(ribbon);
            delete [] col_tmp;
          }else{
            double *col_tmp = RGBReps::GetColourP((*secstr_iter)[1]+1);
            /* This is cylinder rep of alpha helices !! */
            std::vector<Cartesian> lead_in(cartesians.begin(),cartesians.begin()+1*spline_accu);
            std::vector<Cartesian> n1_lead_in(n1_cartesians.begin(),n1_cartesians.begin()+1*spline_accu);
            std::vector<Cartesian> n2_lead_in(n2_cartesians.begin(),n2_cartesians.begin()+1*spline_accu);
            std::vector<Cartesian> cols_lead_in(colours_vec.begin(),colours_vec.begin()+1*spline_accu);
	    Worm *ribbon = new Worm(lead_in,n1_lead_in,n2_lead_in,col_tmp,cartesians[0],cols_lead_in,worm_width,worm_width,worm_width,1.0,ribbon_accu*2,spline_accu);
	    obj.add_primitive(ribbon);
            std::vector<Cartesian> main(cartesians.begin()+1*spline_accu-1,cartesians.end()-2*spline_accu+1);
            std::vector<Cartesian> n1_main(n1_cartesians.begin()+1*spline_accu-1,n1_cartesians.end()-2*spline_accu+1);
            std::vector<Cartesian> n2_main(n2_cartesians.begin()+1*spline_accu-1,n2_cartesians.end()-2*spline_accu+1);
            std::vector<Cartesian> cols_main(colours_vec.begin()+1*spline_accu-1,colours_vec.end()-2*spline_accu+1);
	    ribbon = new Worm(main,n1_main,n2_main,col_tmp,cartesians[0],cols_main,helix_tube_diameter,helix_tube_diameter,helix_tube_diameter,1.0,ribbon_accu*2,spline_accu);
	    obj.add_primitive(ribbon);
            std::vector<Cartesian> lead_out(cartesians.end()-2*spline_accu,cartesians.end());
            std::vector<Cartesian> n1_lead_out(n1_cartesians.end()-2*spline_accu,n1_cartesians.end());
            std::vector<Cartesian> n2_lead_out(n2_cartesians.end()-2*spline_accu,n2_cartesians.end());
            std::vector<Cartesian> cols_lead_out(colours_vec.end()-2*spline_accu,colours_vec.end());
	    ribbon = new Worm(lead_out,n1_lead_out,n2_lead_out,col_tmp,cartesians[0],cols_lead_out,worm_width,worm_width,worm_width,1.0,ribbon_accu*2,spline_accu);
	    obj.add_primitive(ribbon);
            delete [] col_tmp;
          }
        }else if((*secstr_iter)[1]==BETA&&cartesians.size()>4){
          double *col_tmp = RGBReps::GetColourP((*secstr_iter)[1]+1);
	  Ribbon *ribbon;
          if (mode == SPLINE|| mode == FATWORM) {
            if(ribbon_style==0)
	      ribbon = new ArrowHeadRibbon(cartesians,n1_cartesians,n2_cartesians,col_tmp,cartesians[0],colours_vec,worm_width,beta_sheet_width,worm_width,arrow_length,arrow_width,1.0,ribbon_accu*2,spline_accu);
            else
	      ribbon = new ArrowHeadRibbon(cartesians,n1_cartesians,n2_cartesians,col_tmp,cartesians[0],colours_vec,worm_width,beta_sheet_width,worm_width,arrow_length,arrow_width,1.0,-((ribbon_style<<16)|ribbon_accu*2),spline_accu);
            obj.add_primitive(ribbon);
          } else {
	    ribbon = new Ribbon(cartesians,n1_cartesians,n2_cartesians,col_tmp,cartesians[0],colours_vec,worm_width,beta_sheet_width,worm_width,1.0,ribbon_accu*2,spline_accu);
   	   obj.add_primitive(ribbon);
	  }
          delete [] col_tmp;
        }else{
          double *col_tmp = RGBReps::GetColourP((*secstr_iter)[1]+1);
	  Worm *ribbon = new Worm(cartesians,n1_cartesians,n2_cartesians,col_tmp,cartesians[0],colours_vec,worm_width,worm_width,worm_width,1.0,ribbon_accu*2,spline_accu);
	  obj.add_primitive(ribbon);
          delete [] col_tmp;
        }
      }
      cartesians.clear();
      n1_cartesians.clear();
      n2_cartesians.clear();
      colours_vec.clear();
      //std::cout << (*secstr_iter)[0]*params.spline_accu << " " << totalpoints << "\n";
      secstr_iter++;
    }
    splines_iter++;
    n1_splines_iter++;
    n2_splines_iter++;
    if(multicolour) colour_vecs++;
  }

 //std::cout << "done build_spline" << std::endl;
}

ConnectivityDraw::ConnectivityDraw(){
}

void ConnectivityDraw::SetParametersAndCalculate(const Connectivity &connectivity_in, PCMMANManager molhnd_in, Displayobject &obj, int mode, const CParamsManager &params, const CParamsManager &global_params, int nSelAtoms, AtomColourVector *atom_colour_vector, const std::vector<double> &atomRadii, const std::string &texture, const std::string &bumpmap, int stick_colour, int side_to_ribbon, int side_to_worm, int bonds_mode){
  molhnd = molhnd_in;
  connectivity = connectivity_in;
  RedrawPrimitives(obj,mode,params,global_params,nSelAtoms,atom_colour_vector,atomRadii,texture,bumpmap,stick_colour,side_to_ribbon,side_to_worm,bonds_mode);
}

ConnectivityDraw::ConnectivityDraw(const Connectivity &connectivity_in, PCMMANManager molhnd_in, Displayobject &obj, int mode, const CParamsManager &params, const CParamsManager &global_params, int nSelAtoms, AtomColourVector *atom_colour_vector, const std::vector<double> &atomRadii, const std::string &texture, const std::string &bumpmap, int stick_colour, int side_to_ribbon, int side_to_worm,int bonds_mode ){

  molhnd = molhnd_in;
  connectivity = connectivity_in;
  RedrawPrimitives(obj,mode,params,global_params,nSelAtoms,atom_colour_vector,atomRadii,texture,bumpmap,stick_colour,side_to_ribbon,side_to_worm,bonds_mode);
}

void ConnectivityDraw::RedrawPrimitives(Displayobject &obj, int mode, const CParamsManager &params, const CParamsManager &global_params, int nSelAtoms,  AtomColourVector *atom_colour_vector, const std::vector<double> &atomRadii, const std::string &texture, const std::string &bumpmap, int stick_colour , int side_to_ribbon, int side_to_worm, int bonds_mode){

  double width=1.0;
  LineCollection *lines = new LineCollection();
  PolyCollection *polys = new PolyCollection();
  bool warning = false;

  double spheres_size=1.0;
  double cylinders_size=params.GetFloat("cylinder_width");
  int spheres_accu=1,cylinders_accu=4;
  double midpoint_frac0 = 0.5;
  double midpoint_frac1 = 0.5;
  std::vector <Cartesian> pyramid_carts;

  bool dashed = params.GetInt("dashed_bonds"); 
  double dash_length = params.GetFloat("dashed_bond_length"); 

  //std::cout << "RedrawPrimitives " << bonds_mode << " " << side_to_ribbon << std::endl;

  std::vector <int> catom_indices = connectivity.GetCAtomIndex();

  if(mode==SPHERES){
    spheres_size = 1.2;
    spheres_accu = global_params.GetInt("solid_quality"); 
    cylinders_accu = 8;
  }
  if(mode==BALLSTICK){
    spheres_size = 0.6;
    cylinders_size = params.GetFloat("ballstick_stick");
    spheres_accu = global_params.GetInt("solid_quality");
    cylinders_accu = 8+8*global_params.GetInt("solid_quality");
  }
  if(mode==CYLINDERS){
    spheres_size = cylinders_size;
    spheres_accu = global_params.GetInt("solid_quality");
    cylinders_accu = 8+8*global_params.GetInt("solid_quality");
  }
  if(mode==BONDS)
    width = params.GetInt("bond_width");
  else if (mode==FATBONDS)
    width = params.GetInt("fat_bond_width");
  else if (mode==PYRAMIDS) {
    float pyra_size = params.GetFloat("pyramid_size");
    pyramid_carts.push_back(Cartesian(2*pyra_size,pyra_size+0.001,0.001));
    pyramid_carts.push_back(Cartesian(-2*pyra_size+0.001,pyra_size,0.001));
    pyramid_carts.push_back(Cartesian(0.001,-pyra_size,2*pyra_size));
    pyramid_carts.push_back(Cartesian(0.001,-pyra_size,-2*pyra_size+0.001));
    pyramid_carts.push_back(Cartesian(2*pyra_size,pyra_size,0.001));
    pyramid_carts.push_back(Cartesian(0.001,-pyra_size+0.001,2*pyra_size));
  } else
    width = params.GetInt("thin_bond_width");

  std::vector<std::vector<int> > conn_lists = connectivity.GetConnectivityLists();
  std::vector<std::vector<int> > ext_conn_lists = connectivity.GetExternalConnectivityLists();
  PPCAtom atoms = connectivity.GetAtoms();
  int natoms = connectivity.GetNumberOfAtoms();
  std::vector <Cartesian> int_carts = CartesiansFromAtoms(atoms,natoms);
  if(atoms) delete [] atoms;
  std::vector<std::vector <Cartesian> > ext_carts = GetExternalCartesians(molhnd,ext_conn_lists,side_to_ribbon,side_to_worm);

  std::vector<Cartesian> carts(2);

  Cartesian tmp_v;

  double *stick_colour_array=0;
  double *colour_array=0;
  double *colour_array2=0;
  //double colour_array[4];
  //double colour_array2[4];
  if ( stick_colour > 0 ) stick_colour_array = RGBReps::GetColourP(stick_colour);
  if (side_to_ribbon>0 || side_to_worm>0 ) {
    midpoint_frac0 = 0.0;
    midpoint_frac1 = 1.0;
    //std::cout << "Setting midpoint to external cart\n";
  }

  if (bonds_mode != DRAW_INTERNAL_BONDS) {

  for(unsigned i=0;i<ext_conn_lists.size();i++){
    colour_array = atom_colour_vector->GetRGB(i);
    
    for(unsigned j=0;j<ext_conn_lists[i].size();j++){
      if(dashed==true&&(mode==BONDS||mode==FATBONDS||mode==THINBONDS)){
        Cartesian midpoint =midpoint_frac0 *int_carts[i]+midpoint_frac1*ext_carts[i][j];
        carts[1] = int_carts[i];
        carts[0] = midpoint;
        tmp_v = carts[0] -  carts[1];
        if (tmp_v.length() < 0.001 ) {
          if (!warning) {
	    warning = true;
            std::cout << "Error drawing BONDS ext\n" ;
          }
        }
        DashLineElement *line = new DashLineElement(carts,colour_array,carts[1],width,1.0);
        line->SetColourOverride(catom_indices[i]);
        line->SetDashLength(dash_length);
        lines->add_primitive(line);
      }
      if(dashed==false&&(mode==BONDS||mode==FATBONDS||mode==THINBONDS)){
        Cartesian midpoint =midpoint_frac0 *int_carts[i]+midpoint_frac1*ext_carts[i][j];
        carts[0] = int_carts[i];
        carts[1] = midpoint;
        tmp_v = carts[0] -  carts[1];
        if (tmp_v.length() < 0.001 ) {
          if (!warning) {
	    warning = true;
            std::cout << "Error drawing BONDS ext\n" ;
          }
        }
        LineElement *line = new LineElement(carts,colour_array,carts[0],width,1.0);
        line->SetColourOverride(catom_indices[i]);
        lines->add_primitive(line);
      }
      if(mode==CYLINDERS||mode==BALLSTICK){
        Cartesian midpoint =midpoint_frac0 *int_carts[i]+midpoint_frac1*ext_carts[i][j];
        carts[0] = int_carts[i];
        /*
        if(mode==BALLSTICK){
          Cartesian vec = ext_carts[i][j]-int_carts[i];
          vec.normalize();
          double rad = sqrt(atomRadii[i]*atomRadii[i]-cylinders_size*cylinders_size);
          carts[0] = int_carts[i] +rad*vec;
        }
        */
        carts[1] = midpoint;
        tmp_v = carts[0] -  carts[1];
        if (tmp_v.length() < 0.001 ) {
          if ( !warning ) {
            warning = true;
            std::cout << "Error drawing CYLINDER ext\n";
          }
        }
        if (stick_colour > 0 && mode==BALLSTICK && stick_colour_array ) {
          CylinderElement *line = new CylinderElement(carts,stick_colour_array,carts[0],cylinders_size,1.0,cylinders_accu);
          line->SetColourOverride(catom_indices[i]);
          polys->add_primitive(line);
        } else {
          CylinderElement *line = new CylinderElement(carts,colour_array,carts[0],cylinders_size,1.0,cylinders_accu);
          line->SetColourOverride(catom_indices[i]);
          polys->add_primitive(line);
        }
      }
    }
    if(colour_array) delete [] colour_array;
  }
  }

  if (bonds_mode != DRAW_EXTERNAL_BONDS) {

  for(unsigned i=0;i<conn_lists.size();i++){
    colour_array = atom_colour_vector->GetRGB(i);
    if(mode==SPHERES||mode==BALLSTICK){
      if (atomRadii[i]<0.01) {
	//std::cout << "SPHERE radii";
      } else {
      SphereElement *sphere = new SphereElement(int_carts[i],colour_array,int_carts[i],atomRadii[i],1.0,spheres_accu);
      sphere->SetColourOverride(catom_indices[i]);
      polys->add_primitive(sphere);
      }
    }
    if(mode==CYLINDERS){
        //if (spheres_size<0.01) std::cout << "CYLINDER radii";
        SphereElement *sphere = new SphereElement(int_carts[i],colour_array,int_carts[i],spheres_size,1.0,spheres_accu);
        sphere->SetColourOverride(catom_indices[i]);
        polys->add_primitive(sphere);
    }
    if(mode==PYRAMIDS){
      for (unsigned j=0;j<5;j++) {     
        carts[0] = int_carts[i]+pyramid_carts[j];
        carts[1] = int_carts[i]+pyramid_carts[j+1];
        LineElement *line = new LineElement(carts,colour_array,carts[1],width,1.0);
        line->SetColourOverride(catom_indices[i]);
        lines->add_primitive(line);
      }
      carts[0] = int_carts[i]+pyramid_carts[1];
      carts[1] = int_carts[i]+pyramid_carts[3];
      LineElement *line = new LineElement(carts,colour_array,carts[1],width,1.0);
      line->SetColourOverride(catom_indices[i]);
      lines->add_primitive(line);
    }

    for(unsigned j=0;j<conn_lists[i].size();j++){
      colour_array2 = atom_colour_vector->GetRGB(conn_lists[i][j]);
      if(dashed==true&&(mode==BONDS||mode==FATBONDS||mode==THINBONDS)){
        //Cartesian midpoint = Cartesian::MidPoint(int_carts[i],int_carts[conn_lists[i][j]]);
        Cartesian midpoint = 0.50*int_carts[i]+0.50*int_carts[conn_lists[i][j]];
        carts[1] = int_carts[i];
        carts[0] = midpoint;
        tmp_v = carts[0] -  carts[1];
        if (tmp_v.length() < 0.001 ) {
          if (!warning) {
            warning = true;
            std::cout << "Error drawing BONDS \n";
          }
        }
        DashLineElement *line = new DashLineElement(carts,colour_array,carts[1],width,1.0);
        line->SetColourOverride(catom_indices[i]);
        line->SetDashLength(dash_length);
        lines->add_primitive(line);
        /*
        carts[0] = midpoint;
        carts[1] = int_carts[conn_lists[i][j]];
        midpoint = 0.50*int_carts[i]+0.50*int_carts[conn_lists[i][j]];
        line = new DashLineElement(carts,colour_array2,carts[1],width,1.0);
        line->SetColourOverride(catom_indices[conn_lists[i][j]]);
        lines->add_primitive(line);
        */
      }
      if(dashed==false&&(mode==BONDS||mode==FATBONDS||mode==THINBONDS)){
        //Cartesian midpoint = Cartesian::MidPoint(int_carts[i],int_carts[conn_lists[i][j]]);
        Cartesian midpoint = 0.50*int_carts[i]+0.50*int_carts[conn_lists[i][j]];
        carts[0] = int_carts[i];
        carts[1] = midpoint;
        tmp_v = carts[0] -  carts[1];
        if (tmp_v.length() < 0.001 ) {
          if (!warning) {
            warning = true;
            std::cout << "Error drawing BONDS \n";
          }
        }
        LineElement *line = new LineElement(carts,colour_array,carts[0],width,1.0);
        line->SetColourOverride(catom_indices[i]);
        lines->add_primitive(line);
        /*
        carts[0] = midpoint;
        carts[1] = int_carts[conn_lists[i][j]];
        midpoint = 0.50*int_carts[i]+0.50*int_carts[conn_lists[i][j]];
        line = new LineElement(carts,colour_array2,carts[1],width,1.0);
        line->SetColourOverride(catom_indices[conn_lists[i][j]]);
        lines->add_primitive(line);
       */
      }
      if(mode==CYLINDERS||mode==BALLSTICK){
        Cartesian midpoint = Cartesian::MidPoint(int_carts[i],int_carts[conn_lists[i][j]]);
        carts[0] = int_carts[i];
        /*
        if(mode==BALLSTICK){
          Cartesian vec = int_carts[conn_lists[i][j]]-int_carts[i];
          vec.normalize();
          double rad = sqrt(atomRadii[i]*atomRadii[i]-cylinders_size*cylinders_size);
          carts[0] = int_carts[i] +rad*vec;
        }
        */
        carts[1] = midpoint;
        tmp_v = carts[0] -  carts[1];
        if (tmp_v.length() < 0.001 ) {
	  if (!warning) {
	    warning = true;
            std::cout << "Error drawing CYLINDER \n";
          }
        }
        if (stick_colour > 0 && mode==BALLSTICK ) {
          CylinderElement *line = new CylinderElement(carts,stick_colour_array,carts[0],cylinders_size,1.0,cylinders_accu);
          polys->add_primitive(line);
        } else {
          CylinderElement *line = new CylinderElement(carts,colour_array,carts[0],cylinders_size,1.0,cylinders_accu);
          line->SetColourOverride(catom_indices[i]);
          polys->add_primitive(line);
        }
      }
      if(colour_array2) delete [] colour_array2;
    }
    if(conn_lists[i].size()==0 && ext_conn_lists[i].size()==0){
      if(mode==BONDS||mode==FATBONDS||mode==THINBONDS){
        Cartesian p1 = int_carts[i];
        p1.set_x(p1.get_x()-0.2);
        Cartesian p2 = int_carts[i];
        p2.set_x(p2.get_x()+0.2);
        Cartesian origin = int_carts[i];
        carts[0] = p1; 
        carts[1] = p2;
        LineElement* xaxis = new LineElement(carts,colour_array,origin,width,1.0);
        p1 = int_carts[i];
        p1.set_y(p1.get_y()-0.2);
        p2 = int_carts[i];
        p2.set_y(p2.get_y()+0.2);
        carts[0] = p1; 
        carts[1] = p2;
        LineElement* yaxis = new LineElement(carts,colour_array,origin,width,1.0);
        p1 = int_carts[i];
        p1.set_z(p1.get_z()-0.2);
        p2 = int_carts[i];
        p2.set_z(p2.get_z()+0.2);
        carts[0] = p1; 
        carts[1] = p2;
        LineElement* zaxis = new LineElement(carts,colour_array,origin,width,1.0);
        lines->add_primitive(xaxis);
        lines->add_primitive(yaxis);
        lines->add_primitive(zaxis);
      }
    }
    if(colour_array) delete [] colour_array;
  }
  }

  lines->SetSize(width);
  obj.add_primitive(lines);
  obj.add_primitive(polys);
  if ( stick_colour > 0 ) delete [] stick_colour_array;

}

std::vector<Cartesian> GetLineThroughBasePairs(PCResidue res1, PCResidue res2){
  int natoms1;
  PPCAtom atoms1=0;
  res1->GetAtomTable1(atoms1,natoms1);
  int natoms2;
  PPCAtom atoms2=0;
  res2->GetAtomTable1(atoms2,natoms2);

  std::vector<Cartesian> N1;
  std::vector<Cartesian> O1;
  std::vector<Cartesian> N2;
  std::vector<Cartesian> O2;

  for(int i=0;i<natoms1;i++){
    if(!strncmp(atoms1[i]->element," N",2)) N1.push_back(Cartesian(atoms1[i]->x,atoms1[i]->y,atoms1[i]->z));
    if(!strncmp(atoms1[i]->element," O",2)) O1.push_back(Cartesian(atoms1[i]->x,atoms1[i]->y,atoms1[i]->z));
  }
  for(int i=0;i<natoms2;i++){
    if(!strncmp(atoms2[i]->element," N",2)) N2.push_back(Cartesian(atoms2[i]->x,atoms2[i]->y,atoms2[i]->z));
    if(!strncmp(atoms2[i]->element," O",2)) O2.push_back(Cartesian(atoms2[i]->x,atoms2[i]->y,atoms2[i]->z));
  }

  //std::cout << N1.size() << " " << O1.size() << "\n";
  //std::cout << N2.size() << " " << O2.size() << "\n";

  std::vector<Cartesian> vecs;
  Cartesian midpoints(0,0,0);
  //std::cout << res1->GetResName() << "\n";
  std::vector<Cartesian>::iterator O1iter = O1.begin();
  while(O1iter!=O1.end()){
    std::vector<Cartesian>::iterator N2iter = N2.begin();
    while(N2iter!=N2.end()){
      double ONlength = ((*O1iter)-(*N2iter)).length();
      if(ONlength>2.6&&ONlength<3.2){
        //std::cout << "ON? " << ONlength << "\n";
        vecs.push_back((*O1iter)-(*N2iter));
        midpoints += Cartesian::MidPoint((*O1iter),(*N2iter));
      }
      N2iter++;
    }
    O1iter++;
  }
  std::vector<Cartesian>::iterator N1iter = N1.begin();
  while(N1iter!=N1.end()){
    std::vector<Cartesian>::iterator N2iter = N2.begin();
    while(N2iter!=N2.end()){
      double NNlength = ((*N1iter)-(*N2iter)).length();
      if(NNlength>2.6&&NNlength<3.2){
        //std::cout << "NN? " << NNlength << "\n";
        vecs.push_back((*N1iter)-(*N2iter));
	midpoints += Cartesian::MidPoint((*N1iter),(*N2iter));
      }
      N2iter++;
    }
    std::vector<Cartesian>::iterator O2iter = O2.begin();
    while(O2iter!=O2.end()){
      double NOlength = ((*N1iter)-(*O2iter)).length();
      if(NOlength>2.6&&NOlength<3.2){
        //std::cout << "NO? " << NOlength << "\n";
        vecs.push_back((*N1iter)-(*O2iter));
	midpoints += Cartesian::MidPoint((*N1iter),(*O2iter));
      }
      O2iter++;
    }
    N1iter++;
  }

  Cartesian v(0,0,0); 
  Cartesian m(0,0,0); 
  std::vector<Cartesian> vm;
  if(vecs.size()>1){
    for(unsigned ii=0;ii<vecs.size();ii++){
      vecs[ii].normalize();
      v += vecs[ii];
      //std::cout << vecs[ii] << "\n";
    }
    v /= vecs.size();
    m = midpoints/vecs.size();
    vm.push_back(v);
    vm.push_back(m);
  }
  //std::cout << v << ", " << m << "\n";
  return vm;
}

Cartesian GetClosestSplinePoint(const std::vector<Cartesian> &carts, const SplineInfo &splineinfo);
Cartesian GetClosestSplinePoint(const Cartesian &cart, const SplineInfo &splineinfo);

std::vector<Cartesian> GetBasePairEnds(PCResidue res1, PCResidue res2, const SplineInfo &splineinfo){
  std::vector<Cartesian> carts(2);
  PCAtom c11 = res1->GetAtom("C1\'");
  PCAtom c21 = res1->GetAtom("C2\'");
  PCAtom c31 = res1->GetAtom("C3\'");
  PCAtom c41 = res1->GetAtom("C4\'");
  PCAtom o41 = res1->GetAtom("O4\'");

  if(!c11) c11 = res1->GetAtom("C1*");
  if(!c21) c21 = res1->GetAtom("C2*");
  if(!c31) c31 = res1->GetAtom("C3*");
  if(!c41) c41 = res1->GetAtom("C4*");
  if(!o41) o41 = res1->GetAtom("O4*");

  PCAtom c12 = res2->GetAtom("C1\'");
  PCAtom c22 = res2->GetAtom("C2\'");
  PCAtom c32 = res2->GetAtom("C3\'");
  PCAtom c42 = res2->GetAtom("C4\'");
  PCAtom o42 = res2->GetAtom("O4\'");

  if(!c11) c12 = res2->GetAtom("C1*");
  if(!c21) c22 = res2->GetAtom("C2*");
  if(!c31) c32 = res2->GetAtom("C3*");
  if(!c41) c42 = res2->GetAtom("C4*");
  if(!o41) o42 = res2->GetAtom("O4*");

  PCAtom c51 = res1->GetAtom("C5\'");
  if(!c51) c51 = res1->GetAtom("C5*");
  PCAtom c52 = res2->GetAtom("C5\'");
  if(!c52) c52 = res2->GetAtom("C5*");

  if(c11&&c21&&c31&&c41&&o41&&c12&&c22&&c32&&c42&&o42){
    std::vector<Cartesian> carts1;
    std::vector<Cartesian> carts2;
    carts1.push_back(Cartesian(c11->x,c11->y,c11->z));
    carts1.push_back(Cartesian(c21->x,c21->y,c21->z));
    carts1.push_back(Cartesian(c31->x,c31->y,c31->z));
    carts1.push_back(Cartesian(c41->x,c41->y,c41->z));
    carts1.push_back(Cartesian(o41->x,o41->y,o41->z));
    carts2.push_back(Cartesian(c12->x,c12->y,c12->z));
    carts2.push_back(Cartesian(c22->x,c22->y,c22->z));
    carts2.push_back(Cartesian(c32->x,c32->y,c32->z));
    carts2.push_back(Cartesian(c42->x,c42->y,c42->z));
    carts2.push_back(Cartesian(o42->x,o42->y,o42->z));
    carts[0] = Cartesian::MidPoint(carts1);
    carts[1] = Cartesian::MidPoint(carts2);

    Cartesian cart1 = carts[0];
    Cartesian cart2 = carts[1];
    if((cart1-cart2).length()>9.0){
       std::vector<Cartesian> carts_tmp(2);
       Cartesian M = Cartesian::MidPoint(cart1,cart2);
       Cartesian c51c(c51->x,c51->y,c51->z);
       Cartesian c52c(c52->x,c52->y,c52->z);

       Cartesian MP = cart1 - M;
       Cartesian MC = c51c - M;

       Cartesian MPnorm = MP;
       MPnorm.normalize();
       double l = Cartesian::DotProduct(MP,MC)/MP.length();
        
       carts_tmp[0] = M;
       carts_tmp[1] = M+l*MPnorm;
       //carts[0] = GetClosestSplinePoint(carts_tmp,splineinfo);
       carts[0] = GetClosestSplinePoint(carts_tmp[1],splineinfo);

       MP = cart2 - M;
       MC = c52c - M;

       MPnorm = MP;
       MPnorm.normalize();
       l = Cartesian::DotProduct(MP,MC)/MP.length();

       carts_tmp[1] = M+l*MPnorm;
       //carts[1] = GetClosestSplinePoint(carts_tmp,splineinfo);
       carts[1] = GetClosestSplinePoint(carts_tmp[1],splineinfo);

    }
  } else {
    carts.clear();
  }

  return carts;
}

void DrawBaseBlock(PolyCollection *polys, PCResidue res1, double *col1, const CParamsManager &params ){

    //float thickness = params.GetFloat("cylinder_width")-0.02;
    float thickness = params.GetFloat("base_block_thickness")-0.02;
    if (thickness<0.02)thickness=0.1;
 
    PCAtom n1 = res1->GetAtom("N1");
    PCAtom c2 = res1->GetAtom("C2");
    PCAtom n3 = res1->GetAtom("N3");
    PCAtom c4 = res1->GetAtom("C4");
    PCAtom c5 = res1->GetAtom("C5");
    PCAtom c6 = res1->GetAtom("C6");

    if(n1&&c2&&n3&&c4&&c5&&c6){
      Cartesian n1cart(n1->x,n1->y,n1->z);
      Cartesian c2cart(c2->x,c2->y,c2->z);
      Cartesian n3cart(n3->x,n3->y,n3->z);
      Cartesian c4cart(c4->x,c4->y,c4->z);
      Cartesian c5cart(c5->x,c5->y,c5->z);
      Cartesian c6cart(c6->x,c6->y,c6->z);

      std::vector <Cartesian> carts;
      Cartesian up = Cartesian::CrossProduct(n1cart-c2cart,c2cart-n3cart);
      up.normalize(thickness);
      carts.push_back(n1cart+up);
      carts.push_back(c2cart+up);
      carts.push_back(n3cart+up);
      carts.push_back(c4cart+up);
      carts.push_back(c5cart+up);
      carts.push_back(c6cart+up);
      carts.push_back(n1cart+up);
      Cartesian midpoint = Cartesian::MidPoint(carts);
      carts.insert(carts.begin(),midpoint);
      TriangleFanElement* polygon = new TriangleFanElement(carts,col1,midpoint,1.0);
      polys->add_primitive(polygon);
      carts.clear();
      carts.push_back(c6cart-up);
      carts.push_back(c5cart-up);
      carts.push_back(c4cart-up);
      carts.push_back(n3cart-up);
      carts.push_back(c2cart-up);
      carts.push_back(n1cart-up);
      carts.push_back(c6cart-up);
      midpoint = Cartesian::MidPoint(carts);
      carts.insert(carts.begin(),midpoint);
      polygon = new TriangleFanElement(carts,col1,midpoint,1.0);
      polys->add_primitive(polygon);
      carts.clear();
      carts.push_back(n1cart-up);
      carts.push_back(c2cart-up);
      carts.push_back(c2cart+up);
      carts.push_back(n1cart+up);
      QuadElement* quad = new QuadElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(quad);
      carts.clear();
      carts.push_back(c2cart-up);
      carts.push_back(n3cart-up);
      carts.push_back(n3cart+up);
      carts.push_back(c2cart+up);
      quad = new QuadElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(quad);
      carts.clear();
      carts.push_back(c5cart-up);
      carts.push_back(c6cart-up);
      carts.push_back(c6cart+up);
      carts.push_back(c5cart+up);
      quad = new QuadElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(quad);
      carts.clear();
      carts.push_back(c6cart-up);
      carts.push_back(n1cart-up);
      carts.push_back(n1cart+up);
      carts.push_back(c6cart+up);
      quad = new QuadElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(quad);
      carts.clear();
      carts.push_back(n3cart-up);
      carts.push_back(c4cart-up);
      carts.push_back(c4cart+up);
      carts.push_back(n3cart+up);
      quad = new QuadElement(carts,col1,Cartesian::MidPoint(carts),1.0);
      polys->add_primitive(quad);

      PCAtom n9 = res1->GetAtom("N9");
      PCAtom c8 = res1->GetAtom("C8");
      PCAtom n7 = res1->GetAtom("N7");
      if(n9&&c8&&n7){
        Cartesian n9cart(n9->x,n9->y,n9->z);
        Cartesian c8cart(c8->x,c8->y,c8->z);
        Cartesian n7cart(n7->x,n7->y,n7->z);
        carts.clear();
        carts.push_back(c5cart+up);
        carts.push_back(c4cart+up);
        carts.push_back(n9cart+up);
        carts.push_back(c8cart+up);
        carts.push_back(n7cart+up);
        carts.push_back(c5cart+up);
        midpoint = Cartesian::MidPoint(carts);
        carts.insert(carts.begin(),midpoint);
        polygon = new TriangleFanElement(carts,col1,midpoint,1.0);
        polys->add_primitive(polygon);
        carts.clear();
        carts.push_back(n7cart-up);
        carts.push_back(c8cart-up);
        carts.push_back(n9cart-up);
        carts.push_back(c4cart-up);
        carts.push_back(c5cart-up);
        carts.push_back(n7cart-up);
        midpoint = Cartesian::MidPoint(carts);
        carts.insert(carts.begin(),midpoint);
        polygon = new TriangleFanElement(carts,col1,midpoint,1.0);
        polys->add_primitive(polygon);
        carts.clear();
        carts.push_back(c4cart-up);
        carts.push_back(n9cart-up);
        carts.push_back(n9cart+up);
        carts.push_back(c4cart+up);
        quad = new QuadElement(carts,col1,Cartesian::MidPoint(carts),1.0);
        polys->add_primitive(quad);
        carts.clear();
        carts.push_back(n9cart-up);
        carts.push_back(c8cart-up);
        carts.push_back(c8cart+up);
        carts.push_back(n9cart+up);
        quad = new QuadElement(carts,col1,Cartesian::MidPoint(carts),1.0);
        polys->add_primitive(quad);
        carts.clear();
        carts.push_back(c8cart-up);
        carts.push_back(n7cart-up);
        carts.push_back(n7cart+up);
        carts.push_back(c8cart+up);
        quad = new QuadElement(carts,col1,Cartesian::MidPoint(carts),1.0);
        polys->add_primitive(quad);
        carts.clear();
        carts.push_back(n7cart-up);
        carts.push_back(c5cart-up);
        carts.push_back(c5cart+up);
        carts.push_back(n7cart+up);
        quad = new QuadElement(carts,col1,Cartesian::MidPoint(carts),1.0);
        polys->add_primitive(quad);
      }else{
        carts.clear();
        carts.push_back(c4cart-up);
        carts.push_back(c5cart-up);
        carts.push_back(c5cart+up);
        carts.push_back(c4cart+up);
        quad = new QuadElement(carts,col1,Cartesian::MidPoint(carts),1.0);
        polys->add_primitive(quad);
      }
    }
}

void DrawBaseBlocks(Displayobject &obj, CMMANManager *molHnd, int selHnd, PPCAtom selAtoms, int nSelAtoms, AtomColourVector *atom_colour_vector,const CParamsManager &params ){
  PolyCollection *polys = new PolyCollection();
  int C5sel = molHnd->NewSelection();
  molHnd->Select(C5sel,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C5*","C","*",SKEY_NEW);
  molHnd->Select(C5sel,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C5'","C","*",SKEY_OR);
  for(int ii=0;ii<nSelAtoms;ii++){
    if(selAtoms[ii]->isInSelection(C5sel)){
      PCResidue res = selAtoms[ii]->GetResidue();
      if(res){
        int restype = molHnd->GetRestypeCode(res);
	//std::cout << res->name << " " << restype << std::endl;
        if(restype==RESTYPE_NUCL||restype==RESTYPE_DNA||restype==RESTYPE_RNA){
           double *col = atom_colour_vector->GetRGB(ii);
           DrawBaseBlock(polys,res,col, params);
           delete [] col;
        }
      }
    }
  }
  polys->SetAlpha(obj.GetAlpha());
  obj.add_primitive(polys);
}


//-----------------------------------------------------------------
void DrawAnisoU(Displayobject &obj,
     int selHndin, PPCAtom selAtoms, int nSelAtoms, 
     AtomColourVector *atom_colour_vector, const std::vector<double> &atomRadii,
     int style, double scale, 
     const CParamsManager &global_params) {
//-----------------------------------------------------------------

  std::string label;
  Spheroid *ellipse;
  double a,b,c;
  Cartesian a_axis,b_axis,c_axis;

  int accu = global_params.GetInt("solid_quality");
  int show_axes=0, show_solid=0;
  if (style == SPHEROID_AXES) show_axes=1;
  if (style == SPHEROID_SOLID) show_solid=1;
  if (style == SPHEROID_SOLIDAXES) {
    show_axes=1;
    show_solid=1;
  }

  PolyCollection *polys = new PolyCollection();
  for ( int i = 0; i < nSelAtoms; i++ ) {
    Cartesian vertex = Cartesian(selAtoms[i]->x, 
                             selAtoms[i]->y,selAtoms[i]->z);
    double *col = atom_colour_vector->GetRGB(i);
    if (selAtoms[i]->WhatIsSet&ASET_Anis_tFac) {
      double array[9] = { 
	  selAtoms[i]->u11,selAtoms[i]->u12, selAtoms[i]->u13, 
	  selAtoms[i]->u12,selAtoms[i]->u22, selAtoms[i]->u23, 
	  selAtoms[i]->u13,selAtoms[i]->u23, selAtoms[i]->u33 };
      matrix U(3,3,array);
      matrix A = atomRadii[i]*scale*U.Cholesky();
      matrix B(4,4);
      B(0,0) = A(0,0);
      B(1,0) = A(1,0);
      B(1,1) = A(1,1);
      B(2,0) = A(2,0);
      B(2,1) = A(2,1);
      B(2,2) = A(2,2);
      B(3,3) = 1.0;
      ellipse = new Spheroid(vertex, col, vertex,
		    B,
         	    1.0, accu , show_axes, show_solid );
      polys->add_primitive(ellipse);
    } else {
      matrix B(4,4);
      double Bfac = atomRadii[i]*scale * selAtoms[i]->tempFactor/(8*M_PI*M_PI) ;
      B(0,0) = B(1,1) = B(2,2) = Bfac;
      B(3,3) = 1.0;
      ellipse = new Spheroid(vertex, col, vertex,
		    B,
         	    1.0, accu , show_axes, show_solid );
      polys->add_primitive(ellipse);
      
    }
    delete [] col;
  }
  obj.add_primitive(polys);

}


Cartesian GetClosestSplinePoint(const Cartesian &cart, const SplineInfo &splineinfo){

  //double min_t = -1.0;
  double min_dist = 1e8;
  int idx = 0;
  int chain = 0;

  for(unsigned j=0;j<splineinfo.nasplines.size();j++){
    for(unsigned i=0;i<splineinfo.nasplines[j].size();i++){
      double dist = (cart-splineinfo.nasplines[j][i]).length();
      if(dist<min_dist&&dist>0.0){
        min_dist = dist;
        idx = i;
        chain = j;
      }
    }
  }
  return splineinfo.nasplines[chain][idx];
}

Cartesian GetClosestSplinePoint(const std::vector<Cartesian> &carts, const SplineInfo &splineinfo){
  Cartesian closest;

  Cartesian ls = carts[0];
  Cartesian le = carts[1];

  int idx = 0;
  int chain = 0;
  int step = 4;

  double min_t = -1.0;
  double min_dist = 1e8;

  for(unsigned j=0;j<splineinfo.nasplines.size();j++){
    for(unsigned i=0;i<splineinfo.nasplines[j].size();i+=step){
      std::vector<double> dist = DistanceBetweenPointAndLine(ls,le,splineinfo.nasplines[j][i]);
      if(dist[0]<min_dist&&dist[1]>0.0){
        min_dist = dist[0];
        idx = i;
        chain = j;
	min_t = dist[1];
      }
    }
  }

  if(idx<step) idx = step;
  if(idx>(int)(splineinfo.nasplines[chain].size())-step-1) idx = splineinfo.nasplines[chain].size()-step-1;

  /*
  if(idx<step){
    std::vector<double> dist = DistanceBetweenPointAndLine(ls,le,splineinfo.nasplines[chain][0]);
    min_dist = dist[0];
    min_t = dist[1];
    idx = 0;
    std::cout << "Beginning\n";
  } else if(idx>splineinfo.nasplines[chain].size()-step-1) {
    std::vector<double> dist = DistanceBetweenPointAndLine(ls,le,splineinfo.nasplines[chain].back());
    min_dist = dist[0];
    min_t = dist[1];
    idx = splineinfo.nasplines[chain].size()-1;
    std::cout << "End (" << splineinfo.nasplines[chain].size() << "\n";
  } else {
  */
    int idx2 = idx;
    for(int i=idx-step;i<idx+step;i++){
       std::vector<double> dist = DistanceBetweenPointAndLine(ls,le,splineinfo.nasplines[chain][i]);
       if(dist[0]<min_dist){
          min_dist = dist[0];
	  idx2 = i;
	  min_t = dist[1];
       }
    }
  //}
  idx = idx2;

  Cartesian MP=le-ls;
  double length = MP.length();
  MP.normalize();
  Cartesian n1 = splineinfo.n1_nasplines[chain][idx];
  n1.normalize(); // Already normalized ?
  double l = min_dist*(Cartesian::DotProduct(n1,MP));
  l = l/length;

  Cartesian sp = splineinfo.nasplines[chain][idx]; 
  Cartesian new_le = ls+min_t*(le-ls);
  Cartesian l1 = new_le - sp;
  double ang_test = Cartesian::DotProduct(l1,n1);

  if(ang_test>0&&(idx<step||idx>(int)(splineinfo.nasplines[chain].size())-step-1)){
    min_t += l;
  }else{
    min_t -= l;
  }

  if(min_t>0.0)
    return ls+min_t*(le-ls);

  return le;
}

void DrawBasePairs(Displayobject &obj, const CNABasePairs &bp, const SplineInfo &splineinfo, const CParamsManager &params){
  std::vector<std::pair<PCResidue,PCResidue> > base_pairs = bp.GetPairs();
  std::vector<std::pair<double*,double*> > colours = bp.GetColours();
  PolyCollection *polys = new PolyCollection();
  double cylinders_size = params.GetFloat("nucleic_stick_width");
  if (cylinders_size<0.02)cylinders_size=0.1;

  int cylinders_accu = 8;

  for(unsigned ii=0;ii<base_pairs.size();ii++){
    PCResidue res1 = base_pairs[ii].first;
    PCResidue res2 = base_pairs[ii].second;
    double *col1 = colours[ii].first;
    double *col2 = colours[ii].second;

    std::vector<Cartesian> cartsend = GetBasePairEnds(res1,res2,splineinfo);
    if(cartsend.size()>0){
      std::vector<Cartesian> carts(2);
      std::vector<Cartesian> carts1(2);
      std::vector<Cartesian> carts2(2);
      Cartesian Mnew = Cartesian::MidPoint(cartsend[0],cartsend[1]);
      carts[0] = Mnew;
      carts[1] = cartsend[0];
      carts2[0] = Mnew;
      carts2[1] = cartsend[1];

      CylinderElement *line = new CylinderElement(carts,col1,carts[1],cylinders_size,1.0,cylinders_accu);
      polys->add_primitive(line);

      line = new CylinderElement(carts2,col2,carts2[1],cylinders_size,1.0,cylinders_accu);
      polys->add_primitive(line);

    }
  }
  polys->SetAlpha(obj.GetAlpha());
  obj.add_primitive(polys);
}

void DrawLipids(Displayobject &obj, const std::vector<MMUTLipid> &lipids, CMMANManager *molHnd, AtomColourVector *atom_colour_vector, const CParamsManager &params, const CParamsManager &global_params){

  clock_t t1 = clock();
  /* Need colours, etc. */

  if(lipids.size()<1) return;

  PolyCollection *polys = new PolyCollection();
  int spheres_accu=1;
  spheres_accu = global_params.GetInt("solid_quality");
  int spline_accu = 1+global_params.GetInt("solid_quality");
  int ribbon_accu = 4+global_params.GetInt("solid_quality")*4;

  //std::cout << "In routine to DRAWLIPIDS\n";

  //std::cout << "ribbon_accu: " << ribbon_accu << "\n";
  //std::cout << "spline_accu: " << spline_accu << "\n";
  for(unsigned ilipid=0;ilipid<lipids.size();ilipid++){
  MMUTLipid lipid = lipids[ilipid];
  std::vector<std::vector<Cartesian> > sorted_tail_carts = lipid.GetTailCartesians();
  std::vector<std::vector<int> > sorted_tail_serNums = lipid.GetTailSerNums();
  std::vector<std::vector<Cartesian> > all_head_carts = lipid.GetHeadCartesians();
  std::vector<std::vector<int> > all_head_serNums = lipid.GetHeadSerNums();

  int selHnd = lipid.GetMainSelectionHandle();
  PPCAtom atomTable_main_selection=0;
  int nAtoms_main_selection;
  molHnd->GetSelIndex ( selHnd, atomTable_main_selection, nAtoms_main_selection);
  //std::cout << nAtoms_main_selection << " atoms in main selection\n";
  if(nAtoms_main_selection==0) return;

  std::vector<int> serNums;
  for(int jj=0;jj<nAtoms_main_selection;jj++){
    serNums.push_back(atomTable_main_selection[jj]->serNum);
  }

  //std::cout << "Heads:\n";
  std::vector<Cartesian> head_centres;
  for(unsigned ii=0;ii<all_head_carts.size();ii++){
     std::vector<int> head_serNums = all_head_serNums[ii];
     std::vector<Cartesian> head_carts = all_head_carts[ii];
     if(head_carts.size()>0){ //Hopefully
       Cartesian centre = Cartesian::MidPoint(head_carts);
       head_centres.push_back(centre);
       int atom_map=0;
       for(int j=0;j<nAtoms_main_selection;j++) {
         if(head_serNums[0]==serNums[j]){
           atom_map = j;
	   //std::cout << atoms[0]->name << " " << atoms[0]->GetResidueNo() << " " << j << "\n";
           break;
         }
       }
       double *col = atom_colour_vector->GetRGB(atom_map);
       if(head_carts.size()>3){
       std::vector<Cartesian> pca = Cartesian::PrincipalComponentAnalysis(head_carts);
         SpheroidElement *spheroid = new SpheroidElement(centre,col,centre,2.0*sqrt(pca[4].get_x())+.7,2.0*sqrt(pca[4].get_y())+.7,2.0*sqrt(pca[4].get_z())+.7,pca[0],pca[1],pca[2],1.0,spheres_accu);
         polys->add_primitive(spheroid);
       }else if(head_carts.size()>1){
         double max_dist = 0.0;
         for(unsigned ih=0;ih<head_carts.size();ih++){
           if(LineLength(head_carts[ih],centre)>max_dist) max_dist = LineLength(head_carts[ih],centre);
         }
         SphereElement *spheroid = new SphereElement(centre,col,centre,max_dist+.7,1.0,spheres_accu);
         polys->add_primitive(spheroid);
       }else{
         SphereElement *spheroid = new SphereElement(centre,col,centre,1.0,1.0,spheres_accu);
         polys->add_primitive(spheroid);
       }
       delete [] col;
     //delete [] atoms;
     }
  }

  //std::cout << "Tails:\n";
  for(unsigned ii=0;ii<sorted_tail_serNums.size();ii++){
     std::vector<Cartesian> tail_carts = sorted_tail_carts[ii];
     std::vector<int> tail_serNums = sorted_tail_serNums[ii];
     std::vector<int> atom_map;
     // This looks slow, any better way? (But it isn't slow compared to mmdb stuff above).
     for(unsigned i=0;i<tail_serNums.size();i++) {
       int serNum = tail_serNums[i];
       for(int j=0;j<nAtoms_main_selection;j++) {
	  if(serNum==serNums[j]){
            atom_map.push_back(j);
	  }
       }
     }

     if(tail_carts.size()>3) {
       // Drawing from hereonin about ??% time
       int head_offset = 0;
       if(head_centres.size()>0){
	  head_offset = -1;
	  Cartesian posn = tail_carts[0];
	  Cartesian extra_posn = head_centres[0];
	  int first_t = 1;
	  double min_dist = LineLength(extra_posn,tail_carts[0]);
	  for(unsigned ic=1;ic<head_centres.size();ic++){
            double dist = LineLength(head_centres[ic],tail_carts[0]);
            if(dist<min_dist){
              dist = min_dist;             
	      extra_posn = head_centres[ic];
            }
	  }
	  for(unsigned ic=0;ic<head_centres.size();ic++){
            double dist = LineLength(head_centres[ic],tail_carts.back());
            if(dist<min_dist){
              dist = min_dist;             
	      extra_posn = head_centres[ic];
	      posn = tail_carts.back();
	      first_t = 0;
            }
	  }
	  Cartesian dirn = extra_posn-posn;
	  dirn.normalize();
	  if(first_t){
	    double wanted_length = LineLength(tail_carts[1],posn);
            tail_carts.insert(tail_carts.begin(),tail_carts[0]+wanted_length*dirn);
	  }else{
	    head_offset = 1;
	    double wanted_length = LineLength(tail_carts[tail_carts.size()-2],posn);
            tail_carts.insert(tail_carts.end(),tail_carts.back()+wanted_length*dirn);
	  }
       }

       std::vector<Cartesian> tail_n1;
       std::vector<Cartesian>::iterator t_iter = tail_carts.begin()+1;
       while(t_iter!=tail_carts.end()-1){
         tail_n1.push_back(Cartesian::CrossProduct((*(t_iter)-*(t_iter-1)),(*(t_iter+1)-*(t_iter))));
	 tail_n1.back().normalize();
         t_iter++;
       }
       tail_n1.push_back(tail_n1.back());
       tail_n1.insert(tail_n1.begin(),tail_n1[0]);

       t_iter = tail_n1.begin()+1;
       while(t_iter!=tail_n1.end()){
         if(Cartesian::DotProduct(*t_iter,*(t_iter-1))<0.0){
           *(t_iter) = -(*t_iter);
         }
         t_iter++;
       }
       //std::cout << "Size of spline inputs: " << tail_carts.size() << " " << tail_n1.size() << "\n";
       std::cout.flush();
       // And this is what we draw.
       std::vector<Cartesian> tail_spline = BezierCurve(tail_carts,spline_accu);
       std::vector<Cartesian> tail_n1_spline = BezierCurve(tail_n1,spline_accu);
       std::vector<Cartesian> tail_n2_spline;
       std::vector<Cartesian> colour_vector;
       double *col = atom_colour_vector->GetRGB(serNums[0]-1);
       if(head_offset>0){
         for(int iat=0;iat<(int)tail_carts.size()-1;iat++){
           double *atcol = atom_colour_vector->GetRGB(atom_map[iat]);
	   for(int iacc=0;iacc<spline_accu;iacc++)
             colour_vector.push_back(Cartesian(atcol));
           delete [] atcol;
         }
         double *atcol = atom_colour_vector->GetRGB(atom_map[tail_serNums.size()-1]);
	 for(int iacc=0;iacc<spline_accu;iacc++)
           colour_vector.push_back(Cartesian(atcol));
         delete [] atcol;
       } else {
         if(head_offset<0){
           double *atcol = atom_colour_vector->GetRGB(atom_map[0]);
	   for(int iacc=0;iacc<spline_accu;iacc++)
             colour_vector.push_back(Cartesian(atcol));
           delete [] atcol;
         }
         for(int iat=-head_offset;iat<(int)tail_carts.size();iat++){
           double *atcol = atom_colour_vector->GetRGB(atom_map[iat+head_offset]);
	   for(int iacc=0;iacc<spline_accu;iacc++)
             colour_vector.push_back(Cartesian(atcol));
           delete [] atcol;
         }
       }
       delete [] col;
       for(int iacc=0;iacc<spline_accu;iacc++)
         colour_vector.push_back(colour_vector.back());
       for(unsigned ispl=0;ispl<tail_spline.size()-1;ispl++){
	 Cartesian ts = tail_spline[ispl]-tail_spline[ispl+1];
	 Cartesian tn1s = tail_n1_spline[ispl];
	 //std::cout << Cartesian::DotProduct(ts,tn1s) << "\n";
	 ts.normalize();
	 tn1s.normalize();
         tail_n2_spline.push_back(Cartesian::CrossProduct(ts,tn1s));
	 tail_n2_spline.back().normalize();
       }
       tail_n2_spline.push_back(tail_n2_spline.back());
       float worm_width = params.GetFloat("worm_width");
       Ribbon *ribbon = new Ribbon(tail_spline,tail_n1_spline,tail_n2_spline,col,tail_spline[0],colour_vector,worm_width,worm_width,worm_width,1.0,ribbon_accu*2,spline_accu);
       polys->add_primitive(ribbon);
     }
  }
  }
  obj.add_primitive(polys);
  clock_t t2 = clock();
  std::cout << "Time for lipid draw " << ((t2-t1)*1000.0/CLOCKS_PER_SEC)/1000.0<< "\n";
}

