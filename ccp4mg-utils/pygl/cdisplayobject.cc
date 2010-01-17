/*
     pygl/cdisplayobject.cc: CCP4MG Molecular Graphics Program
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
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include <ctype.h>
#include <math.h>

#include "quat.h"
#include "cprimitive.h"
#include "cdisplayobject.h"
#include <stdio.h>
#include "cartesian.h"
#include "help.h"
#include "volume.h"
#include "matrix.h"
#include <vector>

#include "rgbreps.h"

#ifdef __APPLE_CC__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#ifndef M_PI
#define M_PI 3.141592653589793238462643
#define PIBY2 (M_PI * 2)
#endif

struct ZSortPrimitive {
  Primitive *prim;
  double z;
} ZSortPrimitive;

class sort_primitives{
  public:
    int operator()(const struct ZSortPrimitive &p1, const struct ZSortPrimitive &p2) const { 
      return p1.z < p2.z;
    }
};

const std::vector<Primitive*> &Displayobject::GetPrimitives() const {
  return prims;
}
const std::vector<Primitive*> &Displayobject::GetSurfacePrimitives() const {
  return surf_prims;
}
const std::vector<BillBoard*> &Displayobject::GetImagePrimitives() const {
  return image_prims;
}
const std::vector<SimpleText*> &Displayobject::GetTextPrimitives() const {
  return text_prims;
}

void DrawSortedTransparentPrimitives(const std::vector<Displayobject> &objs, int acsize, double xoff, double yoff, std::vector<std::vector<double> > jarray, std::vector<Cartesian> axes, bool antialias, bool rebuilt){

  clock_t t1 = clock();
  Volume clip_vol = GetFrontAndBackClippingPlanes();

  glEnable(GL_LIGHTING);
  glEnable(GL_BLEND);
  glDepthMask(GL_FALSE);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  GLdouble modelMatrix[16];
  glGetDoublev(GL_MODELVIEW_MATRIX,modelMatrix);

  matrix modelmat = matrix(4,4,modelMatrix);

  modelmat = modelmat.Transpose();

  std::vector<struct ZSortPrimitive> transformed_prims;

  std::vector<Displayobject>::const_iterator obj_iter = objs.begin();
  static std::vector<std::vector<Primitive*> > simple_prims;
  int i=0;

  //clock_t t2 = clock();
  //std::cout << "Time for setup " << ((t2-t1)*1000.0/CLOCKS_PER_SEC)/1000.0<< "\n";

  clock_t t3 = clock();
  if(rebuilt){
  std::vector<std::vector<Primitive*> >::iterator prims_iter = simple_prims.begin();
  while(prims_iter!=simple_prims.end()){
    std::vector<Primitive*>::iterator prim_iter = prims_iter->begin();
    while(prim_iter!=prims_iter->end()){
      delete *prim_iter;
      prim_iter++;
    }
    //delete *prims_iter;
    prims_iter++;
  }
  simple_prims.clear();
  t3 = clock();
  //std::cout << "Time for delete " << ((t3-t2)*1000.0/CLOCKS_PER_SEC)/1000.0<< "\n";
  while(obj_iter!=objs.end()){
    simple_prims.push_back(std::vector<Primitive*>(0));
    if(obj_iter->visible){
      std::vector<Primitive*> obj_prims=obj_iter->GetPrimitives();
      obj_prims.insert(obj_prims.begin(),obj_iter->GetSurfacePrimitives().begin(),obj_iter->GetSurfacePrimitives().end());
      matrix objrotmat = (*obj_iter).quat.getMatrix();
      std::vector<Primitive*>::const_iterator prim_iter = obj_prims.begin();

      while(prim_iter!=obj_prims.end()){
        std::vector<Primitive*> this_simple_prims = (*prim_iter)->GetSimplePrimitives(clip_vol,objrotmat,obj_iter->origin);
	if(this_simple_prims.size()>0){
          simple_prims[i].insert(simple_prims[i].end(),this_simple_prims.begin(),this_simple_prims.end());
        }
        prim_iter++;
      }
    }
    obj_iter++; i++;
  }
  }
  clock_t t4 = clock();
  std::cout << "Time for GetPrimitives" << ((t4-t3)*1000.0/CLOCKS_PER_SEC)/1000.0<< "\n";
  i=0;
  obj_iter = objs.begin();
  while(obj_iter!=objs.end()){
    if(obj_iter->visible){
      matrix objrotmat = (*obj_iter).quat.getMatrix();
      std::vector<Primitive*>::const_iterator prim_iter = simple_prims[i].begin();
      Cartesian obj_origin = Cartesian(obj_iter->get_origin()) ;
      while(prim_iter!=simple_prims[i].end()){
        Cartesian transrotprim = objrotmat*((*prim_iter)->get_origin()) + obj_origin;
        double z = (modelmat*transrotprim).get_z();
        struct ZSortPrimitive sort_prim = { *prim_iter, z};
        transformed_prims.push_back(sort_prim);
        prim_iter++;
      }
    }
    obj_iter++; i++;
  }

  clock_t t5 = clock();
  std::cout << "Time for setup sort" << ((t5-t4)*1000.0/CLOCKS_PER_SEC)/1000.0<< "\n";
  std::sort(transformed_prims.begin(),transformed_prims.end(),sort_primitives());
  clock_t t6 = clock();
  std::cout << "Time for sort" << ((t6-t5)*1000.0/CLOCKS_PER_SEC)/1000.0<< "\n";

  std::vector<struct ZSortPrimitive>::iterator k = transformed_prims.begin();
  glPolygonMode(GL_FRONT,GL_FILL);

  GLfloat white[4] = { 1.0, 1.0, 1.0, 1.0 };
  glMaterialfv(GL_FRONT, GL_SPECULAR, white);

  k = transformed_prims.begin();
  while(k!=transformed_prims.end()){
      // All this state change might be(!) expensive, but what can I do?
      if(k->prim->isLine()) {
        glDisable(GL_LIGHTING);
        glLineWidth(k->prim->GetSize());
        k->prim->set_draw_colour();
        k->prim->draw();
      } else {
        glEnable(GL_LIGHTING);
        k->prim->set_draw_colour();
        k->prim->draw();
      }
      k++;
  }
  clock_t t7 = clock();
  std::cout << "Time for draw" << ((t7-t6)*1000.0/CLOCKS_PER_SEC)/1000.0<< "\n";
  std::cout << "Time for transparency" << ((t7-t1)*1000.0/CLOCKS_PER_SEC)/1000.0<< "\n";
  glDepthMask(GL_TRUE);

}

void Displayobject::ApplySymmetryMatrix(int i){
  if(i>int(GetNumSymmetryMatrices())){
     return;
  }
  double *glrotmat = (GetSymmetryMatrix(i).Transpose()).to_dp();
  glMultMatrixd(glrotmat);
  delete [] glrotmat;
}

void Displayobject::ApplySymmetryMatrix_RotationOnly(int i){
  matrix symmat = GetSymmetryMatrix(i).Transpose();
  symmat(3,0) = 0.0;
  symmat(3,1) = 0.0;
  symmat(3,2) = 0.0; 
  double *glrotmat = symmat.to_dp();
  glMultMatrixd(glrotmat);
  delete [] glrotmat;
}

void Displayobject::ApplyTranslation(){
  glTranslatef(origin.get_x(),origin.get_y(),origin.get_z());
}

void Displayobject::ApplyRotation(){
  double *glrotmat = get_rotation_matrix().to_dp();
  glMultMatrixd(glrotmat);
  delete [] glrotmat;
}

Displayobject::~Displayobject(){
  clear_prims();
  clear_images();
  clear_labels();
}


int Displayobject::reInitializeTextPrims(){
  std::vector<SimpleText*>::const_iterator k = text_prims.begin();
  while(k!=text_prims.end()){
    (*k)->initialize();
    k++;
  }
  return 0;
}

Displayobject::Displayobject(){
  
  visible = 1;
  rot.push_back(0.0);
  rot.push_back(0.0);
  rot.push_back(0.0);
  drot.push_back(0.0);
  drot.push_back(0.0);
  drot.push_back(0.0);
  dx = 0.0;
  dy = 0.0;
  dz = 0.0;
  do_rebuild = 0;
  anchored = 0;
  build_disp_list = 1;
  draw_symm = 0;
  draw_unit_cell = 0;
  symm_diff_colour = 1;
  transparent = 0;
  alpha = 1.0;
}

void Displayobject::SetBuildDisplayList(int i){
  build_disp_list = i;
}

int Displayobject::BuildDisplayList(void){
  return build_disp_list;
}

void Displayobject::changevis(void){
  if(visible==1)
    visible = 0;
  else
    visible = 1;
}

int Displayobject::GetNumberOfTextIDS(void) const{
  return text_prims.size();
}

int* Displayobject::GetTextIDS(void) const{
  int *text_ids = new int[text_prims.size()];
  int i = 0;

  std::vector<SimpleText*>::const_iterator k = text_prims.begin();
  while(k!=text_prims.end()){
    text_ids[i] = (*k)->GetID();
    k++; i++;
  }
  return text_ids;
}

void Displayobject::SetTextFont(
  const std::string family,  const std::string weight, 
  const std::string slant, const int size, const int underline ) {
  std::vector<SimpleText*>::const_iterator k = text_prims.begin();
  while(k!=text_prims.end()){
    (*k)->SetFontFamily(family);
    (*k)->SetFontWeight(weight);
    (*k)->SetFontSlant(slant);
    (*k)->SetFontSize(size);
    (*k)->SetFontUnderline(underline);
    //(*k)->initialize();
    (*k)->renderStringToPixmap(); 
    k++;
  }
}


void Displayobject::DeleteText(void) {
  /* Since Text's destructor doesn't seem to do anything, this all looks
     a little pointless?
   */
  //std::vector<SimpleText*>::iterator text_prim_iter = text_prims.begin();
  //while(text_prim_iter!=text_prims.end()){
    //if(*text_prim_iter) delete *text_prim_iter;
    //text_prim_iter++;
  //}
  text_prims.clear();
}

Quat Displayobject::GetUnitCellAlignmentRotation(const std::string &axis){
  Quat q;
  Cartesian x,y,z;
  bool reverse_rot = false;

  if(unit_cell.size()<3)
    return q;

  if(axis=="a"||axis=="A"){
    x = unit_cell[1];
    y = unit_cell[2];
  }
  else if(axis=="b"||axis=="B"){
    x = unit_cell[2];
    y = unit_cell[0];
    reverse_rot = true;
  }
  else if(axis=="c"||axis=="C"){
    x = unit_cell[0];
    y = unit_cell[1];
  }
  else {
    std::cout << "Unknown Crystallographic axis\n";
    return q;
  }
  x.normalize();
  y.normalize();
  z = x.CrossProduct(x,y);

  Cartesian Z(0,0,1);

  double theta = acos(Cartesian::DotProduct(z,Z))*180.0/M_PI;

  Cartesian rot_axis = Cartesian::CrossProduct(z,Z);

  //std::cout << rot_axis << " " << theta << "\n"; std::cout.flush();

  if(fabs(theta)>1e-7)
    q = Quat(rot_axis,1,theta);

  Cartesian xprim = q.getMatrix() * x;
  Cartesian X(1,0,0);
  
  double theta2 = acos(Cartesian::DotProduct(xprim,X))*180.0/M_PI;
  Cartesian rot_axis2 = Cartesian::CrossProduct(xprim,X);

  if(fabs(theta2)>1e-7){
    if(reverse_rot) theta2 *= -1;
    Quat q2(rot_axis2,1,theta2);
    q.postMult(q2);
  }

  //std::cout << rot_axis2 << " " << theta2 << "\n"; std::cout.flush();

  return q;

}

void Displayobject::DrawUnitCell(){
  if(unit_cell.size()<3||!GetDrawUnitCell())
    return;
  glDisable(GL_LIGHTING);
  GLfloat *params = new GLfloat[4];
  glGetFloatv(GL_COLOR_CLEAR_VALUE,params);
  GLfloat colour[4] = {1.0-params[0], 1.0-params[1], 1.0-params[2], 1.0};
  glColor4fv(colour);
  //glEnable(GL_LINE_SMOOTH);
  //glEnable(GL_BLEND);
  //glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
  glBegin(GL_LINES);
  glVertex3d(0,0,0);
  glVertex3d(unit_cell[0].get_x(), unit_cell[0].get_y(), unit_cell[0].get_z());
  glVertex3d(0,0,0);
  glVertex3d(unit_cell[1].get_x(), unit_cell[1].get_y(), unit_cell[1].get_z());
  glVertex3d(0,0,0);
  glVertex3d(unit_cell[2].get_x(), unit_cell[2].get_y(), unit_cell[2].get_z());
  glVertex3d(unit_cell[1].get_x(), unit_cell[1].get_y(), unit_cell[1].get_z());
  glVertex3d((unit_cell[1]+unit_cell[0]).get_x(), (unit_cell[1]+unit_cell[0]).get_y(), (unit_cell[1]+unit_cell[0]).get_z());
  glVertex3d(unit_cell[0].get_x(), unit_cell[0].get_y(), unit_cell[0].get_z());
  glVertex3d((unit_cell[1]+unit_cell[0]).get_x(), (unit_cell[1]+unit_cell[0]).get_y(), (unit_cell[1]+unit_cell[0]).get_z());
  glVertex3d(unit_cell[2].get_x(), unit_cell[2].get_y(), unit_cell[2].get_z());
  glVertex3d((unit_cell[1]+unit_cell[2]).get_x(), (unit_cell[1]+unit_cell[2]).get_y(), (unit_cell[1]+unit_cell[2]).get_z());
  glVertex3d(unit_cell[1].get_x(), unit_cell[1].get_y(), unit_cell[1].get_z());
  glVertex3d((unit_cell[1]+unit_cell[2]).get_x(), (unit_cell[1]+unit_cell[2]).get_y(), (unit_cell[1]+unit_cell[2]).get_z());
  glVertex3d((unit_cell[1]+unit_cell[2]).get_x(), (unit_cell[1]+unit_cell[2]).get_y(), (unit_cell[1]+unit_cell[2]).get_z());
  glVertex3d((unit_cell[1]+unit_cell[2]+unit_cell[0]).get_x(), (unit_cell[1]+unit_cell[2]+unit_cell[0]).get_y(), (unit_cell[1]+unit_cell[2]+unit_cell[0]).get_z());
  glVertex3d((unit_cell[1]+unit_cell[0]).get_x(), (unit_cell[1]+unit_cell[0]).get_y(), (unit_cell[1]+unit_cell[0]).get_z());
  glVertex3d((unit_cell[1]+unit_cell[2]+unit_cell[0]).get_x(), (unit_cell[1]+unit_cell[2]+unit_cell[0]).get_y(), (unit_cell[1]+unit_cell[2]+unit_cell[0]).get_z());
  glVertex3d((unit_cell[0]+unit_cell[2]).get_x(), (unit_cell[0]+unit_cell[2]).get_y(), (unit_cell[0]+unit_cell[2]).get_z());
  glVertex3d((unit_cell[1]+unit_cell[2]+unit_cell[0]).get_x(), (unit_cell[1]+unit_cell[2]+unit_cell[0]).get_y(), (unit_cell[1]+unit_cell[2]+unit_cell[0]).get_z());
  glVertex3d(unit_cell[0].get_x(), unit_cell[0].get_y(), unit_cell[0].get_z());
  glVertex3d((unit_cell[0]+unit_cell[2]).get_x(), (unit_cell[0]+unit_cell[2]).get_y(), (unit_cell[0]+unit_cell[2]).get_z());
  glVertex3d(unit_cell[2].get_x(), unit_cell[2].get_y(), unit_cell[2].get_z());
  glVertex3d((unit_cell[0]+unit_cell[2]).get_x(), (unit_cell[0]+unit_cell[2]).get_y(), (unit_cell[0]+unit_cell[2]).get_z());

  glEnd();
  //glBlendFunc (GL_SRC_ALPHA, GL_ZERO);
  //glDisable(GL_BLEND);
}

void Displayobject::SetUnitCell(const std::vector<Cartesian> &cell_params){
  unit_cell = cell_params;
}

void Displayobject::SetSymmetryMatrixNumbers(const std::vector<int> &symm_nos_in){
  symm_nos = symm_nos_in;
}

void Displayobject::SetSymmetryMatrices(const std::vector<matrix> &symm_mat_in){
  symm_mat = symm_mat_in;
}

unsigned int Displayobject::GetNumSymmetryMatrices() const{
  return symm_mat.size();
}

matrix Displayobject::GetSymmetryMatrix(int i) const{
  if(i>int(GetNumSymmetryMatrices())){
    return matrix(4,kdelta); // return a unit matrix.
  }
  return symm_mat[i];
}

int Displayobject::GetSymmetryMatrixNumber(int i) const {
  if(i>int(GetNumSymmetryMatrices()))
     return 0;
  return symm_nos[i];
}

void Displayobject::SetTextString(int text_id, const char* new_string){
  std::vector<SimpleText*>::iterator k = text_prims.begin();
  while(k!=text_prims.end()){
    if(text_id==(*k)->GetID()){
       (*k)->SetText(std::string(new_string));
       return;
    }
    k++;
  }
}

void Displayobject::DeleteTextPrimitive(int text_id){
  std::vector<SimpleText*> new_text_prims;
  std::vector<SimpleText*>::iterator k = text_prims.begin();
  while(k!=text_prims.end()){
    if(text_id!=(*k)->GetID()){
       new_text_prims.push_back(*k);
    }
    k++;
  }
  text_prims = new_text_prims;
}

void Displayobject::SetTextString(int text_id, std::string new_string){
  std::vector<SimpleText*>::iterator k = text_prims.begin();
  while(k!=text_prims.end()){
    if(text_id==(*k)->GetID()){
       (*k)->SetText(new_string);
       return;
    }
    k++;
  }
}

const char* Displayobject::GetTextString(int text_id) const{
  std::vector<SimpleText*>::const_iterator k = text_prims.begin();
  while(k!=text_prims.end()){
    if(text_id==(*k)->GetID()){
       return (*k)->GetText().c_str();
    }
    k++;
  }
  return 0;
}

/*
void Displayobject::add_primitive(Primitive &prim){
  Primitive *p = &prim;
  prims.push_back(p);
}

void Displayobject::add_text_primitive(SimpleText &prim){
  SimpleText *p = &prim;
  text_prims.push_back(p);
}
*/

void Displayobject::add_surf_primitive(Primitive *prim){
  surf_prims.push_back(prim);
}

void Displayobject::add_primitive(Primitive *prim){
  prims.push_back(prim);
}

void Displayobject::add_text_primitive(SimpleText *prim){
  text_prims.push_back(prim);
}

void Displayobject::add_image_primitive(BillBoard *prim){
  image_prims.push_back(prim);
}

void Displayobject::clear_images(void){
  /*
  std::vector<SimpleBillBoard*>::iterator image_prim_iter = image_prims.begin();
  while(image_prim_iter!=image_prims.end()){
    if(*image_prim_iter) delete *image_prim_iter;
    image_prim_iter++;
  }
  */
  image_prims.clear();
}

void Displayobject::clear_prims(void){
  std::vector<Primitive*>::iterator prim_iter = prims.begin();
  while(prim_iter!=prims.end()){
    //if(*prim_iter) delete *prim_iter;
    prim_iter++;
  }
  prims.clear();
  /*
    //Oh dear, we do have some object ownership confusion. 
    //We want prims deleted from SurfaceDispobj.py and then we
    //want to use same prim again. Very confusing.
  prim_iter = surf_prims.begin();
  while(prim_iter!=surf_prims.end()){
    if(*prim_iter) delete *prim_iter;
    prim_iter++;
  }
  */
  surf_prims.clear();
}

void Displayobject::increase_shininess(double shininess){
  
  lighting.shininess += shininess;
  if(lighting.shininess < 0.0)
    lighting.shininess = 0.0;
  if(lighting.shininess > 128.0)
    lighting.shininess = 128.0;
  rebuild();
}


void Displayobject::increase_specular(std::vector<double>specular){

  lighting.specular[0] = lighting.specular[0] + specular[0];
  lighting.specular[1] = lighting.specular[1] + specular[1];
  lighting.specular[2] = lighting.specular[2] + specular[2];

  if(lighting.specular[0] < 0.0)
    lighting.specular[0] = 0.0;
  if(lighting.specular[0] > 1.0)
    lighting.specular[0] = 1.0;
  if(lighting.specular[1] < 0.0)
    lighting.specular[1] = 0.0;
  if(lighting.specular[1] > 1.0)
    lighting.specular[1] = 1.0;
  if(lighting.specular[2] < 0.0)
    lighting.specular[2] = 0.0;
  if(lighting.specular[2] > 1.0)
    lighting.specular[2] = 1.0;
  rebuild();
}

void Displayobject::increase_ambient(std::vector<double>ambient){

  lighting.ambient[0] = lighting.ambient[0] + ambient[0];
  lighting.ambient[1] = lighting.ambient[1] + ambient[1];
  lighting.ambient[2] = lighting.ambient[2] + ambient[2];

  if(lighting.ambient[0] < 0.0)
    lighting.ambient[0] = 0.0;
  if(lighting.ambient[0] > 1.0)
    lighting.ambient[0] = 1.0;
  if(lighting.ambient[1] < 0.0)
    lighting.ambient[1] = 0.0;
  if(lighting.ambient[1] > 1.0)
    lighting.ambient[1] = 1.0;
  if(lighting.ambient[2] < 0.0)
    lighting.ambient[2] = 0.0;
  if(lighting.ambient[2] > 1.0)
    lighting.ambient[2] = 1.0;
  rebuild();
}

void Displayobject::rebuild(int doit){
  do_rebuild = doit;
}

int Displayobject::get_rebuild(void) const{
  return do_rebuild;
}

void Displayobject::increase_diffuse(std::vector<double>diffuse){

  lighting.diffuse[0] = lighting.diffuse[0] + diffuse[0];
  lighting.diffuse[1] = lighting.diffuse[1] + diffuse[1];
  lighting.diffuse[2] = lighting.diffuse[2] + diffuse[2];

  if(lighting.diffuse[0] < 0.0)
    lighting.diffuse[0] = 0.0;
  if(lighting.diffuse[0] > 1.0)
    lighting.diffuse[0] = 1.0;
  if(lighting.diffuse[1] < 0.0)
    lighting.diffuse[1] = 0.0;
  if(lighting.diffuse[1] > 1.0)
    lighting.diffuse[1] = 1.0;
  if(lighting.diffuse[2] < 0.0)
    lighting.diffuse[2] = 0.0;
  if(lighting.diffuse[2] > 1.0)
    lighting.diffuse[2] = 1.0;
  rebuild();
}

void Displayobject::increase_emission(std::vector<double>emission){

  lighting.emission[0] = lighting.emission[0] + emission[0];
  lighting.emission[1] = lighting.emission[1] + emission[1];
  lighting.emission[2] = lighting.emission[2] + emission[2];

  if(lighting.emission[0] < 0.0)
    lighting.emission[0] = 0.0;
  if(lighting.emission[0] > 1.0)
    lighting.emission[0] = 1.0;
  if(lighting.emission[1] < 0.0)
    lighting.emission[1] = 0.0;
  if(lighting.emission[1] > 1.0)
    lighting.emission[1] = 1.0;
  if(lighting.emission[2] < 0.0)
    lighting.emission[2] = 0.0;
  if(lighting.emission[2] > 1.0)
    lighting.emission[2] = 1.0;
  rebuild();
}

void Displayobject::move(double *d_in){
  move(d_in[0],d_in[1],d_in[2]);
}

void Displayobject::move(const std::vector<double>& d_in){
  move(d_in[0],d_in[1],d_in[2]);
}

void Displayobject::move(const Cartesian &d_in){
  move(d_in.get_x(),d_in.get_y(),d_in.get_z());
}

void Displayobject::move(double dx_in, double dy_in,double  dz_in){
  dx += dx_in;
  dy += dy_in;
  dz += dz_in;
}

void Displayobject::rotate(double dphi, double dchi, double dpsi){
  drot[0] = dphi;
  drot[1] = dchi;
  drot[2] = dpsi;
}

std::vector<double> Displayobject::get_origin(void) const{
  return origin.getxyza_vec();
}

std::vector<double> Displayobject::get_drot(void) const{
  return drot;
}

void Displayobject::apply_rotation_about_axes(double *xaxis, double *yaxis, double *zaxis){
  double *drotloc = new double[3];
  drotloc[0] = drot[0];
  drotloc[1] = drot[1];
  drotloc[2] = drot[2];
  quat.rotate_about_axes(xaxis,yaxis,zaxis,drotloc);
  delete [] drotloc;
  drot[0] = drot[1] = drot[2] = 0.0;
}

matrix Displayobject::get_rotation_matrix() const{
  return quat.getMatrix();
}

void Displayobject::set_origin(double *o_in){
  origin = Cartesian(o_in);
}

void Displayobject::set_origin(const std::vector<double> &o_in){
  origin = Cartesian(o_in);
}

void Displayobject::set_origin(const Cartesian &o_in){
  origin = o_in;
}
   
void Displayobject::set_origin(double o1,double o2, double o3){
  origin.set_x(o1);
  origin.set_y(o2);
  origin.set_z(o3);
}

void Displayobject::move_origin(){

  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT,viewport);
  std::vector<Cartesian> dum = getxyzc(float(viewport[2]/2),float(viewport[3]/2));
  Cartesian xyzo = dum[1];
  Cartesian xyzf = dum[0];
  std::vector<Cartesian> dum2 = getxyzc(float(viewport[2]),float(viewport[3]/2));
  Cartesian xyzr = dum2[1];
  std::vector<Cartesian> dum3 = getxyzc(float(viewport[2]/2),float(viewport[3]));
  Cartesian xyzu = dum3[1];
  Cartesian xaxis = xyzr - xyzo;
  Cartesian yaxis = xyzu - xyzo;
  Cartesian zaxis = xyzf - xyzo;

  xaxis.normalize();
  yaxis.normalize();
  zaxis.normalize();

  std::vector<double> dr;
  dr.push_back(dx*xaxis.get_x());
  dr.push_back(dx*xaxis.get_y());
  dr.push_back(dx*xaxis.get_z());

  std::vector<double> du;
  du.push_back(dy*yaxis.get_x());
  du.push_back(dy*yaxis.get_y());
  du.push_back(dy*yaxis.get_z());

  std::vector<double> df;
  df.push_back(dz*zaxis.get_x());
  df.push_back(dz*zaxis.get_y());
  df.push_back(dz*zaxis.get_z());

  dx = dr[0] + du[0] + df[0];
  dy = dr[1] + du[1] + df[1];
  dz = dr[2] + du[2] + df[2];
  move_origin(dx,dy,dz);
}

void Displayobject::move_origin(double* origin_in){
  std::vector<double> o;
  o.push_back(origin_in[0] + origin.get_x());
  o.push_back(origin_in[1] + origin.get_y());
  o.push_back(origin_in[2] + origin.get_z());
  o.push_back(origin.get_a());
  origin.setxyza_vec(o);
  dx = dy = dz = 0.0;
}

void Displayobject::move_origin(double o1, double o2, double o3){
  std::vector<double> o;
  o.push_back(o1 + origin.get_x());
  o.push_back(o2 + origin.get_y());
  o.push_back(o3 + origin.get_z());
  o.push_back(origin.get_a());
  origin.setxyza_vec(o);
  dx = dy = dz = 0.0;
}

void Displayobject::move_origin(const std::vector<double> &origin_in){
  std::vector<double> o;
  o.push_back(origin_in[0] + origin.get_x());
  o.push_back(origin_in[1] + origin.get_y());
  o.push_back(origin_in[2] + origin.get_z());
  o.push_back(origin.get_a());
  origin.setxyza_vec(o);
  dx = dy = dz = 0.0;
}

void Displayobject::move_origin(const Cartesian &origin_in){
  std::vector<double> o;
  o.push_back(origin_in.get_x() + origin.get_x());
  o.push_back(origin_in.get_y() + origin.get_y());
  o.push_back(origin_in.get_z() + origin.get_z());
  o.push_back(origin.get_a());
  origin.setxyza_vec(o);
  dx = dy = dz = 0.0;
}


void Displayobject::draw_text(const Quat &quat_in, double radius, double ox, double oy, double oz){
  GLboolean clip_test = glIsEnabled(GL_CLIP_PLANE0);
  Volume v = GetClippingPlanes();

  GLdouble projMatrix[16];
  glGetDoublev(GL_PROJECTION_MATRIX,projMatrix);
  GLdouble l = (projMatrix[12] - 1.0)/projMatrix[0]; 
  GLdouble r = (projMatrix[12] + 1.0)/projMatrix[0];
  GLdouble t = (projMatrix[13] - 1.0)/projMatrix[5];
  GLdouble b = (projMatrix[13] + 1.0)/projMatrix[5];
  GLdouble n = 1;
  if(fabs(projMatrix[10])>1e-7){
    n = (projMatrix[14] + 1)/projMatrix[10];
  }
  GLdouble win_h = ((projMatrix[13] + 1.0)/projMatrix[5] - (projMatrix[13] - 1.0)/projMatrix[5])*0.5;
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT,viewport);
  GLdouble view_h = (viewport[3] - viewport[1]);

  GLdouble* matinv = (GLdouble*)quat_in.getInvMatrix().to_dp();
  std::vector<SimpleText*>::const_iterator k = text_prims.begin();

  GLfloat params[4];
  glGetFloatv(GL_COLOR_CLEAR_VALUE,params);
  GLdouble y = params[0]*0.299 + params[1]*0.587 + params[2]*0.114;

  glColor4f(1.0,1.0,1.0,1.0);
  //if(k!=text_prims.end()) (*k)->SetDefaultColour();
  if(k!=text_prims.end()){
     if(clip_test){
        if((*k)->IsBillBoard()){
           glDisable(GL_CLIP_PLANE0);
           glDisable(GL_CLIP_PLANE1);
         } else {
           glEnable(GL_CLIP_PLANE0);
           glEnable(GL_CLIP_PLANE1);
         }
      }
  }
  while(k!=text_prims.end()){
    if(!(*k)->isMultiColoured()){
      glPushMatrix();
      Cartesian v = (*k)->GetVertices()[0]; 
      glTranslatef(v.get_x(),v.get_y(),v.get_z());
      glMultMatrixd(matinv);
      if( (*k)->GetTextureID()==0) (*k)->initialize();
      if(y<0.5)
        glBindTexture( GL_TEXTURE_2D, (*k)->GetTextureID() );
      else
        glBindTexture( GL_TEXTURE_2D, (*k)->GetTextureID_B() );
      //glBindTexture( GL_TEXTURE_2D, (*k)->GetTextureID() );
      double pix_w= (*k)->GetTexture().get_width();
      double pix_h= (*k)->GetTexture().get_height();
      double size = fabs(win_h*pix_h/view_h*radius/60.0);
      double size_offset = -2*fabs(win_h*((*k)->GetFontDescent())/view_h*radius/60.0);
      double ratio = double(pix_w)/pix_h;
      //glDisable(GL_TEXTURE_2D); glColor4f(0,0,0,1);
      glBegin(GL_QUADS);
      glTexCoord2f(0,0);
      glVertex3f(0,size_offset,0);
      glTexCoord2f(1,0);
      glVertex3f(2*size*ratio,size_offset,0);
      glTexCoord2f(1,1);
      glVertex3f(2*size*ratio,2*size+size_offset,0);
      glTexCoord2f(0,1);
      glVertex3f(0,2*size+size_offset,0);
      //(*k)->draw();
      glEnd();
      glPopMatrix();
    }
#ifdef __APPLE_CC__
    //}
#endif
    k++;
  }

  glColor4f(1.0,1.0,1.0,1.0);
  k = text_prims.begin();
  while(k!=text_prims.end()){
    if((*k)->isMultiColoured()){
      if( (*k)->GetTextureID()==0) (*k)->initialize();
      if(y<0.5)
        glBindTexture( GL_TEXTURE_2D, (*k)->GetTextureID() );
      else
        glBindTexture( GL_TEXTURE_2D, (*k)->GetTextureID_B() );
      glPushMatrix();
      Cartesian v = (*k)->GetVertices()[0]; 
      double x_cent_off = 0;
      double y_cent_off = 0;
      if((*k)->GetCentered()){
         x_cent_off = -fabs(win_h*(*k)->GetTextWidth()/view_h*radius/60.0);
         y_cent_off = -fabs(win_h*(*k)->GetTextHeight()/2.0/view_h*radius/60.0);
      }
      double pix_w= (*k)->GetTexture().get_width();
      double pix_h= (*k)->GetTexture().get_height();
      double size = fabs(win_h*pix_h/view_h*radius/60.0);
      double size_offset = -2*fabs(win_h*((*k)->GetFontDescent())/view_h*radius/60.0);
      double ratio = double(pix_w)/pix_h;
      if((*k)->IsBillBoard()){
        glLoadIdentity();
        GLdouble x = (v.get_x()-.5)*(r-l);
        GLdouble y = (v.get_y()-.5)*(b-t);
        GLdouble z =  -n-3;
	size_offset = 0;
        size = fabs(win_h*pix_h/view_h);
        glTranslatef(x,y,z);
      } else {
        glTranslatef(v.get_x(),v.get_y(),v.get_z());
        glMultMatrixd(matinv);
      }
      //glDisable(GL_TEXTURE_2D); glColor4f(0,0,0,1);
      glBegin(GL_QUADS);
      glTexCoord2f(0,0);
      glVertex3f(x_cent_off,y_cent_off+size_offset,0);
      glTexCoord2f(1,0);
      glVertex3f(x_cent_off+2*size*ratio,y_cent_off+size_offset,0);
      glTexCoord2f(1,1);
      glVertex3f(x_cent_off+2*size*ratio,y_cent_off+2*size+size_offset,0);
      glTexCoord2f(0,1);
      glVertex3f(x_cent_off,y_cent_off+2*size+size_offset,0);
      //(*k)->draw();
      glEnd();
      glPopMatrix();
    }
    k++;
  }
  delete [] matinv;
  if(clip_test){
    glEnable(GL_CLIP_PLANE0);
    glEnable(GL_CLIP_PLANE1);
  }
}

void Displayobject::draw_images(void){

  glDepthFunc(GL_ALWAYS);
  std::vector<BillBoard*>::const_iterator k = image_prims.begin();

  while(k!=image_prims.end()){
    (*k)->draw();
    k++;
  }
  glDepthFunc(GL_LESS);
}


std::vector<Primitive *> Displayobject::GetTransparentPrimitives(){
	/* This is dead and useless !!! */
  std::vector<Primitive *> trans_prims;
  std::vector<Primitive*>::iterator k = prims.begin();
  while(k!=prims.end()){
    if((*k)->get_transparent()>0){
       trans_prims.push_back(*k);
    }
    k++;
  }
  return trans_prims;

}

void Displayobject::draw_lines(double *override_colour, int transparent_in, int selective_override) const {
  std::vector<Primitive*>::const_iterator k = prims.begin();

  GLfloat col[4];
  GLfloat *colp=0;

  while(k!=prims.end()){
    if((*k)->isLine()){
      if(override_colour){
        if(selective_override==0||(selective_override==1&&(*k)->GetColourOverride())){
	  col[0] = override_colour[0];
	  col[1] = override_colour[1];
	  col[2] = override_colour[2];
	  col[3] = override_colour[3];
	  colp = col;
        }
      }
      (*k)->set_draw_colour(colp);
      (*k)->draw(override_colour, selective_override);
    }
    k++;
  }
  k = surf_prims.begin();
  while(k!=surf_prims.end()){
    if((*k)->get_transparent()==transparent_in&&(*k)->isLine()){
      if(override_colour){
        if(selective_override==0||(selective_override==1&&(*k)->GetColourOverride())){
	  col[0] = override_colour[0];
	  col[1] = override_colour[1];
	  col[2] = override_colour[2];
	  col[3] = override_colour[3];
	  colp = col;
        }
      }
      (*k)->set_draw_colour(colp);
      (*k)->draw(override_colour, selective_override);
    }
    k++;
  }

}

void Displayobject::SetAlpha(double alpha_in){
  alpha = alpha_in;
  std::vector<Primitive*>::iterator k = prims.begin();
  while(k!=prims.end()){
    (*k)->set_transparent(transparent);
    (*k)->SetAlpha(alpha);
    k++;
  }
  k = surf_prims.begin();
  while(k!=surf_prims.end()){
    (*k)->set_transparent(transparent);
    (*k)->SetAlpha(alpha);
    k++;
  }
}
void Displayobject::set_transparent(int trans){
  transparent = trans;
  std::vector<Primitive*>::iterator k = prims.begin();
  while(k!=prims.end()){
    (*k)->set_transparent(trans);
    (*k)->SetAlpha(alpha);
    k++;
  }
  k = surf_prims.begin();
  while(k!=surf_prims.end()){
    (*k)->set_transparent(trans);
    (*k)->SetAlpha(alpha);
    k++;
  }
}

void Displayobject::draw_prims(double *override_colour, int transparent_in, int selective_override) const {
  GLfloat col[4];
  GLfloat *colp=0;

  glMaterialfv(GL_FRONT, GL_SPECULAR,  lighting.specular);
  glMaterialfv(GL_FRONT, GL_EMISSION,  lighting.emission);
  glMaterialf(GL_FRONT, GL_SHININESS, lighting.shininess);

  std::vector<Primitive*>::const_iterator k;
  k = prims.begin();
  while(k!=prims.end()){
    if((*k)->get_transparent()==transparent_in&&!(*k)->isLine()){
      if(override_colour){
        if(selective_override==0||(selective_override==1&&(*k)->GetColourOverride())){
	  col[0] = override_colour[0];
	  col[1] = override_colour[1];
	  col[2] = override_colour[2];
	  col[3] = override_colour[3];
	  colp = col;
        }
      }
      (*k)->set_draw_colour(colp);
      (*k)->draw(override_colour, selective_override);
    }
    k++;
  }

}

void Displayobject::draw_surf_prims(double *override_colour, int transparent_in, int selective_override) const {
  GLfloat col[4];
  GLfloat *colp=0;

  glMaterialfv(GL_FRONT, GL_SPECULAR,  lighting.specular);
  glMaterialfv(GL_FRONT, GL_EMISSION,  lighting.emission);
  glMaterialf(GL_FRONT, GL_SHININESS, lighting.shininess);

  std::vector<Primitive*>::const_iterator k;
  k = surf_prims.begin();
  while(k!=surf_prims.end()){
    if((*k)->get_transparent()==transparent_in&&!(*k)->isLine()){
      if(override_colour){
        if(selective_override==0||(selective_override==1&&(*k)->GetColourOverride())){
	  col[0] = override_colour[0];
	  col[1] = override_colour[1];
	  col[2] = override_colour[2];
	  col[3] = override_colour[3];
	  colp = col;
        }
      }
      (*k)->set_draw_colour(colp);
      (*k)->draw(override_colour, selective_override);
    }
    k++;
  }

}

void Displayobject::draw_solids(double *override_colour, int transparent_in, int selective_override) const {
  std::vector<Primitive*>::const_iterator k = prims.begin();

  GLfloat col[4];
  GLfloat *colp=0;

  glMaterialfv(GL_FRONT, GL_SPECULAR,  lighting.specular);
  glMaterialfv(GL_FRONT, GL_EMISSION,  lighting.emission);
  glMaterialf(GL_FRONT, GL_SHININESS, lighting.shininess);

  while(k!=prims.end()){
    if((*k)->get_transparent()==transparent_in&&!(*k)->isLine()){
      if(override_colour){
        if(selective_override==0||(selective_override==1&&(*k)->GetColourOverride())){
	  col[0] = override_colour[0];
	  col[1] = override_colour[1];
	  col[2] = override_colour[2];
	  col[3] = override_colour[3];
	  colp = col;
        }
      }
      (*k)->set_draw_colour(colp);
      (*k)->draw(override_colour, selective_override);
    }
    k++;
  }

  k = surf_prims.begin();
  while(k!=surf_prims.end()){
    if((*k)->get_transparent()==transparent_in&&!(*k)->isLine()){
      if(override_colour){
        if(selective_override==0||(selective_override==1&&(*k)->GetColourOverride())){
	  col[0] = override_colour[0];
	  col[1] = override_colour[1];
	  col[2] = override_colour[2];
	  col[3] = override_colour[3];
	  colp = col;
        }
      }
      (*k)->set_draw_colour(colp);
      (*k)->draw(override_colour, selective_override);
    }
    k++;
  }
}

void Displayobject::draw(double *override_colour, int transparent_in, int selective_override) const {
  std::vector<Primitive*>::const_iterator k = prims.begin();

  GLfloat col[4];
  GLfloat *colp=0;

  glMaterialfv(GL_FRONT, GL_SPECULAR,  lighting.specular);
  glMaterialfv(GL_FRONT, GL_EMISSION,  lighting.emission);
  glMaterialf(GL_FRONT, GL_SHININESS, lighting.shininess);

  while(k!=prims.end()){
    if((*k)->get_transparent()==transparent_in){
      if(override_colour){
        if(selective_override==0||(selective_override==1&&(*k)->GetColourOverride())){
	  col[0] = override_colour[0];
	  col[1] = override_colour[1];
	  col[2] = override_colour[2];
	  col[3] = override_colour[3];
	  colp = col;
        }
      }
      (*k)->set_draw_colour(colp);
      (*k)->draw(override_colour, selective_override);
    }
    k++;
  }
}

Primitive* Displayobject::get_primitive(int i) const{
  return prims[i];
}

double* Displayobject::get_primorigin(int i) const{
  Cartesian cart = prims[i]-> get_origin();
  double *pv = new double[4];
  pv[0] = cart.get_x();
  pv[1] = cart.get_y();
  pv[2] = cart.get_z();
  pv[3] = cart.get_a();
  return pv;
}

double* Displayobject::rotate_point(double x, double y, double z){

  Cartesian v = quat.getMatrix()*Cartesian(x,y,z);

  return v.to_dp();
}

double* Displayobject::get_primoriginrot(int i) const{

  Cartesian cart = prims[i]-> get_origin();

  Cartesian v = quat.getMatrix()*cart;

  return v.to_dp();

}

void Displayobject::MoveTextPrimitiveInWindowCoords(SimpleText *text_primitive, double x, double y, double z, double *world_quat_dvals){
  
  std::vector<double> wqd;
  wqd.push_back(world_quat_dvals[0]);
  wqd.push_back(world_quat_dvals[1]);
  wqd.push_back(world_quat_dvals[2]);
  wqd.push_back(world_quat_dvals[3]);
  MoveTextPrimitiveInWindowCoords(text_primitive,x,y,z,wqd);
}

void Displayobject::MoveTextPrimitiveInWindowCoords(SimpleText *text_primitive, double x, double y, double z, const std::vector<double> &world_quat_dvals){
  Quat world_quat;
  world_quat.Setdval(world_quat_dvals);
  world_quat.postMult(quat);

  Cartesian result = world_quat.getMatrix()*Cartesian(x,y,z);

  result += text_primitive->get_origin();
  text_primitive->set_origin(result.to_dp());

}

SimpleText* Displayobject::findtextprimitive(const std::vector<Cartesian> &xyzbox){

  SimpleText *dummy = 0;
  std::vector<Cartesian> primorigin;

  for(unsigned int i=0;i<text_prims.size();i++){
    primorigin.push_back(text_prims[i]->get_origin());
  }

  int nearprim = findprimc(xyzbox,primorigin,origin,quat.getInvMatrix());

  if(nearprim>-1)
    return text_prims[nearprim];
  else
    return dummy;

}          

int Displayobject::findprimitive(const std::vector<Cartesian> &xyzbox){

  std::vector<Cartesian> primorigin;

  for(unsigned int i=0;i<prims.size();i++){
    primorigin.push_back(prims[i]->get_origin());
  }

  int nearprim = findprimc(xyzbox,primorigin,origin,quat.getInvMatrix());

  return nearprim;

}          

void Displayobject::clear_labels(void){
  // This seems to cause a nasty double delete.
  /*
  std::vector<SimpleText*>::iterator text_prim_iter = text_prims.begin();
  while(text_prim_iter!=text_prims.end()){
    if(*text_prim_iter) delete *text_prim_iter;
    text_prim_iter++;
  }
  */
  text_prims.clear();
}

std::vector<int> Displayobject::GetPrimitivesInVolume(Volume volume) const{

  std::vector<Primitive *>:: const_iterator prim_iter = prims.begin();
  Cartesian prim_origin;
  int clicked;
  Plane plane;
  std::vector <Cartesian> points;
  Cartesian n;
  Cartesian p2prim;
  std::vector<int> clicked_prims;

  int i=0;
  while(prim_iter!=prims.end()){
     clicked = 1;
     prim_origin = (*prim_iter)->get_origin();
     prim_origin = quat.getInvMatrix()*prim_origin;
     prim_origin += origin;
     for(int ii=0;ii<volume.GetNumberOfPlanes();ii++){
       plane = volume.GetPlane(ii);
       n = plane.get_normal();
       points = plane.find_points_on_plane();
       p2prim = points[0] - prim_origin;
       n.normalize();
       p2prim.normalize();
       if(n.DotProduct(n,p2prim)<0.0)
         clicked = 0;
     }
     if(clicked) clicked_prims.push_back(i);
     prim_iter++;
     i++;
  }

  return clicked_prims;

}

std::string X112PSFontName(const std::string &name_in, const std::string &weight, const std::string &slant){
  std::string name="";

  int space = 1;
  for(unsigned int i=0;i<name_in.size();i++){
    if(space){
      name += toupper(name_in[i]);
      space = 0;
    }else{
      if(name_in[i]==' ')
       space = 1;
      else
       name += name_in[i];
    }
      
  }

  if(weight=="bold")
    name += "Bold";

  if(slant=="i")
    name += "Italic";

  if(slant=="o")
    name += "Oblique";

  return name;
  
}

void Displayobject::OutputTextLabels(std::ofstream &fp, const Quat &quat_in, double radius, double ox, double oy, double oz){

  GLboolean clip_test = glIsEnabled(GL_CLIP_PLANE0);
  Volume v = GetClippingPlanes();
  std::vector<SimpleText*>::const_iterator kk = text_prims.begin();
  while(kk!=text_prims.end()){
    if(clip_test){
      if((*kk)->IsBillBoard()||(v.PointInVolume(get_rotation_matrix()*((*kk)->GetVertices()[0]+origin)))){
        glDisable(GL_CLIP_PLANE0);
        glDisable(GL_CLIP_PLANE1);
        (*kk)->DrawPovray(fp,quat_in,radius,ox,oy,oz,get_rotation_matrix(),origin,v);
        glEnable(GL_CLIP_PLANE0);
        glEnable(GL_CLIP_PLANE1);
      }
    } else {
      (*kk)->DrawPovray(fp,quat_in,radius,ox,oy,oz,get_rotation_matrix(),origin,v);
    }
    kk++;
  }
  std::vector<BillBoard*>::const_iterator ll = image_prims.begin();
  while(ll!=image_prims.end()){
    if(clip_test){
      glDisable(GL_CLIP_PLANE0);
      glDisable(GL_CLIP_PLANE1);
      (*ll)->DrawPovray(fp,quat_in,radius,ox,oy,oz,get_rotation_matrix(),origin,v);
      glEnable(GL_CLIP_PLANE0);
      glEnable(GL_CLIP_PLANE1);
    } else {
      (*ll)->DrawPovray(fp,quat_in,radius,ox,oy,oz,get_rotation_matrix(),origin,v);
    }
    ll++;
  }
  if(clip_test){
    glEnable(GL_CLIP_PLANE0);
    glEnable(GL_CLIP_PLANE1);
  }
}

void Displayobject::DrawPovRay(std::ofstream &fp, const Quat &quat_in, double radius, double ox, double oy, double oz){

  Volume v = GetClippingPlanes();
  std::vector<Primitive*>::const_iterator k = prims.begin();
  while(k!=prims.end()){
    (*k)->DrawPovray(fp,quat_in,radius,ox,oy,oz,get_rotation_matrix(),origin,v);
    k++;
  }
  k = surf_prims.begin();
  while(k!=surf_prims.end()){
    (*k)->DrawPovray(fp,quat_in,radius,ox,oy,oz,get_rotation_matrix(),origin,v);
    k++;
  }
}

void Displayobject::DrawPostscript(std::ofstream &fp, const Quat &quat_in, double radius, double ox, double oy, double oz){

  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT,viewport);

  Cartesian tl = getxyzc(0,0)[0];
  Cartesian br = getxyzc(viewport[2]/2.0,viewport[3]/2.0)[0];
  tl = quat_in.getInvMatrix()*tl;
  br = quat_in.getInvMatrix()*br;

  double xscale = fabs(tl.get_x()-br.get_x());
  double yscale = fabs(tl.get_y()-br.get_y());

  double xscaleps = 421.0; // For A4 !!
  double yscaleps = 421.0; // For A4 !!

  int xoff = 297;
  int yoff = 421;

  if(yscaleps<xscaleps)
    xscaleps = yscaleps;
  if(xscaleps<yscaleps)
    yscaleps = xscaleps;
  if(yscale<xscale)
    xscale = yscale;
  if(xscale<yscale)
    yscale = xscale;

  Cartesian p;

  Volume v = GetClippingPlanes();
  std::vector<Primitive*>::const_iterator k = prims.begin();
  while(k!=prims.end()){
    std::vector<Cartesian> vertices = (*k)->GetVertices();
    (*k)->DrawPostScript(fp,quat_in,radius,ox,oy,oz,get_rotation_matrix(),origin,xoff,yoff,xscale,yscale,xscaleps,v);
    k++;
  }
  k = surf_prims.begin();
  while(k!=surf_prims.end()){
    (*k)->DrawPostScript(fp,quat_in,radius,ox,oy,oz,get_rotation_matrix(),origin,xoff,yoff,xscale,yscale,xscaleps,v);
    k++;
  }

  /*
  * Ok, now text. This works quite nicely except it doesn't handle any
  * of the markup language used to define colour, subscripts, etc. I'll
  * need to make all that more portable and interpret it here somehow.
  * Tricky. */

  std::vector<BillBoard*>::const_iterator ll = image_prims.begin();
  while(ll!=image_prims.end()){
    (*ll)->DrawPostScript(fp,quat_in,radius,ox,oy,oz,get_rotation_matrix(),origin,double(xoff),double(yoff),xscale,yscale,xscaleps,v);
    ll++;
  }
  std::vector<SimpleText*>::const_iterator l = text_prims.begin();
  while(l!=text_prims.end()){
    (*l)->DrawPostScript(fp,quat_in,radius,ox,oy,oz,get_rotation_matrix(),origin,double(xoff),double(yoff),xscale,yscale,xscaleps,v);
    l++;
  }

}

