/*
     pygl/cdisplayobject.h: CCP4MG Molecular Graphics Program
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


#ifndef _CCP4MG_DISPLAYOBJ_
#define _CCP4MG_DISPLAYOBJ_
#include "quat.h"
#include "volume.h"
#include "cprimitive.h"
#include <vector>
#include <string>
#include <fstream>

#include "matrix.h"


class Displayobject{
  int transparent;
  double alpha;
  int build_disp_list;
  int draw_symm;
  int symm_diff_colour;
  int draw_unit_cell;
  int occDataAttrib;
  int useVBO;
  int useVertexArrays;
 protected:
  std::vector<Primitive*> prims;
  std::vector<Primitive*> surf_prims;
  std::vector<BillBoard*> image_prims;
  std::vector<SimpleText*> text_prims;
  std::vector<BillboardPrimitive*> bill_prims;
  std::vector<BillboardPrimitive*> imposter_prims;
 public:
  Quat camera_quat; // These two are for any drawing which needs to
  Cartesian camera_origin;  // know about camera, ie zsorting.
  int anchored;
  virtual const std::vector<Primitive*> &GetPrimitives() const;
  virtual const std::vector<Primitive*> &GetSurfacePrimitives() const;
  virtual const std::vector<BillBoard*> &GetImagePrimitives() const;
  virtual const std::vector<SimpleText*> &GetTextPrimitives() const;
  virtual const std::vector<BillboardPrimitive*> &GetBillboardPrimitives() const;
  virtual const std::vector<BillboardPrimitive*> &GetImposterPrimitives() const;
  virtual const size_t GetNumberOfPrimitives() const { return prims.size();} ;
  virtual const size_t GetNumberOfSurfacePrimitives() const { return surf_prims.size();} ;
  virtual const size_t GetNumberOfImagePrimitives() const { return image_prims.size();} ;
  virtual const size_t GetNumberOfTextPrimitives() const { return text_prims.size();} ;
  virtual const size_t GetNumberOfBillboardPrimitives() const { return bill_prims.size();} ;
  virtual const size_t GetNumberOfImposterPrimitives() const { return imposter_prims.size();} ;
  std::vector<matrix> symm_mat;
  std::vector<int> symm_nos;
  std::vector<Cartesian> unit_cell;
  Cartesian origin;
  Cartesian com;
  int visible;
  std::vector<double> rot;
  std::vector<double> drot;
  double dx;
  double dy;
  double dz;
  Quat quat;
  Lighting lighting;
  int do_rebuild;
  void rebuild(int doit=1);
  Displayobject();
  virtual ~Displayobject();
  void changevis(void);
  void SetPrimitives(std::vector<Primitive*> prims_in){prims = prims_in;};
  //void add_primitive(Primitive &prim);
  //void add_text_primitive(SimpleText &prim);
  void add_surf_primitive(Primitive *prim);
  void add_primitive(Primitive *prim);
  void add_text_primitive(SimpleText *prim);
  void add_image_primitive(BillBoard *prim);
  void add_billboard_primitive(BillboardPrimitive *prim);
  void add_imposter_primitive(BillboardPrimitive *prim);
  void increase_shininess(double shininess);
  void increase_specular(std::vector<double>specular);
  void increase_ambient(std::vector<double>specular);
  void increase_diffuse(std::vector<double>specular);
  void increase_emission(std::vector<double>specular);
  void move(double dx_in, double dy_in,double  dz_in);
  void move(double *d_in);
  void move(const std::vector<double>& d_in);
  void move(const Cartesian &d_in);
  void rotate(double dphi, double dchi, double dpsi);
  void move_origin();
  void move_origin(double *o_in);
  void move_origin(double o1, double o2, double o3);
  void move_origin(const std::vector<double>& o_in);
  void move_origin(const Cartesian &o_in);
  std::vector<double> get_origin(void) const;
  void set_origin(double o1, double o2, double o3);
  void set_origin(double *o_in);
  void set_origin(const std::vector<double>& o_in);
  void set_origin(const Cartesian &o_in);
  std::vector<double> get_drot(void) const;
  void apply_rotation_about_axes(double *xaxis, double *yaxis, double *zaxis);
  matrix get_rotation_matrix() const;
  void draw_lines(double *override_colour=0, int transparent=0,int selective_override=0) const ;
  void draw_solids(double *override_colour=0, int transparent=0,int selective_override=0) const ;
  void draw_imposters(double *override_colour=0, int transparent=0,int selective_override=0) const ;
  void draw_prims(double *override_colour=0, int transparent=0,int selective_override=0) const ;
  void draw_surf_prims(double *override_colour=0, int transparent=0,int selective_override=0) const ;
  void draw(double *override_colour=0, int transparent=0,int selective_override=0) const ;
  std::vector<Primitive *> GetTransparentPrimitives();
  void draw_text(const Quat &quat_in, double radius, double ox, double oy, double oz, double fontScaling=1.0,int nCustomPlanes=0);
  void draw_text_background(const Quat &quat_in, double radius, double ox, double oy, double oz, double fontScaling=1.0);
  void draw_billboards(const Quat &quat_in, double radius, double ox, double oy, double oz);
  void draw_images(void);
  void clear_prims(void);
  void clear_images(void);
  void clear_labels(void);
  void clear_billboards(void);
  int get_rebuild(void) const;
  void MoveTextPrimitiveInWindowCoords(SimpleText *text_primitive, double x, double y, double z, double *world_quat_dvals);
  void MoveTextPrimitiveInWindowCoords(SimpleText *text_primitive, double x, double y, double z, const std::vector<double> &world_quat_dvals);
  int findprimitive(const std::vector<Cartesian> &xyzbox);
  SimpleText* findtextprimitive(const std::vector<Cartesian> &xyzbox);
  Primitive* get_primitive(int i) const;
  double* get_primorigin(int i) const;
  double* get_primoriginrot(int i) const;
  std::vector<int> GetPrimitivesInVolume(Volume volume) const;
  double* rotate_point(double x, double y, double z);

  // Ability to play with Text, this could be done better by exposing
  // the Text class to Python. But were doing it like this now.
  int *GetTextIDS(void) const;
  void SetTextFont(const std::string family,  const std::string weight, 
                   const std::string slant, const int size, const int underline=0);
  void DeleteText(void);
  int GetNumberOfTextIDS(void) const;
  void SetTextString(int text_id, const char* new_string);
  void SetTextString(int text_id, std::string new_string);
  const char* GetTextString(int text_id) const;
  void DeleteTextPrimitive(int text_id);

  int BuildDisplayList(void);
  void SetBuildDisplayList(int build_or_not);
  void SetUnitCell(const std::vector<Cartesian> &cell_params);
  void DrawUnitCell();
  Quat GetUnitCellAlignmentRotation(const std::string &axis);
  void SetSymmetryMatrices(const std::vector<matrix> &symm_mat_in);
  void SetSymmetryMatrixNumbers(const std::vector<int> &symm_nos);
  unsigned int GetNumSymmetryMatrices() const;
  matrix GetSymmetryMatrix(int i) const;
  int GetSymmetryMatrixNumber(int i) const ;
  void ApplyTranslation();
  void ApplyRotation();
  void ApplySymmetryMatrix(int i);
  void ApplySymmetryMatrix_RotationOnly(int i);
  void SetDrawSymmetry(int draw_symm_in) { draw_symm = draw_symm_in; };
  int GetDrawSymmetry() const { return draw_symm; };
  void SetDrawSymmetryColoured(int symm_diff_colour_in) { symm_diff_colour = symm_diff_colour_in; };
  int GetDrawSymmetryColoured() const { return symm_diff_colour; };
  void SetDrawUnitCell(int draw_unit_cell_in) { draw_unit_cell = draw_unit_cell_in; };
  int GetDrawUnitCell() const { return draw_unit_cell; };
  int IsAnchored() const {return anchored;};
  void SetAnchored(bool _anchored) {anchored=_anchored;};
  void ZoomIn() {};
  void ZoomOut() {};
  void DrawPostscript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz);
  void DrawPovRay(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz);
  void OutputTextLabels(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz);
  void SetAlpha(double);
  double GetAlpha() const {return alpha;};
  void set_transparent(int trans_in);
  int get_transparent() const {return transparent;} ;
  int reInitializeTextPrims();
  const std::vector<Cartesian> &GetUnitCell() const {return unit_cell;};
  const std::vector<matrix> &GetSymmetryMatrices() const {return symm_mat;};
  void SetOccDataAttrib(int occDataAttrib_in) {occDataAttrib=occDataAttrib_in;} ;
  void forceRegenerateSurfaceArrays();
  int GetUseVBO() const {return useVBO;} ; 
  void SetUseVBO(const int _useVBO); 
  int GetUseVertexArrays() const {return useVertexArrays;} ; 
  void SetUseVertexArrays(const int _useVertexArrays); 
};

void DrawSortedTransparentPrimitives(const std::vector<Displayobject> &objs, int acsize, double xoff, double yoff, std::vector<std::vector<double> > jarray, std::vector<Cartesian> axes, bool antialias, bool rebuilt);

#endif
