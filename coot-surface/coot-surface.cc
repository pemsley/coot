/* coot-surface/coot-surface.cc
 * 
 * Copyright 2005 The University of Oxford
 * Copyright 2005 The University of York
 * Author: Martin Noble, Jan Gruber, Paul Emsley
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */


#if defined _MSC_VER
#include <windows.h>
#endif

#if __APPLE__
#   include <OpenGL/gl.h>
#else
#   include <GL/gl.h>
#endif

#include "coot-surface.hh"
#include "surface/CXXCreator.h"
#include "surface/CXXUtils.h"

#include "rgbreps.h"
#include "coot-utils/coot-coord-utils.hh"

// assign_charges is an instruction to CXX surface code to generate
// the charges - sometime we want to do that (e.g. surface on a
// protein), sometimes not (e.g. surface on a ligand).
// 
void
coot::surface::fill_from(mmdb::Manager *mol, int selHnd, float col_scale, bool assign_charges) {


   // debug::
   bool debug = 0;

   mmdb::PPAtom atoms;
   int n_atoms;
   mol->GetSelIndex(selHnd, atoms, n_atoms);
   std::vector<float> save_charges(n_atoms, -100); // all unset

   if (! assign_charges) 
      // save our charges (restore later).
      for (int iat=0; iat<n_atoms; iat++)
	 if (atoms[iat]->charge > CXX_UNSET_CHARGE)
	    save_charges[iat] = atoms[iat]->charge;
   

   if (debug) { 
      std::cout << "========pre =============" << std::endl;
      for (int iat=0; iat<n_atoms; iat++) {
	 std::cout << "  atom: " << atoms[iat]->GetChainID() << " " <<  atoms[iat]->GetSeqNum() << " "
		   << atoms[iat]->name << " charge: " << atoms[iat]->charge << std::endl;
      }
   }

   
   theSurface = new CXX_mot::CXXSurface;
   CXX_mot::CXXChargeTable theChargeTable;
   CXXUtils_old::assignCharge(mol, selHnd, &theChargeTable); // magic happens

   if (! assign_charges) 
      // restore charges 
      for (int iat=0; iat<n_atoms; iat++)
	 if (save_charges[iat] > CXX_UNSET_CHARGE)
	    atoms[iat]->charge = save_charges[iat];
      

   if (debug) {
      std::cout << "========post=============" << std::endl;
      mmdb::PPAtom atoms;
      int n_atoms;
      mol->GetSelIndex(selHnd, atoms, n_atoms);
      for (int iat=0; iat<n_atoms; iat++) {
	 std::cout << "  atom: " << atoms[iat]->GetChainID() << " " <<  atoms[iat]->GetSeqNum() << " "
		   << atoms[iat]->name << " charge: " << atoms[iat]->charge << std::endl;
      }
   }

   
   theSurface->calculateFromAtoms(mol, selHnd, selHnd, 1.4, 0.5, false); // 0.785 is 45
							  // degrees used by
							  // Martin
   evaluateElectrostaticPotential(mol, selHnd, col_scale);
}

void
coot::surface::fill_surface(mmdb::Manager *mol, int SelHnd_selection, int SelHnd_all,
			    float col_scale, bool assign_charges) {

   // e.g. selection handles: residues_of_the_active_site, residues_of_the_chain
   // the first is a subset of the other.
   theSurface = new CXX_mot::CXXSurface;
   if (assign_charges) { 
      CXX_mot::CXXChargeTable theChargeTable;
      CXXUtils_old::assignCharge(mol, SelHnd_all, &theChargeTable);
   } 
   theSurface->calculateFromAtoms(mol, SelHnd_selection, SelHnd_all, 1.4, 0.5, false);
   evaluateElectrostaticPotential(mol, SelHnd_all, col_scale);
} 
			    

// 
void
coot::surface::draw(double *override_colour, int selective_override) {

  double coords[4];
  float whiteColour[] = {1.,1.,1.,1.};
  float blackColor[] = {0.,0.,0.,0.};
  float mat_specular[] = {0.8, 0.8, 0.8, 0.8};
  
  //Save the aspects of gl state that I am going to change
  glPushAttrib(GL_COLOR_BUFFER_BIT | GL_CURRENT_BIT | GL_ENABLE_BIT | GL_LIGHTING_BIT);
  
  glEnable(GL_LIGHTING);
  glEnable(GL_NORMALIZE);
  
  double colour[] = {0.5,0.5,0.5,1.0};
  
  double *vertices = (double *) malloc (3*sizeof(double)*theSurface->numberOfVertices());
  double *normals  = (double *) malloc (3*sizeof(double)*theSurface->numberOfVertices());
  double *colours  = (double *) malloc (4*sizeof(double)*theSurface->numberOfVertices());
  
  for (int i=0; i< theSurface->numberOfVertices(); i++){
    //Use the colour if it has been assigned
    if (!theSurface->getCoord("colour", i, coords)){
      for (int k=0; k<3; k++) {
        colour[k] = coords[k];
      }
      // set alpha to 1 (rather than default 0) needed for screenshot
      colour[3] = 1.0;
    }
    for (int k=0; k<4; k++) colours [4*i+k] = colour[k];
    
    //Copy normal into normals array	
    if (!theSurface->getCoord("normals", i, coords)){
      for (int k=0; k<3; k++) normals[3*i+k] = coords[k];
    }
    
    //Copy vertex into vertices array	
    if (!theSurface->getCoord("vertices", i, coords)){
      for (int k=0; k<3; k++) vertices[3*i +k] = coords[k];
    }
  }
  
  //Tell OpenGL to use these arrays for drawing purposes
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(3, GL_DOUBLE, 0, vertices);
  glEnableClientState(GL_NORMAL_ARRAY);
  glNormalPointer(GL_DOUBLE, 0, normals);
  glEnableClientState(GL_COLOR_ARRAY);
  glColorPointer(4, GL_DOUBLE, 0, colours);
  
  //Let the returned colour dictate: note obligatory order of these calls
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);
  
  //And begin the actual drawing
  glBegin(GL_TRIANGLES);
  
  //Set material properties that are not per-vertex
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
  glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, blackColor);
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 120); 
  
  for (int i=0; i< theSurface->numberOfTriangles(); i++){
     glArrayElement(theSurface->vertex(i, 0));
     glArrayElement(theSurface->vertex(i, 1));
     glArrayElement(theSurface->vertex(i, 2));
  }
  glEnd();
  
  glPopAttrib();
  
  free (colours);
  free (normals);
  free (vertices);
}

void
coot::surface::transparent_draw(float opacity) {

   GLubyte stippleMask[128] =
   {0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55,
    0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55,
    0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55,
    0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55,
    0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55,
    0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55,
    0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55,
    0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55};

   glEnable(GL_POLYGON_STIPPLE);
   glPolygonStipple(stippleMask);
   draw(0,0);
   glDisable(GL_POLYGON_STIPPLE);
   
} 



// colscale is by default 0.2.
void
coot::surface::evaluateElectrostaticPotential(mmdb::Manager *theManager, int selHnd, float col_scale) {
  iEval = 0;
  evaluatePhiAndColourWithDefaultScheme(theManager, selHnd, col_scale);
}

int
coot::surface::evaluatePhiAndColourWithDefaultScheme(mmdb::Manager *theManager, const int selHnd,
						     float col_scale){
   std::cout << "In evaluatePhiAndColourWithDefaultScheme with col_scale: "
	     << col_scale << std::endl;
   CColourScheme defaultScheme;
   std::vector<float> typ;
   typ.push_back(-col_scale);
   typ.push_back(-col_scale);
   typ.push_back(0.0);
   typ.push_back(col_scale);
   typ.push_back(col_scale);
   std::vector<std::string> cols;
   cols.push_back("red");
   cols.push_back("red");
   cols.push_back("white");
   cols.push_back("blue");
   cols.push_back("blue");
   defaultScheme.SetSchemeFloat(typ, cols);
   return evaluatePhiAndColourWithScheme(theManager, selHnd, defaultScheme);
}

int coot::surface::evaluatePhiAndColourWithScheme(mmdb::Manager *theManager, const int selHnd, CColourScheme &colourScheme){
  cout << "In evaluatePhiAndColourWithScheme\n"; cout.flush();
  //Instantiate and calculate electrostatic potential
  CXX_mot::CXXCreator theCreator(theManager, selHnd);
  theCreator.calculate();	
  
  mmdb::realtype a, b, c;
  mmdb::realtype al, be, ga, vol;
  int ncode; 
  theManager->GetCell(a, b, c, al, be, ga, vol, ncode);
  clipper::Cell cell(clipper::Cell_descr(a, b, c, 
 					 clipper::Util::d2rad(al), 
 					 clipper::Util::d2rad(be), 
 					 clipper::Util::d2rad(ga)));

  //Coerce map into clipper NXmap
  clipper::NXmap<double> thePhiMap(theCreator.coerceToClipperMap(cell));
  std::cout << " ====================== xmap created ================ " << std::endl;

  //Interpolate into this map at ssurface vertices
  if (interpolateIntoMap("vertices", "potential", thePhiMap)) {
     return 1;
  }

  //Use these values to apply a colour scheme
  if (colourByScalarValue("potential", colourScheme)) return 1;
  return 0;
}

int coot::surface::colourByScalarValue(const std::string &scalarType, CColourScheme &colourScheme){
  cout << "In colourByScalarValue\n"; cout.flush();
  int scalarHandle = theSurface->getReadScalarHandle(scalarType);
  if (scalarHandle<0) return 1;
  for (int i=0; i<theSurface->numberOfVertices(); i++){
    double scalar;
    int scalarRead = theSurface->getScalar(scalarHandle, i, scalar);
    if (!scalarRead){
       // double *newColour = colourScheme.GetRGB(scalar);
       std::vector<double>  newColour = colourScheme.GetRGB(scalar);
      CXX_mot::CXXCoord colour (newColour[0]/255., newColour[1]/255., newColour[2]/255.);
      theSurface->setCoord("colour", i, colour);
    }
  }
  return 0;
}




coot::CColourScheme::CColourScheme() {
  defColour = 0;
  nTypes = 0;
}
coot::CColourScheme::~CColourScheme() {
}


int
coot::surface::interpolateIntoMap(const std::string &coordinateType, const std::string &scalarType, const clipper::NXmap<double> &aMap) {
   
   // cout << "In interpolateIntoMap\n"; cout.flush();
   int coordHandle = theSurface->getReadVectorHandle(coordinateType);
   if (coordHandle<0) return 1;
   else {
      double *scalarBuffer = new double[theSurface->numberOfVertices()];
      double coords[4];
      for (int i=0; i<theSurface->numberOfVertices(); i++){
	 //Use the colour if it has been assigned
	 if (theSurface->getCoord(coordHandle, i, coords)){
	    cout << "Bizarely no vertices for coordinate "<< i << endl;
	 }
	 clipper::Coord_orth orthogonals(coords[0], coords[1], coords[2]);
	 const clipper::Coord_map mapUnits(aMap.coord_map(orthogonals));
	 scalarBuffer[i] = aMap.interp<clipper::Interp_cubic>( mapUnits );
      }
      theSurface->addPerVertexScalar (scalarType, scalarBuffer);
      delete [] scalarBuffer;
   }
   return 0;
}


// ------------------------------------------------------------------------
//               
// ------------------------------------------------------------------------

int coot::CColourScheme::SetSchemeFloat ( const std::vector<float>& typ , const std::vector<std::string>& cols ) {
  if (typ.size() != cols.size() ) return 1;
  codes = RGBReps::GetColourNumbers (cols);
  nTypes =  cols.size();
  ranges = typ;
  colours = cols;
  mode = "float";
  return 0;
}

//-----------------------------------------------------------------------
std::vector<double> coot::CColourScheme::GetRGB(double value) {
//-----------------------------------------------------------------------
  int n;
  double frac;
  std::vector<double> col(4);
 
  if ( value < ranges[1] ) { 
     return RGBReps::GetColour(codes[0]);
  } else if ( value > ranges[nTypes-2] ) {
     return RGBReps::GetColour(codes[nTypes-1]);
  } else {
    for ( n = 1; n<nTypes-2;n++) {
      if ( value < ranges[n+1] ) {
        frac = (value - ranges[n])/ (ranges[n+1] - ranges[n]);
	if ( blend_mode == HSVCOLOURREP ) {
	    std::vector<double> col1 = RGBReps::GetColourHSV(codes[n]);        
	    std::vector<double> col2 = RGBReps::GetColourHSV(codes[n+1]);
	              for (int i=0;i<3;i++) 
	                col1[i] = (frac * col2[i]) + ((1.0-frac)*col1[i]);
	              col1[3] = 255.0;
           col = RGBReps::hsvtorgb(col1);
        } else {
	             std::vector<double> col1 = RGBReps::GetColour(codes[n]);        
	             std::vector<double> col2 = RGBReps::GetColour(codes[n+1]);
		     //printf("col1 %f %f %f\n",col1[0],col1[1],col1[2]);
		     // printf("col2 %f %f %f\n",col2[0],col2[1],col2[2]);
		     for (int i=0;i<3;i++) 
			  col[i] = (frac * col2[i]) + ((1.0-frac)*col1[i]);
		      col[3] = 255.0;
        }
      break;
      } 
    }
  }
  return col;
}

