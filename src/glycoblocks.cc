/*
   glycoblocks.cc: Moorhen (WebCoot) Molecular Graphics Program
   Derived from From pygl/build_tree_primitives.cc: CCP4MG Molecular Graphics Program
   Copyright (C) 2001-2008 University of York, CCLRC
   Copyright (C) 2009-2011 University of York
   Copyright (C) 2012-2023 STFC

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

// Ported to Coot from Moorhen (wasm_src/glycoblocks.cc, Stuart McNicholas) - the geometry code
// is unchanged, apart from being wrapped in the coot::glycoblocks namespace (so that its
// helper functions do not clash with others in Coot), the include paths, and an optional
// 3D-SNFG style (marked "Coot:") with cubes and spheres rather than flat blocks.

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#undef AddAtom
#endif

#include <iostream>
#include <set>
#include <stdlib.h>
#include "coords/Cartesian.hh"
#include "ccp4mg-utils/matrix.h"
#include <vector>
#include <utility>
#include <algorithm>
#include <map>
#include <math.h>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <mmdb2/mmdb_manager.h>
#include <sstream>
#include <iostream>
#include <iomanip>

#include "api/instancing.hh"
#include "glycoblocks.hh"
#include "coot-utils/oct.hh"
/*
//FIXME - only makes sense in context of molecules_container_t
#include "molecules_container.hh"
#include "coot-utils/coot-h-bonds.hh"
*/

#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif
#define PIBY2 (M_PI * 2)

// TODO - Generate mesh

/*
   cylinders - hbonds, covalent linkages
   cylinders - blocks GLC, BGC, GLA, GAL, MAN, BMA 
   square prism - blocks NDG, A2G, NGA, BM3, KDM, KDN, SIA, SLB
   triangular prism - blocks FUC, FUL
   star prism - blocks XYS, XYP
   square-edged triangular prism - blocks PA1, GCS, GCU, BDP, GTR, ADA, MAV, BEM, IDR

 */

namespace coot {
namespace glycoblocks {

using namespace mmdb;

std::vector<std::string> sugarResidues {"KDN","KDM","A2G", "NDG", "GAL","GLA","GLC","MAN","NGA","NAG","RAM","RIB","RIP","XYS","SUC","RER","BGC","LAT","LBT","BNG","BOG","AGC","GCD","IDT","BM3","BMA","GLO","DMU","LMT","TRE","FFC","XYP","SIA","SLB","FUC","FCA","FCB","FUL","GCS","GCU","BDP","GTR", "ADA", "MAV", "BEM", "IDR","PA1","TDX","CXY","RBY","XYL"};

bool isSugar(const std::string &resName){
    if (std::find(std::begin(sugarResidues), std::end(sugarResidues), resName) != std::end(sugarResidues))
        return true;
    return false;
}

bool isSugar(const char *resName){
    return isSugar(std::string(resName));
}

glm::mat4 get_orientation_matrix(const Cartesian &moving, const Cartesian &target){
    
    Cartesian v = Cartesian::CrossProduct(moving,target);
    double s = v.length();
    double c = dot_product(moving,target);
    //std::cout << "v: " << v << std::endl;
    //std::cout << "s: " << s << std::endl;

    glm::mat4 R(1.0f);

    if(c+1>1e-5){
        matrix vx(3,3);
        vx(0,0) = 0.0;
        vx(0,1) = -v.get_z();
        vx(0,2) =  v.get_y();

        vx(1,0) =  v.get_z();
        vx(1,1) = 0.0;
        vx(1,2) = -v.get_x();

        vx(2,0) = -v.get_y();
        vx(2,1) =  v.get_x();
        vx(2,2) = 0.0;

        matrix vxsq = vx * vx;
        vxsq = 1/(1+c) * vxsq;

        //std::cout << vx << std::endl;
        //std::cout << vxsq << std::endl;

        R[0][0] += vx(0,0) + vxsq(0,0);
        R[0][1] += vx(0,1) + vxsq(0,1);
        R[0][2] += vx(0,2) + vxsq(0,2);
        R[1][0] += vx(1,0) + vxsq(1,0);
        R[1][1] += vx(1,1) + vxsq(1,1);
        R[1][2] += vx(1,2) + vxsq(1,2);
        R[2][0] += vx(2,0) + vxsq(2,0);
        R[2][1] += vx(2,1) + vxsq(2,1);
        R[2][2] += vx(2,2) + vxsq(2,2);

    } else {
        R[1][1] = -1;
    }

    return R;
}

enum { NOLINE, LINE, DASHLINE, ARROW, DASHARROW, CYLINDER, CYLINDERARROW, CONE, DASHCYLINDER, DASHCYLINDERARROW };

std::string FloatToString(float n, int precision){
    std::stringstream stream;
    stream << std::fixed << std::setprecision(precision) << n;
    std::string s = stream.str();
    return s;
}

Cartesian MidPoint(const std::vector<Cartesian> &v1){
    Cartesian v;
    if(v1.size()==0) return v;
    for(unsigned i=0;i<v1.size();i++)
        v += v1[i];
    double frac = 1.0 / v1.size();
    return v.by_scalar(frac);
}

void SelectAminoNotHet (Manager *molHnd, int selHnd, mmdb::SELECTION_TYPE selType, int iModel, mmdb::cpstr Chains, int ResNo1, mmdb::cpstr Ins1,
        int ResNo2, mmdb::cpstr Ins2, mmdb::cpstr RNames, mmdb::cpstr ANames, mmdb::cpstr Elements, mmdb::cpstr altLocs, mmdb::SELECTION_KEY selKey ){
    int aminoSel = molHnd->NewSelection();

    molHnd->Select(aminoSel,selType,iModel,Chains,ResNo1,Ins1,ResNo2,Ins2,RNames,ANames,Elements,altLocs,mmdb::SKEY_NEW);

    if(selType==mmdb::STYPE_RESIDUE){
        int nSelRes;
        mmdb::Residue** selRes;
        molHnd->GetSelIndex(aminoSel,selRes,nSelRes);
        for(int ires=0;ires<nSelRes;ires++){
            if(selRes[ires]&&selRes[ires]->nAtoms>0&&(strncmp(selRes[ires]->name,"MSE",3)==0)){
                continue;
            }
            if(selRes[ires]&&selRes[ires]->nAtoms>0&&selRes[ires]->GetAtom(0)){
                if(selRes[ires]->GetAtom(0)->Het){
                    molHnd->SelectResidue(aminoSel,selRes[ires],mmdb::STYPE_RESIDUE,mmdb::SKEY_CLR,false);
                }
            }
        }
        molHnd->MakeSelIndex(aminoSel);
    }else if(selType==mmdb::STYPE_ATOM){
        int nSelAtoms;
        mmdb::Atom** selAtoms;
        molHnd->GetSelIndex(aminoSel,selAtoms,nSelAtoms);
        for(int iat=0;iat<nSelAtoms;iat++){
            if(selAtoms[iat]){
                if(selAtoms[iat]->Het&&(strncmp(selAtoms[iat]->residue->name,"MSE",3)!=0)){
                    molHnd->SelectAtom(aminoSel,selAtoms[iat],mmdb::SKEY_CLR,false);
                }
            }
        }
        molHnd->MakeSelIndex(aminoSel);
    }
    // else do nothing ... Should we do STYPE_CHAIN and STYPE_MODEL?

    molHnd->Select(selHnd,selType,aminoSel,selKey);

    molHnd->DeleteSelection(aminoSel);
}

class glyco_shapes_t {
  public:
    glyco_shapes_t(){
        // SQUARE
        std::vector<coot::api::vn_vertex> square_vertices;
        square_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1, -1,  1),glm::vec3( 0,  0,  1)));
        square_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1, -1,  1),glm::vec3( 0,  0,  1)));
        square_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1,  1,  1),glm::vec3( 0,  0,  1)));
        square_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1,  1,  1),glm::vec3( 0,  0,  1)));
        square_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1, -1, -1),glm::vec3( 0,  0, -1)));
        square_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1, -1, -1),glm::vec3( 0,  0, -1)));
        square_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1,  1, -1),glm::vec3( 0,  0, -1)));
        square_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1,  1, -1),glm::vec3( 0,  0, -1)));

        square_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1,  1,  1),glm::vec3( 0,  1,  0)));
        square_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1,  1,  1),glm::vec3( 0,  1,  0)));
        square_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1,  1, -1),glm::vec3( 0,  1,  0)));
        square_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1,  1, -1),glm::vec3( 0,  1,  0)));
        square_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1, -1,  1),glm::vec3( 0, -1,  0)));
        square_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1, -1,  1),glm::vec3( 0, -1,  0)));
        square_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1, -1, -1),glm::vec3( 0, -1,  0)));
        square_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1, -1, -1),glm::vec3( 0, -1,  0)));

        square_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1, -1,  1),glm::vec3( 1,  0,  0)));
        square_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1, -1, -1),glm::vec3( 1,  0,  0)));
        square_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1,  1, -1),glm::vec3( 1,  0,  0)));
        square_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1,  1,  1),glm::vec3( 1,  0,  0)));
        square_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1, -1,  1),glm::vec3(-1,  0,  0)));
        square_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1, -1, -1),glm::vec3(-1,  0,  0)));
        square_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1,  1, -1),glm::vec3(-1,  0,  0)));
        square_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1,  1,  1),glm::vec3(-1,  0,  0)));

        std::vector<g_triangle> square_triangles;
        square_triangles.push_back(g_triangle(  0,  1,  2));
        square_triangles.push_back(g_triangle(  0,  2,  3));
        square_triangles.push_back(g_triangle(  5,  4,  6));
        square_triangles.push_back(g_triangle(  6,  4,  7));
        square_triangles.push_back(g_triangle(  8,  9, 10));
        square_triangles.push_back(g_triangle(  8, 10, 11));
        square_triangles.push_back(g_triangle( 13, 12, 14));
        square_triangles.push_back(g_triangle( 14, 12, 15));
        square_triangles.push_back(g_triangle( 16, 17, 18));
        square_triangles.push_back(g_triangle( 16, 18, 19));
        square_triangles.push_back(g_triangle( 21, 20, 22));
        square_triangles.push_back(g_triangle( 22, 20, 23));

        square_geom.vertices = square_vertices;
        square_geom.triangles = square_triangles;
        square_geom.name = "square_mesh";

        // CIRCLE
        std::vector<coot::api::vn_vertex> circle_vertices;
        std::vector<g_triangle> circle_triangles;
        int k = 0;
        int accu = 36;
        for(int j=0;j<360;j=j+360/accu,k+=10){
            double theta = (double)j/360.0 * PIBY2;
            double theta2 = (double)(j+360/accu)/360.0 * PIBY2;
            double x1 = cos(theta);
            double y1 = sin(theta);
            double x2 = cos(theta2);
            double y2 = sin(theta2);
            circle_vertices.push_back(coot::api::vn_vertex(glm::vec3(  0,   0,  1),glm::vec3( 0,  0,  1)));
            circle_vertices.push_back(coot::api::vn_vertex(glm::vec3( x1,  y1,  1),glm::vec3( 0,  0,  1)));
            circle_vertices.push_back(coot::api::vn_vertex(glm::vec3( x2,  y2,  1),glm::vec3( 0,  0,  1)));
            circle_triangles.push_back(g_triangle( k,  k+1,  k+2));
            circle_vertices.push_back(coot::api::vn_vertex(glm::vec3(  0,   0, -1),glm::vec3( 0,  0, -1)));
            circle_vertices.push_back(coot::api::vn_vertex(glm::vec3( x1,  y1, -1),glm::vec3( 0,  0, -1)));
            circle_vertices.push_back(coot::api::vn_vertex(glm::vec3( x2,  y2, -1),glm::vec3( 0,  0, -1)));
            circle_triangles.push_back(g_triangle( k+4,  k+3,  k+5));
            circle_vertices.push_back(coot::api::vn_vertex(glm::vec3( x1,  y1, -1),glm::vec3( x1, y1,  0)));
            circle_vertices.push_back(coot::api::vn_vertex(glm::vec3( x2,  y2, -1),glm::vec3( x2, y2,  0)));
            circle_vertices.push_back(coot::api::vn_vertex(glm::vec3( x2,  y2,  1),glm::vec3( x2, y2,  0)));
            circle_vertices.push_back(coot::api::vn_vertex(glm::vec3( x1,  y1,  1),glm::vec3( x1, y1,  0)));
            circle_triangles.push_back(g_triangle( k+6,  k+7,  k+8));
            circle_triangles.push_back(g_triangle( k+6,  k+8,  k+9));
        }

        circle_geom.vertices = circle_vertices;
        circle_geom.triangles = circle_triangles;
        circle_geom.name = "circle_mesh";

        // TRIANGLE
        std::vector<coot::api::vn_vertex> triangle_vertices;
        std::vector<g_triangle> triangle_triangles;
        triangle_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1, -0.8660,  1),glm::vec3( 0,  0,  1)));
        triangle_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1, -0.8660,  1),glm::vec3( 0,  0,  1)));
        triangle_vertices.push_back(coot::api::vn_vertex(glm::vec3( 0,  0.8660,  1),glm::vec3( 0,  0,  1)));
        triangle_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1, -0.8660, -1),glm::vec3( 0,  0, -1)));
        triangle_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1, -0.8660, -1),glm::vec3( 0,  0, -1)));
        triangle_vertices.push_back(coot::api::vn_vertex(glm::vec3( 0,  0.8660, -1),glm::vec3( 0,  0, -1)));

        triangle_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1, -0.8660, -1),glm::vec3( 0, -1,  0)));
        triangle_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1, -0.8660, -1),glm::vec3( 0, -1,  0)));
        triangle_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1, -0.8660,  1),glm::vec3( 0, -1,  0)));
        triangle_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1, -0.8660,  1),glm::vec3( 0, -1,  0)));

        triangle_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1, -0.8660, -1),glm::vec3( 0.8660,  0.5,  0)));
        triangle_vertices.push_back(coot::api::vn_vertex(glm::vec3( 0,  0.8660, -1),glm::vec3( 0.8660,  0.5,  0)));
        triangle_vertices.push_back(coot::api::vn_vertex(glm::vec3( 0,  0.8660,  1),glm::vec3( 0.8660,  0.5,  0)));
        triangle_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1, -0.8660,  1),glm::vec3( 0.8660,  0.5,  0)));

        triangle_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1, -0.8660, -1),glm::vec3(-0.8660,  0.5,  0)));
        triangle_vertices.push_back(coot::api::vn_vertex(glm::vec3( 0,  0.8660, -1),glm::vec3(-0.8660,  0.5,  0)));
        triangle_vertices.push_back(coot::api::vn_vertex(glm::vec3( 0,  0.8660,  1),glm::vec3(-0.8660,  0.5,  0)));
        triangle_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1, -0.8660,  1),glm::vec3(-0.8660,  0.5,  0)));

        triangle_triangles.push_back(g_triangle(  0,  1,  2));
        triangle_triangles.push_back(g_triangle(  4,  3,  5));
        triangle_triangles.push_back(g_triangle(  6,  7,  8));
        triangle_triangles.push_back(g_triangle(  6,  8,  9));
        triangle_triangles.push_back(g_triangle( 10, 11, 12));
        triangle_triangles.push_back(g_triangle( 10, 12, 13));
        triangle_triangles.push_back(g_triangle( 15, 14, 16));
        triangle_triangles.push_back(g_triangle( 16, 14, 17));

        triangle_geom.vertices = triangle_vertices;
        triangle_geom.triangles = triangle_triangles;
        triangle_geom.name = "triangle_mesh";

        // PENTAGRAM
        // http://mathworld.wolfram.com/RegularPentagon.html
        // https://mathworld.wolfram.com/Pentagram.html R / rho
        double c1 = cos(2.0*M_PI/5.);
        double c2 = cos(M_PI/5.);
        double s1 = sin(2.0*M_PI/5.);
        double s2 = sin(4.0*M_PI/5.);
        double RoverRho = 0.3819660112501046; // sqrt(25 - 11*sqrt(5)) / sqrt(5 - sqrt(5))

        std::vector<coot::api::vn_vertex> pentagram_vertices;
        std::vector<g_triangle> pentagram_triangles;

        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(           s1,           c1,  1),glm::vec3( 0,  0,  1)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(          -s1,           c1,  1),glm::vec3( 0,  0,  1)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(            0,    -RoverRho,  1),glm::vec3( 0,  0,  1)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(            0,            1,  1),glm::vec3( 0,  0,  1)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3( -s2*RoverRho,  c2*RoverRho,  1),glm::vec3( 0,  0,  1)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(  s2*RoverRho,  c2*RoverRho,  1),glm::vec3( 0,  0,  1)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(          -s2,          -c2,  1),glm::vec3( 0,  0,  1)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(            0,    -RoverRho,  1),glm::vec3( 0,  0,  1)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3( -s1*RoverRho, -c1*RoverRho,  1),glm::vec3( 0,  0,  1)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(           s2,          -c2,  1),glm::vec3( 0,  0,  1)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(  s1*RoverRho, -c1*RoverRho,  1),glm::vec3( 0,  0,  1)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(            0,    -RoverRho,  1),glm::vec3( 0,  0,  1)));

        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(           s1,           c1, -1),glm::vec3( 0,  0, -1)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(          -s1,           c1, -1),glm::vec3( 0,  0, -1)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(            0,    -RoverRho, -1),glm::vec3( 0,  0, -1)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(            0,            1, -1),glm::vec3( 0,  0, -1)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3( -s2*RoverRho,  c2*RoverRho, -1),glm::vec3( 0,  0, -1)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(  s2*RoverRho,  c2*RoverRho, -1),glm::vec3( 0,  0, -1)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(          -s2,          -c2, -1),glm::vec3( 0,  0, -1)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(            0,    -RoverRho, -1),glm::vec3( 0,  0, -1)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3( -s1*RoverRho, -c1*RoverRho, -1),glm::vec3( 0,  0, -1)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(           s2,          -c2, -1),glm::vec3( 0,  0, -1)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(  s1*RoverRho, -c1*RoverRho, -1),glm::vec3( 0,  0, -1)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(            0,    -RoverRho, -1),glm::vec3( 0,  0, -1)));

        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(          -s2,          -c2, -1),glm::vec3(  s2, -c2,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(            0,    -RoverRho, -1),glm::vec3(  s2, -c2,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(            0,    -RoverRho,  1),glm::vec3(  s2, -c2,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(          -s2,          -c2,  1),glm::vec3(  s2, -c2,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3( -s1*RoverRho, -c1*RoverRho, -1),glm::vec3( -s1,  c1,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(          -s2,          -c2, -1),glm::vec3( -s1,  c1,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(          -s2,          -c2,  1),glm::vec3( -s1,  c1,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3( -s1*RoverRho, -c1*RoverRho,  1),glm::vec3( -s1,  c1,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(          -s1,           c1, -1),glm::vec3( -s2, -c2,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3( -s1*RoverRho, -c1*RoverRho, -1),glm::vec3( -s2, -c2,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3( -s1*RoverRho, -c1*RoverRho,  1),glm::vec3( -s2, -c2,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(          -s1,           c1,  1),glm::vec3( -s2, -c2,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3( -s2*RoverRho,  c2*RoverRho, -1),glm::vec3( 0,  1,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(          -s1,           c1, -1),glm::vec3( 0,  1,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(          -s1,           c1,  1),glm::vec3( 0,  1,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3( -s2*RoverRho,  c2*RoverRho,  1),glm::vec3( 0,  1,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(            0,            1, -1),glm::vec3( -s1,  c1,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3( -s2*RoverRho,  c2*RoverRho, -1),glm::vec3( -s1,  c1,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3( -s2*RoverRho,  c2*RoverRho,  1),glm::vec3( -s1,  c1,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(            0,            1,  1),glm::vec3( -s1,  c1,  0)));

        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(           s2,          -c2, -1),glm::vec3( -s2, -c2,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(            0,    -RoverRho, -1),glm::vec3( -s2, -c2,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(            0,    -RoverRho,  1),glm::vec3( -s2, -c2,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(           s2,          -c2,  1),glm::vec3( -s2, -c2,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(  s1*RoverRho, -c1*RoverRho, -1),glm::vec3(  s1,  c1,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(           s2,          -c2, -1),glm::vec3(  s1,  c1,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(           s2,          -c2,  1),glm::vec3(  s1,  c1,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(  s1*RoverRho, -c1*RoverRho,  1),glm::vec3(  s1,  c1,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(           s1,           c1, -1),glm::vec3(  s2, -c2,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(  s1*RoverRho, -c1*RoverRho, -1),glm::vec3(  s2, -c2,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(  s1*RoverRho, -c1*RoverRho,  1),glm::vec3(  s2, -c2,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(           s1,           c1,  1),glm::vec3(  s2, -c2,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(  s2*RoverRho,  c2*RoverRho, -1),glm::vec3( 0,  1,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(           s1,           c1, -1),glm::vec3( 0,  1,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(           s1,           c1,  1),glm::vec3( 0,  1,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(  s2*RoverRho,  c2*RoverRho,  1),glm::vec3( 0,  1,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(            0,            1, -1),glm::vec3(  s1,  c1,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(  s2*RoverRho,  c2*RoverRho, -1),glm::vec3(  s1,  c1,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(  s2*RoverRho,  c2*RoverRho,  1),glm::vec3(  s1,  c1,  0)));
        pentagram_vertices.push_back(coot::api::vn_vertex(glm::vec3(            0,            1,  1),glm::vec3(  s1,  c1,  0)));

        pentagram_triangles.push_back(g_triangle(  0,  1,  2));
        pentagram_triangles.push_back(g_triangle(  3,  4,  5));
        pentagram_triangles.push_back(g_triangle(  6,  7,  8));
        pentagram_triangles.push_back(g_triangle(  9, 10, 11));

        pentagram_triangles.push_back(g_triangle(  13,  12,  14));
        pentagram_triangles.push_back(g_triangle(  16,  15,  17));
        pentagram_triangles.push_back(g_triangle(  19,  18,  20));
        pentagram_triangles.push_back(g_triangle(  22,  21,  23));

        pentagram_triangles.push_back(g_triangle(  24,  25,  26));
        pentagram_triangles.push_back(g_triangle(  24,  26,  27));
        pentagram_triangles.push_back(g_triangle(  28,  29,  30));
        pentagram_triangles.push_back(g_triangle(  28,  30,  31));
        pentagram_triangles.push_back(g_triangle(  32,  33,  34));
        pentagram_triangles.push_back(g_triangle(  32,  34,  35));
        pentagram_triangles.push_back(g_triangle(  36,  37,  38));
        pentagram_triangles.push_back(g_triangle(  36,  38,  39));
        pentagram_triangles.push_back(g_triangle(  40,  41,  42));
        pentagram_triangles.push_back(g_triangle(  40,  42,  43));

        pentagram_triangles.push_back(g_triangle(  45,  44,  46));
        pentagram_triangles.push_back(g_triangle(  46,  44,  47));
        pentagram_triangles.push_back(g_triangle(  49,  48,  50));
        pentagram_triangles.push_back(g_triangle(  50,  48,  51));
        pentagram_triangles.push_back(g_triangle(  53,  52,  54));
        pentagram_triangles.push_back(g_triangle(  54,  52,  55));
        pentagram_triangles.push_back(g_triangle(  57,  56,  58));
        pentagram_triangles.push_back(g_triangle(  58,  56,  59));
        pentagram_triangles.push_back(g_triangle(  61,  60,  62));
        pentagram_triangles.push_back(g_triangle(  62,  60,  63));

        pentagram_geom.vertices = pentagram_vertices;
        pentagram_geom.triangles = pentagram_triangles;
        pentagram_geom.name = "pentagram_mesh";

        // BI-SQUARE-1
        std::vector<coot::api::vn_vertex> bi_square_1_vertices;
        bi_square_1_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1, -1,  1),glm::vec3( 0,  0,  1)));
        bi_square_1_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1, -1,  1),glm::vec3( 0,  0,  1)));
        bi_square_1_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1,  1,  1),glm::vec3( 0,  0,  1)));

        bi_square_1_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1, -1, -1),glm::vec3( 0,  0, -1)));
        bi_square_1_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1, -1, -1),glm::vec3( 0,  0, -1)));
        bi_square_1_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1,  1, -1),glm::vec3( 0,  0, -1)));

        bi_square_1_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1, -1,  1),glm::vec3( 0, -1,  0)));
        bi_square_1_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1, -1,  1),glm::vec3( 0, -1,  0)));
        bi_square_1_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1, -1, -1),glm::vec3( 0, -1,  0)));
        bi_square_1_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1, -1, -1),glm::vec3( 0, -1,  0)));

        bi_square_1_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1, -1,  1),glm::vec3( 1,  0,  0)));
        bi_square_1_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1, -1, -1),glm::vec3( 1,  0,  0)));
        bi_square_1_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1,  1, -1),glm::vec3( 1,  0,  0)));
        bi_square_1_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1,  1,  1),glm::vec3( 1,  0,  0)));

        std::vector<g_triangle> bi_square_1_triangles;
        bi_square_1_triangles.push_back(g_triangle(  0,  1,  2));
        bi_square_1_triangles.push_back(g_triangle(  4,  3,  5));

        bi_square_1_triangles.push_back(g_triangle(  7,  6,  8));
        bi_square_1_triangles.push_back(g_triangle(  8,  6,  9));

        bi_square_1_triangles.push_back(g_triangle(  10,  11,  12));
        bi_square_1_triangles.push_back(g_triangle(  10,  12,  13));

        bi_square_1_geom.vertices = bi_square_1_vertices;
        bi_square_1_geom.triangles = bi_square_1_triangles;
        bi_square_1_geom.name = "bi_square_1_mesh";

        // BI-SQUARE-2
        std::vector<coot::api::vn_vertex> bi_square_2_vertices;
        bi_square_2_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1,  1,  1),glm::vec3( 0,  0,  1)));
        bi_square_2_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1,  1,  1),glm::vec3( 0,  0,  1)));
        bi_square_2_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1, -1,  1),glm::vec3( 0,  0,  1)));

        bi_square_2_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1,  1, -1),glm::vec3( 0,  0, -1)));
        bi_square_2_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1,  1, -1),glm::vec3( 0,  0, -1)));
        bi_square_2_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1, -1, -1),glm::vec3( 0,  0, -1)));

        bi_square_2_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1,  1,  1),glm::vec3( 0,  1,  0)));
        bi_square_2_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1,  1,  1),glm::vec3( 0,  1,  0)));
        bi_square_2_vertices.push_back(coot::api::vn_vertex(glm::vec3( 1,  1, -1),glm::vec3( 0,  1,  0)));
        bi_square_2_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1,  1, -1),glm::vec3( 0,  1,  0)));

        bi_square_2_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1, -1,  1),glm::vec3(-1,  0,  0)));
        bi_square_2_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1, -1, -1),glm::vec3(-1,  0,  0)));
        bi_square_2_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1,  1, -1),glm::vec3(-1,  0,  0)));
        bi_square_2_vertices.push_back(coot::api::vn_vertex(glm::vec3(-1,  1,  1),glm::vec3(-1,  0,  0)));

        std::vector<g_triangle> bi_square_2_triangles;
        bi_square_2_triangles.push_back(g_triangle(  0,  1,  2));
        bi_square_2_triangles.push_back(g_triangle(  4,  3,  5));

        bi_square_2_triangles.push_back(g_triangle(  6,  7,  8));
        bi_square_2_triangles.push_back(g_triangle(  6,  8,  9));

        bi_square_2_triangles.push_back(g_triangle(  11,  10,  12));
        bi_square_2_triangles.push_back(g_triangle(  12,  10,  13));

        bi_square_2_geom.vertices = bi_square_2_vertices;
        bi_square_2_geom.triangles = bi_square_2_triangles;
        bi_square_2_geom.name = "bi_square_2_mesh";

        std::vector<coot::api::vn_vertex> sphere_vertices;
        std::vector<g_triangle> sphere_triangles;
        // Coot: a real unit sphere. (In Moorhen this was a dummy, as its renderer draws the
        // "spherical-atoms" geometry as impostor spheres, but Coot draws the triangles.)
        std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > octasphere = tessellate_octasphere(3);
        for (const auto &p : octasphere.first)
            sphere_vertices.push_back(coot::api::vn_vertex(p, p)); // unit sphere: the normal is the position
        sphere_triangles = octasphere.second;

        sphere_geom.vertices = sphere_vertices;
        sphere_geom.triangles = sphere_triangles;
        sphere_geom.name = "spherical-atoms";

    }

    coot::instanced_geometry_t square_geom;
    coot::instanced_geometry_t circle_geom;
    coot::instanced_geometry_t triangle_geom;
    coot::instanced_geometry_t pentagram_geom;
    coot::instanced_geometry_t bi_square_1_geom;
    coot::instanced_geometry_t bi_square_2_geom;
    coot::instanced_geometry_t sphere_geom;
    // Coot: 3D-SNFG style - draw the circle shapes (hexoses) as spheres rather than discs
    bool spheres_for_circles = false;
};

std::vector<double> DistanceBetweenPointAndLine(const Cartesian &ls, const Cartesian &le, const Cartesian &p){

    std::vector<double> ret(3);
    ret[0] = -1.0;
    ret[1] = -1.0;
    ret[2] = -1.0;

    double linesize = (le-ls).length();
    if(fabs(linesize)<1e-6){
        printf("Zero length line in DistanceBetweenPointAndLine\n");
        return ret;
    }
    /* t is value in line equation p = p1 + t(p2-p1) */
    double t = dot_product(p-ls,le-ls) / (linesize*linesize);

    Cartesian pt = ls + (le-ls).by_scalar(t);

    ret[0] = (pt-p).length();
    ret[1] = t;
    return ret;

}

Cartesian MidPoint(const Cartesian &v1, const Cartesian &v2){
    Cartesian v;
    double x = v1.x()+(v2.x()-v1.x())/2;
    double y = v1.y()+(v2.y()-v1.y())/2;
    double z = v1.z()+(v2.z()-v1.z())/2;
    v.set_them(x,y,z);
    return v;
}


std::vector<double> DistanceBetweenTwoLines(const Cartesian &p1, const Cartesian &p2, const Cartesian &p3, const Cartesian &p4){

    Cartesian a1 = p2-p1;
    Cartesian a2 = p4-p3;

    std::vector<double> ret(3);
    ret[0] = -1.0;
    ret[1] = -1.0;
    ret[2] = -1.0;

    if(a1.length()==0)
        return ret;

    Cartesian dum;

    Cartesian n = dum.CrossProduct(a1,a2);
    if(fabs(n.length())<1e-6){
        return DistanceBetweenPointAndLine(p1,p2,p3);
    }

    if(a2.length()==0)
        return ret;

    n.normalize();

    double dist = fabs(dot_product(p4-p1,n));

    double a1sq = a1.length()*a1.length();
    double a2sq = a2.length()*a2.length();

    double u = (a1sq*(dot_product(a2,p3) - dot_product(a2,p1)) + 
            dot_product(a1,a2)*(dot_product(a1,p1) - dot_product(a1,p3)))/
        (dot_product(a1,a2)*dot_product(a1,a2) - a1sq*a2sq);

    double t = (dot_product(a1,p3) + u * dot_product(a1,a2) - dot_product(a1,p1)) / a1sq;

    ret[0] = dist;
    ret[1] = t;
    ret[2] = u;

    return ret;

}

std::vector<Cartesian> DrawSugarBlockInt(mmdb::Residue* res1, int selHnd, mmdb::Atom* at1,mmdb::Atom* at2,mmdb::Atom* at3,mmdb::Atom* at4,mmdb::Atom* at5,mmdb::Atom* at6, mmdb::Residue* res, mmdb::Manager *molHnd, bool two_colour, float sugar_block_thickness, float sugar_block_scale, glyco_shapes_t &glyco_shapes){

    // FIXME, likely broken for F6P,RIB and friends.
    //std::cout << "Possibly drawing " << res->GetResName() << std::endl;
    std::vector<Cartesian> retval;

    const auto blue_col = glm::vec4(0,0.56,0.73,1);
    const auto green_col = glm::vec4(0,0.65,0.31,1);
    const auto yellow_col = glm::vec4(0.8,0.63,0,1);
    const auto red_col = glm::vec4(0.92,0.1,0.14,1);
    const auto orange_col = glm::vec4(1,0.49,0,1);

    // Should not be able to get here with these not already registered.
    int udd_C1X = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C1X" );
    int udd_C1Y = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C1Y" );
    int udd_C1Z = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C1Z" );
    int udd_C3X = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C3X" );
    int udd_C3Y = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C3Y" );
    int udd_C3Z = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C3Z" );
    int udd_C4X = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C4X" );
    int udd_C4Y = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C4Y" );
    int udd_C4Z = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C4Z" );
    int udd_C5X = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C5X" );
    int udd_C5Y = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C5Y" );
    int udd_C5Z = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C5Z" );
    int udd_C2X = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C2X" );
    int udd_C2Y = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C2Y" );
    int udd_C2Z = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C2Z" );

    if(at1&&at2&&at3&&at4&&at5){
        //std::cout << "DrawSugarBlockInt have at least 5 atoms.\n";
        Cartesian xaxis(1,0,0);
        Cartesian yaxis(0,1,0);
        Cartesian zaxis(0,0,1);
        std::vector<Cartesian> norm_vec;
        std::vector<Cartesian> mid_vec;
        Cartesian cat1 = Cartesian(at1->x,at1->y,at1->z);
        Cartesian cat2 = Cartesian(at2->x,at2->y,at2->z);
        Cartesian cat3 = Cartesian(at3->x,at3->y,at3->z);
        Cartesian cat4 = Cartesian(at4->x,at4->y,at4->z);
        Cartesian cat5 = Cartesian(at5->x,at5->y,at5->z);
        res->PutUDData(udd_C1X,cat1.get_x()); res->PutUDData(udd_C1Y,cat1.get_y()); res->PutUDData(udd_C1Z,cat1.get_z());
        res->PutUDData(udd_C2X,cat2.get_x()); res->PutUDData(udd_C2Y,cat2.get_y()); res->PutUDData(udd_C2Z,cat2.get_z());
        res->PutUDData(udd_C3X,cat3.get_x()); res->PutUDData(udd_C3Y,cat3.get_y()); res->PutUDData(udd_C3Z,cat3.get_z());
        res->PutUDData(udd_C4X,cat4.get_x()); res->PutUDData(udd_C4Y,cat4.get_y()); res->PutUDData(udd_C4Z,cat4.get_z());
        res->PutUDData(udd_C5X,cat5.get_x()); res->PutUDData(udd_C5Y,cat5.get_y()); res->PutUDData(udd_C5Z,cat5.get_z());

        double d12 = (cat2-cat1).length();
        double d23 = (cat3-cat2).length();
        double d34 = (cat4-cat3).length();
        double d45 = (cat5-cat4).length();
        Cartesian c123 = Cartesian::CrossProduct((cat2-cat1),(cat3-cat2));
        Cartesian c234 = Cartesian::CrossProduct((cat3-cat2),(cat4-cat3));
        Cartesian c345 = Cartesian::CrossProduct((cat4-cat3),(cat5-cat4));
        if((d12>1.0&&d12<3.2)&&(d23>1.0&&d23<3.2)&&(d34>1.0&&d34<3.2)&&(d45>1.0&&d45<3.2)){
            if(at6){
                //std::cout << "Six membered test.\n";
                Cartesian cat6 = Cartesian(at6->x,at6->y,at6->z);
                double d56 = (cat6-cat5).length();
                double d61 = (cat6-cat1).length();
                if((d61>1.0&&d61<3.2)&&(d56>1.0&&d56<3.2)){
                    norm_vec.push_back(c123);
                    norm_vec.push_back(c234);
                    norm_vec.push_back(c345);
                    mid_vec.push_back(cat1);
                    mid_vec.push_back(cat2);
                    mid_vec.push_back(cat3);
                    mid_vec.push_back(cat4);
                    mid_vec.push_back(cat5);
                    Cartesian c456 = Cartesian::CrossProduct((cat5-cat4),(cat6-cat5));
                    Cartesian c561 = Cartesian::CrossProduct((cat6-cat5),(cat1-cat6));
                    norm_vec.push_back(c456);
                    norm_vec.push_back(c561);
                    mid_vec.push_back(cat6);
                    //std::cout << "We have 6 membered ring\n";
                } else if(d61>4.0&&d61<7.0&&strncmp(res->GetResName(),"XLS",3)==0) {
                    norm_vec.push_back(c123);
                    norm_vec.push_back(c234);
                    norm_vec.push_back(c345);
                    mid_vec.push_back(cat1);
                    mid_vec.push_back(cat2);
                    mid_vec.push_back(cat3);
                    mid_vec.push_back(cat4);
                    mid_vec.push_back(cat5);
                    Cartesian c456 = Cartesian::CrossProduct((cat5-cat4),(cat6-cat5));
                    Cartesian c561 = Cartesian::CrossProduct((cat6-cat5),(cat1-cat6));
                    norm_vec.push_back(c456);
                    norm_vec.push_back(c561);
                    mid_vec.push_back(cat6);
                    //std::cout << "Linear?\n";
                }
            } else {
                //std::cout << "Five membered test.\n";
                double d51 = (cat5-cat1).length();
                if((d51>1.0&&d51<3.2)){
                    norm_vec.push_back(c123);
                    norm_vec.push_back(c234);
                    norm_vec.push_back(c345);
                    mid_vec.push_back(cat1);
                    mid_vec.push_back(cat2);
                    mid_vec.push_back(cat3);
                    mid_vec.push_back(cat4);
                    mid_vec.push_back(cat5);
                    Cartesian c451 = Cartesian::CrossProduct((cat5-cat4),(cat1-cat5));
                    norm_vec.push_back(c451);
                    //std::cout << "We have 5 membered ring\n";
                }
            }
        } else {
            //std::cout << "Fail first test\n";
        }

        float sugar_block_radius = 1.4 * sugar_block_scale;

        if(norm_vec.size()>0){
            //std::cout << "We have some normals" << std::endl;
            Cartesian normal = MidPoint(norm_vec);
            normal.normalize();
            //std::cout << "We have some coords" << std::endl;
            //for(int i=0;i<mid_vec.size();i++){
                //std::cout << mid_vec[i] << std::endl;
            //}
            Cartesian centre = MidPoint(mid_vec);
            retval.push_back(centre);
            retval.push_back(normal);
            //std::cout << "res name tests " << res->GetResName() << std::endl;
            if(strncmp(res->GetResName(),"PA1",3)==0||strncmp(res->GetResName(),"BM3",3)==0||strncmp(res->GetResName(),"NAC",3)==0||strncmp(res->GetResName(),"NDG",3)==0||strncmp(res->GetResName(),"NAG",3)==0||strncmp(res->GetResName(),"NBG",3)==0||strncmp(res->GetResName(),"NGA",3)==0||strncmp(res->GetResName(),"NG1",3)==0||strncmp(res->GetResName(),"NG6",3)==0||strncmp(res->GetResName(),"GCS",3)==0||strncmp(res->GetResName(),"6MN",3)==0||strncmp(res->GetResName(),"GLP",3)==0||strncmp(res->GetResName(),"GP4",3)==0||strncmp(res->GetResName(),"A2G",3)==0){

                bool diagonal = false;
                if(strncmp(res->GetResName(),"GCS",3)==0||strncmp(res->GetResName(),"PA1",3)==0) {
                    diagonal = true;
                }
                two_colour = two_colour && diagonal;

                //std::cout << "should draw a glycoblock!! " << res->GetResName() << std::endl;
                Cartesian c1c4 = cat1-cat4;
                c1c4.normalize();

                glm::vec4 block_col;
                Cartesian c1c4p = Cartesian::CrossProduct(c1c4,normal);
                c1c4p.normalize();
                Cartesian normalp = Cartesian::CrossProduct(c1c4p,c1c4);

                glm::mat4 R = get_orientation_matrix(normalp,Cartesian(0.0,0.0,1.0));
                glm::vec4 xprime = R * glm::vec4(1.0,0.0,0.0,1.0);
                glm::mat4 R2 = get_orientation_matrix(c1c4,Cartesian(-xprime[0],-xprime[1],-xprime[2]));
                glm::mat4 R3 = R2*R;

                glm::vec3 block_pos(centre.get_x(),centre.get_y(),centre.get_z());
                if(strncmp(res->GetResName(),"GCS",3)==0||strncmp(res->GetResName(),"PA1",3)==0||strncmp(res->GetResName(),"NDG",3)==0||strncmp(res->GetResName(),"NAG",3)==0||strncmp(res->GetResName(),"NAC",3)==0||strncmp(res->GetResName(),"NBG",3)==0)
                    block_col = blue_col;
                if(strncmp(res->GetResName(),"A2G",3)==0||strncmp(res->GetResName(),"NGA",3)==0||strncmp(res->GetResName(),"NG1",3)==0||strncmp(res->GetResName(),"NG6",3)==0)
                    block_col = yellow_col;
                if(strncmp(res->GetResName(),"BM3",3)==0)
                    block_col = green_col;

                glm::vec3 block_size(sugar_block_scale,sugar_block_scale,sugar_block_thickness);
                coot::instancing_data_type_B_t block(block_pos, block_col, block_size, R3);
                if(diagonal){
                    glyco_shapes.bi_square_1_geom.instancing_data_B.push_back(block);
                    block_col = glm::vec4(1,1,1,1);
                    coot::instancing_data_type_B_t block2(block_pos, block_col, block_size, R3);
                    glyco_shapes.bi_square_2_geom.instancing_data_B.push_back(block2);
                } else {
                    glyco_shapes.square_geom.instancing_data_B.push_back(block);
                }

            } else if(strncmp(res->GetResName(),"BMA",3)==0||strncmp(res->GetResName(),"MAN",3)==0||strncmp(res->GetResName(),"GLA",3)==0||strncmp(res->GetResName(),"GAL",3)==0||strncmp(res->GetResName(),"GLC",3)==0||strncmp(res->GetResName(),"BGC",3)==0){

                glm::vec4 circle_col;

                Cartesian c1c4 = cat1-cat4;
                c1c4.normalize();

                Cartesian c1c4p = Cartesian::CrossProduct(c1c4,normal);
                c1c4p.normalize();
                Cartesian normalp = Cartesian::CrossProduct(c1c4p,c1c4);

                glm::mat4 R = get_orientation_matrix(normalp,Cartesian(0.0,0.0,1.0));

                glm::vec3 circle_pos(centre.get_x(),centre.get_y(),centre.get_z());
                if(strncmp(res->GetResName(),"GLC",3)==0||strncmp(res->GetResName(),"BGC",3)==0)
                    circle_col = blue_col;
                if(strncmp(res->GetResName(),"GAL",3)==0||strncmp(res->GetResName(),"GLA",3)==0)
                    circle_col = yellow_col;
                if(strncmp(res->GetResName(),"MAN",3)==0||strncmp(res->GetResName(),"BMA",3)==0)
                    circle_col = green_col;

                if (glyco_shapes.spheres_for_circles) {
                    // Coot: 3D-SNFG style - a sphere with the radius of the disc
                    glm::vec3 sphere_size(sugar_block_scale,sugar_block_scale,sugar_block_scale);
                    coot::instancing_data_type_A_t sphere(circle_pos, circle_col, sphere_size);
                    glyco_shapes.sphere_geom.instancing_data_A.push_back(sphere);
                } else {
                    glm::vec3 circle_size(sugar_block_scale,sugar_block_scale,sugar_block_thickness);
                    coot::instancing_data_type_B_t test_circle(circle_pos, circle_col, circle_size, R);
                    glyco_shapes.circle_geom.instancing_data_B.push_back(test_circle);
                }

                //Put c2,c3,c5 on the circle.
                Cartesian oc2 = cat2-centre;
                oc2.normalize();
                Cartesian p = Cartesian::CrossProduct(oc2,normalp);
                Cartesian r = Cartesian::CrossProduct(normalp,p).by_scalar(sugar_block_scale);
                Cartesian c2p = centre + r;
                res->PutUDData(udd_C2X,c2p.get_x()); res->PutUDData(udd_C2Y,c2p.get_y()); res->PutUDData(udd_C2Z,c2p.get_z());

                Cartesian oc3 = cat3-centre;
                oc3.normalize();
                p = Cartesian::CrossProduct(oc3,normalp);
                r = Cartesian::CrossProduct(normalp,p).by_scalar(sugar_block_scale);
                Cartesian c3p = centre + r;
                res->PutUDData(udd_C3X,c3p.get_x()); res->PutUDData(udd_C3Y,c3p.get_y()); res->PutUDData(udd_C3Z,c3p.get_z());

                Cartesian oc4 = cat4-centre;
                oc4.normalize();
                p = Cartesian::CrossProduct(oc4,normalp);
                r = Cartesian::CrossProduct(normalp,p).by_scalar(sugar_block_scale);
                Cartesian c4p = centre + r;
                res->PutUDData(udd_C4X,c4p.get_x()); res->PutUDData(udd_C4Y,c4p.get_y()); res->PutUDData(udd_C4Z,c4p.get_z());

                Cartesian oc5 = cat5-centre;
                oc5.normalize();
                p = Cartesian::CrossProduct(oc5,normalp);
                r = Cartesian::CrossProduct(normalp,p).by_scalar(sugar_block_scale);
                Cartesian c5p = centre + r;
                res->PutUDData(udd_C5X,c5p.get_x()); res->PutUDData(udd_C5Y,c5p.get_y()); res->PutUDData(udd_C5Z,c5p.get_z());

            } else if(strncmp(res->GetResName(),"FCA",3)==0
                    ||strncmp(res->GetResName(),"FCB",3)==0
                    ||strncmp(res->GetResName(),"FUC",3)==0
                    ||strncmp(res->GetResName(),"FUL",3)==0
                    ){
                /* Draw a triangle */
                //std::cout << "Draw a triangle" << std::endl;
                Cartesian c1c4 = cat1-cat4;
                c1c4.normalize();

                // ADD TO MESH - FUC and friends

                glm::vec4 triangle_col;
                triangle_col = red_col;
                c1c4.normalize();

                Cartesian c1c4p = Cartesian::CrossProduct(c1c4,normal);
                c1c4p.normalize();
                Cartesian normalp = Cartesian::CrossProduct(c1c4p,c1c4);

                glm::mat4 R = get_orientation_matrix(normalp,Cartesian(0.0,0.0,1.0));
                glm::vec4 xprime = R * glm::vec4(1.0,0.0,0.0,1.0);
                glm::mat4 R2 = get_orientation_matrix(c1c4,Cartesian(-xprime[0],-xprime[1],-xprime[2]));
                glm::mat4 R3 = R2*R;

                glm::vec3 triangle_pos(centre.get_x(),centre.get_y(),centre.get_z());
                glm::vec3 triangle_size(sugar_block_scale,sugar_block_scale,sugar_block_thickness);
                coot::instancing_data_type_B_t test_triangle(triangle_pos, triangle_col, triangle_size, R3);
                glyco_shapes.triangle_geom.instancing_data_B.push_back(test_triangle);

            } else if(strncmp(res->GetResName(),"BEM",3)==0
                    ||strncmp(res->GetResName(),"GTR",3)==0
                    ||strncmp(res->GetResName(),"ADA",3)==0
                    ||strncmp(res->GetResName(),"DGU",3)==0
                    ||strncmp(res->GetResName(),"KDN",3)==0
                    ||strncmp(res->GetResName(),"SI3",3)==0
                    ||strncmp(res->GetResName(),"NCC",3)==0

                    ||strncmp(res->GetResName(),"IDR",3)==0
                    ||strncmp(res->GetResName(),"GC4",3)==0
                    ||strncmp(res->GetResName(),"GCD",3)==0
                    ||strncmp(res->GetResName(),"GCU",3)==0
                    ||strncmp(res->GetResName(),"GCV",3)==0
                    ||strncmp(res->GetResName(),"GCW",3)==0
                    ||strncmp(res->GetResName(),"IDS",3)==0
                    ||strncmp(res->GetResName(),"REL",3)==0

                    ||strncmp(res->GetResName(),"SIA",3)==0

                    ){
                bool horizontal = false;
                bool vertical = false;
                bool invert_colour = false;
                if(strncmp(res->GetResName(),"IDR",3)==0
                        ||strncmp(res->GetResName(),"GC4",3)==0
                        ||strncmp(res->GetResName(),"GCD",3)==0
                        ||strncmp(res->GetResName(),"GCU",3)==0
                        ||strncmp(res->GetResName(),"GCV",3)==0
                        ||strncmp(res->GetResName(),"GCW",3)==0
                        ||strncmp(res->GetResName(),"IDS",3)==0
                        ||strncmp(res->GetResName(),"REL",3)==0) {
                    horizontal = true;
                }
                if(strncmp(res->GetResName(),"GTR",3)==0
                        ||strncmp(res->GetResName(),"ADA",3)==0
                        ||strncmp(res->GetResName(),"DGU",3)==0
                        ||strncmp(res->GetResName(),"BEM",3)==0) {
                    vertical = true;
                }
                if(strncmp(res->GetResName(),"IDR",3)==0
                        ||strncmp(res->GetResName(),"BEM",3)==0) {
                    invert_colour = true;
                }

                two_colour = two_colour && (horizontal||vertical);

                /* Draw a diamond */
                //std::cout << "Draw a diamond" << std::endl;
                Cartesian c1c4 = cat1-cat4;
                c1c4.normalize();
                // ADD TO MESH TriangleElement *tri;

                // ADD TO MESH - BEM and friends, horizontal/vertical ?

            } else if(strncmp(res->GetResName(),"XLS",3)==0
                    ||strncmp(res->GetResName(),"CXY",3)==0
                    ||strncmp(res->GetResName(),"RBY",3)==0
                    ||strncmp(res->GetResName(),"TDX",3)==0
                    ||strncmp(res->GetResName(),"XYL",3)==0 /* Should all be orange */
                    ||strncmp(res->GetResName(),"XYS",3)==0
                    ||strncmp(res->GetResName(),"XYP",3)==0
                    ){
                /* Draw a star! */
                //std::cout << "Draw a star!" << std::endl;
                Cartesian c1c4 = cat1-cat4;
                c1c4.normalize();

                // ADD TO MESH - XLS and friends

                glm::vec4 pentagram_col = orange_col;
                c1c4.normalize();

                Cartesian c1c4p = Cartesian::CrossProduct(c1c4,normal);
                c1c4p.normalize();
                Cartesian normalp = Cartesian::CrossProduct(c1c4p,c1c4);

                glm::mat4 R = get_orientation_matrix(normalp,Cartesian(0.0,0.0,1.0));
                glm::vec4 xprime = R * glm::vec4(1.0,0.0,0.0,1.0);
                glm::mat4 R2 = get_orientation_matrix(c1c4,Cartesian(-xprime[0],-xprime[1],-xprime[2]));
                glm::mat4 R3 = R2*R;

                glm::vec3 pentagram_pos(centre.get_x(),centre.get_y(),centre.get_z());
                glm::vec3 pentagram_size(sugar_block_scale,sugar_block_scale,sugar_block_thickness);

                coot::instancing_data_type_B_t block(pentagram_pos, pentagram_col, pentagram_size, R3);
                glyco_shapes.pentagram_geom.instancing_data_B.push_back(block);

            } else {
                //std::cout << "I don't know what to draw for " << res->GetResName() << std::endl;
                // ADD TO MESH - unknown
            }
        }
    }
    return retval;
}

static std::vector<std::vector<std::vector<std::string> > > sugarAtomNames;

void populateSugarAtomNames(){
    // FIXME - We should almost certainly check that residue names match as well. This requires user being able
    // to provide custom saccharide residues....
    // We also need to be checking so that we know what shape to draw. Maybe Jon's code would be helpful here.

    //Complicated ?
    // RIB     C1',...,C4',O4' (5-membered)
    // F6P     C2,...,C5,O5 (5-membered)
    // GLO     O1,C1 ... C5(O5),C6,O6 (linear! Ignored?)
    // SUC     C1,...,C5,O5;C1',...,C4',O4' (6,5-membered)
    // DMU     C1,...,C5,O5;    C10,C5,C7,C8,C9,01 (6,6-membered)

    //Simpler
    // NAG     C1,...,C5,O5; (6-membered)
    // FFC     C1A,...,C5A,O5A; C1B,...,C5B,O5B; (6,6-membered)
    // LAT,LBT C1,...,C5,O5;    C1',...,C5',O5' (6,6-membered)
    // LMT     C1',...,C5',O5'; C1B,...,C5B,O5B; (6,6-membered)
    // TRE     C1P,...,C5P,O5P; C1B,...,C5B,O5B; (6,6-membered)

    std::vector<std::vector<std::string> > ribAtomNames;
    std::vector<std::vector<std::string> > nagAtomNames;
    std::vector<std::vector<std::string> > ffcAtomNames;
    std::vector<std::vector<std::string> > latAtomNames;
    std::vector<std::vector<std::string> > lmtAtomNames;
    std::vector<std::vector<std::string> > treAtomNames;
    std::vector<std::vector<std::string> > dmuAtomNames;
    std::vector<std::vector<std::string> > siaAtomNames;
    std::vector<std::vector<std::string> > a2gAtomNames;
    std::vector<std::string> ribR1AtomNames;
    std::vector<std::string> nagR1AtomNames;
    std::vector<std::string> ffcR1AtomNames;
    std::vector<std::string> latR1AtomNames;
    std::vector<std::string> lmtR1AtomNames;
    std::vector<std::string> treR1AtomNames;
    std::vector<std::string> dmuR1AtomNames;
    std::vector<std::string> ffcR2AtomNames;
    std::vector<std::string> latR2AtomNames;
    std::vector<std::string> lmtR2AtomNames;
    std::vector<std::string> treR2AtomNames;
    std::vector<std::string> dmuR2AtomNames;
    std::vector<std::string> siaR1AtomNames;
    std::vector<std::string> a2gR1AtomNames;

    ffcR1AtomNames.push_back(std::string("C1A"));
    ffcR1AtomNames.push_back(std::string("C2A"));
    ffcR1AtomNames.push_back(std::string("C3A"));
    ffcR1AtomNames.push_back(std::string("C4A"));
    ffcR1AtomNames.push_back(std::string("C5A"));
    ffcR1AtomNames.push_back(std::string("O5A"));
    ffcR2AtomNames.push_back(std::string("C1B"));
    ffcR2AtomNames.push_back(std::string("C2B"));
    ffcR2AtomNames.push_back(std::string("C3B"));
    ffcR2AtomNames.push_back(std::string("C4B"));
    ffcR2AtomNames.push_back(std::string("C5B"));
    ffcR2AtomNames.push_back(std::string("O5B"));
    ffcAtomNames.push_back(ffcR1AtomNames);
    ffcAtomNames.push_back(ffcR2AtomNames);
    sugarAtomNames.push_back(ffcAtomNames);

    latR1AtomNames.push_back(std::string("C1"));
    latR1AtomNames.push_back(std::string("C2"));
    latR1AtomNames.push_back(std::string("C3"));
    latR1AtomNames.push_back(std::string("C4"));
    latR1AtomNames.push_back(std::string("C5"));
    latR1AtomNames.push_back(std::string("O5"));
    latR2AtomNames.push_back(std::string("C1'"));
    latR2AtomNames.push_back(std::string("C2'"));
    latR2AtomNames.push_back(std::string("C3'"));
    latR2AtomNames.push_back(std::string("C4'"));
    latR2AtomNames.push_back(std::string("C5'"));
    latR2AtomNames.push_back(std::string("O5'"));
    latAtomNames.push_back(latR1AtomNames);
    latAtomNames.push_back(latR2AtomNames);
    sugarAtomNames.push_back(latAtomNames);

    lmtR1AtomNames.push_back(std::string("C1'"));
    lmtR1AtomNames.push_back(std::string("C2'"));
    lmtR1AtomNames.push_back(std::string("C3'"));
    lmtR1AtomNames.push_back(std::string("C4'"));
    lmtR1AtomNames.push_back(std::string("C5'"));
    lmtR1AtomNames.push_back(std::string("O5'"));
    lmtR2AtomNames.push_back(std::string("C1B"));
    lmtR2AtomNames.push_back(std::string("C2B"));
    lmtR2AtomNames.push_back(std::string("C3B"));
    lmtR2AtomNames.push_back(std::string("C4B"));
    lmtR2AtomNames.push_back(std::string("C5B"));
    lmtR2AtomNames.push_back(std::string("O5B"));
    lmtAtomNames.push_back(lmtR1AtomNames);
    lmtAtomNames.push_back(lmtR2AtomNames);
    sugarAtomNames.push_back(lmtAtomNames);

    treR1AtomNames.push_back(std::string("C1P"));
    treR1AtomNames.push_back(std::string("C2P"));
    treR1AtomNames.push_back(std::string("C3P"));
    treR1AtomNames.push_back(std::string("C4P"));
    treR1AtomNames.push_back(std::string("C5P"));
    treR1AtomNames.push_back(std::string("O5P"));
    treR2AtomNames.push_back(std::string("C1B"));
    treR2AtomNames.push_back(std::string("C2B"));
    treR2AtomNames.push_back(std::string("C3B"));
    treR2AtomNames.push_back(std::string("C4B"));
    treR2AtomNames.push_back(std::string("C5B"));
    treR2AtomNames.push_back(std::string("O5B"));
    treAtomNames.push_back(treR1AtomNames);
    treAtomNames.push_back(treR2AtomNames);
    sugarAtomNames.push_back(treAtomNames);

    ribR1AtomNames.push_back(std::string("C1'"));
    ribR1AtomNames.push_back(std::string("C2'"));
    ribR1AtomNames.push_back(std::string("C3'"));
    ribR1AtomNames.push_back(std::string("C4'"));
    ribR1AtomNames.push_back(std::string("O4'"));
    ribAtomNames.push_back(ribR1AtomNames);
    sugarAtomNames.push_back(ribAtomNames);

    nagR1AtomNames.push_back(std::string("C1"));
    nagR1AtomNames.push_back(std::string("C2"));
    nagR1AtomNames.push_back(std::string("C3"));
    nagR1AtomNames.push_back(std::string("C4"));
    nagR1AtomNames.push_back(std::string("C5"));
    nagR1AtomNames.push_back(std::string("O5"));
    nagAtomNames.push_back(nagR1AtomNames);
    sugarAtomNames.push_back(nagAtomNames);

    siaR1AtomNames.push_back(std::string("C2"));
    siaR1AtomNames.push_back(std::string("C3"));
    siaR1AtomNames.push_back(std::string("C4"));
    siaR1AtomNames.push_back(std::string("C5"));
    siaR1AtomNames.push_back(std::string("C6"));
    siaR1AtomNames.push_back(std::string("O6"));
    siaAtomNames.push_back(siaR1AtomNames);
    sugarAtomNames.push_back(siaAtomNames);

    a2gR1AtomNames.push_back(std::string("C1"));
    a2gR1AtomNames.push_back(std::string("C2"));
    a2gR1AtomNames.push_back(std::string("C3"));
    a2gR1AtomNames.push_back(std::string("C4"));
    a2gR1AtomNames.push_back(std::string("C5"));
    a2gR1AtomNames.push_back(std::string("O"));
    a2gAtomNames.push_back(a2gR1AtomNames);
    sugarAtomNames.push_back(a2gAtomNames);

}

std::vector<Cartesian> DrawSugarBlock(mmdb::Residue* res1, int selHnd, mmdb::Manager *molHnd, bool two_colour, float sugar_block_thickness, float sugar_block_scale, glyco_shapes_t &glyco_shapes){
    //std::cout << "DrawSugarBlock " << res1->GetResName() << "\n";

    if(sugarAtomNames.size()==0){
        populateSugarAtomNames();
    }

    std::set<std::string> inRingMatch;
    std::vector<Cartesian> first_ring; // Coot: centre and normal of the first ring drawn

    for(unsigned isugartype=0;isugartype<sugarAtomNames.size();isugartype++){
        for(unsigned iring=0;iring<sugarAtomNames[isugartype].size();iring++){
            if(sugarAtomNames[isugartype][iring].size()>4){

                int nring = 0;

                auto pos = inRingMatch.find(sugarAtomNames[isugartype][iring][0]);
                if(pos!=inRingMatch.end()) nring++;
                pos = inRingMatch.find(sugarAtomNames[isugartype][iring][1]);
                if(pos!=inRingMatch.end()) nring++;
                pos = inRingMatch.find(sugarAtomNames[isugartype][iring][2]);
                if(pos!=inRingMatch.end()) nring++;
                pos = inRingMatch.find(sugarAtomNames[isugartype][iring][3]);
                if(pos!=inRingMatch.end()) nring++;
                pos = inRingMatch.find(sugarAtomNames[isugartype][iring][4]);
                if(pos!=inRingMatch.end()) nring++;

                if(nring>0) continue; // A sugar atom shouldn't be in more than 1 ring?

                const char* c1name = sugarAtomNames[isugartype][iring][0].c_str();
                const char* c2name = sugarAtomNames[isugartype][iring][1].c_str();
                const char* c3name = sugarAtomNames[isugartype][iring][2].c_str();
                const char* c4name = sugarAtomNames[isugartype][iring][3].c_str();
                const char* c5name = sugarAtomNames[isugartype][iring][4].c_str();
                const char* o5name = NULL;
                mmdb::Atom* c1 = res1->GetAtom(c1name);
                mmdb::Atom* c2 = res1->GetAtom(c2name);
                mmdb::Atom* c3 = res1->GetAtom(c3name);
                mmdb::Atom* c4 = res1->GetAtom(c4name);
                mmdb::Atom* c5 = res1->GetAtom(c5name);
                mmdb::Atom* o5 = NULL;
                if(sugarAtomNames[isugartype][iring].size()==6){
                    o5name = sugarAtomNames[isugartype][iring][5].c_str();
                    o5 = res1->GetAtom(o5name);
                    //std::cout << "Looking for " << c1name << " " << c2name << " " << c3name << " " << c4name << " " << c5name << " " << o5name << "\n";
                } else {
                    //std::cout << "Looking for " << c1name << " " << c2name << " " << c3name << " " << c4name << " " << c5name << "\n";
                }
                if(c1&&c2&&c3&&c4&&c5&&(o5||sugarAtomNames[isugartype][iring].size()!=6)){
                    //std::cout << "Have them all (1)\n";
                    inRingMatch.insert(sugarAtomNames[isugartype][iring][0]);
                    inRingMatch.insert(sugarAtomNames[isugartype][iring][1]);
                    inRingMatch.insert(sugarAtomNames[isugartype][iring][2]);
                    inRingMatch.insert(sugarAtomNames[isugartype][iring][3]);
                    inRingMatch.insert(sugarAtomNames[isugartype][iring][4]);
                    std::vector<Cartesian> res = DrawSugarBlockInt(res1, selHnd, c1,c2,c3,c4,c5,o5,res1, molHnd, two_colour, sugar_block_thickness, sugar_block_scale, glyco_shapes );
                    //std::cout << "retval size: " << res.size() << std::endl;
                    if(iring==sugarAtomNames[isugartype].size()-1&&res.size()==2) return res;
                    if(first_ring.empty()&&res.size()==2) first_ring = res; // Coot
                }

                c1 = res1->GetAtom(c1name,NULL,"A");
                c2 = res1->GetAtom(c2name,NULL,"A");
                c3 = res1->GetAtom(c3name,NULL,"A");
                c4 = res1->GetAtom(c4name,NULL,"A");
                c5 = res1->GetAtom(c5name,NULL,"A");
                if(sugarAtomNames[isugartype][iring].size()==6){
                    o5 = res1->GetAtom(o5name,NULL,"A");
                }
                if(c1&&c2&&c3&&c4&&c5&&(o5||sugarAtomNames[isugartype][iring].size()!=6)){
                    if(c1->isInSelection(selHnd)&&c2->isInSelection(selHnd)&&c3->isInSelection(selHnd)&&c4->isInSelection(selHnd)&&c5->isInSelection(selHnd)&&o5->isInSelection(selHnd)){
                        //std::cout << "Have them all (2)\n";
                        inRingMatch.insert(sugarAtomNames[isugartype][iring][0]);
                        inRingMatch.insert(sugarAtomNames[isugartype][iring][1]);
                        inRingMatch.insert(sugarAtomNames[isugartype][iring][2]);
                        inRingMatch.insert(sugarAtomNames[isugartype][iring][3]);
                        inRingMatch.insert(sugarAtomNames[isugartype][iring][4]);
                        std::vector<Cartesian> res = DrawSugarBlockInt(res1, selHnd, c1,c2,c3,c4,c5,o5,res1, molHnd, two_colour, sugar_block_thickness, sugar_block_scale, glyco_shapes );
                        if(iring==sugarAtomNames[isugartype].size()-1&&res.size()==2) return res;
                        if(first_ring.empty()&&res.size()==2) first_ring = res; // Coot
                    }
                }
                c1 = res1->GetAtom(c1name,NULL,"B");
                c2 = res1->GetAtom(c2name,NULL,"B");
                c3 = res1->GetAtom(c3name,NULL,"B");
                c4 = res1->GetAtom(c4name,NULL,"B");
                c5 = res1->GetAtom(c5name,NULL,"B");
                if(sugarAtomNames[isugartype][iring].size()==6){
                    o5 = res1->GetAtom(o5name,NULL,"B");
                }
                if(c1&&c2&&c3&&c4&&c5&&(o5||sugarAtomNames[isugartype][iring].size()!=6)){
                    if(c1->isInSelection(selHnd)&&c2->isInSelection(selHnd)&&c3->isInSelection(selHnd)&&c4->isInSelection(selHnd)&&c5->isInSelection(selHnd)&&o5->isInSelection(selHnd)){
                        //std::cout << "Have them all (3)\n";
                        inRingMatch.insert(sugarAtomNames[isugartype][iring][0]);
                        inRingMatch.insert(sugarAtomNames[isugartype][iring][1]);
                        inRingMatch.insert(sugarAtomNames[isugartype][iring][2]);
                        inRingMatch.insert(sugarAtomNames[isugartype][iring][3]);
                        inRingMatch.insert(sugarAtomNames[isugartype][iring][4]);
                        std::vector<Cartesian> res = DrawSugarBlockInt(res1, selHnd, c1,c2,c3,c4,c5,o5,res1, molHnd, two_colour, sugar_block_thickness, sugar_block_scale, glyco_shapes );
                        if(iring==sugarAtomNames[isugartype].size()-1&&res.size()==2) return res;
                        if(first_ring.empty()&&res.size()==2) first_ring = res; // Coot
                    }
                }
                c1 = res1->GetAtom(c1name,NULL,"C");
                c2 = res1->GetAtom(c2name,NULL,"C");
                c3 = res1->GetAtom(c3name,NULL,"C");
                c4 = res1->GetAtom(c4name,NULL,"C");
                c5 = res1->GetAtom(c5name,NULL,"C");
                if(sugarAtomNames[isugartype][iring].size()==6){
                    o5 = res1->GetAtom(o5name,NULL,"C");
                }
                if(c1&&c2&&c3&&c4&&c5&&(o5||sugarAtomNames[isugartype][iring].size()!=6)){
                    if(c1->isInSelection(selHnd)&&c2->isInSelection(selHnd)&&c3->isInSelection(selHnd)&&c4->isInSelection(selHnd)&&c5->isInSelection(selHnd)&&o5->isInSelection(selHnd)){
                        //std::cout << "Have them all (4)\n";
                        inRingMatch.insert(sugarAtomNames[isugartype][iring][0]);
                        inRingMatch.insert(sugarAtomNames[isugartype][iring][1]);
                        inRingMatch.insert(sugarAtomNames[isugartype][iring][2]);
                        inRingMatch.insert(sugarAtomNames[isugartype][iring][3]);
                        inRingMatch.insert(sugarAtomNames[isugartype][iring][4]);
                        std::vector<Cartesian> res = DrawSugarBlockInt(res1, selHnd, c1,c2,c3,c4,c5,o5,res1, molHnd, two_colour, sugar_block_thickness, sugar_block_scale, glyco_shapes );
                        if(iring==sugarAtomNames[isugartype].size()-1&&res.size()==2) return res;
                        if(first_ring.empty()&&res.size()==2) first_ring = res; // Coot
                    }
                }
            }
        }
    }
    return first_ring;
}

//coot::instanced_mesh_t molecules_container_t::DrawSugarBlocks(mmdb::Manager *molHnd, const std::string &cid_str) const {
coot::instanced_mesh_t DrawSugarBlocks(mmdb::Manager *molHnd, const std::string &cid_str, bool three_d, float scale_factor, const glm::vec4 &block_colour) {

    const char * cid = cid_str.c_str();
    int selHnd = molHnd->NewSelection();
    molHnd->Select(selHnd, mmdb::STYPE_RESIDUE, cid, mmdb::SKEY_NEW);

    coot::instanced_mesh_t mesh;

    mmdb::PPResidue selResidues;
    int nSelResidues;
    molHnd->GetSelIndex(selHnd, selResidues, nSelResidues);

    //std::cout << "DrawSugarBlocks, nSelResidues " << nSelResidues << std::endl;

    int two_colour = 0;
    //std::cout << "DrawSugarBlocks, two colour: " << two_colour << "\n";

    const float width = 0.4;

    coot::instanced_geometry_t linkages;

    glyco_shapes_t glyco_shapes;
    glyco_shapes.spheres_for_circles = three_d;

    double stick_col[4];
    bool have_stick_col = false;

    int first_nmr_model=-1;
    int n_nmr_models=molHnd->GetNumberOfModels();
    if(n_nmr_models>0){
        for(int i=1;i<=n_nmr_models;i++){
            if(molHnd->GetModel(i)){
                first_nmr_model = i;
                break;
            }
        }
    } else {
        first_nmr_model=1;
    }
    const char *amino_acid = "GLY,ALA,VAL,PRO,SER,THR,LEU,ILE,CYS,ASP,GLU,ASN,GLN,ARG,LYS,MET,MSE,HIS,PHE,TYR,TRP,HCS,ALO,PDD";
    int proteinSelHnd = molHnd->NewSelection();
    SelectAminoNotHet(molHnd, proteinSelHnd,mmdb::STYPE_ATOM,0,"*",mmdb::ANY_RES,"*",mmdb::ANY_RES,"*",amino_acid,"*","*","*",mmdb::SKEY_NEW);
    int nProteinAtoms;
    mmdb::Atom** ProteinAtoms=0;
    molHnd->GetSelIndex(proteinSelHnd,ProteinAtoms,nProteinAtoms);
    //std::cout << "first_nmr_model " << first_nmr_model << std::endl;
    //std::cout << "Selected " << nProteinAtoms << " atoms" << std::endl;

    bool drawProteinInteractions = false;
    bool drawCovalentProteinInteractions = false;
    bool labelProteinInteractions = false;
    double interaction_cylinder_radius = 0.1;
    int interaction_line_width = 1;
    int interaction_style = DASHCYLINDER;

    double cylinders_accu = 20;
    int quality = 1;
    if(quality==0){
        cylinders_accu = 20;
    } else if(quality==1){
        cylinders_accu = 40;
    } else if(quality==2){
        cylinders_accu = 60;
    }

    double interaction_colour[4] = {0.0,0.0,0.0,1.0};
    float interaction_dash_length = 0.36;
    float sugar_block_thickness = 0.1;
    float sugar_block_scale = 1.0;

    /*
       if(params.GetString("glycoblock_interaction_colour")!=std::string("")&&params.GetString("glycoblock_interaction_colour")!=std::string("default")){
       std::vector<double> stick_col_v = RGBReps::GetColour(params.GetString("glycoblock_interaction_colour"));
       interaction_colour[0] = stick_col_v[0];
       interaction_colour[1] = stick_col_v[1];
       interaction_colour[2] = stick_col_v[2];
       interaction_colour[3] = stick_col_v[3];
       }
     */
    sugar_block_thickness = 0.4;//params.GetFloat("glycoblock_thickness");
    sugar_block_scale = 1.4;//params.GetFloat("glycoblock_scale");
    float stick_thickness = sugar_block_thickness * 0.5;
    float ball_thickness = sugar_block_thickness * 0.9;
    // Coot: user scale for the blocks (the linkage sticks and balls keep their size)
    sugar_block_thickness *= scale_factor;
    sugar_block_scale *= scale_factor;
    interaction_cylinder_radius = 0.1;//params.GetFloat("glycoblock_interaction_cylinder_radius");
    interaction_line_width = 2;//params.GetInt("glycoblock_interaction_line_width");
    std::string interaction_style_pref = "Dashed cylinder";//params.GetString("glycoblock_interaction_style");
    if(interaction_style_pref=="Dashed line"){
        interaction_style = DASHLINE;
    } else {
        interaction_style = DASHCYLINDER;
    }
    interaction_dash_length = 0.36;//params.GetFloat("glycoblock_interaction_dash_length");
    drawProteinInteractions = 1;//params.GetInt("glycoblock_draw_protein_interactions");
    drawCovalentProteinInteractions = 1;//params.GetInt("glycoblock_draw_cov_protein_interactions");
    labelProteinInteractions = 0;//params.GetInt("glycoblock_label_protein_interactions");


    //FIXME - Not cylinders in general!
    // ADD TO MESH SphereCollection *spheres = new SphereCollection();
    // ADD TO MESH CylinderCollection *cyls = new CylinderCollection();
    // ADD TO MESH TriangleCollection *tris = new TriangleCollection();
    // ADD TO MESH DashCylinderCollection *dash_cylinders = new DashCylinderCollection();
    // ADD TO MESH DashLinesCollection *dash_lines_collection = new DashLinesCollection();
    char ResidueID1[30];
    char ResidueID2[30];
    int C1sel = molHnd->NewSelection();
    int ND2sel = molHnd->NewSelection();
    // FFC     C1A,...,C5A,O5A;C1B,...,C5B,O5B; (6,6-membered)
    // RIB     C1',...,C4',O4' (5-membered)
    // F6P     C2,...,C5,O5 (5-membered)
    // SUC     C1,...,C5,O5;C1',...,C4',O4' (6,5-membered)
    // LAT,LBT C1,...,C5,O5;C1',...,C5',O5' (6,6-membered)
    // GLO     O1,C1 ... C5(O5),C6,O6 (linear!)
    // DMU     C1,...,C5,O5;C10,C5,C7,C8,C9,01 (6,6-membered)
    // LMT     C1',...,C5',O5';C1B,...,C5B,O5B; (6,6-membered)
    // TRE     C1P,...,C5P,O5P;C1B,...,C5B,O5B; (6,6-membered)

    molHnd->Select(C1sel,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C1","C","*",SKEY_NEW);
    molHnd->Select(C1sel,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C2","C","*",SKEY_OR);
    molHnd->Select(C1sel,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C1A","C","*",SKEY_OR);
    molHnd->Select(C1sel,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C1B","C","*",SKEY_OR);
    molHnd->Select(C1sel,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C1P","C","*",SKEY_OR);
    molHnd->Select(C1sel,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C1'","C","*",SKEY_OR);

    molHnd->Select(ND2sel,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","ND2","N","*",SKEY_OR);
    // Coot: O-linked glycans too
    molHnd->Select(ND2sel,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","OG","O","*",SKEY_OR);
    molHnd->Select(ND2sel,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","OG1","O","*",SKEY_OR);
    // Coot: and C-mannosylation (C1 to Trp CD1)
    molHnd->Select(ND2sel,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","TRP","CD1","C","*",SKEY_OR);

    mmdb::Atom** nd2Atoms;
    int nNd2Atoms;
    molHnd->GetSelIndex(ND2sel,nd2Atoms,nNd2Atoms);

    mmdb::Atom** c1Atoms;
    int nC1Atoms;
    molHnd->GetSelIndex(C1sel,c1Atoms,nC1Atoms);

    std::cout << "Selected " << nNd2Atoms << " ND2 atoms" << std::endl;
    std::cout << "Selected " << nC1Atoms << " C1 atoms" << std::endl;

    int udd_C1X = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C1X" );
    if (udd_C1X<=0) udd_C1X = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C1X");

    int udd_C1Y = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C1Y" );
    if (udd_C1Y<=0) udd_C1Y = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C1Y");

    int udd_C1Z = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C1Z" );
    if (udd_C1Z<=0) udd_C1Z = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C1Z");

    int udd_C3X = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C3X" );
    if (udd_C3X<=0) udd_C3X = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C3X");

    int udd_C3Y = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C3Y" );
    if (udd_C3Y<=0) udd_C3Y = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C3Y");

    int udd_C3Z = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C3Z" );
    if (udd_C3Z<=0) udd_C3Z = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C3Z");

    int udd_C4X = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C4X" );
    if (udd_C4X<=0) udd_C4X = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C4X");

    int udd_C4Y = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C4Y" );
    if (udd_C4Y<=0) udd_C4Y = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C4Y");

    int udd_C4Z = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C4Z" );
    if (udd_C4Z<=0) udd_C4Z = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C4Z");

    int udd_C5X = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C5X" );
    if (udd_C5X<=0) udd_C5X = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C5X");

    int udd_C5Y = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C5Y" );
    if (udd_C5Y<=0) udd_C5Y = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C5Y");

    int udd_C5Z = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C5Z" );
    if (udd_C5Z<=0) udd_C5Z = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C5Z");

    int udd_C2X = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C2X" );
    if (udd_C2X<=0) udd_C2X = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C2X");

    int udd_C2Y = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C2Y" );
    if (udd_C2Y<=0) udd_C2Y = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C2Y");

    int udd_C2Z = molHnd->GetUDDHandle ( UDR_RESIDUE,"GLYCO_BLOCK_C2Z" );
    if (udd_C2Z<=0) udd_C2Z = molHnd->RegisterUDReal(UDR_RESIDUE,"GLYCO_BLOCK_C2Z");

    /*
    //FIXME - Need geom
    std::vector<coot::h_bond> protein_sugar_hbonds;
    std::vector<coot::h_bond> sugar_sugar_hbonds;
    coot::h_bonds hb;
    protein_sugar_hbonds = hb.get_mcdonald_and_thornton(selHnd, proteinSelHnd, molHnd, geom, 3.8);
    sugar_sugar_hbonds = hb.get_mcdonald_and_thornton(selHnd, selHnd, molHnd, geom, 3.8);
    */

    //std::cout << "Start drawing loop" << std::endl;
    // Coot: ring centres, so that 3D-SNFG linkages can join the blocks centre-to-centre
    std::map<mmdb::Residue *, Cartesian> block_centres;
    for(int ii=0;ii<nSelResidues;ii++){
        mmdb::Residue* res = selResidues[ii];
        if(res&&isSugar(res->name)){
            // Coot: in the 3D-SNFG style the prisms are as thick as they are wide (squares -> cubes)
            float block_thickness = three_d ? sugar_block_scale : sugar_block_thickness;
            std::vector<Cartesian> cn = DrawSugarBlock(res,selHnd, molHnd, two_colour, block_thickness, sugar_block_scale, glyco_shapes);
            if (cn.size() == 2) block_centres[res] = cn[0];
        }
    }
    // Coot: one colour for all the blocks (e.g. to show a reference model), rather than SNFG.
    // Done before the linkages are added, so they stay grey.
    if (block_colour.r >= 0.0f) {
        for (auto geom : {&glyco_shapes.square_geom, &glyco_shapes.circle_geom, &glyco_shapes.triangle_geom,
                          &glyco_shapes.pentagram_geom, &glyco_shapes.bi_square_1_geom,
                          &glyco_shapes.bi_square_2_geom, &glyco_shapes.sphere_geom}) {
            for (auto &inst : geom->instancing_data_A) inst.colour = block_colour;
            for (auto &inst : geom->instancing_data_B) inst.colour = block_colour;
        }
    }
    for(int ii=0;ii<nSelResidues;ii++){
        mmdb::Residue* res = selResidues[ii];
        if(res&&isSugar(res->name)){
            mmdb::realtype UDDX, UDDY, UDDZ;
            res->GetUDData(udd_C1X,UDDX); res->GetUDData(udd_C1Y,UDDY); res->GetUDData(udd_C1Z,UDDZ);
            Cartesian c1(UDDX,UDDY,UDDZ);
            for(int jj=0;jj<ii;jj++){
                mmdb::realtype UDDX2, UDDY2, UDDZ2;
                mmdb::Residue* res2 = selResidues[jj];
                if(res2&&isSugar(res2->name)){
                    double min_length = 1e6;
                    int i_c2 = -1;
                    res2->GetUDData(udd_C2X,UDDX2); res2->GetUDData(udd_C2Y,UDDY2); res2->GetUDData(udd_C2Z,UDDZ2);
                    Cartesian c2_2(UDDX2,UDDY2,UDDZ2);
                    Cartesian b = c2_2-c1;
                    Cartesian final_b;
                    if(b.length()<min_length){
                        min_length = b.length();
                        final_b = b;
                        i_c2 = 2;
                    }
                    res2->GetUDData(udd_C3X,UDDX2); res2->GetUDData(udd_C3Y,UDDY2); res2->GetUDData(udd_C3Z,UDDZ2);
                    Cartesian c2_3(UDDX2,UDDY2,UDDZ2);
                    b = c2_3-c1;
                    if(b.length()<min_length){
                        min_length = b.length();
                        final_b = b;
                        i_c2 = 3;
                    }
                    res2->GetUDData(udd_C4X,UDDX2); res2->GetUDData(udd_C4Y,UDDY2); res2->GetUDData(udd_C4Z,UDDZ2);
                    Cartesian c2_4(UDDX2,UDDY2,UDDZ2);
                    b = c2_4-c1;
                    if(b.length()<min_length){
                        min_length = b.length();
                        final_b = b;
                        i_c2 = 4;
                    }
                    res2->GetUDData(udd_C5X,UDDX2); res2->GetUDData(udd_C5Y,UDDY2); res2->GetUDData(udd_C5Z,UDDZ2);
                    Cartesian c2_5(UDDX2,UDDY2,UDDZ2);
                    b = c2_5-c1;
                    if(b.length()<min_length){
                        min_length = b.length();
                        final_b = b;
                        i_c2 = 5;
                    }
                    if(min_length<4.0){
                        // Coot: the linkage atoms are buried at the faces of solid blocks, so in
                        // the 3D-SNFG style join the ring centres instead, with no end balls
                        if (three_d && block_centres.count(res) && block_centres.count(res2)) {
                            Cartesian cc = block_centres[res2] - block_centres[res];
                            float cc_length = cc.length();
                            Cartesian mid = MidPoint(block_centres[res], block_centres[res2]);
                            cc.normalize();
                            glm::mat4 R = get_orientation_matrix(cc,Cartesian(0.0,0.0,1.0));
                            glm::vec4 bond_col(0.5,0.5,0.5,1);
                            glm::vec3 bond_size(stick_thickness, stick_thickness, cc_length*.5);
                            coot::instancing_data_type_B_t bond_cylinder(glm::vec3(mid.get_x(),mid.get_y(),mid.get_z()), bond_col, bond_size, R);
                            glyco_shapes.circle_geom.instancing_data_B.push_back(bond_cylinder);
                            continue;
                        }
                        glm::vec3 bond_pos(c1.get_x()+0.5*final_b.get_x(),c1.get_y()+0.5*final_b.get_y(),c1.get_z()+0.5*final_b.get_z());
                        glm::vec3 s1_pos(c1.get_x(),c1.get_y(),c1.get_z());
                        glm::vec3 s2_pos(c1.get_x()+final_b.get_x(),c1.get_y()+final_b.get_y(),c1.get_z()+final_b.get_z());
                        final_b.normalize();
                        glm::mat4 R = get_orientation_matrix(final_b,Cartesian(0.0,0.0,1.0));
                        glm::vec4 bond_col(0.5,0.5,0.5,1);
                        glm::vec3 bond_size(stick_thickness, stick_thickness,min_length*.5);
                        coot::instancing_data_type_B_t bond_cyclinder(bond_pos, bond_col, bond_size, R);
                        glyco_shapes.circle_geom.instancing_data_B.push_back(bond_cyclinder);

                        glm::vec3 sphere_size(ball_thickness, ball_thickness,sugar_block_thickness);
                        coot::instancing_data_type_A_t bond_s1(s1_pos, bond_col, sphere_size);
                        glyco_shapes.sphere_geom.instancing_data_A.push_back(bond_s1);
                        coot::instancing_data_type_A_t bond_s2(s2_pos, bond_col, sphere_size);
                        glyco_shapes.sphere_geom.instancing_data_A.push_back(bond_s2);
                    }
                }
            }
        }
    }

    // Coot: link each glycan to its protein residue - a stick from the CA (on the ribbon) to
    // the block centre (3D) or to C1 (flat), for an ND2/OG/OG1/Trp CD1 that is bonded to a C1
    for(int ii=0;ii<nSelResidues;ii++){
        mmdb::Residue* res = selResidues[ii];
        if(!res||!isSugar(res->name)) continue;
        mmdb::Atom* c1_at = res->GetAtom("C1");
        if(!c1_at) continue;
        Cartesian c1_pos(c1_at->x,c1_at->y,c1_at->z);
        for(int ia=0;ia<nNd2Atoms;ia++){
            mmdb::Atom* link_at = nd2Atoms[ia];
            if(!link_at||link_at->isTer()) continue;
            mmdb::Residue* prot_res = link_at->GetResidue();
            if(!prot_res||isSugar(prot_res->name)) continue;
            Cartesian link_pos(link_at->x,link_at->y,link_at->z);
            if((c1_pos-link_pos).length()>2.0) continue;
            mmdb::Atom* ca_at = prot_res->GetAtom("CA");
            if(!ca_at) continue;
            Cartesian start(ca_at->x,ca_at->y,ca_at->z);
            Cartesian end = c1_pos;
            if(three_d && block_centres.count(res)) end = block_centres[res];
            Cartesian b = end - start;
            float b_length = b.length();
            Cartesian mid = MidPoint(start,end);
            b.normalize();
            glm::mat4 R = get_orientation_matrix(b,Cartesian(0.0,0.0,1.0));
            glm::vec4 bond_col(0.5,0.5,0.5,1);
            glm::vec3 bond_size(stick_thickness, stick_thickness, b_length*.5);
            coot::instancing_data_type_B_t bond_cylinder(glm::vec3(mid.get_x(),mid.get_y(),mid.get_z()), bond_col, bond_size, R);
            glyco_shapes.circle_geom.instancing_data_B.push_back(bond_cylinder);
            glm::vec3 sphere_size(ball_thickness, ball_thickness, ball_thickness);
            glyco_shapes.sphere_geom.instancing_data_A.push_back(coot::instancing_data_type_A_t(glm::vec3(start.get_x(),start.get_y(),start.get_z()), bond_col, sphere_size));
            if(!three_d)
                glyco_shapes.sphere_geom.instancing_data_A.push_back(coot::instancing_data_type_A_t(glm::vec3(end.get_x(),end.get_y(),end.get_z()), bond_col, sphere_size));
            break;
        }
    }

    mesh.add(glyco_shapes.square_geom);
    mesh.add(glyco_shapes.circle_geom);
    mesh.add(glyco_shapes.triangle_geom);
    mesh.add(glyco_shapes.pentagram_geom);
    mesh.add(glyco_shapes.bi_square_1_geom);
    mesh.add(glyco_shapes.bi_square_2_geom);
    mesh.add(glyco_shapes.sphere_geom);

    molHnd->DeleteSelection(C1sel);
    molHnd->DeleteSelection(ND2sel);
    molHnd->DeleteSelection(selHnd);
    molHnd->DeleteSelection(proteinSelHnd);

    return mesh;

}

} // namespace glycoblocks
} // namespace coot
