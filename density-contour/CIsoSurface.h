/* density-contour/CIsoSurface.h
 *
 * Copyright 2000 Paul Bourke
 * Copyright 2000 Cory Gene Bloyd
 * Copyright 2005 The University of York
 *
 * Author: Raghavendra Chandrashekara, Paul Bourke and Cory Gene Bloyd
 *         Paul Emsley and Kevin Cowtan
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
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

#ifndef CISOSURFACE_H
#define CISOSURFACE_H
// File Name: CIsoSurface.h
// Last Modified: 5/8/2000
// Author: Raghavendra Chandrashekara (basesd on source code
// provided by Paul Bourke and Cory Gene Bloyd)
// Email: rc99@doc.ic.ac.uk, rchandrashekara@hotmail.com
//
// Description: This is the interface file for the CIsoSurface class.
// CIsoSurface can be used to construct an isosurface from a scalar
// field.

#include <map>
#include <vector>
#include "Vectors.h"

#include <string>

#include "coords/Cartesian.h"
// Clipper stuff
#include "clipper/core/xmap.h"
#include "clipper/core/nxmap.h"

#include "density-contour-triangles.hh"

typedef std::vector<TRIANGLE> TRIANGLEVECTOR;


template <class T> class CIsoSurface {
public:
	// Constructor and destructor.
	CIsoSurface();
	~CIsoSurface();

	// Generates the isosurface from the scalar field contained in the
	// buffer ptScalarField[].
	void GenerateSurface(const T* ptScalarField, T tIsoLevel, unsigned int nCellsX, unsigned int nCellsY,  unsigned int nCellsZ, float fCellLengthX, float fCellLengthY, float fCellLengthZ);

	// Called with an Xmap.
	//
	// I suppose we could overload the function name.... Nah, let's not.
	//
	//vector<CartesianPair> GenerateSurface_from_Xmap(const clipper::Xmap<T>& crystal_map,
	//			       const  T tIsoLevel);

	// We overload the function name this time.
	//
	// vector<CartesianPair>
	coot::CartesianPairInfo
	  GenerateSurface_from_Xmap(const clipper::Xmap<T>& crystal_map,
				    T tIsoLevel,
				    float box_radius, // half length
				    coot::Cartesian centre_point,
				    int isample_step,
				    int iream_start, int n_reams,
				    bool is_em_map);

	coot::CartesianPairInfo
	  GenerateSurface_from_NXmap(const clipper::NXmap<T>& nx_map,
				    T tIsoLevel,
				    float box_radius, // half length
				    coot::Cartesian centre_point,
				    int isample_step); // is EM map

	coot::density_contour_triangles_container_t
	  GenerateTriangles_from_Xmap(const clipper::Xmap<T>& crystal_map,
				      T tIsoLevel,
				      float box_radius, // half length
				      coot::Cartesian centre_point,
				      int isample_step, int iream_start, int n_reams, bool is_em_map);

	std::pair<int, int> rangeify(const clipper::Grid_map &grid, int isample_step, int isection_start,
				     int n_sections) const;

	// Returns true if a valid surface has been generated.
	bool IsSurfaceValid();

	// Deletes the isosurface.
	void DeleteSurface();

	// Returns the length, width, and height of the volume in which the
	// isosurface in enclosed in.  Returns -1 if the surface is not
	// valid.
	int GetVolumeLengths(float& fVolLengthX, float& fVolLengthY, float& fVolLengthZ);

	// PE adds
	unsigned int nTriangles(void);

	// PE adds
	void morphVertices(void);

	// PE adds
	void writeTriangles(std::string);

	// PE adds
	coot::CartesianPairInfo
	   returnTriangles(const clipper::Xmap<T>& xmap,
			   const clipper::Coord_frac& base,
			   float radius,
			   coot::Cartesian centre,
			   bool is_em_map) const;

	// PE adds
	coot::CartesianPairInfo
	  returnTriangles(const clipper::NXmap<T>& nx_map,
			  const clipper::Coord_frac& base,
			  float radius,
			  coot::Cartesian centre) const; // certainly is EM map


	// PE adds
	void check_max_min_vertex_index_from_triangles(void);

	// PE adds
	void check_max_min_vertices(void);

protected:
	// The number of vertices which make up the isosurface.
	unsigned int m_nVertices;

	// The vertices which make up the isosurface.
	POINT3D* m_ppt3dVertices;

	// The number of triangles which make up the isosurface.
	unsigned int m_nTriangles;

	// The indices of the vertices which make up the triangles.
	unsigned int* m_piTriangleIndices;

	// The number of normals.
	unsigned int m_nNormals;

	// The normals.
	VECTOR3D* m_pvec3dNormals;

	// List of POINT3Ds which form the isosurface.
	ID2POINT3DID m_i2pt3idVertices;

	// List of TRIANGLES which form the triangulation of the isosurface.
	TRIANGLEVECTOR m_trivecTriangles;

	// Returns the edge ID.
	unsigned int GetEdgeID(unsigned int nX, unsigned int nY, unsigned int nZ, unsigned int nEdgeNo);

	// Returns the vertex ID.
	unsigned int GetVertexID(unsigned int nX, unsigned int nY, unsigned int nZ);

	// Calculates the intersection point of the isosurface with an
	// edge.
	POINT3DID CalculateIntersection(unsigned int nX, unsigned int nY, unsigned int nZ, unsigned int nEdgeNo);

	// Interpolates between two grid points to produce the point at which
	// the isosurface intersects an edge.
	POINT3DID Interpolate(float fX1, float fY1, float fZ1, float fX2, float fY2, float fZ2, T tVal1, T tVal2);

	// Renames vertices and triangles so that they can be accessed more
	// efficiently.
	void RenameVerticesAndTriangles();

	// used by above
	static void
	  rename_tris_in_thread(const std::pair<unsigned int, unsigned int> &idx_range,
				TRIANGLEVECTOR &tv, const ID2POINT3DID &point_map);


	// Calculates the normals.
	void CalculateNormals();

	// No. of cells in x, y, and z directions.
	unsigned int m_nCellsX, m_nCellsY, m_nCellsZ;

	// Cell length in x, y, and z directions.
	float m_fCellLengthX, m_fCellLengthY, m_fCellLengthZ;

	// The buffer holding the scalar field.
	const T* m_ptScalarField;

	// The isosurface value.
	T m_tIsoLevel;

	// Indicates whether a valid surface is present.
	bool m_bValidSurface;

	// Lookup tables used in the construction of the isosurface.
	static const unsigned int m_edgeTable[256];
	// PE changes to int (some values are negative)
	static const int m_triTable[256][16];

	// PE adds
	bool isSmallTriangle(unsigned int i);

	// PE adds
	void adjustVertices(unsigned int i);


};


// This is a list of vertices (basically, indices)
//
class to_vertex_list_t {

   //vector<bool> vertex_list;
   int *vertex_list;
   int vertex_list_size;  // the size of the array
   int n_vertices;        // the maximum index filled so far.

 public:
   to_vertex_list_t();
   to_vertex_list_t(const to_vertex_list_t &a);
   void Copy(const to_vertex_list_t &a);
   ~to_vertex_list_t();

   const to_vertex_list_t& operator=(const to_vertex_list_t &a);

   void add(int i);
   bool contains(int i);
   // bool operator[](unsigned int) const;
};

// This is a list of vertices to which there may be connections to
// other vertices.
//
// It is a container class.
//
// I loathe this type of programming. I loathe it, I loathe it, I
// loathe it, I loathe it, I loathe it.  I've spent two days on this now
// and it still doesn't work.  Grrrr.  Waaagh.... and it would be so
// simple in scheme...
//
class done_line_list_t {

   to_vertex_list_t *from_vertices;
   void resize_and_copy(int j);

 public:

   done_line_list_t();

   ~done_line_list_t();

   int from_vertices_size;  // the size of the array
   int max_from_vertex;     // the maximum vertex encountered so far.

   // to_vertex_list_t operator[](unsigned int) const;
   to_vertex_list_t getVertex(unsigned int i) const;

   void mark_as_done(int i, int j);
   bool done_before(int i, int j) ;  // question and manipulation of class
};


#endif // CISOSURFACE_H
