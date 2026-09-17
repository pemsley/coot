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

#include "coords/Cartesian.hh"
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
				      int isample_step, int iream_start, int n_reams,
                                      bool is_em_map,
                                      bool use_vertex_gradients_for_map_normals_flag);

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

	// De-duplication of vertices during marching cubes. Each grid edge that is
	// cut by the isosurface produces exactly one vertex, shared by the triangles
	// of all cells that touch that edge. The edge id space is dense and bounded
	// (see GetEdgeID()/GetVertexID()), so instead of a std::map keyed by edge id
	// we directly address a flat array: m_edge_to_vertex_index[edge_id] holds the
	// compacted vertex index (-1 if this edge has no vertex yet), and
	// m_edge_vertices holds the unique intersection points in first-encounter order.
	std::vector<int> m_edge_to_vertex_index;
	std::vector<POINT3DID> m_edge_vertices;

	// De-duplicating store of the intersection point on the given edge (insert if
	// absent - the first point stored for an edge id wins, as with the old map).
	void store_edge_vertex(unsigned int edge_id, const POINT3DID &pt);

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


#endif // CISOSURFACE_H
