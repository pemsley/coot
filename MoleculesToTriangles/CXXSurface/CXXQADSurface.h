#include <iostream>
#include <vector>
#include "clipper/clipper.h"
#include "mmdb2/mmdb_manager.h"
#include "mmdb2/mmdb_tables.h"
#include "mmdb2/mmdb_uddata.h"
#include <math.h>
#include "CXXCircle.h"
#include "CXXCoord.h"
#include "CXXNewHood.h"

typedef struct {
	double x, y, z, t;
}Point3D;
class CXXQADSurface {
	
	typedef struct {
		char *longname;         /* long name of object */
		char *shortname;        /* short name of object */
		char *dual;             /* long name of dual */
		int numverts;           /* number of vertices */
		int numedges;           /* number of edges */
		int numfaces;           /* number of faces */
		Point3D v[162];    /* the vertices */
		int f[320*4];     /* the faces */
	} Polyinfo;
	
	typedef struct {
		clipper::Coord_orth p[8];
		double val[8];
	} GRIDCELL;
	
private: 
		double sample;
	double probeRadius;
	double atomRadius;
	mmdb::Atom**	selectedAtoms;
	int nSelectedAtoms;
	clipper::Xmap<double> theDoubleMap;	
	clipper::Xmap<int> theFlagMap;	
	clipper::Spacegroup clipperSpacegroup; 
	clipper::Cell clipperCell;
	clipper::Grid_sampling clipperGridSampling;
	clipper::Grid_range clipperGridRange;
	int prepareGrids();
	int makeDistanceSqMap();
	int allowProbesToEat();
	int allowProbeToEatWithinGridRange(clipper::Coord_orth probeCoordOrth, clipper::Grid_range theRange);
	void dump(clipper::Grid_sampling theObject);
	void dump(clipper::Coord_grid theObject);
	void dump(clipper::Grid_range theObject);
	int selHndl;
	class mmdb::Manager *theMMDBManager;
	int coordIsBuriedByNeighbours(clipper::Coord_orth &point, int iAtom);
	void addProbe(clipper::Coord_orth);
	void copyFlagToDouble();
	
	static int edgeTable[];
	static int triTable[][16];
	static int triTablePrime[][16];
	static int nTriangles[];
	
	std::vector<clipper::Coord_orth> vertices;
	std::vector<int> triangles;
	std::vector<clipper::Coord_orth>vertexNormals;
	std::vector<clipper::Coord_orth>probePositions;
	std::vector<std::vector<int> > neighbourhoods;
	std::vector<double> atomRadii;
	int nVdwProbePositions;
	
	clipper::Coord_orth coordOrthInterp(double isolevel, clipper::Coord_orth p1, clipper::Coord_orth p2, 
										double v1, double v2);
	int contourPixel(clipper::Xmap_base::Map_reference_coord index, double isoLevel, 
					 clipper::Coord_orth vertlist[], 
					 std::vector<clipper::Coord_orth> &vertVector,
					 int triangles[][3]);
	int contourMap(double isoLevel);
	int compareVertices(clipper::Coord_orth &v1, clipper::Coord_orth &v2);
	int calculateAveragedNormals();
	double maxAtomRadius;
	double fastGetAtomRadius(int iAtom) {return atomRadii[iAtom];}
	double getAtomRadius(mmdb::Atom*);
	int transformTriTable();	
	int setInaccessibleDistanceSq();
	int sqrtDistanceSq();
	int addProbesFromVdwSurface();
	
	int toruses();
	clipper::Grid_range &gdIntersection(clipper::Grid_range &g0, clipper::Grid_range &g1, clipper::Grid_range &g2);
public:
		enum Flag {
			Solvent = 0, vdW = 1, Inaccessible = 2, Accessible = 3
		};
	
	CXXQADSurface(mmdb::Manager* theMMDBManager, int selHndl, 
				  double probeRadius_in, double sample_in);
	~CXXQADSurface();	
	clipper::Cell &getCell();
	clipper::Xmap<double> &getDoubleMap();
	std::vector<clipper::Coord_orth> &getVertices();
	std::vector<clipper::Coord_orth> &getNormals();
	std::vector<int> &getTriangles();
};
