
#include "clipper/clipper.h"
#include "clipper/clipper-ccp4.h"

using namespace clipper;

#include <iostream>
#include <string.h>
#include "CXXSurfaceMaker.h"
#include "CXXCreator.h"
#ifndef  __MMDB_Manager__
#include "mmdb2/mmdb_manager.h"
#endif
#include "CXXSurfaceVertex.h"
#include "CXXSphereTriangleEdge.h"
#include "CXXCoord.h"
#include "CXXUtils.h"
#include <math.h>

int usage (string arg);

int main (int argc, char * const argv[]) {
	mmdb::Manager*  theMMDBManager;
	int RC;
	string insrf, outsrf;
	
	mmdb::InitMatType();
	
	std::string inputCoordinateName = " ";
	std::string inputSurfaceName = " ";
	std::string outputSurfaceName = " ";
	std::string inputMapName = " ";
	std::string outputMapName = " ";
	std::string vertexName = "vertices";
	
	double angle = 30.;
	double probeRadius = 1.5;
	double sampling = 0.7;
	int analytical = 1;
	bool evaluateElectrostatics = true;
	
	if (argc<2){
		usage(argv[0]);
		exit (0);
	}
	
	for (int iArg = 1; iArg < argc; iArg++){
		string arg(argv[iArg]);
		if (arg == "-usage" || arg == "-u" || arg == "-help" || arg == "-h"){
			usage(argv[0]);
			exit (0);
		}
		else if (arg == "-probe"){
			iArg++;
			probeRadius = atof(argv[iArg]);
		}
		else if (arg == "-angle"){
			analytical = 1;
			iArg++;
			angle = atof(argv[iArg]);
		}
		else if (arg == "-vertex"){
			iArg++;
			vertexName = argv[iArg];
		}
		else if (arg == "-mapout"){
			iArg++;
			outputMapName = argv[iArg];
		}
		else if (arg == "-mapin"){
			iArg++;
			inputMapName = argv[iArg];
		}
		else if (arg == "-surfacein"){
			iArg++;
			inputSurfaceName = argv[iArg];
		}
		else if (arg == "-surfaceout"){
			iArg++;
			outputSurfaceName = argv[iArg];
		}
		else if (arg == "-pdbin" || arg == "xyzin"){
			iArg++;
			inputCoordinateName = argv[iArg];
		}
		else if (arg.substr(0,1) != "-"){
			if (inputCoordinateName == " ") inputCoordinateName = arg;
			else outputSurfaceName = arg;
		}
		else if (arg == "-no"){
			evaluateElectrostatics= false;
		}
	}
	
	//Read input coordinates if appropriate
	int selHnd;
	if (inputCoordinateName != " "){
		char *fn = new char[(strlen(inputCoordinateName.c_str())+1)];
		strcpy (fn, inputCoordinateName.c_str());
		theMMDBManager = new mmdb::Manager();
		theMMDBManager->SetFlag( mmdb::MMDBF_PrintCIFWarnings );
		RC = theMMDBManager->ReadCoorFile (fn);
		// if RC > 0 reading in file has failed - say so and quit
		if (RC) {
			std::cout << "error could not read file ...[" << fn << "]\n";
			return 2;
		}
		else {
			cout << "Reading coordinate file " << inputCoordinateName << endl;
		}
		
		//make new selection containing all atoms for now
		selHnd = theMMDBManager->NewSelection();
		theMMDBManager->SelectAtoms( selHnd, 0,"*",mmdb::ANY_RES,"*",mmdb::ANY_RES,"*","*","*","*","*" );	
		
		//Add united Atom Radius
		CXXUtils::assignUnitedAtomRadius(theMMDBManager, selHnd);
	}
	
	//
	//Get hold of a surface: May be readin, or aclculated
	CXXSurfaceMaker *calculatedSurface = new CXXSurfaceMaker();
	if (inputSurfaceName != " "){
		calculatedSurface->readGraspFile(inputCoordinateName);
		cout << "Input surface :\n";
		cout << calculatedSurface->report();
	}
	else {
		if (analytical){
			calculatedSurface->calculateFromAtoms (theMMDBManager, selHnd, selHnd, probeRadius, angle*2.*M_PI / 360., false);
		}
	}	
	
	clipper::Cell cell;
	clipper::NXmap<double> theClipperNXMap;

	if (evaluateElectrostatics){
		//Generate or read a potential map		
		if (inputMapName != " "){
			CCP4MAPfile inputMap;
			inputMap.open_read(inputMapName);
			inputMap.import_nxmap(theClipperNXMap);
		}
		else {
			//Instantiate an electrostatics map and cause it to calculate itself
			CXXChargeTable theChargeTable;
			CXXUtils::assignCharge(theMMDBManager, selHnd, &theChargeTable);
			CXXCreator *theCreator = new CXXCreator(theMMDBManager, selHnd);
			theCreator->calculate();
			theClipperNXMap = theCreator->coerceToClipperMap(cell);
		}
		
		// Now bring the surface and the map together
		double coords[4];
        for (std::vector<CXXSurface>::iterator subSurfaceIter= calculatedSurface->getChildSurfaces().begin();
             subSurfaceIter != calculatedSurface->getChildSurfaces().end();
             subSurfaceIter++){
            int potentialHandle = subSurfaceIter->getScalarHandle("potential");
            
            for (int i=0; i< subSurfaceIter->numberOfVertices(); i++){
                //Use the colour if it has been assigned
                if (subSurfaceIter->getCoord(vertexName, i, coords)){
                    cout << "Bizarely no vertices for coordinate "<< i << endl;
                }
                Coord_orth orthogonals(coords[0], coords[1], coords[2]);
                double potential;
                const Coord_map mapUnits(theClipperNXMap.coord_map(orthogonals));
                potential = theClipperNXMap.interp<Interp_cubic>( mapUnits );
                subSurfaceIter->setScalar(potentialHandle, i, potential);
            }
        }
	}
	
	//Output surface file if desired
	if (outputSurfaceName != " "){
		calculatedSurface->writeAsGrasp(outputSurfaceName);
		cout << calculatedSurface->report();
	}
	
	//Output potential map if so requested
	if (outputMapName != " "){
		CCP4MAPfile outputMapFile;
		outputMapFile.open_write(outputMapName);
		outputMapFile.set_cell(cell);
		outputMapFile.export_nxmap(theClipperNXMap);
		outputMapFile.close_write();
	}
	
	delete calculatedSurface;
	return 0;
}

int usage(std::string argv0){
	cout << "Usage: " << argv0 << "[-probe radius  <1.5>] [-angle triangle_size_in_degrees <30.0>]\n"
	<<"\t [-vertex name <vertices>] [-grid surface_grid_in_angstroms]\n"
	<<"\t [-mapout Output_CCP4_map_file_name] [-mapin Input_CCP4_map_file_name]\n"
	<<"\t [-surfacein Input_Grasp_surface_file_name] [Input_PDB_file_name] Output_Grasp_surface_file_name\n\n";
	cout << "The program assigns a potential labelled \"potential\" to the vertices of a grasp format surface\n\n";
	cout << "The potential comes either from electrostatic potential calculated from the provided PDB (default)\n";
	cout << "\t or from interpolation into a ccp4 format map provided with the -mapin specifier\n\n";
	cout << "The surface is either calculated from the provided PDB (default)\n";
	cout << "\t or from the surface provided with the -surfacein specifier\n\n";
	cout << "Optionally, the potential calculated can be written to a CCP4 format map (-mapout specifier)\n\n";
	cout << "By default surfaces are calculated on a rather fine angular tolerance (-angle specifier: smaller=smoother)\n";
	cout << "\t alternatively a less accurate surface algorithm can be specified with the -grid specifier with an arguement of about 0.5\n";
	return 0;
}
