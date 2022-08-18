
/* definiton of memberfunctions of SolventMap*/

#include "CXXFFTSolventMap.h"
//#include "SolventMapBindings.h"

SolventMap::SolventMap () { //dummy constructor ...
}

SolventMap::SolventMap (fftw_real dGrid, fftw_real rProbe, fftw_real xmin, fftw_real xmax, 
                        fftw_real ymin, fftw_real ymax, fftw_real zmin, fftw_real zmax) {
	
	
	/* Copy information about location of protein and grid/probe parameters */
	
	angstroemMinMax[0]=xmin;
	angstroemMinMax[1]=xmax;
	angstroemMinMax[2]=ymin;
	angstroemMinMax[3]=ymax;
	angstroemMinMax[4]=zmin;
	angstroemMinMax[5]=zmax;
	
	// these shoud denote a box with positive dimensions - check this WARNING
	
	
	probeRadius = rProbe;
	gridSpacing = dGrid;
	countRapid = 1;
	
	/* check for obvious error in parameters*/
	
	if (probeRadius < 0 | gridSpacing < 0) {
		CXXException theException = CXXException("ERROR: SolventMap, negative probeRadius or gridSpacing - check parameter list?\n");
		throw theException;
	}
    
	/* reserve memory on heap and get grid to point to first position, set all gridpoints to zero*/
	
	positionGrid(); 
	optimizeGridDim();
	
	/* reserve memory on heap and get grid to point to first position, set all gridpoints to zero*/
	
	if ((dim[0] <= 0) | (dim[1] <= 0) | (dim[2] <= 0)) {
		
		CXXException theException = CXXException("ERROR in: SolventMap::SolventMap(...) - zero or negative map dimension");	
		throw theException;
		
	}
	
	Grid = new fftw_real[dim[0]*dim[1]*dim[2]];
	SolidGrid = new fftw_real[dim[0]*dim[1]*dim[2]];
	
	distanceGrid = new fftw_real[dim[0]*dim[1]*dim[2]];
	
	if (Grid == 0 || SolidGrid == 0 || distanceGrid == 0) {
		CXXException theException = CXXException(" ERROR: in: SolventMap::SolventMap(...) - Could not reserve suffiecent memory for grid");
		throw theException;
	}
	
	
	
	/* set all values in parentGrid to 1 ( bulk solvent, all points allowed probe positions)*/
	for (int x = 0; x < dim[0]; x++) {
		for (int y  = 0; y < dim[1]; y++) {
			for (int z  = 0; z < dim[2]; z++) {
				setGrid3i (x, y, z, 1.);
				setDistanceGrid3i (x, y, z, 0.0);
				setSolidGrid3i (x, y, z, 0.);
				
				
			}
		}
	}
}



SolventMap::~SolventMap () {
	delete [] Grid;
	delete [] SolidGrid;
	delete [] distanceGrid;
}

void SolventMap::optimizeGridDim () {
	
	int Flag = 0;
	int temp[3] = {1,1,1};
	int i;
	
	for (i = 0; i < 3; i++) {
		Flag = 0;
		
		/* See if dim[i] can be expressed as product of prime numbers smaller than 11 if not add, two (to keep even) and try again*/
		while (Flag == 0) {
			while(dim[i]%11 == 0) {
				dim[i] = dim[i]/11;
				temp[i] = temp[i]*11;
			}
			while(dim[i]%7 == 0) {
				dim[i] = dim[i]/7;
				temp[i] = temp[i]*7;
			}
			while(dim[i]%5 == 0) {
				dim[i] = dim[i]/5;
				temp[i] = temp[i]*5;
			}	  
			while(dim[i]%3 == 0) {
				dim[i] = dim[i]/3;
				temp[i] = temp[i]*3;
			}
			while(dim[i]%2 == 0) {
				dim[i] = dim[i]/2;
				temp[i] = temp[i]*2;
			}
			if (dim[i] != 1) {
				dim[i]=dim[i]*temp[i]+2;
				intOriginOfGrid[i] -= 1;
				originOfGrid[i] = intOriginOfGrid[i] * gridSpacing;
				temp[i]=1;
				Flag = 0;
			}
			else {
				dim[i]=temp[i];
				Flag = 1;
			}	
		}
	}
}


void SolventMap::positionGrid () {
	float xyzMin[3], xyzMax[3];
	float uvwMin[3], uvwMax[3];
	int intLimitOfGrid[3];
	int i;
	
	/* Find position of box given by users dimension  plus two probeRadii padding*/
	
	for (i = 0; i < 3; i++) {
		
		xyzMin[i] = angstroemMinMax[2*i] - 2*probeRadius;
		xyzMax[i] = angstroemMinMax[2*i + 1] + 2*probeRadius;
		
		uvwMin[i] = (xyzMin[i] / gridSpacing) - 1;
		uvwMax[i] = (xyzMax[i] / gridSpacing) + 1;
		
		if (xyzMin[i] < 0)  intOriginOfGrid[i] = (int) (uvwMin[i]-1.0);
		else intOriginOfGrid[i] = (int) uvwMin[i];
		
		if (xyzMax[i] < 0)  intLimitOfGrid[i] = (int) uvwMax[i];
		else intLimitOfGrid[i] = (int) (uvwMax[i]+1.0);
		
		originOfGrid[i] = intOriginOfGrid[i] * gridSpacing;
		dim[i] = (intLimitOfGrid[i] - intOriginOfGrid[i]) + 1;
		dim[i] = ( dim[i]%2 ? dim[i]+1 : dim[i]);  
	}  
}


fftw_real SolventMap::getSpacing() {
	return gridSpacing; 
}

fftw_real SolventMap::getProbeRadius() {
	return probeRadius; 
}

int SolventMap::getFloatOrigin(fftw_real *theFloatOrigin) {
	int i;
	for (i = 0; i < 3; i++) {
		theFloatOrigin[i] = originOfGrid[i]; 
	}
	return 0;
}


int SolventMap::getIntOrigin(int *theIntOrigin) {  
	int i;
	for (i = 0; i < 3; i++) {
		theIntOrigin[i] = intOriginOfGrid[i]; 
	}
	return 0;
}

int SolventMap::getExtent(int *theExtent) {  
	int i;
	for (i = 0; i < 3; i++) {
		theExtent[i] = dim[i]; 
	}
	return 0;
}

void SolventMap::setGrid3i (int i, int j, int k, float value){
    
	if((i*dim[1]*dim[2] + j*dim[2] + k) >= dim[0]*dim[1]*dim[2]){
		CXXException theException = CXXException("ERROR in: SolventMap::setGrid3i(...) - index error");
		throw theException;
	}
	Grid[i*dim[1]*dim[2] + j*dim[2] + k] = value;
}


void SolventMap::setGrid3iv (int *ijk, float value){
	setGrid3i(ijk[0], ijk[1], ijk[2], value);    
}

void SolventMap::setDistanceGrid3iv (int *ijk, float value){
	setDistanceGrid3i(ijk[0], ijk[1], ijk[2], value);
}

void SolventMap::setDistanceGrid3i (int i, int j, int k, float value){
	
	if((i*dim[1]*dim[2] + j*dim[2] + k) >= dim[0]*dim[1]*dim[2]){
		CXXException theException = CXXException("ERROR in: SolventMap::setDistanceGrid3i(...) - index error");
		throw theException;
	}
	distanceGrid[i*dim[1]*dim[2] + j*dim[2] + k] = value;	
}


void SolventMap::setSolidGrid3iv (int *ijk, float value){
	setSolidGrid3i(ijk[0], ijk[1], ijk[2], value);
}

void SolventMap::setSolidGrid3i (int i, int j, int k, float value){
	
	if((i*dim[1]*dim[2] + j*dim[2] + k) >= dim[0]*dim[1]*dim[2]){
		CXXException theException = CXXException("ERROR in: SolventMap::setGridSolid3i(...) - index error");
		throw theException;
	}
	SolidGrid[i*dim[1]*dim[2] + j*dim[2] + k] = value;
}



fftw_real SolventMap::getDistanceAt3i (int i, int j, int k){
	if((i*dim[1]*dim[2] + j*dim[2] + k) >= dim[0]*dim[1]*dim[2]){
		CXXException theException = CXXException("ERROR in: SolventMap::getDistanceAtd3i(...) - index error");
		throw theException;
	}
	return (distanceGrid[i*dim[1]*dim[2] + j*dim[2] + k]);
}


fftw_real SolventMap::getDistanceAt3iv (int *ijk){
	return (getDistanceAt3i(ijk[0], ijk[1], ijk[2]));
}

fftw_real SolventMap::getSolidAt3i (int i, int j, int k){
	if((i*dim[1]*dim[2] + j*dim[2] + k) >= dim[0]*dim[1]*dim[2]){
		CXXException theException = CXXException("ERROR in: SolventMap::setSolidAt3i(...) - index error");
		throw theException;
	}
	return (SolidGrid[i*dim[1]*dim[2] + j*dim[2] + k]);
}

fftw_real SolventMap::getSolidAt3iv (int *ijk){
	return (getSolidAt3i(ijk[0], ijk[1], ijk[2]));
}

void SolventMap::prepareForRapidFFT() {
	countRapid = 0;
}

void SolventMap::makeDistMap (int sampleNrPerGrid) {
	int i,j;
	fftw_real deltaR;
	deltaR = gridSpacing/sampleNrPerGrid;
	fftw_real* testProbeRadius = new fftw_real[2*sampleNrPerGrid+1];
	fftw_real** contactNrGrid = new fftw_real*[2*sampleNrPerGrid+1];
	fftw_real* GridPointer;
	
	prepareForRapidFFT();
	
	/* Null Grid for reference*/
	
	GridPointer = new fftw_real[dim[0]*dim[1]*dim[2]];
	for (j =0; j < dim[0]*dim[1]*dim[2]; j++) {
		GridPointer[j] = 0;
	}
	contactNrGrid[0] = GridPointer;
	
	
	/* make contact Nr grids for radii between probeRadius - 2*gridSpacing to probeRadius + 2*gridSpacing...*/
	
	for (i = 0; i < 2*sampleNrPerGrid+1; i++) {
		testProbeRadius[i+1] = ((probeRadius - gridSpacing) + i*deltaR);
		cout << i+1 << " FFT: Now testing: rProbe = " << testProbeRadius[i+1] << "\n";
		convoluteSolidProbe (testProbeRadius[i+1], 1,0,0 );
		GridPointer = new fftw_real[dim[0]*dim[1]*dim[2]];
		
		for (j =0; j < dim[0]*dim[1]*dim[2]; j++) {
			GridPointer[j] = SolidGrid[j];
		}
		contactNrGrid[i+1] = GridPointer;
	}
	
	
	/* run through successive contactNr grids and note for whichtestProbeRadius*/
	/* gridpoints "light up" for the first time (distance from closest probe..)*/
	fftw_real* theGrid;
	fftw_real* lastGrid;    
	
	for (i = 0; i < dim[0]*dim[1]*dim[2]; i++) {
		distanceGrid[i] = 0;
	}
	
	for (j = 0; j <  2*sampleNrPerGrid+1; j++) {
		for (i = 0; i < (dim[0]*dim[1]*dim[2]); i++) {
			
			lastGrid = contactNrGrid[j];
			theGrid = contactNrGrid[j+1];
			if (theGrid[i] > 0 && lastGrid[i] == 0) {
				distanceGrid[i] = testProbeRadius[j+1];
			}
		}
	}
    
	/* these are the ones that can not be accessed at all   */
	
	for (i = 0; i < dim[0]*dim[1]*dim[2]; i++) {
		if (distanceGrid[i] == 0)
			distanceGrid[i] = probeRadius+3*gridSpacing;
	}
	
	/* clean up and log off */
	for (i = 0; i < 2*sampleNrPerGrid + 1; i++) {
		theGrid = contactNrGrid[i];
		delete [] theGrid;
	}
    
	delete [] testProbeRadius;
	delete [] contactNrGrid;
	cleanUpForRapidFFT();
    
}  



int SolventMap::convoluteSolidProbe (fftw_real choosenProbeRadius, int rapidFlag, int smoothFlag, float smoothRadius) {
	
	/* Dynamic allocation of storage for probe array and complex arrays for fourier FFTs*/
	
	fftw_real *Probe;
	fftw_complex *FFTProbe;
	fftw_complex *FourierScratch;
	int i, j, k, x, y, z;
	
	Probe = new fftw_real[dim[0]*dim[1]*dim[2]];
	FFTProbe = new fftw_complex[dim[0]*dim[1]*(dim[2]/2+1)];
	FourierScratch = new fftw_complex[dim[0]*dim[1]*(dim[2]/2+1)];  
	
	if ((rapidFlag == 0) | (rapidFlag == 1 && countRapid == 0)) {
		FFTGrid =  new fftw_complex[dim[0]*dim[1]*(dim[2]/2+1)];
		cout << "Making new FFTGrid \n";
	}
	
	/* exeption stuff here....*/
	
	if (FFTProbe == 0 | FFTGrid == 0 | FourierScratch == 0 | Probe == 0) {
		CXXException theException(" ERROR (SolventMap::convoluteSphere() ): could not reserve suffiecent memory !");
		throw theException;
	}
	
	
	/* open or create file: ESPfile. If old ESPfile present import wisdom from file */
	
        /*
	FILE *ESPfile = NULL;
	if (rapidFlag == 0 | (rapidFlag == 1 && countRapid == 0)) {
		ESPfile = fopen("ESPfile", "r+");
		if (ESPfile == NULL) {
			ESPfile=fopen("ESPfile", "w+");
			cout << "No old wisdom -  new ESPfile created\n";
		}
		else {
			cout << "Found old wisdom\n";
			fftw_import_wisdom_from_file(ESPfile);
		}
	}
        */
	
	/* make plans for FFT in both directions*/
    
	if ((rapidFlag == 0) | (rapidFlag == 1 && countRapid == 0)) {
		complexToRealPlan = rfftw3d_create_plan 
		(dim[0], dim[1], dim[2], FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_USE_WISDOM);
		realToComplexPlan = rfftw3d_create_plan 
			(dim[0], dim[1], dim[2], FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_USE_WISDOM);  
	}
	cout << "FFTW plans created " << realToComplexPlan << "\n"; cout.flush();
	
	/* export new wisdom to file and close file*/
	
	//int status; /* not used yet... ERROR handling yet to come....*/
	
        /*
	if (rapidFlag == 0 | (rapidFlag == 1 && countRapid == 0)) {
		fftw_export_wisdom_to_file(ESPfile);
		status=fclose(ESPfile);
	}
	cout << "FFTW plans exported\n"; cout.flush();
        */
	
	/* Make a real space image of the surfacing probe in row major..*/
	
	fftw_real  probeRadiusSq, DistSq, Scale;
	probeRadiusSq = choosenProbeRadius*choosenProbeRadius;
	fftw_real GrSpSq;
	GrSpSq = gridSpacing*gridSpacing;
	
	for (x = 0; x < dim[0]; x++) {
		for (y = 0; y < dim[1]; y++) {
			for (z = 0; z < dim[2]; z++) {
				Probe[x*dim[1]*dim[2]+y*dim[2]+z] = 0;
			}
		}
	}
	
	for (x = 0; x < dim[0]; x++) {
		DistSq =  x*x*GrSpSq;
		if ( DistSq < probeRadiusSq) {
			for (y = 0; y < dim[1]; y++) {
				DistSq =  x*x*GrSpSq + y*y*GrSpSq;
				if (DistSq < probeRadiusSq) {
					for (z = 0; z < dim[2]; z++) {
						DistSq =  x*x*GrSpSq + y*y*GrSpSq +  z*z*GrSpSq;
						if (DistSq < probeRadiusSq) {
							float value;
							int x1, y1, z1;
							value = 1.;
							x1 = (dim[0]-x)%dim[0];
							y1 = (dim[1]-y)%dim[1];
							z1 = (dim[2]-z)%dim[2];
							Probe[x*dim[1]*dim[2]+  y*dim[2]+ z] = value;
							Probe[x1*dim[1]*dim[2]+ y*dim[2]+ z] = value;
							Probe[x*dim[1]*dim[2]+ y1*dim[2]+ z] = value;
							Probe[x1*dim[1]*dim[2]+y1*dim[2]+ z] = value;
							Probe[x*dim[1]*dim[2]+  y*dim[2]+z1] = value;
							Probe[x1*dim[1]*dim[2]+ y*dim[2]+z1] = value;
							Probe[x*dim[1]*dim[2]+ y1*dim[2]+z1] = value;
							Probe[x1*dim[1]*dim[2]+y1*dim[2]+z1] = value;
						}
					}
				}
			}      
		}
	}
	cout << "FFTW Probe map generated\n"; cout.flush();
	
	if ((rapidFlag == 0) | (rapidFlag == 1 && countRapid == 0)) {
		int i;
		for (i = 0; i < dim[0]*dim[1]*(dim[2]/2+1); i++) { 
			FFTGrid[i].re = 0;        
			FFTGrid[i].im = 0;        
		}
	}
	
	for (i = 0; i < dim[0]*dim[1]*(dim[2]/2+1); i++) {
		FFTProbe[i].re = 0;
		FFTProbe[i].im = 0;   
		FourierScratch[i].re = 0;    
		FourierScratch[i].im = 0;    
	}
	cout << "FFTW intermediate and target arrays emptied\n";cout.flush();
	
	/* do the fourier transforms by passing start address of row major arrays to FFTW	  */
	rfftwnd_one_real_to_complex (realToComplexPlan, Probe, FFTProbe);
	cout << "FFTW Probe fourier transformed\n";cout.flush();
	
	
	/* FFT Protein Grid only if old one is not recycled*/
	if ((rapidFlag == 0) | (rapidFlag == 1 && countRapid == 0)) {
		rfftwnd_one_real_to_complex (realToComplexPlan, Grid, FFTGrid);
		countRapid = 1;
	}
	cout << "FFTW Protein fourier transformed\n";cout.flush();
	
	fftw_real * RealScratch;
	/* Convolute  FFT of Probe with FFT of Grid*/
	Scale = (dim[0]*dim[1]*dim[2]);       
	for (i = 0; i < dim[0]*dim[1]*(dim[2]/2+1); i++ ) {
		FourierScratch[i].re = (FFTGrid[i].re * FFTProbe[i].re - FFTGrid [i].im * FFTProbe[i].im)/Scale;
		FourierScratch[i].im = (FFTGrid[i].re * FFTProbe[i].im + FFTGrid [i].im * FFTProbe[i].re)/Scale;
	}
	cout << "FFTW Convolution effected\n";cout.flush();
	/* free some memory and clean up*/
	delete [] FFTProbe;
	delete [] Probe;
	
	/*smooth stuff*/
	cout << "Probe space freed\n";cout.flush();
	
	if (smoothFlag) {
		fftw_real *smoothProbe;
		fftw_complex *FFTSmoothProbe;
		
		smoothProbe = new fftw_real[dim[0]*dim[1]*dim[2]];
		FFTSmoothProbe = new fftw_complex[dim[0]*dim[1]*(dim[2]/2+1)];
		
		for (i = 0; i < dim[0]*dim[1]*dim[2]; i++) {
			SolidGrid[i] = 0.;
		}
		/* transform map back into real space... */
		rfftwnd_one_complex_to_real (complexToRealPlan, FourierScratch, SolidGrid);
		/*
		 *	Set all non zero points to a value of 0, otherwise 1
		 */
		for (i = 0; i <dim[0] ; i++) {
			for (j = 0; j <dim[1] ; j++) {
				for (k = 0; k <dim[2] ; k++) {
					smoothProbe[i*dim[1]*dim[2]+j*dim[2]+k] = 0;
					if (SolidGrid[i*dim[1]*dim[2]+j*dim[2]+k] < 0.001){ 
						SolidGrid[i*dim[1]*dim[2]+j*dim[2]+k] = 1;
					}
					else { 
						SolidGrid[i*dim[1]*dim[2]+j*dim[2]+k] = 0;
					}
				}
			}
		}  
		
		/*  Pop back into reciprocal space*/
		for (i = 0; i < dim[0]*dim[1]*(dim[2]/2+1); i++) {
			FourierScratch[i].re = FFTSmoothProbe[i].re = 0;    
			FourierScratch[i].im = FFTSmoothProbe[i].im = 0;    
		}
		
		rfftwnd_one_real_to_complex (realToComplexPlan, SolidGrid, FourierScratch);
		
		//Generate a real space smoothing probe, which has a value of 1 at the origin, 
		// and a value of max(0, 1-distance_from_origin/smoothing_radius) elsewhere
		float smoothRadiusSq;
		smoothRadiusSq = smoothRadius*smoothRadius;
		cout << "Smoothing by " << smoothRadius << " " << smoothRadiusSq << "\n";
		smoothRadiusSq = smoothRadius*smoothRadius;
		for (x = 0; x < dim[0]; x++) {
			DistSq =  x*x*GrSpSq;
			if ( DistSq < smoothRadiusSq) {
				for (y = 0; y < dim[1]; y++) {
					DistSq =  x*x*GrSpSq + y*y*GrSpSq;
					if (DistSq < smoothRadiusSq) {
						for (z = 0; z < dim[2]; z++) {
							DistSq =  x*x*GrSpSq + y*y*GrSpSq +  z*z*GrSpSq;
							if (DistSq < smoothRadiusSq) {
								float value;
								int x1, y1, z1;
								value = (1.0-(pow(DistSq,0.5f)/smoothRadius));
								x1 = (dim[0]-x)%dim[0];
								y1 = (dim[1]-y)%dim[1];
								z1 = (dim[2]-z)%dim[2];
								smoothProbe[x*dim[1]*dim[2]+  y*dim[2]+ z] = value;
								smoothProbe[x1*dim[1]*dim[2]+ y*dim[2]+ z] = value;
								smoothProbe[x*dim[1]*dim[2]+ y1*dim[2]+ z] = value;
								smoothProbe[x1*dim[1]*dim[2]+y1*dim[2]+ z] = value;
								smoothProbe[x*dim[1]*dim[2]+  y*dim[2]+z1] = value;
								smoothProbe[x1*dim[1]*dim[2]+ y*dim[2]+z1] = value;
								smoothProbe[x*dim[1]*dim[2]+ y1*dim[2]+z1] = value;
								smoothProbe[x1*dim[1]*dim[2]+y1*dim[2]+z1] = value;
							}
						}
					}
				}
			}      
		}
		
		//Prepare to Fourier transform the smoothing probe, and transform it into reciprocal space
		for ( i = 0; i < dim[0]*dim[1]*(dim[2]/2+1); i++) {
			FFTSmoothProbe[i].re = 0;
			FFTSmoothProbe[i].im = 0;
		}
		rfftwnd_one_real_to_complex (realToComplexPlan, smoothProbe, FFTSmoothProbe);
		
		/* convolute smoothing probe with the modified mask*/
		for (i = 0; i < dim[0]*dim[1]*(dim[2]/2+1); i++ ) {
			FourierScratch[i].re = 
			(FourierScratch[i].re * FFTSmoothProbe[i].re - FourierScratch[i].im * FFTSmoothProbe[i].im)
			/Scale;
			FourierScratch[i].im = 
				(FourierScratch[i].re * FFTSmoothProbe[i].im + FourierScratch[i].im * FFTSmoothProbe[i].re)
				/Scale; 
		}  
	}
	
	/* free some more memory*/
	
	if (rapidFlag == 0){
		delete [] FFTGrid;
		rfftwnd_destroy_plan(realToComplexPlan);
	}
	
	/* transform map back into real space... */
	RealScratch = new fftw_real[dim[0]*dim[1]*dim[2]];
	for ( i = 0; i < dim[0]*dim[1]*dim[2]; i++) {
		RealScratch[i] = 0.;
	}
	rfftwnd_one_complex_to_real (complexToRealPlan, FourierScratch, RealScratch);
	
	
	/* finish cleaning up  */
	
	if (rapidFlag == 0) {
		rfftwnd_destroy_plan(complexToRealPlan);
	}
	
	delete [] FourierScratch;  
	
	/* transformations could happen here, truncate */
	
	for ( i = 0; i <dim[0] ; i++) {
		for (j = 0; j <dim[1] ; j++) {
			for (k = 0; k <dim[2] ; k++) {
				float value;
				value = getDistanceAt3i(i, j, k);
				if (value<0.1) setDistanceGrid3i(i, j, k, 100.);
				value = getDistanceAt3i(i, j, k);
				setDistanceGrid3i(i,j,k,100.-value);
			}
		}
	}
	
	for ( i = 0; i <dim[0] ; i++) {
		for (j = 0; j <dim[1] ; j++) {
			for (k = 0; k <dim[2] ; k++) {
				float value;
				if (smoothFlag) {
					if (/* DISABLES CODE */ (0)&&(value = getDistanceAt3i(i, j, k)) > 0.){
					}
					else {
						value = RealScratch[(i*dim[1]*dim[2]+j*dim[2]+k)];
						value = value * (probeRadius/2.87);
					}
					setSolidGrid3i (i, j, k, value);
				}
				else {
					value = floor(RealScratch[(i*dim[1]*dim[2]+j*dim[2]+k)] + 0.5);
					setSolidGrid3i (i, j, k, value);
				}
			}
		}  
	}
	delete [] RealScratch;
	
	return 0;
}




void SolventMap::makeDualContactMap () {
	int i;
	convoluteSolidProbe(probeRadius, 0,0,0);
	for (i = 0; i < dim[0]*dim[1]*dim[2]; i++) {    
		if ( SolidGrid[i] > 0 )
			SolidGrid[i] = 1;
	} 
}

void SolventMap::cleanUpForRapidFFT() {
	countRapid = 0;
	rfftwnd_destroy_plan(realToComplexPlan);
	rfftwnd_destroy_plan(complexToRealPlan);
	delete [] FFTGrid;
	cout << "FFTGrid obliterated, Plan gone -  SNAFU\n";
}

int SolventMap::atomAt3f (float x, float y, float z, float AtomRadius) {
	return(atomAt3d ((fftw_real)x, (fftw_real)y, (fftw_real)z, AtomRadius));
}

int SolventMap::atomAt3fv (float xyz[3], float AtomRadius) {
	return (atomAt3d ((fftw_real)xyz[0], (fftw_real)xyz[1], (fftw_real)xyz[2], AtomRadius));
}

int SolventMap::atomAt3dv (fftw_real xyz[3], fftw_real AtomRadius) {
	return (atomAt3d (xyz[0], xyz[1], xyz[2], AtomRadius));
}

int SolventMap::atomAt3d (fftw_real x, fftw_real y, fftw_real z, fftw_real AtomRadius) {
	float xyz[3], uvw[3];
	int i, j, k, iUvwMin[3], iUvwMax[3];
	float uvwRadius, uvwRadiusSq, distSq, xyzRadius;
    
	xyz[0] = x; xyz[1] = y; xyz[2] = z;
	xyzRadius = AtomRadius + probeRadius;
	uvwRadius = xyzRadius / gridSpacing;
	uvwRadiusSq = uvwRadius * uvwRadius;
	for (i=0; i<3; i++){
		uvw[i] = (xyz[i] - originOfGrid[i]) / gridSpacing;
		iUvwMin[i] = (int) (uvw[i] - uvwRadius);
		iUvwMax[i] = (int) (uvw[i] + uvwRadius);
		iUvwMin[i] -= 1;
		iUvwMax[i] += 1;
		if (iUvwMin[i] < 0) iUvwMin[i] = 0;
		if (iUvwMin[i] > dim[i]-1) iUvwMin[i] = dim[i]-1;
		if (iUvwMax[i] < 0) iUvwMax[i] = 0;
		if (iUvwMax[i] > dim[i]-1) iUvwMax[i] = dim[i]-1;
	}
    
	for (i=iUvwMin[0]; i<iUvwMax[0]; i++){
		for (j=iUvwMin[1]; j<iUvwMax[1]; j++){
			for (k=iUvwMin[2]; k<iUvwMax[2]; k++){
				distSq = ((float)i-uvw[0])*((float)i-uvw[0]) +
				((float)j-uvw[1])*((float)j-uvw[1]) +
				((float)k-uvw[2])*((float)k-uvw[2]);
				if (distSq < uvwRadiusSq) {
					setGrid3i (i, j, k, 0.);
					float dist = sqrt(distSq);
					dist *= gridSpacing;
					float value = 100. - (xyzRadius-dist);
					float oldValue = getDistanceAt3i(i, j, k);
					if (oldValue > 0.1){
						if (oldValue < 100. - probeRadius){
							if (value < oldValue){
								setDistanceGrid3i(i, j, k, value);
							}
						}
						else if (value < 100.-probeRadius) {
							setDistanceGrid3i(i, j, k, value);
						}  
						else {
							setDistanceGrid3i(i,  j, k, 200.);
						}
					}
					else {
						setDistanceGrid3i(i, j, k, value);
					}
				}
			}
		}
	}
	
	return 0;
}



int SolventMap::dumpXSlice (int type, int SliceNr, int broad) {
	int i;
	fstream dump ("dumpfile",ios::out|ios::app );
	
	if (SliceNr > dim[0] | SliceNr < 1) {
		dump << "ERROR: Slice Nr outside of range\n"; 
		return 1;
	}
	
	
	dump << "\nAtomGrid: Slicing through 3D grid at x= " << SliceNr <<". Directions: (horizontal/vetical)  <-> (z/y): \n\n";  
	for (i = (dim[2]*dim[1]*(SliceNr-1)) ; i < (dim[1]*dim[2]*SliceNr); i++) {
		if (i%dim[2]==0 && i!= (dim[2]*dim[1]*(SliceNr-1)))  
		{
			dump << "\n\n";
		}
		
		dump.width(broad);
		dump.precision(2);
		/* only for test case truncate output:    if (i%dim[2]>6 && i%dim[2]<33)*/
		if (type == 2) dump << distanceGrid[i];
		if (type == 1) dump << Grid[i];
	}
	dump << "\n\n";
	return 0;
}

/*===============================================================================*/
/* C- bindings*/

#define CSymbol "C"

extern CSymbol void *cSolventMap (fftw_real dGrid, fftw_real rProbe, fftw_real xmin, 
								  fftw_real xmax, fftw_real ymin, fftw_real ymax, fftw_real zmin, fftw_real zmax){
	return ((void *) new SolventMap( dGrid, rProbe, xmin, xmax, ymin, ymax, zmin, zmax));
}

extern CSymbol fftw_real cGetProbeRadius (void *theMap){
	return (((SolventMap *)theMap)->getProbeRadius());
}

extern CSymbol fftw_real cGetSpacing (void *theMap){
	return (((SolventMap *)theMap)->getSpacing());
}

extern CSymbol int cGetFloatOrigin (void *theMap, fftw_real *theFloatOrigin){
	return (((SolventMap *)theMap)->getFloatOrigin(theFloatOrigin));
}

extern CSymbol int cGetIntOrigin (void *theMap, int *theIntOrigin){
	return (((SolventMap *)theMap)->getIntOrigin(theIntOrigin));
}

extern CSymbol int cGetExtent (void *theMap, int *theExtent){
	return (((SolventMap *)theMap)->getExtent(theExtent));
}

extern CSymbol fftw_real cGetSolidAt3i (void *theMap, int i, int j, int k){
	return (((SolventMap *)theMap)->getSolidAt3i(i, j, k));
}

extern CSymbol fftw_real cGetSolidAt3iv (void *theMap, int *ijk){
	return (((SolventMap *)theMap)->getSolidAt3iv(ijk));
}

extern CSymbol fftw_real cGetDistanceAt3i (void *theMap, int i, int j, int k){
	return (((SolventMap *)theMap)->getDistanceAt3i(i, j, k));
}

extern CSymbol fftw_real cGetDistanceAt3iv (void *theMap, int *ijk){
	return (((SolventMap *)theMap)->getDistanceAt3iv(ijk));
}

extern CSymbol int const AtomAt3fv (void *theMap, float xyz[3], float AtomRadius){
	return (((SolventMap *)theMap)->atomAt3fv(xyz, AtomRadius));
}
extern CSymbol int const AtomAt3f (void *theMap, float x, float y, float z, float AtomRadius){
	return (((SolventMap *)theMap)->atomAt3f(x, y, z, AtomRadius));
}
extern CSymbol int const AtomAt3dv (void *theMap, fftw_real xyz[3], fftw_real AtomRadius){
	return (((SolventMap *)theMap)->atomAt3dv(xyz, AtomRadius));
}
extern CSymbol int const AtomAt3d (void *theMap, fftw_real x, fftw_real y, fftw_real z, fftw_real AtomRadius){
	return (((SolventMap *)theMap)->atomAt3d(x, y, z, AtomRadius));
}
extern CSymbol int cConvoluteSolidProbe (void *theMap, fftw_real ProbeRadius, int flag, int smoothFlag, float smoothRadius){
	return (((SolventMap *)theMap)->convoluteSolidProbe (ProbeRadius, flag, smoothFlag, smoothRadius));
}

extern CSymbol void cMakeDistMap (void *theMap, int sampleNr){
	((SolventMap *)theMap)->makeDistMap (sampleNr);
	return;
}

















