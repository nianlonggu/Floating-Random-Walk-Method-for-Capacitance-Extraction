#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <fstream>
#include "FRWcapext.h"
#include <unistd.h>
using namespace std;



bool FRWControl::updateProgress(int tempWalkNum, int masterCondIndex,  double temp_VGSArea ){
	currentNumOfWalkOfCond[masterCondIndex]+=tempWalkNum;

	for(int index=0; index< FRWControl::capMatrix[masterCondIndex].size(); index++ ){
		double average_cap=FRWControl::capMatrix[masterCondIndex][index].sumOfCap/currentNumOfWalkOfCond[masterCondIndex];
		double variance_cap=FRWControl::capMatrix[masterCondIndex][index].sumOfCapSquared/currentNumOfWalkOfCond[masterCondIndex] - pow(average_cap,2.0);
		double std_var_cap=pow( variance_cap/currentNumOfWalkOfCond[masterCondIndex], 0.5  );
		double cap_err=100.0*3*std_var_cap/average_cap;
		std_var_cap*=temp_VGSArea;
		FRWControl::capMatrix[masterCondIndex][index].estimateErr=cap_err;
		FRWControl::capMatrix[masterCondIndex][index].std_var=std_var_cap;
	}


	int temp=currentNumOfWalkOfCond[masterCondIndex]*100 / totalNumOfWalkOfCond[masterCondIndex] ;
	if( temp> currentProgressOfCond[masterCondIndex]  ){
		currentProgressOfCond[masterCondIndex] =temp;
		cout<<"\rConductor No "<<masterCondIndex<<"/"<<FRWControl::conductorList.size()-1<<":\tFRW Path number"<<temp<<"%"<<flush;  //<<"%\tstd_err: "<<std_var_cap<<"\t percent: "<<cap_err<<"%"<<endl;
	}

	if( currentProgressOfCond[masterCondIndex]>=100 ){
		return true;
	}else{

		for(int index=0; index< FRWControl::capMatrix[masterCondIndex].size(); index++ ){
			if( FRWControl::capMatrix[masterCondIndex][index].estimateErr>2 ){
				return false;
			}
		}
		return true;
	}

}


Cube::Cube(){
	x1=x2=y1=y2=z1=z2=0;
}

Cube::Cube(double x_1, double y_1 , double z_1, double size ){
	x1=x_1;
	y1=y_1;
	z1=z_1;
	x2=x_1+size;
	y2=y_1+size;
	z2=z_1+size;
}

Cube::Cube( FPoint &minPoint, double size ){
	x1=minPoint.x1;
	y1=minPoint.y1;
	z1=minPoint.z1;
	x2=x1+size;
	y2=y1+size;
	z2=z1+size;
}


Cube::Cube( double x_1, double x_2, double y_1, double y_2, double z_1, double z_2 ){
	x1=x_1;
	x2=x_2;
	y1=y_1;
	y2=y_2;
	z1=z_1;
	z2=z_2;
}

double Cube::size(){
	return x2-x1;
}


GaussianSurface::GaussianSurface(){
	x1=x2=y1=y2=z1=z2=0.0;
	setSurfaceArea();

}

GaussianSurface::GaussianSurface( FPoint &minPoint, FPoint &maxPoint ,std::vector<double> extensionDistance   ){
	x1=minPoint.x1;
	y1=minPoint.y1;
	z1=minPoint.z1;
	x2=maxPoint.x1;
	y2=maxPoint.y1;
	z2=maxPoint.z1;
	extensionDis=extensionDistance;
	setSurfaceArea();
}

GaussianSurface::GaussianSurface(double x_1, double x_2, double y_1, double y_2, double z_1, double z_2 ,std::vector<double> extensionDistance ){
	x1=x_1;
	x2=x_2;
	y1=y_1;
	y2=y_2;
	z1=z_1;
	z2=z_2;
	extensionDis=extensionDistance;
	setSurfaceArea();
}


//enum NormalDirection{  PX, PY, PZ, NX, NY, NZ };
void GaussianSurface::setSurfaceArea(){
	surfaceArea=std::vector<double> (6);
	surfaceArea[PX]=surfaceArea[NX]= ( y2-y1 ) *(z2-z1);
	surfaceArea[PY]=surfaceArea[NY]= (x2-x1)*(z2-z1);
	surfaceArea[PZ]=surfaceArea[NZ]= (x2-x1)*(y2-y1);
}

double GaussianSurface::area(){
	double len=x2-x1;
	double width=y2-y1;
	double height=z2-z1;

	return 2*( len*width + len*height + width*height  );
}

bool GaussianSurface::isPointEnclosedBySurface( FPoint &point ){
	return  point.x1>x1 && point.x1<x2 && point.y1>y1 && point.y1 <y2 && point.z1 >z1 && point.z1 <z2;
}

bool GaussianSurface::isPointOnSurface(FPoint &point){
	if(  point.x1 >=x1 && point.x1 <=x2 && point.y1>=y1 &&point.y1 <= y2  ){
		if( fabs( point.z1-z1)<1E-6  || fabs( point.z1-z2 )<1E-6   ){
			return true;
		}
	}
	if( point.x1 >= x1 && point.x1<=x2 && point.z1 >=z1 && point.z1 <=z2   ){
		if(  fabs(point.y1-y1 )< 1E-6 || fabs( point.y1-y2 )<1E-6  ){
			return true;
		}
	}
	if(  point.y1 >=y1 && point.y1 <=y2 && point.z1 >=z1 && point.z1 <=z2   ){
		if(  fabs(point.x1-x1 )< 1E-6 || fabs( point.x1-x2 )< 1E-6  ){
			return true;
		}
	}
	return false;

}

//@author:Dragon
// for the 6 surfaces of GaussianSurface, if point is on that surface, the retuened values will be:
/*
PX   1
NX  -1
PY    2
NY  -2           
PZ    3
NZ   -3

 if  point is not on the surface, the return value will be 0
*/
int GaussianSurface::getNormalDirOfPoint(FPoint  &point){
	if(  point.y1 >=y1 && point.y1 <=y2 && point.z1 >=z1 && point.z1 <=z2   ){
		if(  fabs(point.x1-x1 )< 1E-7 ){
			return -1;
		}else if(fabs( point.x1-x2 )<1E-7 ){
			return 1;
		}
	}
	if( point.x1 >= x1 && point.x1<=x2 && point.z1 >=z1 && point.z1 <=z2   ){
		if(  fabs(point.y1-y1 )< 1E-7 ){
			return -2;
		}else if(fabs( point.y1-y2 )<1E-7  ){
			return 2;
		}
	}
	if(  point.x1 >=x1 && point.x1 <=x2 && point.y1>=y1 &&point.y1 <= y2  ){
		if( fabs( point.z1-z1)<1E-7 ){
			return -3;
		}else if(fabs( point.z1-z2 )<1E-7   ){
			return 3;
		}
	}
	//point is not on the surface
	return 0;
}

bool GaussianSurface::operator ==( GaussianSurface &rhs  ){
	if(  x1==rhs.x1 &&y1==rhs.y1 &&z1==rhs.z1 &&x2==rhs.x2 &&y2==rhs.y2 &&z2==rhs.z2  ){
		return true;
	}
	return false;
}

bool GaussianSurface::isNeighborWith( GaussianSurface &gaussian ){
	return x1<=gaussian.x2 && x2>=gaussian.x1 && y1<=gaussian.y2 && y2>=gaussian.y1 && z1<=gaussian.z2 && z2 >= gaussian.z1;
}


double myDistance( Cube &cube,  Rectangle &rect ){
	double xDis=max (   cube.x1-rect.x2, rect.x1-cube.x2   );
	double yDis=max(   cube.y1-rect.y2, rect.y1-cube.y2   );
	double zDis=max(  cube.z1-rect.z2, rect.z1-cube.z2   );
	double temp= max( xDis, yDis  );
	return max(temp, zDis);

}

void *FRW( void *condIn ){


//设置线程亲和度
	// cpu_set_t mask;
	// CPU_ZERO(&mask);
	// pthread_mutex_lock(&processorIndex_mutex);
	// CPU_SET(FRWControl::currentProcessorIndex, &mask);
	// if(pthread_setaffinity_np(pthread_self(), sizeof(mask), &mask) < 0) {
 //        cout<<"set thread affinity failed"<<endl;
 //    }
 //    FRWControl::currentProcessorIndex= (FRWControl::currentProcessorIndex+1)%FRWControl::totalProcessorNum;
	// pthread_mutex_unlock(&processorIndex_mutex);


	long randSeed_uniform_01=0;
	long randSeed_int_normalDir=0;
	long randSeed_int_GVT=0;  //green function value table

	pthread_mutex_lock(&getRandSeed_mutex);

	int threadIndex=FRWControl::currentThreadIndex++;
	if(FRWControl::randSeedVec.size()==0 ){
		cout<<"Err: no randSeed in randSeedVec!"<<endl;
		exit(1);
	}
	randSeed_uniform_01=FRWControl::randSeedVec.back();
	FRWControl::randSeedVec.pop_back();

	if(FRWControl::randSeedVec.size()==0 ){
		cout<<"Err: no randSeed in randSeedVec!"<<endl;
		exit(1);
	}
	randSeed_int_normalDir = FRWControl::randSeedVec.back();
	FRWControl::randSeedVec.pop_back();	

	if(FRWControl::randSeedVec.size()==0 ){
		cout<<"Err: no randSeed in randSeedVec!"<<endl;
		exit(1);
	}
	randSeed_int_GVT =FRWControl::randSeedVec.back();
	FRWControl::randSeedVec.pop_back();
	pthread_mutex_unlock(&getRandSeed_mutex);

	//@author:Dragon
	//initialize boost random generator, use randSeed!
	boost::mt19937 gen_uni_01(randSeed_uniform_01);
	boost::mt19937 gen_int_normalDir(randSeed_int_normalDir );
	boost::mt19937 gen_int_GVT(randSeed_int_GVT );
	//here we use a group of rand seed generator!
	boost::uniform_01<> dist_uni_01;
	boost::uniform_int<> dist_int_normalDir(0,5);
	boost::uniform_int<> dist_int_GVT(0, FRWControl::GreenVT.size()-1);

	boost::variate_generator<boost::mt19937&,boost::uniform_01<> > rand_01(gen_uni_01, dist_uni_01 );
	boost::variate_generator<boost::mt19937&,boost::uniform_int<> >rand_normalDir(gen_int_normalDir, dist_int_normalDir );
	boost::variate_generator<boost::mt19937&,boost::uniform_int<> > rand_GVT(gen_int_GVT, dist_int_GVT );


	FPoint pointOnVGS, nextPoint, point2d;
	Cube cube;
	Cell *rootCell=NULL, *leafCell=NULL;

	bool isFinish=false;
	int loopTimer=0;
	int Nt=0, Ng=0;
	//生成本地电容矩阵
	CapacitanceMatrix capMatrix=generateCapMatrix( FRWControl::conductorList );

	int condIndex=0;
	ConductorList::iterator condIt;
	for(condIt=FRWControl::conductorList.begin(); condIt!=FRWControl::conductorList.end();condIt++, condIndex++ ){
		while( true  )  { 
					
        	FloatingRandomWalk(  pointOnVGS, cube, rootCell, leafCell, nextPoint,  point2d,  *condIt, FRWControl::conductorList ,  
                          	FRWControl::gridOctree, capMatrix ,FRWControl::GreenVT,FRWControl::GExVT,FRWControl::GEyVT,FRWControl::GEzVT ,rand_01, rand_normalDir, Nt, Ng);

        	if(++loopTimer>=10000 ){
        		isFinish=false;
        		pthread_mutex_lock(&updateFRWData_mutex);
        		condIt->gaussianSurfaceList.Nt+=Nt;
        		condIt->gaussianSurfaceList.Ng+=Ng;
        		condIt->gaussianSurfaceList.VGSArea= double(condIt->gaussianSurfaceList.Ng)/double(condIt->gaussianSurfaceList.Nt)* condIt->gaussianSurfaceList.SumOfBGSArea;  

        		for(int j=0; j< capMatrix.size(); j++ ){
        			FRWControl::capMatrix[condIndex][j].update(capMatrix[condIndex][j]);
        			capMatrix[condIndex][j].reset();
        		}
        		isFinish=FRWControl::updateProgress(loopTimer,condIndex, condIt->gaussianSurfaceList.VGSArea);
        		pthread_mutex_unlock(&updateFRWData_mutex);
        		//note: betwwen lock and unlock there should be on jump out!!!
        		loopTimer=0;
        		if(isFinish){
        			break;
        		}
      	  	}
    
		}


	}
	pthread_exit(NULL);

}



void FloatingRandomWalk(  FPoint &pointOnVGS, Cube &cube, Cell* &rootCell, Cell* &leafCell, FPoint &nextPoint,  FPoint &point2d, Conductor  &cond,  ConductorList &condList,  GridOctree &gridOctree,  CapacitanceMatrix &capMatrix 
	,std::vector<std::vector<double> > &GreenVT, vector<vector<double> > &GExVT, vector<vector<double> > &GEyVT, vector<vector<double> > &GEzVT ,
	boost::variate_generator<boost::mt19937&,boost::uniform_01<> > &rand_01, 
    boost::variate_generator<boost::mt19937&,boost::uniform_int<> > &rand_normalDir , int &Nt, int &Ng ){

	// cout<<"next FRW started"<<endl;
	int numOfStep=1;
	double PDFOnVGS;
	int normalDirOfVGS=0;
	int normalDirOfCube=0;
	// FPoint 
	// int Nt, Ng;
	double PDFIntegralElement;
	// Nt=Ng=0;
	PDFIntegralElement=0;
	sampleOnVGSS( pointOnVGS, cond.gaussianSurfaceList , PDFOnVGS, normalDirOfVGS ,  Nt, Ng, PDFIntegralElement, rand_01 ); ///--------------------------------------
	// Cube 
	generateMaxTransitionCube(cube,  rootCell, leafCell ,pointOnVGS,  gridOctree, normalDirOfCube );        ///-------------------
	
	generatePointOnCubeSurface(nextPoint, point2d ,cube, normalDirOfCube, GreenVT ,rand_01, rand_normalDir );   ///----------------

	// double 
	double greenFuncValue=surfaceGreenFunction(  cube, nextPoint, normalDirOfCube, GreenVT  );
	// double 
	double greenFuncElectricField=greenFunctionForElectricField(  cube, nextPoint,  normalDirOfVGS, normalDirOfCube, GExVT, GEyVT, GEzVT );
	
	int unitNormalOfVGS;
	switch(normalDirOfVGS){
		case NX:
		case NY:
		case NZ: unitNormalOfVGS=-1;break;
		case PX:
		case PY:
		case PZ: unitNormalOfVGS=1;break;
		default:cout<<"wrong normalDirOfVGS!"<<endl; exit(1);
	}
	
	double capacitance= FRWControl::dielectricConstant  * greenFuncElectricField*unitNormalOfVGS / greenFuncValue;

	

	while(true){
		//if point is out of GDS::zone or on the border of GDS::zone, return.
		//Note that point on the GDS::zone is not permitted, or generateGridOctree will break down
		if( nextPoint.x1 <= GDS::zone.x1 || nextPoint.x1 >=GDS::zone.x2 || nextPoint.y1 <= GDS::zone.y1 || nextPoint.y1 >= GDS::zone.y2
			|| nextPoint.z1<=GDS::zone.z1 || nextPoint.z1 >=GDS::zone.z2     ){

			double dis2GDSZone=myDistance( nextPoint, GDS::zone);

			if( dis2GDSZone > 1000*GDS::zone.size() ){
				// cout<<"point is too far away from the GDS::zone! END"<<endl;
				return;
			}
			dis2GDSZone+=0.5*GridOctree::gridCellSize;
			cube.setCoordinate( nextPoint.x1- dis2GDSZone , nextPoint.x1+dis2GDSZone, nextPoint.y1- dis2GDSZone, nextPoint.y1+dis2GDSZone, 
					nextPoint.z1- dis2GDSZone, nextPoint.z1+dis2GDSZone   );

			generatePointOnCubeSurface( nextPoint, point2d, cube , normalDirOfCube, GreenVT, rand_01, rand_normalDir  );

			continue;
		}

		for( SubConductorList::iterator neighborCondIt= cube.neighborCondList.begin(); neighborCondIt != cube.neighborCondList.end(); neighborCondIt++   ){
			//@author:Dragon
			//next Point is on the surface of a subCond
			
			if(  fabs(  myDistance( nextPoint,  *neighborCondIt ) )< 1E-7  ){
				int masterCondID= cond.ID;
				int slaveCondID = neighborCondIt->fatherConductorID;
				int masterCondIndex= getCondIndexByID( condList,  masterCondID  );
				int slaveCondIndex=getCondIndexByID( condList, slaveCondID );
				// cout<<masterCondIndex<<"   ,   "<<slaveCondIndex<<endl;

				capMatrix[ masterCondIndex ][slaveCondIndex].update( capacitance );
				return;
			}
		}
		generateMaxTransitionCube(  cube,  rootCell, leafCell ,nextPoint, gridOctree, normalDirOfCube );
		generatePointOnCubeSurface( nextPoint, point2d, cube, normalDirOfCube, GreenVT ,rand_01, rand_normalDir );	
		numOfStep++;
	}

}


//@author:Dragon
//before selecting point from cube, we should use srand( unsigned( time(NULL))) to initiate the seed (just once!)
void generatePointOnCubeSurface(FPoint &pointOnCube, FPoint &point2d,  Cube &cube  ,int &normalDirOfCube , vector<vector<double> > &GreenVT ,
	boost::variate_generator<boost::mt19937&,boost::uniform_01<> > &rand_01, 
    boost::variate_generator<boost::mt19937&,boost::uniform_int<> > &rand_normalDir ){

	// gettimeofday(&t_p_start, NULL );

	int rand_surface=rand_normalDir();
	normalDirOfCube =rand_surface;
	switch(rand_surface ){
		case NY: {

			selectPointFromPlane( point2d,  cube.x1, cube.x2, cube.z1,cube.z2,  cube.size(), greenFuncProgression_NX, greenFuncProgression_NY, GreenVT , rand_01 );
			pointOnCube.setCoordinate( point2d.x1, cube.y1,  point2d.y1 );
			break ;
		}
		case PY: {

			selectPointFromPlane(point2d, cube.x1, cube.x2, cube.z1,cube.z2,  cube.size(), greenFuncProgression_NX, greenFuncProgression_NY ,GreenVT, rand_01 );
			pointOnCube.setCoordinate(  point2d.x1,  cube.y2,   point2d.y1 );
			break ;
		}
		case NX: {

			selectPointFromPlane(point2d, cube.y1, cube.y2, cube.z1, cube.z2,  cube.size(), greenFuncProgression_NX,greenFuncProgression_NY ,GreenVT ,rand_01 );
			pointOnCube.setCoordinate(  cube.x1,   point2d.x1,   point2d.y1  );
			break ;
		}
		case PX: {

			selectPointFromPlane(point2d, cube.y1, cube.y2, cube.z1, cube.z2,  cube.size(), greenFuncProgression_NX,greenFuncProgression_NY ,GreenVT , rand_01 );
			pointOnCube.setCoordinate( cube.x2,  point2d.x1, point2d.y1  );
			break ;
		}
		case PZ: {

			selectPointFromPlane( point2d, cube.x1, cube.x2, cube.y1, cube.y2 , cube.size(), greenFuncProgression_NX, greenFuncProgression_NY ,GreenVT, rand_01 );
			pointOnCube.setCoordinate(  point2d.x1,  point2d.y1, cube.z2  );
			break ;
		}
		case NZ:{

			selectPointFromPlane( point2d, cube.x1, cube.x2, cube.y1, cube.y2 , cube.size(), greenFuncProgression_NX, greenFuncProgression_NY ,GreenVT, rand_01 );
			pointOnCube.setCoordinate(  point2d.x1, point2d.y1, cube.z1  );
			break ;
		}
		default:{
			cout<<"Err: wrong hop direction!"<<endl;
			exit(1);
		}
	}

	// gettimeofday(&t_p_end, NULL);
	// generatePointTime+= t_p_end.tv_sec*1000000+t_p_end.tv_usec-(t_p_start.tv_sec*1000000+t_p_start.tv_usec  );

}

void generatePointOnCubeSurface( FPoint& fp,  FPoint &point2d,  Cube &cube , vector<vector<double> > &GreenVT ,
	boost::variate_generator<boost::mt19937&,boost::uniform_01<> > &rand_01, 
    boost::variate_generator<boost::mt19937&,boost::uniform_int<> > &rand_normalDir  ){

	int temp;
	generatePointOnCubeSurface( fp, point2d ,cube, temp, GreenVT, rand_01, rand_normalDir );

	return ;
}


//use acception-rejection sampling method to create the desired distribution.
void selectPointFromPlane( FPoint &fp,   double x1, double x2, double y1, double y2 , double size, int nX, int nY , vector<vector<double> > &GreenVT, boost::variate_generator<boost::mt19937&,boost::uniform_01<> > &rand_01  ){
	double x;
	double y;
	double vertical_value=0;
	double greenFunc_value=0;
	int nx,ny;
	int GreenVTSize=GreenVT.size();
	// double nz;
	////rand
	while( true ){
		//vertical_value ranges from 0 to 1. That's enough because the max value of greenFunc_value is less than 0.5
		x=rand_01()*size;
		y=rand_01()*size;
		vertical_value=rand_01()*0.5;    
		nx= int( x/(size)/(1.0/GreenVTSize));
		ny= int( y/(size)/(1.0/GreenVTSize));
		if(nx<0||nx> GreenVTSize || ny<0 || ny >GreenVTSize ){
			cout<<"Wrong green surface locate!!"<<endl;
			exit(1);
		}
  	 	nx=	(nx==GreenVTSize? nx-1: nx);
   	 	ny= (ny==GreenVTSize? ny-1:ny);

   	 	greenFunc_value=GreenVT[ny][nx];
		if(vertical_value <= greenFunc_value){
				fp.setCoordinate( x1+ x  , y1+ y,double(0.0)); 
				return;
		}

	}
}



// enum DirOfRectDistance{  PX, PY, PZ, NX, NY, NZ };
//scale_factor is not less than 1!
GaussianSurface generateBGS( SubConductor &subCond, ConductorList  &condList   , double scale_factor  ){


	std::vector<double> tempExtensionDis(6);
	// int rectDirs[6]={ PX, PY, PZ, NX, NY, NZ  };
	int dirOfRect=0;
	int disMax;


	for(  int dirOfRect  =0; dirOfRect <6; dirOfRect++  ){
		//disMax is infinite!
		//@author:Dragon
		//I think disMax should be the distance between the conductor and the gds border!
		// disMax= unsigned( -1)>>2 ;  
		disMax=distanceToGdsBorder( subCond, dirOfRect  );
		for(   ConductorList::iterator condIt = condList.begin(); condIt != condList.end(); condIt ++    ){
		//@author:Dragon
		//very important! when generate block Gaussian surface, do not consider the conductor which the current subConductor belongs to, or there will be 
		// conditions where extensionDis==0, this will bring unexpected influence. Besides, it is meaningless  to consider the sibling subConductors because their
		// voltage is the same.
		
			// if( condIt->ID ==subCond.fatherConductorID ){
			// 	continue; 
			// }
			for( SubConductorList::iterator subCondIt = condIt->subConductorList.begin()  ;   
				subCondIt != condIt->subConductorList.end(); subCondIt ++    ){
				//@author:Dragon
				// if subCond and *subCondIt are actually the same one, then omit it !
				if( subCond== *subCondIt  ){
					continue;
				}
				int distanceAB= myDistance( subCond, *subCondIt );
				// if subCond and *subCondIt are directed connected, then omit it !
				if( condIt->ID == subCond.fatherConductorID && distanceAB==0  ){
					continue;
				}
				if(  myDistance( subCond, *subCondIt , dirOfRect) == distanceAB  ){
					disMax=min(  disMax,  distanceAB  );
				}
			}
		}
		//@author:Dragon
		//disMax can be 0 when two subConds are adjecent. This is permitted.
		tempExtensionDis[dirOfRect]=disMax/2.0;
		if( fabs(tempExtensionDis[dirOfRect])< 1E-5  ){
			cout<<"The extensionSize of BGS cannot be 0!!"<<endl;
			exit(1);
		}
	}

	// cout<<"------------------"<<extensionDis[NX]<<","<<extensionDis[PX]<<","<<extensionDis[NY]<<","<<extensionDis[PY]<<","<<extensionDis[NZ]<<","<<extensionDis[PZ]<<endl;
//@author:Dragon
//if two subConductors are adjacent, then on this side the extension distance of Gaussian Surface is 0, this side does not need to be considered because it will not be sampled
// double extensionDistanceLimit= min( extensionDis[PX] ,extensionDis[PY] ,extensionDis[PZ] ,extensionDis[NX] ,extensionDis[NY] ,extensionDis[NZ] );
	double extensionDistanceLimit=0;
	for( dirOfRect=0; dirOfRect <6; dirOfRect++  ){
		
			if( extensionDistanceLimit==0  ){
				extensionDistanceLimit=tempExtensionDis[dirOfRect];
			}else{
				if( tempExtensionDis[dirOfRect]  < extensionDistanceLimit ){
					extensionDistanceLimit=tempExtensionDis[dirOfRect];
				}
			}
		
	}

	extensionDistanceLimit= extensionDistanceLimit*scale_factor;

	for(dirOfRect=0; dirOfRect<6; dirOfRect++ ){
		tempExtensionDis[dirOfRect]=min( extensionDistanceLimit,  tempExtensionDis[dirOfRect]  );
	}


	// if(subCond.isVia){
	// 	tempExtensionDis[NX]=tempExtensionDis[PX]=tempExtensionDis[PY]=tempExtensionDis[NY]=0.6*subCond.size();
	// }

	return  GaussianSurface( subCond.x1-tempExtensionDis[NX], subCond.x2+tempExtensionDis[PX], subCond.y1-tempExtensionDis[NY], subCond.y2+tempExtensionDis[PY] ,
					subCond.z1- tempExtensionDis[NZ], subCond.z2+ tempExtensionDis[PZ]  , tempExtensionDis  );

}

void generateGaussianSurfaceOfConductor( ConductorList &condList, double scale_factor  ){
	//@author:Dragon
	//Here we need to update the VGSArea of the conductor simultaeously!
	for( ConductorList::iterator condIt=condList.begin(); condIt !=condList.end(); condIt++   ){
		for( SubConductorList::iterator  subCondIt=condIt->subConductorList.begin(); 
			subCondIt != condIt->subConductorList.end();  subCondIt ++    ){
			GaussianSurface blockGaussianSurface=generateBGS( *subCondIt, condList,scale_factor );
			condIt->gaussianSurfaceList.push_back( blockGaussianSurface  );
			double minExtensionDisOfBGS=min(  blockGaussianSurface.extensionDis[PX]  ,   blockGaussianSurface.extensionDis[PY]  , 
						         blockGaussianSurface.extensionDis[PZ]  ,   blockGaussianSurface.extensionDis[NX]  , 
						         blockGaussianSurface.extensionDis[NY]  ,   blockGaussianSurface.extensionDis[NZ]  );

			if(condIt->gaussianSurfaceList.minExtensionDis==0){
				condIt->gaussianSurfaceList.minExtensionDis= minExtensionDisOfBGS;
			}else{
				if(   condIt->gaussianSurfaceList.minExtensionDis > minExtensionDisOfBGS   ){
					condIt->gaussianSurfaceList.minExtensionDis = minExtensionDisOfBGS;
				}
			}
	//@author:Dragon
	//Normally for each BGS in the BGS list, there are about 2 surfaces out of 6 faces are enclosed by or directedly connected to other Gaussian surfaces
	//which means only 4 faces contribute to the area of the virtual Gaussian Surface, so we set the initial VGSArea as 2/3 of the sum of the area of all the BG
			condIt->gaussianSurfaceList.VGSArea+=2.0/3.0*blockGaussianSurface.area();
			condIt->gaussianSurfaceList.SumOfBGSArea+=blockGaussianSurface.area();

		}
	}

	//检查与一个子高斯面相邻的所有子高斯面
	for( ConductorList::iterator condIt= condList.begin(); condIt != condList.end(); condIt++  ){	
		neighborGaussianSurfaceCheck(  condIt->gaussianSurfaceList );
	}
}


//@author:Dragon
//return the distance from rectA (subConductor) to the border of GDS zone. The direction needs to be taken into consideration
int distanceToGdsBorder( Rectangle &rectA,  int normalDirect  ){
	switch(normalDirect){
		case PX:{
			return  GDS::zone.x2-rectA.x2;
		}
		case NX:{
			return	rectA.x1-GDS::zone.x1;
		}
		case PY:{
			return	GDS::zone.y2-rectA.y2;
		}
		case NY:{
			return  rectA.y1-GDS::zone.y1;
		}
		case PZ:{
			return GDS::zone.z2-rectA.z2;
		}
		case NZ:{
			return rectA.z1-GDS::zone.z1;
		}
		default:{
			cout<<"Error: wrong normalDirect!"<<endl;
			exit(1);
		}

	}

}
 
int myDistance( Rectangle &rectA, Rectangle &rectB, int dirOfRectA ){
	switch(dirOfRectA){
		case PX: return rectB.x1-rectA.x2;
		case PY: return rectB.y1-rectA.y2;
		case PZ: return rectB.z1-rectA.z2;
		case NX: return rectA.x1-rectB.x2;
		case NY: return rectA.y1-rectB.y2;
		case NZ: return rectA.z1-rectB.z2;
		default:{
			cout<<"Error: wrong direction of Rect!"<<endl;
			exit(1);
		}

	}
}



void neighborGaussianSurfaceCheck( GaussianSurfaceList &BGSList ){

	GaussianSurfaceList tempList=BGSList;
	for( GaussianSurfaceList::iterator BGSIt = BGSList.begin(); BGSIt != BGSList.end(); BGSIt ++  ){
		for( GaussianSurfaceList::iterator otherBGSIt=tempList.begin(); otherBGSIt != tempList.end() ; otherBGSIt++ ){
			if(  *BGSIt == *otherBGSIt   ){
				continue;
			}
			if( BGSIt->isNeighborWith(  *otherBGSIt )  ){
				 BGSIt->neighborList.push_back( *otherBGSIt );
			}

		}
	}


}



void sampleOnVGSS( FPoint &pointOnVGS, GaussianSurfaceList &BGSList ,double &PDFOnVGS  ,int &normalDirOfVGS , int &Nt, int &Ng , double &PDFIntegralElement, boost::variate_generator<boost::mt19937&,boost::uniform_01<> > &rand_01 ){
	// if( FRWControl::isRandSeedSet==false  ){
	// 	srand( unsigned ( time( NULL ) )  );
	// 	FRWControl::isRandSeedSet=true;
	// }

	//randomly select a BGS from the BGSlist according to their areas
	// double randArea=rand()/( double(RAND_MAX)/BGSList.SumOfBGSArea );
	double randArea=rand_01()*double(BGSList.SumOfBGSArea);
	GaussianSurfaceList::iterator BGSIt;
	while(true){
		for(BGSIt=BGSList.begin(); BGSIt != BGSList.end(); BGSIt ++ ){
			randArea-=BGSIt->area();
			if( randArea<=0 ){
				break;
			}
		}
		//@author:Dragon
		//because this is the double computation, in case some possible error, we add while loop to make sure BGSIt is a meaningful BGS iterator
		if( BGSIt != BGSList.end() ){
			break;
		}
	}
	//randomly select a point from the BGS uniformly.
	//enum NormalDirection{  PX, PY, PZ, NX, NY, NZ };
	randArea =rand_01()*BGSIt->area();
	int normalDir;
	while(true){
		for(normalDir=0; normalDir <6 ; normalDir++  ){
			randArea-=BGSIt->surfaceArea[normalDir];
			if(randArea<=0 ){
				break;
			}

		}
		if(normalDir !=6 ){
			break;
		}
	}
	double xLen=BGSIt->x2-BGSIt->x1;
	double yLen=BGSIt->y2-BGSIt->y1;
	double zLen=BGSIt->z2-BGSIt->z1;
	double dX,dY,dZ;
	dX=rand_01()*xLen;
	dY=rand_01()*yLen;
	dZ=rand_01()*zLen;
	switch(normalDir){
		case PX: pointOnVGS.setCoordinate( BGSIt->x2, BGSIt->y1+dY, BGSIt->z1+dZ ); break;
		case PY: pointOnVGS.setCoordinate( BGSIt->x1+dX,  BGSIt->y2, BGSIt->z1+ dZ  ); break;
		case PZ: pointOnVGS.setCoordinate( BGSIt->x1+dX , BGSIt->y1+dY, BGSIt->z2  ); break;
		case NX: pointOnVGS.setCoordinate( BGSIt->x1, BGSIt->y1+dY, BGSIt->z1+dZ ); break;
		case NY: pointOnVGS.setCoordinate( BGSIt->x1+dX,  BGSIt->y1, BGSIt->z1+ dZ  ); break;
		case NZ: pointOnVGS.setCoordinate( BGSIt->x1+dX , BGSIt->y1+dY, BGSIt->z1  ); break;
		default:{cout<<"Err: wrong normal direction"<<endl; exit(1);}
	}
	//Here we generate one point., so Nt should add 1;
//这个地方可以加速
	Nt++;

	int nc=1;
	int normalDirOfPointOnBGS=BGSIt->getNormalDirOfPoint(  pointOnVGS );
	int normalDirOfPointOnOtherBGS;

	for(  std::list<GaussianSurface>::iterator neighborBGSIt= BGSIt->neighborList.begin(); neighborBGSIt !=BGSIt->neighborList.end();  neighborBGSIt++  ){
	
		if(   neighborBGSIt->isPointEnclosedBySurface( pointOnVGS )  ){
			sampleOnVGSS( pointOnVGS, BGSList, PDFOnVGS  ,  normalDirOfVGS , Nt, Ng, PDFIntegralElement, rand_01 );
			return;
		}

		//@author:Dragon
		//if p is alos on some other BGS whose outer normal direction is opposite to that of curreeent BGS then we should resample
		normalDirOfPointOnOtherBGS=neighborBGSIt->getNormalDirOfPoint(pointOnVGS  );
		if( normalDirOfPointOnOtherBGS ){
			if( ( normalDirOfPointOnBGS+normalDirOfPointOnOtherBGS )==0 ){
				sampleOnVGSS(pointOnVGS, BGSList ,PDFOnVGS, normalDirOfVGS  , Nt, Ng, PDFIntegralElement, rand_01 );
				return;
			}
			nc++;
		}
	}
	double r1= rand_01();
	if( r1 > 1.0/nc ){
		sampleOnVGSS(pointOnVGS, BGSList ,PDFOnVGS,normalDirOfVGS   , Nt, Ng, PDFIntegralElement, rand_01  );
		return;
	}
	//@author:Dragon
	// up to here we sample a point from the virtual gaussian surface uniformly. so at this point Ng should add 1
	Ng++;
	//update the value of VGSArea;
	//PDF on the Gaussian Surface
	// PDFOnVGS=  1.0/ BGSIt->extensionDis[normalDir ];
	// //update the value of the BGSList.PDFIntegralOnVGS
	// PDFIntegralElement+=BGSList.VGSArea*PDFOnVGS  ;
	// double r2= rand()/( double(RAND_MAX) );
	// if( r2 > PDFOnVGS/( 1.0/ BGSList.minExtensionDis )  ){
	// 	sampleOnVGSS(pointOnVGS, BGSList ,PDFOnVGS,normalDirOfVGS   , Nt, Ng, PDFIntegralElement);
	// 	return;
	// }
	// normalDirOfVGS=normalDir;
	// PDFOnVGS=  1.0/ BGSIt->extensionDis[normalDir ];
	//--------------------------
	normalDirOfVGS=normalDir;

}


double surfaceGreenFunction(  Cube &cube, FPoint &p , int normalDirOfCube, vector<vector<double> > &GreenVT ){
	double x,y;
	switch( normalDirOfCube ){
		case NY:
		case PY: x= p.x1-cube.x1; y=p.z1-cube.z1; break;
		case NX:
		case PX: x= p.y1-cube.y1; y= p.z1-cube.z1; break;
		case PZ:
		case NZ: x=p.x1-cube.x1; y= p.y1-cube.y1; break;
		default:{
			cout<<"Err: wrong cube dir!"<<endl;
			exit(1);
		}
	}
	double cubeSize=cube.size();

/**
	This version is to conduct real-time computation
*/
	// int nx, ny;
	// double nz;
	// double greenFuncValue=0;
	// for( nx=1; nx<=10; nx++    ){
	// 	for( ny=1; ny<=10; ny++  ){
	// 		nz=pow(   pow( double(nx), 2.0 )+pow(double(ny), 2), 0.5 );
	// 		greenFuncValue+= 2.0*sin( M_PI * nx / 2 )*sin( M_PI * ny / 2)*sin( M_PI* nx * x / cubeSize )*sin( M_PI * ny *y / cubeSize ) / cosh(M_PI* nz/2 );

	// 	}
	// }
	// return  greenFuncValue/pow( double(cubeSize),2.0);

/**
	This version is to use offline table
*/

	// int GreenVTSize=GreenVT.size();
	// x=x/cubeSize;
	// y=y/cubeSize;
	// int nx= int( x /(1.0/GreenVTSize));
	// int ny= int( y/ (1.0/GreenVTSize));
	// if(nx<0||nx> GreenVTSize || ny<0 || ny >GreenVTSize ){
	// 	cout<<"Wrong green surface locate!!"<<endl;
	// 	exit(1);
	// }
 //    nx=	(nx==GreenVTSize? nx-1: nx);
 //    ny= (ny==GreenVTSize? ny-1:ny);
 //    //@author:Dragon
 //    /*根据输出的矩阵的特点： a[i][j]，i表示第i行，对应的是y的坐标， j是第j列对应的是x的坐标
	// 	------------> x+
	// 	|
	// 	|
	// 	|
	// 	V
	// 	 y+
	// 	 所以这里的ny是第一维度，nx是第二维度
 //    */
 //    return GreenVT[ny][nx]/pow(double(cubeSize),2.0);


    int GreenVTSize=GreenVT.size();
	x=x/cubeSize;
	y=y/cubeSize;
	int nx= int( x /(1.0/GreenVTSize));
	int ny= int( y/ (1.0/GreenVTSize));
	if(nx<0||nx> GreenVTSize || ny<0 || ny >GreenVTSize ){
		cout<<"Wrong green surface locate!!"<<endl;
		exit(1);
	}
    nx=	(nx==GreenVTSize? nx-1: nx);
    ny= (ny==GreenVTSize? ny-1:ny);

    return GreenVT[ny][nx]/pow(double(cubeSize),2.0);

}




//enum NormalDirection{  PX, PY, PZ, NX, NY, NZ };
/*
conduct this conversion to the normalDirection
PX 	1
PY 	2
PZ 	3
NX 	-1
NY 	-2
NZ 	-3
*/
//note:  x1, y1,  is the relative coordinates with respect to the min coordinates on the cube
// for symmetry, nX should be equal to nY
// int normalDirUnsigned2Signed[6]={ 1, 2 , 3 , -1, -2, -3   };
double GEx( double x, double y ,double L, vector<vector<double> > &GExVT  ){

/**
	This version is to conduct real-time computation
*/
	// int nx, ny;
	// double nz;
	// double GE=0;
	// for(nx=1;nx<=10; nx++  ){
	// 	for( ny=1; ny<=10; ny++ ){
	// 		nz= pow( pow(double(nx),2.0)+pow( double(ny),2.0 ), 0.5      );
	// 		GE+= (-1.0)*2*M_PI * nx * cos(M_PI*nx/2.0)*sin(M_PI*ny/2.0)/cosh(M_PI*nz/2.0) *sin(M_PI*nx*x/L)*sin(M_PI*ny*y/L);
	// 	}

	// }
	// return GE/pow(L, 3.0);

/****
	This version is to use the offline data table
*/
	x=x/L;
	y=y/L;
	int GExVTSize=GExVT.size();
	int nx= int( x /(1.0/GExVTSize));
	int ny= int( y/ (1.0/GExVTSize));
	if(nx<0||nx> GExVTSize || ny<0 || ny >GExVTSize ){
		cout<<"Wrong GEx surface locate!!"<<endl;
		exit(1);
	}
    nx=	(nx==GExVTSize? nx-1: nx);
    ny= (ny==GExVTSize? ny-1:ny);
    return GExVT[ny][nx]/pow(L, 3.0);




}

double GEy(double x, double y ,double L,vector<vector<double> > &GEyVT  ){
/**
	This version is to conduct real-time computation
*/
	// int nx, ny;
	// double nz;
	// double GE=0;
	// for(nx=1;nx<=10; nx++  ){
	// 	for( ny=1; ny<=10; ny++ ){
	// 		nz= pow( pow(double(nx),2.0)+pow( double(ny),2.0 ), 0.5      );
	// 		GE+= (-1.0)*2*M_PI * ny * sin(M_PI*nx/2.0)*cos(M_PI*ny/2.0)/cosh(M_PI*nz/2.0) *sin(M_PI*nx*x/L)*sin(M_PI*ny*y/L);
	// 	}

	// }
	// return GE/pow(L, 3.0);

/****
	This version is to use the offline data table
*/
	x=x/L;
	y=y/L;
	int GEyVTSize=GEyVT.size();
	int nx= int( x /(1.0/GEyVTSize));
	int ny= int( y/ (1.0/GEyVTSize));
	if(nx<0||nx> GEyVTSize || ny<0 || ny >GEyVTSize ){
		cout<<"Wrong GEy surface locate!!"<<endl;
		exit(1);
	}
    nx=	(nx==GEyVTSize? nx-1: nx);
    ny= (ny== GEyVTSize? ny-1:ny);
    return GEyVT[ny][nx]/pow(L, 3.0);

}

double GEz(double x, double y ,double L ,vector<vector<double> > &GEzVT  ){
/**
	This version is to conduct real-time computation
*/
	// int nx, ny;
	// double nz;
	// double GE=0;
	// for(nx=1;nx<=10; nx++  ){
	// 	for( ny=1; ny<=10; ny++ ){
	// 		nz= pow( pow(double(nx),2.0)+pow( double(ny),2.0 ), 0.5      );
	// 		GE+= (-1.0)* 2*M_PI * nz * sin(M_PI*nx/2.0)*sin(M_PI*ny/2.0)/sinh(M_PI*nz/2.0) *sin(M_PI*nx*x/L)*sin(M_PI*ny*y/L);
	// 	}

	// }		
	// return GE/pow(L, 3.0) ;

/****
	This version is to use the offline data table
*/
	x=x/L;
	y=y/L;
	int GEzVTSize=GEzVT.size();
	int nx= int( x /(1.0/GEzVTSize));
	int ny= int( y/ (1.0/GEzVTSize));
	if(nx<0||nx> GEzVTSize || ny<0 || ny >GEzVTSize ){
		cout<<"Wrong GEz surface locate!!"<<endl;
		exit(1);
	}
    nx=	(nx==GEzVTSize? nx-1: nx);
    ny= (ny== GEzVTSize? ny-1:ny);
    return GEzVT[ny][nx]/pow(L, 3.0);

}

double greenFunctionForElectricField(  Cube &cube,  FPoint &p,  int normalDirOfVGS, int normalDirOfCube ,vector<vector<double> > &GExVT, vector<vector<double> > &GEyVT, vector<vector<double> > &GEzVT  ){
	
	double GE=0;
	double L=cube.size();
	int sign;
	double x,y;
	switch(normalDirOfCube ){
		case NX:
		case PX: x=p.y1-cube.y1; y=p.z1-cube.z1;break;
		case NY:
		case PY: x=p.x1-cube.x1; y=p.z1-cube.z1;break;
		case NZ:
		case PZ: x=p.x1-cube.x1; y=p.y1-cube.y1;break;
		default:{
			cout<<"wrong normalDirOfCube!"<<endl;
			exit(1);
		}

	}


	if( normalDirOfVGS==PZ || normalDirOfVGS==NZ  ){
		if( normalDirOfCube==PZ ||normalDirOfCube==NZ ){
			//式子3
			sign=normalDirOfCube==PZ?1:-1;
			GE=sign* GEz( x, y, L , GEzVT );
		}else{
			//式子2
			GE=GEy(x,y,L ,GEyVT );
		}

	}else if( normalDirOfVGS==NX || normalDirOfVGS==PX ){
		if(normalDirOfCube==NX||normalDirOfCube==PX ){
			//式子3
			sign= normalDirOfCube==PX?1:-1;
			GE=sign* GEz(x,y,L, GEzVT );

		}else{
			//式子1
			GE=GEx(x,y,L, GExVT );
		}


	}else if( normalDirOfVGS== NY || normalDirOfVGS==PY ){
		if(normalDirOfCube==NY ||normalDirOfCube==PY ){
			sign= normalDirOfCube==PY?1:-1;
			GE=sign* GEz(x,y,L ,GEzVT );
		}else if( normalDirOfCube==NX || normalDirOfCube==PX ){
			//式子1
			GE=GEx(x,y,L,GExVT );
		}else{
			//式子2
			GE=GEy(x,y,L ,GEyVT );
		}

	}else{
		cout<<"wrong normal direction!!"<<endl;
		exit(1);
	}

	return GE;
}


//@author:Dragon
//Here we assume that FPoint must be within the GDS::zone, if p is on the surface of GDS::zone or on the outside of GDS::zone, the random wlak will terminate
//对于Cell的操作改成使用指针传递
Cell* locateOctreeRootCell(  FPoint &p, GridOctree &gridOctree  ){
	int gridCellSize=GridOctree::gridCellSize;
	int xn= int (p.x1- GDS::zone.x1)/gridCellSize;
	int yn= int (p.y1 -GDS::zone.y1)/gridCellSize;
	int zn= int (p.z1 -GDS::zone.z1) /gridCellSize;

	// cout<<xn<<"  , "<<yn<<" , "<<zn<<endl;

	Cell *rootCell;
	rootCell = &(gridOctree[xn][yn][zn]);
	if(p.x1-rootCell->x1<-1E-7||p.x1-rootCell->x2>1E-7||p.y1-rootCell->y1<-1E-7||p.y1-rootCell->y2>1E-7 ||p.z1-rootCell->z1<-1E-7||p.z1-rootCell->z2>1E-7 ){
		cout<<"p is not in root Cell!!!"<<endl;
		cout<<"P:"<<p.x1<<","<<p.y1<<","<<p.z1<<endl;
		cout<<"CELL:"<<rootCell->x1<<","<<rootCell->x2<<","<<rootCell->y1<<","<<rootCell->y2<<","<<rootCell->z1<<","<<rootCell->z2<<endl;
		exit(1);
	}
	return rootCell;
}


Cell* locateCellFromRoot( FPoint &p, Cell *rootCell , Cube &cube ){
	
	if(rootCell->isFilledWithCond ){
		cout<<"FRW error in  locateCellFromRoot, FRWcapext.cpp: FRW walks into the inside of a conductor! "<<endl;
		exit(1);
	}

	if( !rootCell->hasChildCell ){
		return rootCell;
	}

	for( SubCellList::iterator subCellIt = rootCell->subCellList.begin(); subCellIt != rootCell->subCellList.end(); subCellIt++     ){
		if(   p.x1 - subCellIt->x1>-1E-7 && p.x1 -subCellIt->x2 <1E-7 
			&& p.y1 - subCellIt->y1 > -1E-7 && p.y1 - subCellIt->y2  <1E-7 
				&& p.z1 - subCellIt->z1 > -1E-7 && p.z1 - subCellIt->z2 <1E-7 ){
			//@author:Dragon
			//这里放宽了定位的条件，目的是确保一定可以定位到cell里面

			return locateCellFromRoot(  p, &(*subCellIt), cube );
			
		}

	}
	cout<<"locate Cell accoring to Point failed !!!"<<endl;
	exit(1);
}





void generateMaxTransitionCube(Cube &maxCube, FPoint &p, Cell* &cell ){
//@author:Dragon
// find the minimum distance from point to the candidate SubConductors, then compare the value with distance limit the lesser is the size of transition cube
//the initial value of cubeRadius is the distance from point to the border of inflated Cell.	

	// double x_Dis=max( cell->x1-cell->extensionSize-p.x1, p.x1- ( cell->x2 +cell->extensionSize )    );
	// double y_Dis=max( cell->y1-cell->extensionSize-p.y1, p.y1- ( cell->y2 +cell->extensionSize )    );
	// double z_Dis=max( cell->z1-cell->extensionSize-p.z1, p.z1- ( cell->z2 +cell->extensionSize )    );
	double x_Dis=max( cell->x1-cell->distanceLimit-p.x1, p.x1- ( cell->x2 +cell->distanceLimit  )    );
	double y_Dis=max( cell->y1-cell->distanceLimit-p.y1, p.y1- ( cell->y2 +cell->distanceLimit )    );
	double z_Dis=max( cell->z1-cell->distanceLimit-p.z1, p.z1- ( cell->z2 +cell->distanceLimit )    );

	double temp_XY=max( x_Dis,y_Dis );
	double cubeRadius= fabs( max( temp_XY,z_Dis));


	// cout<<"IIIIOOOOOOIIIIII"<<endl;
// the traditional way:
	

	for( SubConductorList::iterator candidateIt= cell->candidateList.begin(); candidateIt != cell->candidateList.end(); candidateIt++  ){
		//the improved way!!
		if(  myDistance ( *candidateIt, *cell )  - cubeRadius >1E-7   ){    //should not be >=
			break;
		}
		double tempDis=myDistance ( p, *candidateIt ) ;
		if(  tempDis -cubeRadius <= -1E-7   ){
			cubeRadius=tempDis ;
			maxCube.neighborCondList.clear();
			maxCube.neighborCondList.push_back( *candidateIt  );
		}else if( fabs(tempDis - cubeRadius)<1E-7  ){
			maxCube.neighborCondList.push_back( *candidateIt );
		}

	}

	maxCube.setCoordinate( p.x1-cubeRadius , p.x1 +cubeRadius, p.y1-cubeRadius, p.y1+cubeRadius, p.z1-cubeRadius, p.z1+cubeRadius );
}


void generateMaxTransitionCube( Cube &maxCube, Cell* &rootCell, Cell* &leafCell,  FPoint &p, GridOctree &gridOctree, int &normalDirOfCube ){

	totalStepNum++;
	// gettimeofday(&t_locatecell_start, NULL);
	//@author:Dragon
	//first check if the point is within the leafCell itself

//------------method1
	if( leafCell!=NULL && p.x1 - leafCell->x1> -1E-7 && p.x1 - leafCell->x2 < 1E-7 &&
		p.y1 - leafCell->y1> -1E-7 && p.y1 - leafCell->y2 < 1E-7 &&
		p.z1 - leafCell->z1> -1E-7 && p.z1 - leafCell->z2 < 1E-7   ){
		hitTime++;
	}	
	else{ 
		rootCell= locateOctreeRootCell( p, gridOctree  );

		int nx=(int) ((p.x1-rootCell->x1)/Cell::subGridSize);
		int ny=(int) ((p.y1-rootCell->y1)/Cell::subGridSize);
		int nz=(int) ((p.z1-rootCell->z1)/Cell::subGridSize);
		int subGridNum=rootCell->subGridCellVec.size();

		nx= (nx==subGridNum? nx--: nx);
		ny= (ny==subGridNum? ny--: ny);
		nz= (nz==subGridNum? nz--: nz);
		leafCell= rootCell->subGridCellVec[nx][ny][nz];

		if(leafCell->hasChildCell){
			leafCell= locateCellFromRoot( p, leafCell, maxCube );
		}

	}

	// gettimeofday(&t_locatecell_end, NULL );
	// locateCellTime+= t_locatecell_end.tv_sec*1000000+t_locatecell_end.tv_usec-(t_locatecell_start.tv_sec*1000000+t_locatecell_start.tv_usec );


	// gettimeofday(&t_cube_start, NULL);
	generateMaxTransitionCube(  maxCube ,p,  leafCell );

	// gettimeofday(&t_cube_end,NULL );
	// generateMaxCubeTime+= t_cube_end.tv_sec*1000000+t_cube_end.tv_usec - ( t_cube_start.tv_sec*1000000+t_cube_start.tv_usec);

	//to check if I walk right!
	// for(ConductorList::iterator condIt=FRWControl::conductorList.begin();condIt !=FRWControl::conductorList.end(); condIt++   ){
	// 	for(SubConductorList::iterator subCondIt=condIt->subConductorList.begin(); subCondIt!= condIt->subConductorList.end(); subCondIt++ ){
	// 		if(myDistance( maxCube, *subCondIt  )<-1E-5 ){
	// 			cout<<"walk wrong!!"<<endl;
	// 			exit(1);
	// 		}

	// 	}
	// }
	
	
	
}

/*
Cij  the capacitance from conductor i to conductor j

C11  C12 ... C1n
C21  C22 ... C2n
...
Cn1  Cn2 ... Cnn
*/

CapacitanceMatrix generateCapMatrix( ConductorList  &condList ){
	int numCond= condList.size();
	if(numCond==0){
		cout<<"Err: the number of conductor is 0!"<<endl;
		exit(1);
	}
	std::vector<Capacitance> v1(numCond);
	std::vector<  std::vector<Capacitance> > capMatrix( numCond, v1  );

	for(int i=0; i<capMatrix.size(); i++ ){
		for(int j=0; j<capMatrix[i].size(); j++){
			capMatrix[i][j].masterCondID=i;
			capMatrix[i][j].slaveCondID=j;
		}
	}
	return capMatrix;
}
//getCondIndexInConductorListByID
int getCondIndexByID( ConductorList &condList ,  int ID){
	int index=-1;
	for( ConductorList::iterator condIt=condList.begin(); condIt != condList.end(); condIt++  ){
		index++;
		if(  condIt->ID==ID  ){
			return index;
		}
	}
	return -1;
}


double myDistance( FPoint &p, Cube &cube  ){
	double xDis= max(  cube.x1-p.x1, p.x1-cube.x2  );
	double yDis= max(  cube.y1-p.y1, p.y1-cube.y2 );
	double zDis=max( cube.z1-p.z1, p.z1-cube.z2 );
	double temp= max(xDis,yDis);
	return max(temp,zDis);

}


vector<string> split( string str, string pattern){
	vector<string> ret;
	if(pattern.empty()){
		return ret;
	}
	size_t start=0,index=str.find_first_of(pattern,0);
	while(index!=string::npos)
	{
		if(start!=index){
			ret.push_back(str.substr(start,index-start));
		}
		start=index+1;
		index=str.find_first_of(pattern,start);
	}
	if(!str.substr(start).empty()){
		ret.push_back(str.substr(start));
	}
	return ret;
}

void loadDataTable(const char *dataFile, vector<vector<double> > &dataTable , int &tableSize ){
	ifstream fin(dataFile);
	string line;
	while(getline(fin, line)){
		std::vector<string> strVec=split(line,"\t");
		std::vector<double> numVec;
		for(int i=0; i<strVec.size(); i++){
			double temp;
			if(sscanf( strVec[i].c_str(),"%lf", &temp ) !=1){
				cout<<"Data Table format error!"<<endl;
				exit(1);
			}
			numVec.push_back(temp);
		}
		dataTable.push_back(numVec);
	}
	fin.close();
	tableSize=dataTable.size();
}

void copy(ofstream &fout, string &filename ){
	ifstream fin(filename.c_str());
	string str;
	while(getline(fin,str)){
		fout<<str<<"\n"<<flush;
	}
	fin.close();
}
