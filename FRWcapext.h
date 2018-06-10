#include <iostream>
#include <string>
#include "spacemanagement.cpp"
#include <sys/time.h>
#include <unistd.h>
#include <boost/random.hpp>
#include <pthread.h>


using namespace std;

const int greenFuncProgression_NX=5;
const int greenFuncProgression_NY=5;


// double testTime=0;
double generateMaxCubeTime=0;
double generatePointTime=0;
double locateCellTime=0;
struct timeval t_cube_start,t_cube_end, t_p_start, t_p_end, t_locatecell_start, t_locatecell_end  ;
struct timeval tv_begin,tv_end;


long hitTime=0;
long totalStepNum=0;

pthread_mutex_t getRandSeed_mutex;
pthread_mutex_t updateFRWData_mutex;
pthread_mutex_t processorIndex_mutex;  //在设定线程亲和力的时候对线程加锁



class Cube{
public:
	double x1,x2,y1,y2,z1,z2;

	//@author:Dragon
	//neighborCondList stores the subconductors which is in connection with the cube's surface;
	SubConductorList neighborCondList;

	Cube();
	Cube( FPoint &minPoint, double size );
	Cube( double x_1, double y_1 , double z_1, double size);
	Cube(double x_1, double x_2, double y_1, double y_2, double z_1, double z_2  );

	double size();
	void setCoordinate(  double x_1, double x_2, double y_1, double y_2, double z_1, double z_2  ){
		x1=x_1;
		x2=x_2;
		y1=y_1;
		y2=y_2;
		z1=z_1;
		z2=z_2;
	}

};

//@author:Dragon
//the returned value is the absolute coordinate of the point
void generatePointOnCubeSurface( FPoint &pointOnCube, FPoint &point2d,  Cube &cube , int &normalDirOfCube , vector<vector<double> > &GreenVT,
	boost::variate_generator<boost::mt19937&,boost::uniform_01<> > &rand_01, 
    boost::variate_generator<boost::mt19937&,boost::uniform_int<> > &rand_normalDir );
void generatePointOnCubeSurface( FPoint &pointOnCube, FPoint &point2d, Cube &cube  , vector<vector<double> > &GreenVT,
	boost::variate_generator<boost::mt19937&,boost::uniform_01<> > &rand_01, 
    boost::variate_generator<boost::mt19937&,boost::uniform_int<> > &rand_normalDir );

void selectPointFromPlane( FPoint &fp,   double x1, double x2, double y1, double y2 , double size, int nX, int nY , vector<vector<double> > &GreenVT ,
		boost::variate_generator<boost::mt19937&,boost::uniform_01<> > &rand_01);

// void testGeneratePoint();

class GaussianSurface{
public:
	double x1,x2,y1,y2,z1,z2;
	std::list<GaussianSurface> neighborList;

	GaussianSurface();
	GaussianSurface( FPoint &minPoint, FPoint &maxPoint , std::vector<double> extensionDistance);
	GaussianSurface(double x_1, double x_2, double y_1, double y_2, double z_1, double z_2, std::vector<double> extensionDistance  );

	double area();
	bool isPointOnSurface(FPoint &point);
	bool isPointEnclosedBySurface( FPoint &point );
	int getNormalDirOfPoint(FPoint  &point);
	std::vector<double> extensionDis;
	std::vector<double> surfaceArea;
	void setSurfaceArea();
	bool operator ==( GaussianSurface &rhs  );
	bool isNeighborWith( GaussianSurface &gaussian );

};



enum NormalDirection{  PX, PY, PZ, NX, NY, NZ };
void generateGaussianSurfaceOfConductor( ConductorList &condList, double scale_factor  );

GaussianSurface generateBGS( SubConductor &subCond, ConductorList  &condList   , double scale_factor   );
int distanceToGdsBorder( Rectangle &rectA,  int normalDirect  );
//@author:Dragon
//myDistance take the relative direction of rectA and rectB into consideration.
int myDistance( Rectangle &rectA, Rectangle &rectB, int dirOfRectA );




void neighborGaussianSurfaceCheck( GaussianSurfaceList &BGSList );

//@author:Dragon
//sample a point from the virtual Gaussian surface of a conductor
void sampleOnVGSS( FPoint &pointOnVGS, GaussianSurfaceList &BGSList , double &PDFOnVGS ,  int &normalDirOfVGS,  int &Nt, int  &Ng, double &PDFIntegralElement, boost::variate_generator<boost::mt19937&,boost::uniform_01<> > &rand_01);

//@author:Dragon
//Here the x1, y1 is the relative coordinate
// double surfaceGreenFunction( double cubeSize, double x1, double y1,  int nX, int nY  );
double surfaceGreenFunction(  Cube &cube,  FPoint &p, int normalDirOfCube, vector<vector<double> > &GreenVT );

double greenFunctionForElectricField(  Cube &cube,  FPoint &p, int normalDirOfVGS, int normalDirOfCube ,vector<vector<double> > &GExVT, vector<vector<double> > &GEyVT, vector<vector<double> > &GEzVT  );

//@author:this fucntion can be improved
Cell *locateOctreeRootCell( FPoint &p, GridOctree &gridOctree );

Cell* locateCellFromRoot(  FPoint &p, Cell *rootCell ,Cube &cube );





void generateMaxTransitionCube(Cube &maxCube, FPoint &p, Cell* &cell );
void generateMaxTransitionCube( Cube &maxCube, Cell* &rootCell, Cell* &leafCell, FPoint &p, GridOctree &gridOctree,int &normalDirOfCube );


class Capacitance{
public:
	int masterCondID;
	int slaveCondID;
	double sumOfCap;
	double sumOfCapSquared;
	double capacitance;
	double estimateErr;
	double std_var;

	Capacitance(){
		masterCondID=slaveCondID=-1;
		sumOfCap=0;
		sumOfCapSquared=0;
		capacitance=0;
		estimateErr=10000;
		std_var=10000;
	}
	void update( double cap ){
		sumOfCap+=cap;
		sumOfCapSquared+=pow(cap,2.0);
	}

	void update(Capacitance &capEntity){
		sumOfCap+= capEntity.sumOfCap;
		sumOfCapSquared+= capEntity.sumOfCapSquared;
	}

	void reset(){
		sumOfCap=0;
		sumOfCapSquared=0;
		capacitance=0;
		estimateErr=10000;
		std_var=10000;
	}

};

typedef std::vector< std::vector<Capacitance> > CapacitanceMatrix;



CapacitanceMatrix generateCapMatrix( ConductorList  &condList );

int getCondIndexByID( ConductorList &condList, int  ID );

void FloatingRandomWalk(  FPoint &pointOnVGS, Cube &cube, Cell* &rootCell, Cell* &leafCell, FPoint &nextPoint,  FPoint &point2d, Conductor  &cond,  ConductorList &condList,  GridOctree &gridOctree,  CapacitanceMatrix &capMatrix , 
		std::vector<std::vector<double> > &GreenVT, vector<vector<double> > &GExVT, vector<vector<double> > &GEyVT, vector<vector<double> > &GEzVT,
		boost::variate_generator<boost::mt19937&,boost::uniform_01<> > &rand_01, 
    boost::variate_generator<boost::mt19937&,boost::uniform_int<> > &rand_normalDir , int &Nt, int &Ng );

double myDistance( Cube &cube,  Rectangle &rect );
double myDistance( FPoint &p, Cube &cube  );


class FRWControl{
public:
	static bool isRandSeedSet;
	static double dielectricConstant;
	static ConductorList conductorList;
	static GridOctree gridOctree;
	static CapacitanceMatrix capMatrix;
	static vector<vector<double> > GreenVT;
	static vector<vector<double> > GExVT;
	static vector<vector<double> > GEyVT;
	static vector<vector<double> > GEzVT;
	static int GreenVTSize;
	static int GExVTSize;
	static int GEyVTSize;
	static int GEzVTSize;

	static int totalProcessorNum;
	static int currentProcessorIndex;

	static int currentThreadIndex;

	static vector<long> totalNumOfWalkOfCond;
	static vector<long> currentNumOfWalkOfCond;
	static vector<long> currentProgressOfCond;
	static bool updateProgress(int tempWalkNum, int masterCondIndex ,double temp_VGSArea );

	static vector<long> randSeedVec;
	static double maxPermittedCapErr;


};
bool FRWControl::isRandSeedSet=false;
double FRWControl::dielectricConstant=3.9*8.85;        //unit: pF/m
double FRWControl::maxPermittedCapErr=2;
vector<long> FRWControl::totalNumOfWalkOfCond;
vector<long> FRWControl::currentNumOfWalkOfCond;
vector<long> FRWControl::currentProgressOfCond;

int FRWControl::GreenVTSize=0;
int FRWControl::GExVTSize=0;
int FRWControl::GEyVTSize=0;
int FRWControl::GEzVTSize=0;
int FRWControl::totalProcessorNum=0;
int FRWControl::currentProcessorIndex=0;
int FRWControl::currentThreadIndex=0;

ConductorList FRWControl::conductorList; //= ConductorList();
GridOctree FRWControl::gridOctree;//= GridOctree();
CapacitanceMatrix FRWControl::capMatrix;//= CapacitanceMatrix();
vector<vector<double> >  FRWControl::GreenVT;
vector<vector<double> >  FRWControl::GExVT;
vector<vector<double> >  FRWControl::GEyVT;
vector<vector<double> >  FRWControl::GEzVT;
vector<long> FRWControl::randSeedVec;


vector<string> split( string str, string pattern);
void loadDataTable(const char *dataFile, vector<vector<double> > &dataTable, int &tableSize  );

void copy(ofstream &fout, string &filename );
