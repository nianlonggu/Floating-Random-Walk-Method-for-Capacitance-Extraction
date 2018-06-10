#include <iostream>
#include "geo2condlist/geoloader.cpp"

using namespace std;



// int maxCandidateListLen=10;
// int minCellSize=180;
// int extensionSize =600;

//@author: NOTE: in the FRW process, the datatype of the coordinates is double ! for the reason of accuracy!
class FPoint{
public:
	double x1,y1,z1;
	FPoint();
	FPoint(double x_1, double y_1, double z_1 );
	void setCoordinate( double x_1, double y_1, double z_1 ){
		x1=x_1;
		y1=y_1;
		z1=z_1;
	}
	
};

//@author:Dragon
//为了精确计算，把Cell的坐标改成double类型
class Cell
{

public:
	static  int maxCandidateListLen;
	static  double minCellSize;                 //2 times of the min width (nm)
	static  double subGridSize;                 //这个变量是用来规定每一个八叉树空间的栅格的尺寸
	//static  int extensionSize;             //about 7 times of the min width (nm)
	double x1, x2, y1 ,y2, z1, z2;


	SubConductorList candidateList;

	//author:Dragon
	//this variable store the total set of the candidates, only useful for the rootCell of the Octree!
	SubConductorList candidateListOfRootCell;

	//@author:Dragon
	//distanceLimit 	represents the max distance from any point in the cell to the surrounding subConductor.
	double extensionSize;
	double distanceLimit;
	bool isFilledWithCond;  //to express if the Cell is totally filled with SubConductor. Note: subConductor, not the whole Condcutor!
	bool isFilledChecked;	//to express if the FilledWithCond is checked.
	list<Cell> subCellList;     //using list is more flexible than vector


	list<Cell*> candidateCellList;
	FPoint centerPoint;
	bool hasChildCell;

	vector<vector<vector<Cell*> > > subGridCellVec;

	Cell();
	Cell( double x_1, double x_2, double y_1, double y_2, double z_1, double z_2  , double extension);
	Cell(double x_1, double x_2, double y_1, double y_2, double z_1, double z_2  , double extension, double disLimit);

	double size();
	void divideIntoSubCell( double extension );
	Cell getInflatedCell( double extension  );  // not used for now
	bool isInsectWith(Rectangle &subCond );
	bool isIncludedIn(Rectangle &subCond );

};

//注意在本例子中一个Cell肯定是标准的立方体

// int Cell::extensionSize=-1;

typedef std::list<Cell> SubCellList;

class GridOctree:public std::vector< std::vector< std::vector<Cell> > > {
public:
	static int gridCellSize;
};

int GridOctree::gridCellSize=320 ;
int Cell::maxCandidateListLen=4;
//这个表示从根节点（即第一级的3D数组中的一个单元）开始看，八叉树的深度
int depthOfOctreeFromRoot=4;
//这个表示从第二级3D数组的一个单元开始看，八叉树的深度
//depthOfOctreeFromSubGrid=depthOfOctreeFromRoot时表示没有第二级grid,当depthOfOctreeFromSubGrid等于0时，意味着在定位Cell的时候使用的
//是两级3D数组，并没有应用八叉树
int depthOfOctreeFromSubGrid=0;

double Cell::minCellSize= GridOctree::gridCellSize /pow(2.0, double(depthOfOctreeFromRoot) ) ;
double Cell::subGridSize= pow(2.0, double(depthOfOctreeFromSubGrid))*Cell::minCellSize ;     



void generateGridOctree(GridOctree &gridOctree,  int gridCellSize, ConductorList &conductorList );


GridOctree partitionDomain(  double cellSize, int &nX, int &nY, int &nZ );


void  insertToOctree( SubConductor &subCond, Cell &cellT);

void candidateCheck ( SubConductor &subCond,  Cell &cellT );

//distance is an overloaded function

double myDistance( FPoint  &p, Rectangle &rect );
double myDistance( Rectangle &rectA, Cell &cellT );
int myDistance(Rectangle &rectA, Rectangle &rectB );
double myDistance( FPoint &p, Cell &cell  );

int max( int a, int b );
double max( double a, double b  );

int min(int a, int b);
int min( int a1, int a2, int a3, int a4, int a5  , int a6  );
double min( double a1, double a2 );
double min( double a1, double a2, double a3, double a4, double a5, double a6  );


//rectA dominate rectB with recpect to rectT
bool dominate( Rectangle &rectA ,Rectangle &rectB, Cell &cellT );
//to check if the surface of the conductor is within the cell----Time-consuming!
bool isCondSurfaceInCell( ConductorList::iterator condIt , Cell &cellT );

Rectangle getGdsZone(ConductorList &conductorList);


void updateSubGridCellVec(Cell &rootCell,  Cell *cellP,  vector<vector<vector<Cell*> > >  &subGridCellVec );
void updateSubGridCellVec(GridOctree &gridOctree);




