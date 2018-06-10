#include <iostream>
#include "spacemanagement.h"

using namespace std;

FPoint::FPoint(){
	x1=y1=z1=0;
}

FPoint::FPoint( double x_1, double y_1, double z_1 ){
	x1=x_1;
	y1=y_1;
	z1=z_1;
}


Cell::Cell()
{	
	//@author:Dragon
	//the initial value of Cell::extensionSize is -1
	//Let the initial value of distanceLimit be set as the Cell::extensionSize can block the subConductors to accomplish the imcomplete candidateList
	x1=x2=y1=y2=z1=z2=0;
	distanceLimit=extensionSize =-1;	 
	isFilledWithCond=false;
	isFilledChecked=false;
	hasChildCell=false;

	
 }

Cell::Cell(double x_1, double x_2, double y_1, double y_2, double z_1, double z_2  , double extension ){
	x1=x_1;
	x2=x_2;
	y1=y_1;
	y2=y_2;
	z1=z_1;
	z2=z_2;
	distanceLimit= extensionSize = extension;
	isFilledWithCond=false;
	isFilledChecked=false;
	hasChildCell=false;

	
}

Cell::Cell(double x_1, double x_2, double y_1, double y_2, double z_1, double z_2  , double extension, double disLimit)  {
	x1=x_1;
	x2=x_2;
	y1=y_1;
	y2=y_2;
	z1=z_1;
	z2=z_2;
	extensionSize=extension;
	distanceLimit=disLimit;
	isFilledWithCond=false;
	isFilledChecked=false;
	hasChildCell=false;

	
}

//This fuction is used to inflate the Cell to construct the incomplete candidate list.
Cell  Cell::getInflatedCell(double extension){
	return Cell( x1-extension, x2+extension, y1-extension, y2+extension, z1-extension, z2+extension , -1 );
}

//@author:Dragon
//这里为了方便起见cell的设置成标准的立方体
double Cell::size(){
	if( x2-x1<1E-5 ||y2-y1<1E-5 ||z2-z1<1E-5 ){
		cout<<"The size of Cell is 0 !";
		exit(1);
	}
	return max(x2-x1, max(y2-y1, z2-z1) ) ;
}

bool Cell::isInsectWith(Rectangle &subCond ){
	return x2-subCond.x1>1E-7 && x1-subCond.x2<-1E-7&& y2-subCond.y1>1E-7 && y1-subCond.y2<-1E-7&&
			 z2-subCond.z1>1E-7 && z1-subCond.z2<-1E-7; 
}


bool Cell::isIncludedIn(Rectangle &subCond ){
	return x2-subCond.x2<-1E-7&&x1-subCond.x1>1E-7&&y2-subCond.y2<-1E-7&&y1-subCond.y1>1E-7&&
			 z2-subCond.z2<-1E-7&&z1-subCond.z1>1E-7 ;
}

Rectangle getGdsZone(ConductorList &conductorList){
	long int count=0;
	Rectangle gdsZone;
	for( ConductorList::iterator condIt=conductorList.begin(); condIt != conductorList.end(); condIt ++ ){
		for( SubConductorList::iterator subCondIt=condIt->subConductorList.begin(); subCondIt != condIt->subConductorList.end(); subCondIt++  ){
			if(count==0 ){
				gdsZone=Rectangle( *subCondIt );	
			}else{
				if( subCondIt->x1 < gdsZone.x1  ){
					gdsZone.x1=subCondIt->x1;
				}
				if( subCondIt->x2 > gdsZone.x2  ){
					gdsZone.x2=subCondIt->x2;
				}

				if( subCondIt->y1 < gdsZone.y1  ){
					gdsZone.y1=subCondIt->y1;
				}

				if( subCondIt->y2 > gdsZone.y2  ){
					gdsZone.y2=subCondIt->y2;
				}
				if( subCondIt->z1 < gdsZone.z1  ){
					gdsZone.z1=subCondIt->z1;
				}
				if( subCondIt->z2 > gdsZone.z2  ){
					gdsZone.z2=subCondIt->z2;
				}
			
			}
			count++;
		}
	}
	return gdsZone;
}




//@author:Dragon
//in generateGridOctree we should inflate the border of GDS::zone by the size of gridCellSize, for security reason
void generateGridOctree( GridOctree &gridOctree,  int gridCellSize, ConductorList  &conductorList ){
	int nX=0,nY=0,nZ=0;
	//@author:Dragon
	//set Cell::extensionSize into gridCellSize. Cell::extensionSize is actually the initial value of distanceLimit of Cell
	// Cell::extensionSize=gridCellSize;
	//update the size of GDS::zone----inflate it by the size of gridCellSize;
	GDS::zone.x1 -=gridCellSize;
	GDS::zone.x2 +=gridCellSize;
	GDS::zone.y1 -=gridCellSize;
	GDS::zone.y2 +=gridCellSize;
	GDS::zone.z1 -=gridCellSize;
	GDS::zone.z2 +=gridCellSize;

	gridOctree=partitionDomain(double(gridCellSize) , nX, nY, nZ  );

	int ix1,ix2,iy1,iy2,iz1,iz2;
	int tempI;
	for( ConductorList::iterator condIt=conductorList.begin();  condIt != conductorList.end(); condIt ++    ){
		for(  SubConductorList::iterator subCondIt= condIt->subConductorList.begin(); 
			subCondIt !=condIt->subConductorList.end(); subCondIt ++    ){
			
			// tempI=(subCondIt->x1-GDS::zone.x1)/gridCellSize ;
			// ix1=max(  ( subCondIt->x1-GDS::zone.x1)%gridCellSize==0? tempI-2:tempI-1 , 0    );

			// tempI=int(ceil (double(subCondIt->x2 -  GDS::zone.x1)/double(gridCellSize) ) );
			// ix2=min(  (subCondIt->x2 - GDS::zone.x1)%gridCellSize==0? tempI+1:tempI ,  nX -1  );

			// tempI=( subCondIt->y1 - GDS::zone.y1 )/gridCellSize;
			// iy1=max(  ( subCondIt->y1 - GDS::zone.y1 )%gridCellSize==0? tempI-2: tempI-1 , 0    );

			// tempI=int(ceil (double(subCondIt->y2 -  GDS::zone.y1)/double(gridCellSize) ) );
			// iy2=min(   (subCondIt->y2 -  GDS::zone.y1)%gridCellSize==0? tempI+1:tempI   ,  nY -1  );

			// tempI=( subCondIt->z1 - GDS::zone.z1 )/gridCellSize;
			// iz1=max(  ( subCondIt->z1 - GDS::zone.z1 )%gridCellSize==0? tempI-2: tempI-1 , 0    );

			// tempI=int(ceil (double(subCondIt->z2 -  GDS::zone.z1)/double(gridCellSize) ) );
			// iz2=min(   (subCondIt->z2 -  GDS::zone.z1)%(gridCellSize)==0? tempI+1:tempI,  nZ -1  );

//扩大检查初始的candidateList的范围，尽可能产生更大的Cube
			tempI=(subCondIt->x1-GDS::zone.x1)/gridCellSize -1;
			ix1=max(  ( subCondIt->x1-GDS::zone.x1)%gridCellSize==0? tempI-2:tempI-1 , 0    );

			tempI=int(ceil (double(subCondIt->x2 -  GDS::zone.x1)/double(gridCellSize) ) ) +1;
			ix2=min(  (subCondIt->x2 - GDS::zone.x1)%gridCellSize==0? tempI+1:tempI ,  nX -1  );

			tempI=( subCondIt->y1 - GDS::zone.y1 )/gridCellSize -1;
			iy1=max(  ( subCondIt->y1 - GDS::zone.y1 )%gridCellSize==0? tempI-2: tempI-1 , 0    );

			tempI=int(ceil (double(subCondIt->y2 -  GDS::zone.y1)/double(gridCellSize) ) ) +1;
			iy2=min(   (subCondIt->y2 -  GDS::zone.y1)%gridCellSize==0? tempI+1:tempI   ,  nY -1  );

			tempI=( subCondIt->z1 - GDS::zone.z1 )/gridCellSize -1;
			iz1=max(  ( subCondIt->z1 - GDS::zone.z1 )%gridCellSize==0? tempI-2: tempI-1 , 0    );

			tempI=int(ceil (double(subCondIt->z2 -  GDS::zone.z1)/double(gridCellSize) ) ) +1;
			iz2=min(   (subCondIt->z2 -  GDS::zone.z1)%(gridCellSize)==0? tempI+1:tempI,  nZ -1  );

			for( int i= ix1;  i<=ix2 ;i++  ){
				for( int j=iy1; j <=iy2; j++  ){
					for(  int k=iz1; k<=iz2; k++   ){
					
						candidateCheck(  *subCondIt,  gridOctree[i][j][k]  );
					}
				}
			}
		}
	}



	for( int i=0;  i< nX; i++   ){
		for( int j=0; j< nY; j++  ){
			for( int k=0; k< nZ; k++ ){
				//@author:Dragon
	//Before insert SubCond to the Cell, we first need to initiate the root Cell. In other words, we need to make a copy of the candidateList of 
	//gridOctree[i][j][k], then clear the candidateList of gridOctree[i][j][k].
	//If we do not initiate the candidateList, we may insert the same candidate twice, because a candidate cannot dominate itself.
				gridOctree[i][j][k].candidateListOfRootCell =gridOctree[i][j][k].candidateList;	
				gridOctree[i][j][k].candidateList.clear();				
				for(  SubConductorList::iterator  subCondIt=  gridOctree[i][j][k].candidateListOfRootCell.begin(); 
					subCondIt != gridOctree[i][j][k].candidateListOfRootCell.end(); subCondIt++   ){
					insertToOctree( *subCondIt,  gridOctree[i][j][k] );
				}


			}
		}
	}

	//@author:Dragon
	//改进方法：为八叉树创建索引数组（也称为2级三维数组）
	updateSubGridCellVec(gridOctree);


}





GridOctree partitionDomain( double cellSize,  int &nX, int &nY, int &nZ ){

	//inflate the original gds zone for security reason
	int x1=GDS::zone.x1;   
	int x2=GDS::zone.x2;
	int y1=GDS::zone.y1;
	int y2=GDS::zone.y2;
	int z1=GDS::zone.z1;
	int z2=GDS::zone.z2;
	
	nX=   int(ceil(double(x2-x1)/double(cellSize)));
	nY=   int( ceil(double(y2-y1)/double(cellSize)));
	nZ=   int( ceil(double(z2-z1)/double(cellSize)));

//@author:Dragon
	//This method can also work:
	//typedef std::vector< std::vector< std::vector<Cell>> >  GridOctree;
	//use nX,nYnZ to decide the size of gridOctree;
	// vector<Cell> vec1(nZ);
	// vector< vector<Cell> > vec2(nY , vec1  );   //don't forget this !
	// vector<  vector< vector<Cell> > > gridOctree(nX, vec2);
	// int i,j,k;
	// for(k=0;k<nZ;k++ ){
	// 	for( j=0;j<nY;j++ ){
	// 		for( i=0;i<nX;i++ ){

	// 			Cell cell=Cell(  x1+cellSize*i, x1+cellSize*(i+1) ,  y1+cellSize*j , 
	// 				y1+cellSize*(j+1) , z1+cellSize*k, z1+cellSize*(k+1) ,cellSize );
	// 			gridOctree[i][j][k]=cell;
	// 			//set the distance limit of each cell in gridOctree to cellSize;
	// 			
	// 		}
	// 	}
	// }
	GridOctree gridOctree;
	int i,j,k;
	for( i=0;	 i<nX; i++ ){
		vector< vector<Cell> > vecYZ;
		for( j=0; j<nY; j++  ){
			vector < Cell > vecZ;
			for( k=0; k<nZ; k++  ){
				Cell cell=Cell(  double(x1+cellSize*i), double(x1+cellSize*(i+1)) , double( y1+cellSize*j ), 
	 				double(y1+cellSize*(j+1)) , double(z1+cellSize*k), double(z1+cellSize*(k+1)) , double(cellSize) );

				vecZ.push_back(cell);

			}
			vecYZ.push_back(vecZ);
		}
		gridOctree.push_back(vecYZ);
	}



	return  gridOctree;
}

double scale_of_extensionSize ;

void  insertToOctree( SubConductor &subCond, Cell &cellT  ){
	if(cellT.size()<1E-5){
		cout<<"Error: cellT is 2-d or 1-d. This program only handles 3-dimensional space"<<endl;
		return;
	}
	//@author:Dragon
	//by seting the initial value of distanceLimit, We can generate incomplete candidatelist
	if( cellT.distanceLimit !=-1 && myDistance( subCond,cellT )- cellT.distanceLimit>1E-7 ){
		return;
	}
	//if cellT is filled with conductor, then it's useless, omit it.
	if(  cellT.isFilledWithCond ){
		return;
	}

	//if cellT is a leaf node
	if(  cellT.subCellList.size()==0  ){
		candidateCheck(  subCond, cellT );
		if( cellT.candidateList.size()>Cell::maxCandidateListLen &&   cellT.size() > Cell::minCellSize   ){

		//@author:Dragon
		//In the future we may make some improvement about the extensionSize
			cellT.divideIntoSubCell(  cellT.extensionSize );
			// cellT.divideIntoSubCell(cellT.distanceLimit*scale_of_extensionSize);
			for( SubConductorList::iterator candidateIt= cellT.candidateList.begin(); candidateIt != cellT.candidateList.end(); candidateIt++  ){
				for( SubCellList::iterator subCellIt=cellT.subCellList.begin(); subCellIt!= cellT.subCellList.end(); subCellIt++ ){
					insertToOctree( *candidateIt, *subCellIt);
				}

			}
		}


	}else{
		for( SubCellList::iterator subCellIt=cellT.subCellList.begin(); subCellIt!= cellT.subCellList.end(); subCellIt++ ){
			insertToOctree( subCond, *subCellIt);
		}
	}

}


//@author:Dragon
/*
I designed a class, SubConductor which inherit Rectangle. SubConductor contains a pointer which points the Conductor that includes this SubConductor.
Need to check if cell is totally included in a SubConductor.

*/


void candidateCheck(SubConductor &subCond,  Cell &cellT){
	//@author:Dragon
	//Here we add a condition, then we will generate incomplete candidate list.
	//Now we drop it because the imcomplete candidateList can be achieved by set initial distanceLimit
	// if(  myDistance(cellT, subCond )> Cell::extensionSize    ){    //Cell::extensionSize is a positive value.
	// 	return;
	// }

	//check if cellT is totally filled with SubConductor.
	if( cellT.isFilledWithCond ){
		cout<<"Error: try to do candidateCheck on a cellT which is filled with subConductor"<<endl;
		return;
	}
	if( cellT.isFilledChecked ==false ){

	//@author:Dragon 
	//Here is the my original method to check if the cell is in the Conductor, but it is time-consuming, so I droped it.
	//I kept it for further comparision
	//Basic idea is is cell is intersecting with a subConductor, then check if any surface is in that cell. If yes, then cell not in conductor, if not, then cell is in conductor


	//@author:Dragon
	//This is my new method, which is more efficient.
		// if( subCond.isIntersecting3d(cellT) ){	
		if( cellT.isInsectWith(subCond) ){
			if(  cellT.isIncludedIn( subCond )  ){
				cellT.isFilledWithCond=true;
			}else{
				cellT.isFilledWithCond=false;
			}
			cellT.isFilledChecked=true;
		}

	}


	double dis=myDistance( subCond, cellT );
	double cellSize=cellT.size();
	//dis==distanceLimit, still need to recondisder rect, cannot omit it simplely
	if(  cellT.distanceLimit!=-1 && dis - cellT.distanceLimit>1E-7  ) {    
		return;
	}
	//@author: should make sure that the candidateList is not empty first!! or there will be NULL pointer error!
	//@author:Dragon
	// if you conduct a delete or add operation on a list, you should break out the loop! or there will be wrong pointer!
	SubConductorList::iterator candidateRectIt;
	bool isContinue=false;
	while(true){
		isContinue=false;
		for(  candidateRectIt=cellT.candidateList.begin(); candidateRectIt != cellT.candidateList.end(); candidateRectIt++   ){
		
			if(  dominate(  *candidateRectIt , subCond, cellT )  ){
				return;
			}
			else if(  dominate( subCond, *candidateRectIt, cellT )   ){
		
				cellT.candidateList.erase( candidateRectIt );
				isContinue=true;
				break;
			}
		}
		if( isContinue){
			continue;
		}
		if( candidateRectIt == cellT.candidateList.end() ){
			break;
		}

	}
	
//sorting the candidate conductors in the ascending order of their distance to the cell
	if(cellT.candidateList.size()==0){
		cellT.candidateList.push_back( subCond );
	}else{
		double disSubCond_CellT=myDistance(subCond, cellT);
		SubConductorList::iterator candidateIt;
		for( candidateIt=cellT.candidateList.begin(); candidateIt != cellT.candidateList.end(); candidateIt++   ){
			if(   disSubCond_CellT - myDistance(  *candidateIt, cellT )<-1E-7  ){
				cellT.candidateList.insert( candidateIt, subCond  );
				break;
			}
		}
		if( candidateIt == cellT.candidateList.end() ){
			cellT.candidateList.push_back( subCond );
		}
	}

	
	if(    cellT.distanceLimit==-1 ||  (dis + cellSize ) - cellT.distanceLimit < -1E-7  ){
		cellT.distanceLimit=dis+cellSize;
	}

}

int max(int a, int b){
	return a>b? a:b;
}

double max( double a, double b  ){
	return a-b>1E-7 ?a:b;
}


int min( int a, int b  ){ 
	return a<b? a:b;
}

int min( int a1,int a2, int a3, int a4, int a5  , int a6  ){
	return min(min(min(min( min(a1,a2), a3), a4), a5),a6);

}

double min( double a1, double  a2){
	return a1-a2<-1E-7 ?a1:a2;
}

double min( double a1, double a2, double a3, double a4, double a5, double a6  ){
	return min(min(min(min( min(a1,a2), a3), a4), a5),a6);
}



//x-dis between point p and cuboid A is the larger one of xmin(A)-p.x and p.x-xmax(A)
double myDistance( FPoint  &p, Rectangle &rect ){
	double x_dis= max( rect.x1-p.x1, p.x1-rect.x2 );
	double y_dis= max(  rect.y1-p.y1, p.y1-rect.y2  );
	double z_dis=max(rect.z1-p.z1, p.z1-rect.z2 );
	double temp=max(x_dis,y_dis);
	return max(temp,z_dis);
}


double myDistance( FPoint &p, Cell &cell  ){

	double x_dis= max( cell.x1-p.x1, p.x1-cell.x2 );
	double y_dis= max(  cell.y1-p.y1, p.y1-cell.y2  );
	double z_dis=max(cell.z1-p.z1, p.z1-cell.z2 );
	double temp=max(x_dis,y_dis);
	return max(temp,z_dis);
}

double myDistance(Rectangle &rectA, Cell &cellT ){
	double x_dis=max(   rectA.x1-cellT.x2, cellT.x1-rectA.x2 );
	double y_dis=max(  rectA.y1-cellT.y2, cellT.y1-rectA.y2  );
	double z_dis=max(   rectA.z1-cellT.z2,  cellT.z1-rectA.z2 );
	double temp=max(x_dis, y_dis );
	return max( temp,z_dis );

}

//@author:Dragon
//这个函数用来计算两个subConductor之间的距离
//distance(rectA,rectB) represents the minimum distance of d( p in rectA, rectB  )
int myDistance(Rectangle &rectA, Rectangle &rectB ){
	int x_dis=max(   rectA.x1- rectB.x2, rectB.x1-rectA.x2 );
	int y_dis=max(  rectA.y1-rectB.y2, rectB.y1-rectA.y2  );
	int z_dis=max(   rectA.z1-rectB.z2,  rectB.z1-rectA.z2 );
	int temp=max(x_dis, y_dis );
	return max( temp,z_dis );
}

//rectA dominate rectB with recpect to rectT
//here we do not consider that rectT is inclued in rectA or rectB, deal with that later.
bool dominate( Rectangle &rectA ,Rectangle &rectB, Cell &cellT){

	int nVertex=8;
	int i=0;
	double disTA=0;
	double disTB=0;
	double maxDistanceTA=0;
	double minDistanceTB=0;
	vector<FPoint> vertex(nVertex);
	//counter-clock

	vertex[0]=FPoint(cellT.x1, cellT.y1, cellT.z1 );   
	vertex[1]=FPoint(cellT.x2, cellT.y1, cellT.z1 );
	vertex[2]=FPoint(cellT.x2, cellT.y2, cellT.z1 );   
	vertex[3]=FPoint(cellT.x1, cellT.y2, cellT.z1 );
	vertex[4]=FPoint(cellT.x1, cellT.y1, cellT.z2 );   
	vertex[5]=FPoint(cellT.x2, cellT.y1, cellT.z2 );
	vertex[6]=FPoint(cellT.x2, cellT.y2, cellT.z2 );   
	vertex[7]=FPoint(cellT.x1, cellT.y2, cellT.z2 );
	//if rectB intersects with rectT, then the minimum distance of d(p of rectB, rectT ) is 0, so rectB cannot be dominated
	if(cellT.isInsectWith(rectB) ){

		return false;                           
	}

	for(i=0;i<nVertex;i++){
		disTA=myDistance(vertex[i], rectA );
		disTA=disTA<0?0:disTA;
		if(disTA - maxDistanceTA>1E-7 ){
			maxDistanceTA=disTA;
		}
	}

	minDistanceTB=myDistance( rectB, cellT );
	//if maxDistanceTA < minDistanceTB , rectA will dominate rectB with respect to rectT
	return (maxDistanceTA-minDistanceTB<-1E-7);

}



void Cell::divideIntoSubCell( double extension){
	int nSubCell=8;


	double centerX= (x1+x2)/2.0;
	double centerY= (y1+y2)/2.0;
	double centerZ= (z1+z2)/2.0;

	//@author:Dragon
	//使用centerPoint来快速定位point所在的Cell
	//这里插入的顺序也要按照对应的编码进行排列
	/*
			2    3
	
			0    1
	*/
	centerPoint.setCoordinate(centerX, centerY, centerZ );
	hasChildCell=true;
	subCellList.push_back(Cell( x1, centerX, y1, centerY, z1, centerZ ,extension  ));
	subCellList.push_back( Cell( centerX, x2, y1, centerY, z1, centerZ ,extension )  );
	subCellList.push_back(Cell( x1, centerX, centerY, y2, z1, centerZ  ,extension));
	subCellList.push_back(Cell( centerX, x2, centerY, y2, z1, centerZ ,extension));
	subCellList.push_back(Cell( x1, centerX, y1, centerY, centerZ ,z2 ,extension));
	subCellList.push_back( Cell( centerX, x2, y1, centerY, centerZ ,z2 ,extension )  );
	subCellList.push_back(Cell( x1, centerX, centerY, y2, centerZ,z2  ,extension));
	subCellList.push_back(Cell( centerX, x2, centerY, y2, centerZ,z2 ,extension ));
}





void updateSubGridCellVec(Cell &rootCell,  Cell *cellP,  vector<vector<vector<Cell*> > >  &subGridCellVec ){
	if(!cellP->hasChildCell || cellP->size()-Cell::subGridSize<1E-5  ){
		int startXn=(int) nearbyint((cellP->x1-rootCell.x1)/Cell::subGridSize);
		int startYn=(int) nearbyint((cellP->y1-rootCell.y1)/Cell::subGridSize);
		int startZn=(int) nearbyint((cellP->z1-rootCell.z1)/Cell::subGridSize);
		int endXn=((int) nearbyint((cellP->x2-rootCell.x1)/Cell::subGridSize))-1 ;
		int endYn=((int) nearbyint((cellP->y2-rootCell.y1)/Cell::subGridSize))-1 ;
		int endZn=((int) nearbyint((cellP->z2-rootCell.z1)/Cell::subGridSize))-1;
		if( fabs((cellP->x1-rootCell.x1)/Cell::subGridSize - startXn) > 1E-5 ||  
			fabs((cellP->y1-rootCell.y1)/Cell::subGridSize - startYn) > 1E-5 ||  
			fabs((cellP->z1-rootCell.z1)/Cell::subGridSize - startZn) > 1E-5 ||  
			fabs((cellP->x2-rootCell.x1)/Cell::subGridSize -1 - endXn) > 1E-5 ||  
			fabs((cellP->y2-rootCell.y1)/Cell::subGridSize -1 - endYn) > 1E-5 ||  
			fabs((cellP->z2-rootCell.z1)/Cell::subGridSize -1 - endZn) > 1E-5   ){
			cout<<"Err: wrong Cell location!"<<endl;
			exit(1);
		}

		for(int i=startXn; i<= endXn; i++ ){
			for(int j=startYn; j<=endYn; j++ ){
				for(int k=startZn; k<=endZn; k++ ){
					subGridCellVec[i][j][k]=cellP;
				}
			}
		}
	}
	else{
		for(SubCellList::iterator subCellIt=cellP->subCellList.begin(); subCellIt!= cellP->subCellList.end(); subCellIt++   ){
			updateSubGridCellVec( rootCell, &(*subCellIt), subGridCellVec );
		}

	}
}


void updateSubGridCellVec(GridOctree &gridOctree){
	for(int i=0; i<gridOctree.size(); i++ ){
		for(int j=0; j<gridOctree[i].size(); j++ ){
			for(int k=0; k<gridOctree[i][j].size(); k++ ){
				int subGridNum=(int) nearbyint( GridOctree::gridCellSize/Cell::subGridSize );
				Cell *cellP=NULL;
				std::vector<Cell*> vz(subGridNum, cellP);
				std::vector<std::vector<Cell*> > vy(subGridNum, vz );
				std::vector<std::vector<std::vector<Cell*> > > cellPVec(subGridNum, vy);
				updateSubGridCellVec(gridOctree[i][j][k], &(gridOctree[i][j][k])  ,  cellPVec );
				gridOctree[i][j][k].subGridCellVec=cellPVec;
			}

		}
	}

}


