#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <pthread.h>
#include <vector>
using namespace std;

int numOfGrid=2000;
int NTHREAD=40;
vector<string> outNameList(4);

double unitSurfaceGreenFunction(  double x, double y,  int nX, int nY  );
void *generateGreenFuncTable(void * ind);
void copy(ofstream &fout, string &filename );
double unitGEx( double x, double y ,int nX, int nY );
double unitGEy( double x, double y ,int nX, int nY );
double unitGEz( double x, double y ,int nX, int nY );

int main(int argc, char *argv[] ){

	numOfGrid=atoi(argv[1]);

	outNameList[0]="Green";
	outNameList[1]="GEx";
	outNameList[2]="GEy";
	outNameList[3]="GEz";

	srand(unsigned(time(NULL) )  );
	int arr_index[NTHREAD];
	pthread_t tids[NTHREAD];
	pthread_attr_t attr;
	pthread_attr_init(&attr );
	pthread_attr_setdetachstate( &attr, PTHREAD_CREATE_JOINABLE );

	for( int i=0; i<NTHREAD; i++ ){
		arr_index[i]=i;
	}	

	for( int i=0; i<NTHREAD; i++ ){

		if( pthread_create(&tids[i], &attr,generateGreenFuncTable , (void*)&(arr_index[i]) ) !=0    ){
			cout<< "ERR!!"<<endl;
			return 1;
		}
	}
	pthread_attr_destroy( &attr );
	void *status;
	for( int i=0; i<NTHREAD; i++ ){
		if(    pthread_join(  tids[i], &status   )  !=0 ){
			cout<<"ETTR"<<endl;
			return 1;
		}

	}

	for(int i=0;i<NTHREAD/10; i++){
		string finalOutName= outNameList[i]+".txt";
		ofstream fout( finalOutName.c_str() );
		for(int j=0; j<10; j++){
			char fileIndex[256];
			sprintf(fileIndex,"%d", j );
			string infileName(outNameList[i]);
			infileName+=fileIndex;
			infileName+=".txt";
			copy(fout, infileName);
			remove(infileName.c_str() );
		}
		fout.close();
	}

}

void *generateGreenFuncTable(void * ind){
	int index= *((int *)ind);
	int localFileIndex= index %10;
	string outName=outNameList[index/10];
	char fileIndex[256];
	sprintf(fileIndex,"%d", localFileIndex );
	string outfileName(outName);
	outfileName+=fileIndex;
	outfileName+=".txt";

	ofstream fout(outfileName.c_str());
	if(!fout.is_open() ){
		cout<<"file open failed!"<<endl;
		exit(1);
	}
	//预先计算格林函数矩阵， 将0-1的区间平均分成2000*2000的矩阵，在每一个小格子里随机取10个点算出格林函数求平均作为该格子范围内的
	//格林函数的估计值
	
	//@author:Dragon
	//注意矩阵的格式
	/*	同一行的元素所对应的y坐标相同
	
	坐标系如下，所以如果将其一行一行的装进2维数组里时， a[i][j]， i表示第i行，对应的是y的坐标， j对应的是x的坐标！！
		------------> x+
		|
		|
		|
		V
		 y+

	*/

	double x,y;

	for(int j= localFileIndex *( numOfGrid/10 ); j<(localFileIndex+1)*(numOfGrid/10) ; j++){
		for(int i=0; i<numOfGrid; i++ ){
			int k=0;
			double green=0;
			while(k++<10){
			 	x=1.0/double(numOfGrid)*i + (rand()/(double(RAND_MAX)))/numOfGrid;
				y=1.0/double(numOfGrid)*j+ (rand()/(double(RAND_MAX)))/numOfGrid;

				if(index/10 ==0){
					green+=unitSurfaceGreenFunction(x,y,10,10 );
				}else if(index/10 ==1 ){
					green+=unitGEx(x,y,10,10);
				}else if(index/10 ==2 ){
					green+=unitGEy(x,y,10,10);
				}else if(index/10==3 ){
					green+=unitGEz(x,y,10,10 );
				}else{
					cout<<"Wrong index!"<<endl;
					exit(1);
				}
			}
			green=green/10.0;
			fout<<setprecision(9)<<green<<"\t";
		}
		fout<<"\n";
		cout<<"\r"<<j<<"/"<<numOfGrid<<endl;
	}
	fout.close();
	pthread_exit(0);
}


double unitSurfaceGreenFunction(  double x, double y,  int nX, int nY  ){
	int nx, ny;
	double nz;
	double greenFuncValue=0;
	for( nx=1; nx<=nX; nx++    ){
		for( ny=1; ny<=nY; ny++  ){
			nz=pow(   pow( double(nx), 2.0 )+pow(double(ny), 2), 0.5 );
			greenFuncValue+=2.0*sin(M_PI*nx/2.0)*sin(M_PI*ny/2.0)*sin(M_PI*nx*x)*sin(M_PI*ny*y)/ cosh(M_PI* nz/2.0 );

		}
	}
	return greenFuncValue;
}

double unitGEx( double x, double y ,int nX, int nY ){
	int nx, ny;
	double nz;
	double GE=0;
	for(nx=1;nx<=nX; nx++  ){
		for( ny=1; ny<=nY; ny++ ){
			nz= pow( pow(double(nx),2.0)+pow( double(ny),2.0 ), 0.5      );
			GE+= (-1.0)*2*M_PI * nx * cos(M_PI*nx/2.0)*sin(M_PI*ny/2.0)/cosh(M_PI*nz/2.0) *sin(M_PI*nx*x)*sin(M_PI*ny*y);
		}

	}
	return GE;
}

double unitGEy(double x, double y ,int nX, int nY ){
	int nx, ny;
	double nz;
	double GE=0;
	for(nx=1;nx<=nX; nx++  ){
		for( ny=1; ny<=nY; ny++ ){
			nz= pow( pow(double(nx),2.0)+pow( double(ny),2.0 ), 0.5      );
			GE+= (-1.0)*2*M_PI * ny * sin(M_PI*nx/2.0)*cos(M_PI*ny/2.0)/cosh(M_PI*nz/2.0) *sin(M_PI*nx*x)*sin(M_PI*ny*y);
		}

	}
	return GE;
}

double unitGEz(  double x, double y ,int nX, int nY ){
	int nx, ny;
	double nz;
	double GE=0;
	for(nx=1;nx<=nX; nx++  ){
		for( ny=1; ny<=nY; ny++ ){
			nz= pow( pow(double(nx),2.0)+pow( double(ny),2.0 ), 0.5      );
			GE+= (-1.0)* 2*M_PI * nz * sin(M_PI*nx/2.0)*sin(M_PI*ny/2.0)/sinh(M_PI*nz/2.0) *sin(M_PI*nx*x)*sin(M_PI*ny*y);
		}

	}		
	return GE;
}


void copy(ofstream &fout, string &filename ){
	ifstream fin(filename.c_str());
	string str;
	while(getline(fin,str)){
		fout<<str<<"\n"<<flush;
	}
	fin.close();
}
