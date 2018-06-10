#include "FRWcapext.cpp"

using namespace std;


int main(int argc, char* argv[]  ){

	// int condIndex=atoi(argv[1]);
	// char* condIndexStr=argv[1];
	int NTHREAD= atoi(argv[2]) ;
	FRWControl::maxPermittedCapErr=atof(argv[3]); //最大容许的电容误差百分比
	long loopNum= atol(argv[4]);
	stringstream ss;
	ss<<argv[5];
	string outFileName=ss.str();
	ss.str("");
	ss<<argv[1];
	string geoFileName=ss.str();

gettimeofday(&tv_begin,NULL );

	GeoLoader* geoLoader=new GeoLoader();	
	geoLoader->loadGeo(geoFileName.c_str());
	// ConductorList  conductorList;
	geoLoader->generateConductorList(FRWControl::conductorList , true);
	GDS::zone=getGdsZone(FRWControl::conductorList);

gettimeofday(&tv_end,NULL);
cout<<"ConductorList generated.----- Time consumption: "<<( (tv_end.tv_sec*1000000+tv_end.tv_usec)-(tv_begin.tv_sec*1000000+tv_begin.tv_usec))/1000<<" ms" <<endl;
	

gettimeofday(&tv_begin,NULL );
	loadDataTable("DataTables/Green.txt", FRWControl::GreenVT ,FRWControl::GreenVTSize );
	loadDataTable("DataTables/GEx.txt", FRWControl::GExVT , FRWControl::GExVTSize );
	loadDataTable("DataTables/GEy.txt", FRWControl::GEyVT , FRWControl::GEyVTSize );
	loadDataTable("DataTables/GEz.txt", FRWControl::GEzVT , FRWControl::GEzVTSize );
gettimeofday(&tv_end,NULL);
cout<<"Data table loaded. ----- Time consumption: "<<( (tv_end.tv_sec*1000000+tv_end.tv_usec)-(tv_begin.tv_sec*1000000+tv_begin.tv_usec))/1000<<" ms" <<endl;

gettimeofday(&tv_begin,NULL );
	generateGridOctree(FRWControl::gridOctree, GridOctree::gridCellSize , FRWControl::conductorList );
gettimeofday(&tv_end,NULL);
cout<<"GridOctree generated.----- Time consumption: "<<( (tv_end.tv_sec*1000000+tv_end.tv_usec)-(tv_begin.tv_sec*1000000+tv_begin.tv_usec))/1000<<" ms" <<endl;

gettimeofday(&tv_begin,NULL );
	double scale_factor=2;
	generateGaussianSurfaceOfConductor(FRWControl::conductorList, scale_factor);
gettimeofday(&tv_end,NULL);
cout<<"VGS generated. ----- Time consumption: "<<( (tv_end.tv_sec*1000000+tv_end.tv_usec)-(tv_begin.tv_sec*1000000+tv_begin.tv_usec))/1000<<" ms" <<endl;
	
	FRWControl::capMatrix= generateCapMatrix( FRWControl::conductorList );
	FRWControl::totalNumOfWalkOfCond = vector<long>(FRWControl::conductorList.size(), loopNum );
	FRWControl::currentNumOfWalkOfCond =vector<long>(FRWControl::conductorList.size(),0 );
	FRWControl::currentProgressOfCond=vector<long>(FRWControl::conductorList.size(),0);
	cout<<"Global variable initialized!"<<endl;


	// //@author:Dragon
	// //获得总的处理器数
	// FRWControl::totalProcessorNum=sysconf(_SC_NPROCESSORS_CONF);

	srand( unsigned (time(NULL)) );
	FRWControl::isRandSeedSet=true;

    //@author:Dragon
    //generate a random seed array
    //利用boost Random库为每一个线程生成一组随机数生成器，randSeedVec为其提供初始化种子
	int maxNumOfThread=500;
	for(int r=0;r<maxNumOfThread; r++ ){
		long rv=rand();
		bool isSame=false;
		for(int k=0; k<FRWControl::randSeedVec.size(); k++){
			if( FRWControl::randSeedVec[k]==rv ){
				isSame=true;
			}
		}
		if(!isSame){
			FRWControl::randSeedVec.push_back(rv);
		}
	}


	cout<<"FRW walk started."<<endl;
	//start counting
gettimeofday(&tv_begin,NULL );

	pthread_t tids[NTHREAD];
	pthread_attr_t attr;
	pthread_attr_init(&attr );
	pthread_attr_setdetachstate( &attr, PTHREAD_CREATE_JOINABLE );

	pthread_mutex_init( &getRandSeed_mutex, NULL );
	pthread_mutex_init( &updateFRWData_mutex, NULL );
	pthread_mutex_init(&processorIndex_mutex, NULL);

	for( int i=0; i<NTHREAD; i++ ){
		if(    pthread_create(&tids[i], &attr, FRW, NULL ) !=0         ){
			cout<< " pthread_create ERROR!!"<<endl;
			exit(1);
		}
		//@author:Dragon
		//使线程错开，尽量避免同时访问锁
		// sleep(1);
	}
	pthread_attr_destroy( &attr );
	void *status;
	for( int i=0; i<NTHREAD; i++ ){
		if(    pthread_join(  tids[i], &status   )  !=0 ){
			cout<<"pthread_join  ERROR!!"<<endl;
			exit(1);
		}
	}

	pthread_mutex_destroy( &getRandSeed_mutex );
	pthread_mutex_destroy( &updateFRWData_mutex);
	pthread_mutex_destroy(&processorIndex_mutex);

	//end counting
gettimeofday(&tv_end,NULL);
cout<<"\nFRW walk ended. ----- Time consumption: "<<(tv_end.tv_sec-tv_begin.tv_sec)<<" s" <<endl;

	cout<<"locateCellTime: "<<locateCellTime/1000000<<endl;
	cout<<"generateMaxCubeTime: "<<generateMaxCubeTime/1000000<<endl;
	cout<<"generatePointTime: "<<generatePointTime/1000000<<endl;

	cout<<"hitTime: "<<hitTime<<endl;
	cout<<"totalStepNum: "<<totalStepNum<<endl;
	cout<<"hit rate: "<<double(hitTime)/totalStepNum<<endl;


	int temp_index=0;
	for(ConductorList::iterator cond=FRWControl::conductorList.begin(); cond!=FRWControl::conductorList.end(); cond++ ){
		for(int column=0; column<FRWControl::capMatrix[temp_index].size(); column++ ){
			FRWControl::capMatrix[temp_index][column].capacitance=FRWControl::capMatrix[temp_index][column].sumOfCap*cond->gaussianSurfaceList.VGSArea ;
			FRWControl::capMatrix[temp_index][column].capacitance/=FRWControl::currentNumOfWalkOfCond[temp_index];
		}
		temp_index++;
	}

	ofstream capFout(outFileName.c_str());

	for(int i=0; i<FRWControl::capMatrix.size();i++ ){
		for(int j=0; j<FRWControl::capMatrix[i].size(); j++ ){
			cout<<FRWControl::capMatrix[i][j].capacitance;
			capFout<<FRWControl::capMatrix[i][j].capacitance;

			if( j<FRWControl::capMatrix[i].size()-1 ){
				cout<<"\t"<<flush;
				capFout<<"\t"<<flush;
			}
		}
		cout<<endl;
		capFout<<"\r\n";
	}
	capFout.close();

	return 0;
	
}

