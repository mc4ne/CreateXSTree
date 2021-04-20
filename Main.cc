#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <thread>
#include <vector>
#include <chrono>


//#define RATES_DEBUG 1

using namespace std;

#include "HMSXSTree.h"
#include "SHMSXSTree.h"
#include "XSTree.h"
#include "ExtractAcceptance.h"
#include "ACCTools.h"

extern double A1NOptics();
extern double* GetOpticsRate(double pBeamCurrent, double pBeamE, double pDetectorAngle, double pDetectorMomentum, string pDetectorName, int pElasOnly=0);
extern double A1NRates();
//Method 1 
extern double GetRate(double pBeamCurrent, double pBeamE, double pDetectorAngle, double pDetectorMomentum, string pDetectorName, int pElasOnly=0);
//Method 2
extern void GetRate(double pBeamCurrent_uA, double pBeamE_GeV, double pDetectorAngle_deg, double pDetectorMomentum_GeV, int Det, int pFullAcc=1);

int getopticsrate_main(int argc, char** argv)
{
  if(argc<6) {
  cout<<" Error: you need to provide at least 5 arguments!\n"
      <<" Calculate rates using method 1.\n"
      <<" Usage: "<<argv[0]<<" <BeamCurrent_uA> <Beam_GeV> <DetectorAngle_deg> <DetectorMomentum_GeV> <DetectorName=HMS|SHMS> [ElasOnly=0]\n"
      <<"        All energies are in GeV unit. All angles are in degree unit.\n"
      <<"        ElasOnly=-1: pure inelastic for full acceptance.\n"
      <<"        ElasOnly=0:  inelastic + elastic for full acceptance.\n"
      <<"        ElasOnly=1:  pure elastic for full acceptance.\n"
      <<"        ElasOnly=-30: pure inelastic for 2-SC-Bar acceptance.\n"
      <<"        ElasOnly=30: inelastic + elastic for 2-SC-Bar acceptance.\n"
      <<"        ElasOnly=31: pure elastic for 2-SC-Bar acceptance.\n"
      <<endl;
    exit(-1);
  }
  const double degree = asin(1.0)/90.0;
  double pBeamCurrent, pBeamE, pDetectorAngle, pDetectorMomentum;
  string pDetectorName;

  pBeamCurrent = atof(argv[1]); 
  pBeamE = atof(argv[2]); 
  pDetectorAngle = atof(argv[3])*degree; 
  pDetectorMomentum = atof(argv[4]); 
  pDetectorName = argv[5];
  int pElasOnly = 0;
  if(argc>6) pElasOnly = atol(argv[6]); 
  
  cout<<"  Beam="<<pBeamE<<"  DetAngle="<<pDetectorAngle/degree<<"  pDetectorMomentum="<<pDetectorMomentum<<endl;

  //extern double* GetOpticsRate(double pBeamCurrent, double pBeamE, double pDetectorAngle, double pDetectorMomentum, string pDetectorName, int pElasOnly=0);
  GetOpticsRate(pBeamCurrent, pBeamE, pDetectorAngle, pDetectorMomentum, pDetectorName, pElasOnly);
  return 0;
}

int getrate1_main(int argc, char** argv)
{
  if(argc<6) {
  cout<<" Error: you need to provide at least 5 arguments!\n"
      <<" Calculate rates using method 1.\n"
      <<" Usage: "<<argv[0]<<" <BeamCurrent_uA> <Beam_GeV> <DetectorAngle_deg> <DetectorMomentum_GeV> <DetectorName=HMS|SHMS> [ElasOnly=0]\n"
      <<"        All energies are in GeV unit. All angles are in degree unit.\n"
      <<"        ElasOnly=-1: pure inelastic for full acceptance.\n"
      <<"        ElasOnly=0:  inelastic + elastic for full acceptance.\n"
      <<"        ElasOnly=1:  pure elastic for full acceptance.\n"
      <<"        ElasOnly=2:  inelastic + elastic for full acceptance, with cut of 1.10<W<1.35.\n"
      <<"        ElasOnly=4:  inelastic + elastic for full acceptance. with cut of 2.00<W<100.0\n"
      <<"        ElasOnly=-30: pure inelastic for 2-SC-Bar acceptance.\n"
      <<"        ElasOnly=30: inelastic + elastic for 2-SC-Bar acceptance.\n"
      <<"        ElasOnly=31: pure elastic for 2-SC-Bar acceptance.\n"
      <<endl;
    exit(-1);
  }
  const double degree = asin(1.0)/90.0;
  double pBeamCurrent, pBeamE, pDetectorAngle, pDetectorMomentum;
  string pDetectorName;

  pBeamCurrent = atof(argv[1]); 
  pBeamE = atof(argv[2]); 
  pDetectorAngle = atof(argv[3])*degree; 
  pDetectorMomentum = atof(argv[4]); 
  pDetectorName = argv[5];
  int pElasOnly = 0;
  if(argc>6) pElasOnly = atol(argv[6]); 
  
  cout<<"  Beam="<<pBeamE<<"  DetAngle="<<pDetectorAngle/degree<<"  pDetectorMomentum="<<pDetectorMomentum<<endl;

  GetRate(pBeamCurrent, pBeamE, pDetectorAngle, pDetectorMomentum, pDetectorName, pElasOnly);
  return 0;
}

int getrate2_main(int argc, char** argv)
{
  if(argc<8) {
  cout<<" Error: you need to provide 7 arguments!\n"
      <<" Calculate rates using method 2.\n"
      <<" Usage: "<<argv[0]<<" <BeamCurrent_uA> <Beam_GeV> <DetectorAngle_deg> <DetectorMomentum_GeV> <Detector=1 HMS|2 SHMS> <FullAcceptance=0|1> <xstree_file>\n"
      <<"        All energies are in GeV unit. All angles are in degree unit.\n"
      <<"        FullAcceptance==0:  will use full acceptance, otherwise only for 2-SC-Bar.\n"
      <<endl;
    exit(-1);
  }
  const double degree = asin(1.0)/90.0;
  double pBeamCurrent, pBeamE, pDetectorAngle, pDetectorMomentum;
  int pDetector, pFullAcceptance;

  pBeamCurrent = atof(argv[1]); 
  pBeamE = atof(argv[2]); 
  pDetectorAngle = atof(argv[3])*degree; 
  pDetectorMomentum = atof(argv[4]); 
  pDetector = atol(argv[5]);
  pFullAcceptance = atol(argv[6]);
  string xstree_file = argv[7];
  
  TFile *file = TFile::Open(xstree_file.c_str());
  
  cout<<"\n Calculate rates for \""<<xstree_file<<"\" using method 2.\n";
  cout<<"  Beam="<<pBeamE<<"  DetAngle="<<pDetectorAngle/degree<<"  pDetectorMomentum="<<pDetectorMomentum<<endl;

  GetRate(pBeamCurrent, pBeamE, pDetectorAngle, pDetectorMomentum, pDetector, pFullAcceptance);
  file->Close();
  
  return 0;
}

int a1nrate_main(int argc, char** argv)
{
  if(argc<2) {
    cout<<" Usage: "<<argv[0]<<" [no_argument_needed]\n"
        <<endl;
  }
  A1NRates();
  return 0;
}

int a1noptics_main(int argc, char** argv)
{
  if(argc<2) {
    cout<<" Usage: "<<argv[0]<<" [no_argument_needed]\n"
        <<endl;
  }
  A1NOptics();
  return 0;
}

int xstree_main(int argc, char** argv)
{
  if(argc<6) {
    cout<<" Error: you need to provide 6 arguments!\n"
        <<" Usage: "<<argv[0]<<" <pElasOnly> <pBeam_GeV> <pDetectorAngle_deg> <pDetectorMomentum_GeV> <pDetector=1,10 HMS|2,20 SHMS> <rootfile>\n"
        <<"        if pDetector==10 or 20, will use XStree to do the job \n"
        <<"        if pElasOnly!=0, will fill the ntuple using elastic XS other than Bosted inelastic XS\n"
        <<"        All energies are in GeV unit. All angles are in degree unit."
        <<endl;
    exit(-1);
  }
  const double degree = asin(1.0)/90.0;
  int pElasOnly=0;
  double pBeamE, pDetectorAngle, pDetectorMomentum;
  int pDetector=1;
  TString infile;

  pElasOnly = atol(argv[1]); 
  pBeamE = atof(argv[2]); 
  pDetectorAngle = atof(argv[3])*degree; 
  pDetectorMomentum = atof(argv[4]); 
  pDetector = atol(argv[5]);
  infile = argv[6];
  
  cout<<"  Beam="<<pBeamE<<"  DetAngle="<<pDetectorAngle/degree<<"  pDetectorMomentum="<<pDetectorMomentum<<endl;

  if(pDetector==1) {
    HMSXSTree* pHMS = new HMSXSTree(infile.Data());
    pHMS->SetPara(pBeamE, pDetectorAngle, pDetectorMomentum);
    pHMS->SetElas(pElasOnly);
    pHMS->Run();
    delete pHMS;
  } else if(pDetector==2) {
    SHMSXSTree* pSHMS = new SHMSXSTree(infile.Data());
    pSHMS->SetPara(pBeamE, pDetectorAngle, pDetectorMomentum);
    pSHMS->SetElas(pElasOnly);
    pSHMS->Run();
    delete pSHMS;
  } else if(pDetector==10 || pDetector==20) {
    XSTree* pXSTree = new XSTree(infile.Data(),pDetector/10);
    pXSTree->SetPara(pBeamE, pDetectorAngle, pDetectorMomentum);
    pXSTree->SetElas(pElasOnly);
    pXSTree->Run();
    delete pXSTree;
  }
  return 0;
}


int extractacc_main(int argc, char** argv)
{
  if(argc<4) {
    cout<<" Error: you need to provide 3 arguments!\n"
        <<" Usage: "<<argv[0]<<"  <pDetector=1 HMS|2 SHMS> <pType=1|2> <infile> [nthread=1]\n"
        <<"        pType=2 means only 2 SC bars are turned on\n"
        <<"        infile is a txt file that list the absolute path of all source root files \n"
        <<"        nthread is number of threads you want to run with. \n"
        <<endl;
    exit(-1);
  }
  
  int pDet = atol(argv[1]);
  int pType = atol(argv[2]);
  TString infile = argv[3];
  int nthread = 1;
  if(argc>4) nthread = atol(argv[4]); 
      
  clock_t start;
  double cpu_time_used;

  start = clock();

  std::vector<ExtractAcceptance*> pAcc;
  for(int i=0;i<nthread;i++) {
    pAcc.push_back(new ExtractAcceptance(infile.Data(),pDet,pType,i,nthread));
  }
  
 //construct thread, it will run immediately
  std::vector<std::thread> threadObj;
  for (int i = 0; i < nthread; i++) {
    threadObj.push_back(std::thread(&ExtractAcceptance::Run, pAcc[i]));
    std::this_thread::sleep_for(std::chrono::milliseconds(200));
  }

  //synchronize threads:
  for(int i=0;i<nthread;i++) {
    threadObj[i].join();
  }

  //now merge all results  
  for(int i=1;i<nthread;i++) {
    pAcc[0]->MergeResult(pAcc[1]);
  }
  
  cpu_time_used = ((double) (clock() - start)) / CLOCKS_PER_SEC;
  cout<<"All threads are done. Time elapsed "<<cpu_time_used<<endl;
  
  //now create output file
  pAcc[0]->EndOfRun();
  
  return 0;
}

int main(int argc, char** argv)
{
  if(argc<2) {
    cout<<" Error: you need to provide at least 2 arguments!\n"
        <<" Usage: "<<argv[0]<<"  <task=a1nrate|extractacc|xstree|getrate> <other_arguments>\n"
        <<"       task==a1nrate:    get all A1N rates for all kinematic points \n"
        <<"       task==extractacc: extract HMS or SHMS acceptance \n"
        <<"       task==xstree:     add xs branches into mc-signle-arm output ntuple \n"
        <<"       task==getrate1:   get rate for one single kinematics point using method 1\n"
        <<"       task==getrate2:   get rate for given xstree_file using method 2\n"
        <<"       task==a1noptics:  get all A1N optics rates for all kinematic points \n"
        <<"       task==getopticsrate: get rates for each C12 foil for given kinematic\n"
        <<endl
        <<"       For details of 'other_arguments' in each task, type 'help' after that task \n"
        <<"       For example: "<<argv[0]<<" extractacc help \n"
        <<endl;
    exit(-1);
  }

  string task=argv[1];
  if(task == "a1nrate") a1nrate_main(argc-1,&argv[1]);
  else if(task == "extractacc") extractacc_main(argc-1,&argv[1]);
  else if(task == "xstree") xstree_main(argc-1,&argv[1]);
  else if(task == "getrate1") getrate1_main(argc-1,&argv[1]);
  else if(task == "getrate2") getrate2_main(argc-1,&argv[1]);
  else if(task == "a1noptics") a1noptics_main(argc-1,&argv[1]);
  else if(task == "getopticsrate") getopticsrate_main(argc-1,&argv[1]);
  else {
    cout<<" this given task '"<<task<<"' is not curruetly support ...\n";
  }
  
  return 0;
}
