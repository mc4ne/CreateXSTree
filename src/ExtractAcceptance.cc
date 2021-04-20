//implement of  ExtractAcceptance class-
#include <iostream>
#include <fstream>
#include <string>

#include "ExtractAcceptance.h"
#include "ACCInc.h"
#include "ACCTools.h"
#include <math.h>

using namespace std;

ExtractAcceptance::ExtractAcceptance(const char* filename, int det, int type, int tid, int nthread):
  mThreadID(tid),mTotalThread(nthread),mDet(det),mType(type)
{
#ifdef ExtractAcceptance_Debug 
  if(ExtractAcceptance_Debug>=4) 
  {
    cout<<">>>>>ExtractAcceptance::ExtractAcceptance() <<<<<"<<endl;
  }
#endif
   
  mDetName = (mDet==1)?"HMS":"SHMS";
  if(mTotalThread==1)  mThreadID=0;
  
  int bin1=ACCEPTANCE::kYtarBinNum; 
  int bin2=ACCEPTANCE::kDeltaBinNum; 
  int bin3=ACCEPTANCE::kThetaBinNum; 
  int bin4=ACCEPTANCE::kPhiBinNum;
  int iDefVal=0;
  double dDefVal=0;

  N_Inc   =ACCEPTANCE::IniDynamicArray(bin1,bin2,bin3,bin4,iDefVal);
  N_Inc_true=ACCEPTANCE::IniDynamicArray(bin1,bin2,bin3,bin4,iDefVal);
  Acc_Inc   =ACCEPTANCE::IniDynamicArray(bin1,bin2,bin3,bin4,dDefVal);

  CreateFileList(filename);
}


ExtractAcceptance::~ExtractAcceptance()
{   
  int bin1=ACCEPTANCE::kYtarBinNum; 
  int bin2=ACCEPTANCE::kDeltaBinNum; 
  int bin3=ACCEPTANCE::kThetaBinNum; 

  ACCEPTANCE::FreeDynArray(N_Inc,bin1,bin2,bin3);
  ACCEPTANCE::FreeDynArray(N_Inc_true,bin1,bin2,bin3);
  ACCEPTANCE::FreeDynArray(Acc_Inc,bin1,bin2,bin3);
}

//read the content of the given file and store it into the vector
void ExtractAcceptance::CreateFileList(const char* filename)
{
#ifdef ExtractAcceptance_Debug 
  if(ExtractAcceptance_Debug>=4) cout<<"ExtractAcceptance::CreateFileList()"<<endl;
#endif

  if(ACCEPTANCE::AccessFilePath(filename) == 0) {
    cout<<"ExtractAcceptance::CreateFileList(): inputfile \""<<filename<<"\" does not exist, I quit...\n";
    exit(-2);
  }
  
  std::ifstream ifs(filename);
  
  string tmpStr0;
  string tmpStr;
  while (!ifs.eof()) {
    tmpStr0.clear();   //empty the string
    ifs>>tmpStr0;

    //remove leading and trailing spaces
    tmpStr=ACCEPTANCE::trim(tmpStr0);

    //ignore this line if it starts with #
    size_t found = tmpStr.find_first_of("#");
    if (found != std::string::npos) {
      //found it, truncate the string
      if(found == 0)  continue;
      else tmpStr = tmpStr.substr(0, found);
    }
    if(tmpStr.length() < 1)  continue;

    if(ACCEPTANCE::AccessFilePath(tmpStr.c_str())) {
      mVFileList.push_back(tmpStr);
    }
  }
  ifs.close();

  if(mVFileList.size()<1) {
    cout<<"ExtractAcceptance::CreateFileList(): no valid input root files, I quit...\n";
    exit(-2);    
  }
  
  //debug: show the list
  cout<<"Here is the source root file list: \n";
  for(int i=0;i<mVFileList.size();i++) {
    cout<<mVFileList[i]<<endl;
  }
  return;
}

void ExtractAcceptance::BeginOfRun()
{
  // Declaration of histograms for Beta Check
#ifdef ExtractAcceptance_Debug 
  if(ExtractAcceptance_Debug>=4) cout<<"ExtractAcceptance::BeginOfRun()"<<endl;
#endif
}

/////////////////////////////////////////////////////////////////////////////////////
void ExtractAcceptance::EndOfRun()
{
#ifdef ExtractAcceptance_Debug 
  if(ExtractAcceptance_Debug>=4) cout<<"ExtractAcceptance::UserEndOfRun() "<<endl; 
#endif
   
  DoAccCal();
  FillAccTree();
}

/////////////////////////////////////////////////////////////////////////////////////
void ExtractAcceptance::Run()
{
  BeginOfRun();
  
#ifdef ExtractAcceptance_Debug 
  if(ExtractAcceptance_Debug>=4) cout<<" ExtractAcceptance::Run() "<<endl;
#endif

  for(int i=mThreadID;i<mVFileList.size();i+=mTotalThread) {
    
    cout<<" Thread \t"<<mThreadID+1<<": ExtractAcceptance::Run()\n\t";
    LoadAFile(mVFileList[i].c_str(),mDet);
    
    if(!fChain) {
      cout<<"Warning! Something wrong in ReadSingleArm::LoadAFile(), fChain is NULL...\n";
      continue; //exit(-1);
    }
    
    Long64_t nentries = fChain->GetEntriesFast();
#ifdef ExtractAcceptance_Debug 
    if(ExtractAcceptance_Debug>=6) nentries=1000;
#endif

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
#ifdef ExtractAcceptance_Debug 
      if(ExtractAcceptance_Debug>=5)
        cout<<" ExtractAcceptance::Run() is processing event "<<std::setw(6)<<jentry+1<<"\n";
#endif

#ifdef ExtractAcceptance_Debug 
      if( ((jentry+1)%10000) == 0)
        cout<<" ExtractAcceptance::Run() is processing event "<<std::setw(6)<<jentry+1<<" ... \r";
#endif
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      fChain->GetEntry(jentry);

      //apply cuts
      //if(stop_id && stop_id<100) continue;
        
      if(ACCEPTANCE::GetIndex(psytari,psdeltai,psxptari,psyptari,iYtarIdx,iDeltaIdx,iThetaIdx,iPhiIdx)>=0) 
      {
        N_Inc_true[iYtarIdx][iDeltaIdx][iThetaIdx][iPhiIdx]++;
        if(stop_id == 0) {
          //apply your extra cuts here ...
          if(mType==1) N_Inc[iYtarIdx][iDeltaIdx][iThetaIdx][iPhiIdx]++;
          else {
            if( mType==2 && fabs(psxfp)<8.0) N_Inc[iYtarIdx][iDeltaIdx][iThetaIdx][iPhiIdx]++;
          }
        }
      }
      
    }

    //empty the memory for next file
    delete fChain->GetCurrentFile();
    
    cout<<endl;
  }

  //EndOfRun();  //since this is multiple thread mode, user need to merge result before calling end of run
  return;
}


/////////////////////////////////////////////////////////////////////////////////////  
void ExtractAcceptance::PrintACCTable(int level, const char *filename, char *headblock, int ****N_det, int ****N_true)
{ 
  ACCEPTANCE::PrintACCTable(ACCEPTANCE::kYtarBin, ACCEPTANCE::kYtarBinNum, 
    ACCEPTANCE::kDeltaBin, ACCEPTANCE::kDeltaBinNum,
    ACCEPTANCE::kThetaBin, ACCEPTANCE::kThetaBinNum,
    ACCEPTANCE::kPhiBin, ACCEPTANCE::kPhiBinNum, 
    level, filename, headblock,N_det,N_true);  
}

/////////////////////////////////////////////////////////////////////////////////////  
void ExtractAcceptance::PrintACCTable(int level, const char *filename, char *headblock, double ****N_det, double ****N_true)
{ 
  ACCEPTANCE::PrintACCTable(ACCEPTANCE::kYtarBin, ACCEPTANCE::kYtarBinNum, 
    ACCEPTANCE::kDeltaBin, ACCEPTANCE::kDeltaBinNum,
    ACCEPTANCE::kThetaBin, ACCEPTANCE::kThetaBinNum,
    ACCEPTANCE::kPhiBin, ACCEPTANCE::kPhiBinNum, 
    level, filename, headblock,N_det,N_true);  
}
/////////////////////////////////////////////////////////////////////////////////////  
void ExtractAcceptance::DoAccCal()
{
#ifdef ExtractAcceptance_Debug 
  if(ExtractAcceptance_Debug>=4) cout<<"ExtractAcceptance::DoAccCal() "<<endl; 
#endif
  int bin1=ACCEPTANCE::kYtarBinNum; 
  int bin2=ACCEPTANCE::kDeltaBinNum; 
  int bin3=ACCEPTANCE::kThetaBinNum; 
  int bin4=ACCEPTANCE::kPhiBinNum;

  char filename[100],str[512];
  int iw,iq,it,ip;

  char tablename[100];
  char *name_exc[]={"ACCEPTANCE_HMS_1","kYtarBinNum","kDeltaNum","kThetaBinNum","kPhiBinNum"};

  //now is time to get the acc table
  for(iw=0;iw<bin1;iw++)
  {  
    for(iq=0;iq<bin2;iq++)
    {  
      for(it=0;it<bin3;it++)
      {
        for( ip=0;ip<bin4;ip++)
        {
          if(N_Inc_true[iw][iq][it][ip]!=0)
          {
            Acc_Inc[iw][iq][it][ip]=(double)(N_Inc[iw][iq][it][ip])/(double)(N_Inc_true[iw][iq][it][ip]);
          }

#ifdef ExtractAcceptance_Debug 
          if(ExtractAcceptance_Debug>=5) 
          {  
            if(N_Inc[iw][iq][it][ip]!=0) 
            cout<<"N_Inc_true= "<< N_Inc_true[iw][iq][it][ip]
                <<"  N_Inc = "<< N_Inc[iw][iq][it][ip]
                <<"  Acc_Inc = "<< Acc_Inc[iw][iq][it][ip]
                <<endl;

          }
#endif
        }
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////
  //write the acc table
  sprintf(tablename,"ACCEPTANCE_%s_%d",mDetName.c_str(),mType);
  name_exc[0]=tablename;
  ACCEPTANCE::CreateIncFile(Acc_Inc,bin1,bin2,bin3,bin4,name_exc,10);

  /*
  //write det count and true count inc files
  sprintf(tablename,"N_det_%s_%d",mDetName.c_str(),mType);
  name_exc[0]=tablename;
  ACCEPTANCE::CreateIncFile(N_Inc,bin1,bin2,bin3,bin4,name_exc,10);
  sprintf(tablename,"N_true_%s_%d",mDetName.c_str(),mType);
  name_exc[0]=tablename;
  ACCEPTANCE::CreateIncFile(N_Inc_true,bin1,bin2,bin3,bin4,name_exc,10);
  */

  ////////////////////////////////////////////////////////////////////////////////////////  
  //another way to create the acc table for only one table each time  
  sprintf(tablename,"ACCEPTANCE_%s_%d_L1.txt",mDetName.c_str(),mType);
  PrintACCTable(1, tablename, 0, N_Inc, N_Inc_true);
  sprintf(tablename,"ACCEPTANCE_%s_%d_L0.txt",mDetName.c_str(),mType);
  PrintACCTable(0, tablename, 0, N_Inc, N_Inc_true);


  ////////////////////////////////////////////////////////////////////////////////////////  
  //I also generate the whole table for other people
  int LEVEL=3;
  ofstream fout;
  sprintf(filename,"ACCEPTANCE_%s_%d.txt",mDetName.c_str(),mType);
  fout.open(filename);
  fout<<"//static const int kYtarBinNum="<<ACCEPTANCE::kYtarBinNum<<";     //"<<endl;
  fout<<"//static const int kDeltaBinNum="<<ACCEPTANCE::kDeltaBinNum<<";    //"<<endl;
  fout<<"//static const int kThetaBinNum="<<ACCEPTANCE::kThetaBinNum<<";  //"<<endl;
  fout<<"//static const int kPhiBinNum="<<ACCEPTANCE::kPhiBinNum<<";   //"<<endl;

  fout<<"//Ytar  Delta  Theta    Phi Thrown#    Det#     Acc"<<endl;
  for(iw=0;iw<ACCEPTANCE::kYtarBinNum;iw++)
  {  
    for(iq=0;iq<ACCEPTANCE::kDeltaBinNum;iq++)
    {  
      for(it=0;it<ACCEPTANCE::kThetaBinNum;it++)
      {
        for(ip=0;ip<ACCEPTANCE::kPhiBinNum;ip++)
        {
          sprintf(str,"%6.2f %6.2f %6.3f %6.3f %7d %7d %7.4f",
            (ACCEPTANCE::kYtarBin[iw]+ACCEPTANCE::kYtarBin[iw+1])/2.,
            (ACCEPTANCE::kDeltaBin[iq]+ACCEPTANCE::kDeltaBin[iq+1])/2.,
            (ACCEPTANCE::kThetaBin[it]+ACCEPTANCE::kThetaBin[it+1])/2.,
            (ACCEPTANCE::kPhiBin[ip]+ACCEPTANCE::kPhiBin[ip+1])/2.,
            N_Inc_true[iw][iq][it][ip],
            N_Inc[iw][iq][it][ip],
            (N_Inc_true[iw][iq][it][ip]==0)? 0.0:(double)N_Inc[iw][iq][it][ip]/(double)N_Inc_true[iw][iq][it][ip]
          );
          //cout<<str<<endl;
          if(LEVEL>=3)  fout<<str<<endl;
          else if(LEVEL==2)
          { //skip those bins thrown# is zero	
            if (N_Inc_true[iw][iq][it][ip]>0)  fout<<str<<endl;
          }
          else //if(LEVEL<=1)
          { //skip those bins N_Exc# is zero
            if (N_Inc[iw][iq][it][ip]>0)  fout<<str<<endl;
          }
        }
      }  
    } 
  }  
  fout.close();

}

/////////////////////////////////////////////////////////////////////////////////////
void ExtractAcceptance::Reset()
{
}

/////////////////////////////////////////////////////////////////////////////////////
void ExtractAcceptance::MergeResult(ExtractAcceptance *pAcc)
{
  int bin1=ACCEPTANCE::kYtarBinNum;
  int bin2=ACCEPTANCE::kDeltaBinNum;
  int bin3=ACCEPTANCE::kThetaBinNum;
  int bin4=ACCEPTANCE::kPhiBinNum;
  
  ACCEPTANCE::AddArrayBToA(N_Inc, pAcc->N_Inc,bin1,bin2,bin3,bin4);
  ACCEPTANCE::AddArrayBToA(N_Inc_true, pAcc->N_Inc_true,bin1,bin2,bin3,bin4);
}

/////////////////////////////////////////////////////////////////////////////////////
void ExtractAcceptance::FillAccTree()
{    
  char filename[255];
  sprintf(filename,"ACCEPTANCE_%s_%d.root",mDetName.c_str(),mType);
  TFile *fAccFile=new TFile(filename,"recreate","acceptance file");
  int bin1=ACCEPTANCE::kYtarBinNum;
  int bin2=ACCEPTANCE::kDeltaBinNum;
  int bin3=ACCEPTANCE::kThetaBinNum;
  int bin4=ACCEPTANCE::kPhiBinNum; 

  //////////////////////////////////////////////////////////////////////////
  //create trees
  int iw,iq,it,ip;
  float pYtar,pDelta,pTheta,pPhi; //mean value of these 4 variable
  float acc;
  int n,n_true;

  TTree *exc=new TTree("acc","accceptance tree");
  exc->Branch("iyt",&iw,"iyt/I");
  exc->Branch("ide",&iq,"ide/I");
  exc->Branch("ith",&it,"ith/I");
  exc->Branch("iph",&ip,"iph/I");
  exc->Branch("Ytar",&pYtar,"Ytar/F");
  exc->Branch("Delta",&pDelta,"Delta/F");
  exc->Branch("Theta",&pTheta,"Theta/F");
  exc->Branch("Phi",&pPhi,"Phi/F");
  exc->Branch("acc",&acc,"acc/F");
  exc->Branch("n",&n,"n/I");
  exc->Branch("n_true",&n_true,"n_true/I");
  
  for(iw=0;iw<bin1;iw++)
  {
    for(iq=0;iq<bin2;iq++)
    {  
      for(it=0;it<bin3;it++)
      {
        for(ip=0;ip<bin4;ip++)  
        {
          pYtar=(ACCEPTANCE::kYtarBin[iw]+ACCEPTANCE::kYtarBin[iw+1])/2.0;
          pDelta=(ACCEPTANCE::kDeltaBin[iq]+ACCEPTANCE::kDeltaBin[iq+1])/2.0;
          pTheta=(ACCEPTANCE::kThetaBin[it]+ACCEPTANCE::kThetaBin[it+1])/2.0;
          pPhi=(ACCEPTANCE::kPhiBin[ip]+ACCEPTANCE::kPhiBin[ip+1])/2.0;

          n=N_Inc[iw][iq][it][ip];
          n_true=N_Inc_true[iw][iq][it][ip];
          acc=(n>0 && n_true>0)?float(n)/float(n_true):0.0;
          exc->Fill();
        }
      }
    }
  }
  exc->Write();
  exc->Delete();

  /////////////////////////////////////////////////////////////////
  TTree *bininfo=new TTree("bininfo","accceptance binning information tree");

  //please note that const type is not accepted by root tree::branch()
  //try the following 3 lines you will understand this

  //static  int kYtarBinNum=19;  //Ytar for exclusive, 120 MeV each bin
  //static  const int kYtarBinNum=19;  //Ytar for exclusive, 120 MeV each bin
  //bininfo->Branch("kYtarBinNum",&kYtarBinNum,"kYtarBinNum/I");
    
  int nbin1=1+ACCEPTANCE::kYtarBinNum;
  int nbin2=1+ACCEPTANCE::kDeltaBinNum;
  int nbin3=1+ACCEPTANCE::kThetaBinNum;
  int nbin4=1+ACCEPTANCE::kPhiBinNum;
  float afloat=0.0;
  float* vbin1 = ACCEPTANCE::CopyToTypeT(nbin1,ACCEPTANCE::kYtarBin,afloat);
  float* vbin2 = ACCEPTANCE::CopyToTypeT(nbin2,ACCEPTANCE::kDeltaBin,afloat);
  float* vbin3 = ACCEPTANCE::CopyToTypeT(nbin3,ACCEPTANCE::kThetaBin,afloat);
  float* vbin4 = ACCEPTANCE::CopyToTypeT(nbin4,ACCEPTANCE::kPhiBin,afloat);
  
  bininfo->Branch("kYtarBinNum",&nbin1,"kYtarBinNum/I");  
  bininfo->Branch("kYtarBin",vbin1,"kYtarBin[kYtarBinNum]/F");  
  bininfo->Branch("kDeltaBinNum",&nbin2,"kDeltaBinNum/I");
  bininfo->Branch("kDeltaBin",vbin2,"kDeltaBin[kDeltaBinNum]/F");
  bininfo->Branch("kThetaBinNum",&nbin3,"kThetaBinNum/I");
  bininfo->Branch("kThetaBin",vbin3,"kThetaBin[kThetaBinNum]/F");  
  bininfo->Branch("kPhiBinNum",&nbin4,"kPhiBinNum/I");
  bininfo->Branch("kPhiBin",vbin4,"kPhiBin[kPhiBinNum]/F");  
 
  bininfo->Fill();
  bininfo->Write();
  bininfo->Delete();

  ACCEPTANCE::FreeDynArray(vbin1);
  ACCEPTANCE::FreeDynArray(vbin2);
  ACCEPTANCE::FreeDynArray(vbin3);
  ACCEPTANCE::FreeDynArray(vbin4);
  /////////////////////////////////////////////////////////////////////////////////////////////

  fAccFile->Close();
  fAccFile->Delete();
}

