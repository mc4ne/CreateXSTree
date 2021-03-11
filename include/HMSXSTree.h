//HMSXSTree.h 
// create my own tree from mc-single-arm
// ---------------------------------------------

#ifndef _HMSXSTree_
#define _HMSXSTree_

#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <exception>
#include <TROOT.h>
#include <TStyle.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TObject.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TProfile.h>
#include <TString.h>
#include "ReadHMS.h"

using namespace std;


#define HMSXSTree_Debug 4
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class HMSXSTree:ReadHMS
{
public:
    HMSXSTree(const char* filename);
    virtual ~HMSXSTree();

    void  BeginOfRun();
    void  Run(int nevents_to_process=-1);
    void  EndOfRun();
    ////////////////////////////////////////////////////////////////////
    //P. Bosted XS model wrapper
    double GetXS(int Z, int N, float Ei, float Ef, float Theta, float Tb=-0.000, float Ta=-0.000);
    double GetXS_GE180(float Ei, float Ef, float Theta, float Tb=-0.000, float Ta=-0.000);
    
public:
    void SetTreeLevel(int pLevel) {mTreeLevel=pLevel;cout<<"Update HMSXSTree::mTreeLevel to "<<mTreeLevel<<endl; };
    int  GetTreeLevel() {return mTreeLevel;};
    void Reset();  //reset the tree variables    
   
    void SetElas(int v) {mElasOnly=v;};
    void SetPara(float pBeam, float pTheta, float pMom) {mBeam=pBeam;mDetAngle=pTheta;mDetMom=pMom;};
    //================================================================
private:
    //mTreeLevel controls the tree output level
    //0 fill all into the tree
    //1 skip 2-particle events, keep only D(e,e'pi-)X or D(e,e'Ps)X events
    //2 keep only D(e,e'pi-)X events and 0.7<Missing Mass<1.2 only
    int mTreeLevel;
    int mElasOnly;

private:
    TFile*   mFile;
    TTree*   mTree;
    TString  mOutFileName;
    
    //-----------------------------block3 start  -----------------------
    //ini
    float vx0,vy0,vz0,p0,theta0,phi0;
    float xtg0,ytg0,ztg0,xptg0,yptg0,delta0;

    //intermedia
    int    ixs,iys;  //sieve hole index
    float  xsieve,ysieve;
    float  xfp,yfp,xpfp,ypfp;
    int    istop;

    //recon
    float  vx,vy,vz,p,theta,phi;
    float  xtg,ytg,ztg,xptg,yptg,delta;

    //kin
    float  nu,Q2,W,xbj;
    //rate

    float p_rate_he3, p_rate_ge180, p_rate_up, p_rate_down, p_rate_c12, p_rate_n2, p_spec, p_accept, th_accept, ph_accept, dens_he3, dens_ge180, dens_c12, dens_n2, tar_len_ge180, tar_len_up, tar_len_down, tar_len_he3, tar_len_c12, tar_len_n2, n_trials;
    //yield
    float  p_yield_he3, p_yield_ge180, p_yield_up, p_yield_down, p_yield_c12, p_yield_n2, beam_curr;

    //-----------------------------block3 end  -------------------------

    //-----------------------------block2 start  -----------------------
    //this block will not be filled if mTreeLevel>=2
    float  xs_1h,xs_3he,xs_4he,xs_12c,xs_14n,xs_27al,xs_ge180; 
    //-----------------------------block2 end  -------------------------
    
    //-----------------------------block1 start  -----------------------
    //this block will not be filled if mTreeLevel>=1
    int    mDet;
    float  mBeam,mDetAngle,mDetMom;
    //-----------------------------block1 end  -------------------------

};
#endif //_HMSXSTree_
