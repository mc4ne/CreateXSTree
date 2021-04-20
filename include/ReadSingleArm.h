//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 11 16:23:38 2019 by ROOT version 5.28/00
// from TTree h1/HUT
// found on file: hms_12.5deg_carbon_2070.root
//////////////////////////////////////////////////////////

#ifndef ReadSingleArm_h
#define ReadSingleArm_h
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
using namespace std;

class ReadSingleArm {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Declaration of leaf types
  Float_t         psxfp;
  Float_t         psyfp;
  Float_t         psxpfp;
  Float_t         psypfp;
  Float_t         psxtari;
  Float_t         psztari;
  Float_t         psytari;
  Float_t         psdeltai;
  Float_t         psyptari;
  Float_t         psxptari;
  Float_t         psztar;
  Float_t         psytar;
  Float_t         psdelta;
  Float_t         psyptar;
  Float_t         psxptar;
  Float_t         fry;
  Float_t         xsnum;
  Float_t         ysnum;
  Float_t         xsieve;
  Float_t         ysieve;
  Float_t         stop_id;
  Float_t         vxi;
  Float_t         vyi;
  Float_t         vzi;

  // List of branches
  TBranch        *b_psxfp;   //!
  TBranch        *b_psyfp;   //!
  TBranch        *b_psxpfp;   //!
  TBranch        *b_psypfp;   //!
  TBranch        *b_psxtari;   //!
  TBranch        *b_psztari;   //!
  TBranch        *b_psytari;   //!
  TBranch        *b_psdeltai;   //!
  TBranch        *b_psyptari;   //!
  TBranch        *b_psxptari;   //!
  TBranch        *b_psztar;   //!
  TBranch        *b_psytar;   //!
  TBranch        *b_psdelta;   //!
  TBranch        *b_psyptar;   //!
  TBranch        *b_psxptar;   //!
  TBranch        *b_fry;   //!
  TBranch        *b_xsnum;   //!
  TBranch        *b_ysnum;   //!
  TBranch        *b_xsieve;   //!
  TBranch        *b_ysieve;   //!
  TBranch        *b_stop_id;   //!
  TBranch        *b_vxi;   //!
  TBranch        *b_vyi;   //!
  TBranch        *b_vzi;   //!

  ReadSingleArm();
  ReadSingleArm(const char* filename, int hms_or_shms);
  virtual ~ReadSingleArm();
  void LoadAFile(const char* filename, int hms_or_shms);
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree, int hms_or_shms);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ReadSingleArm_cxx
ReadSingleArm::ReadSingleArm()
{
  fChain=0;
  //do nothing
}

ReadSingleArm::ReadSingleArm(const char* filename, int hms_or_shms)
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree. 
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
  if (!f) {
    cout<<"ReadSingleArm() is reading file \""<<filename<<"\" ...\n";
    f = new TFile(filename);
  }
  TTree* tree = 0;
  if(hms_or_shms==1) {
    tree = (TTree*)gDirectory->Get("h1");
    if(!tree) cout<<"could not find h1 tree in this file \""<<filename<<"\".\n";
  }
  else {
    tree = (TTree*)gDirectory->Get("h1411");
    if(!tree) cout<<"could not find h1411 tree in this file \""<<filename<<"\".\n";
  }
  
  Init(tree,hms_or_shms);
}

ReadSingleArm::~ReadSingleArm()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

void ReadSingleArm::LoadAFile(const char* filename, int hms_or_shms)
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree. 
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
  if (!f) {
    cout<<"ReadSingleArm::LoadAFile() is reading file \""<<filename<<"\"  ...\n";
    f = new TFile(filename);
  }
  TTree* tree = 0;
  if(hms_or_shms==1) {
    tree = (TTree*)gDirectory->Get("h1");
    if(!tree) cout<<"could not find h1 tree in this file \""<<filename<<"\".\n";
  }
  else {
    tree = (TTree*)gDirectory->Get("h1411");
    if(!tree) cout<<"could not find h1411 tree in this file \""<<filename<<"\".\n";
  }
  
  Init(tree,hms_or_shms);
}

Int_t ReadSingleArm::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t ReadSingleArm::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void ReadSingleArm::Init(TTree *tree, int hms_or_shms)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  if(hms_or_shms==1) {
    fChain->SetBranchAddress("hsxfp", &psxfp, &b_psxfp);
    fChain->SetBranchAddress("hsyfp", &psyfp, &b_psyfp);
    fChain->SetBranchAddress("hsxpfp", &psxpfp, &b_psxpfp);
    fChain->SetBranchAddress("hsypfp", &psypfp, &b_psypfp);
    fChain->SetBranchAddress("hsxtari", &psxtari, &b_psxtari);
    fChain->SetBranchAddress("hsztari", &psztari, &b_psztari);
    fChain->SetBranchAddress("hsytari", &psytari, &b_psytari);
    fChain->SetBranchAddress("hsdeltai", &psdeltai, &b_psdeltai);
    fChain->SetBranchAddress("hsyptari", &psyptari, &b_psyptari);
    fChain->SetBranchAddress("hsxptari", &psxptari, &b_psxptari);
    fChain->SetBranchAddress("hsztar", &psztar, &b_psztar);
    fChain->SetBranchAddress("hsytar", &psytar, &b_psytar);
    fChain->SetBranchAddress("hsdelta", &psdelta, &b_psdelta);
    fChain->SetBranchAddress("hsyptar", &psyptar, &b_psyptar);
    fChain->SetBranchAddress("hsxptar", &psxptar, &b_psxptar);
  } else {
    fChain->SetBranchAddress("psxfp", &psxfp, &b_psxfp);
    fChain->SetBranchAddress("psyfp", &psyfp, &b_psyfp);
    fChain->SetBranchAddress("psxpfp", &psxpfp, &b_psxpfp);
    fChain->SetBranchAddress("psypfp", &psypfp, &b_psypfp);
    fChain->SetBranchAddress("psxtari", &psxtari, &b_psxtari);
    fChain->SetBranchAddress("psztari", &psztari, &b_psztari);
    fChain->SetBranchAddress("psytari", &psytari, &b_psytari);
    fChain->SetBranchAddress("psdeltai", &psdeltai, &b_psdeltai);
    fChain->SetBranchAddress("psyptari", &psyptari, &b_psyptari);
    fChain->SetBranchAddress("psxptari", &psxptari, &b_psxptari);
    fChain->SetBranchAddress("psztar", &psztar, &b_psztar);
    fChain->SetBranchAddress("psytar", &psytar, &b_psytar);
    fChain->SetBranchAddress("psdelta", &psdelta, &b_psdelta);
    fChain->SetBranchAddress("psyptar", &psyptar, &b_psyptar);
    fChain->SetBranchAddress("psxptar", &psxptar, &b_psxptar);
  }
  
  fChain->SetBranchAddress("fry", &fry, &b_fry);
  fChain->SetBranchAddress("xsnum", &xsnum, &b_xsnum);
  fChain->SetBranchAddress("ysnum", &ysnum, &b_ysnum);
  fChain->SetBranchAddress("xsieve", &xsieve, &b_xsieve);
  fChain->SetBranchAddress("ysieve", &ysieve, &b_ysieve);
  fChain->SetBranchAddress("stop_id", &stop_id, &b_stop_id);
  fChain->SetBranchAddress("vxi", &vxi, &b_vxi);
  fChain->SetBranchAddress("vyi", &vyi, &b_vyi);
  fChain->SetBranchAddress("vzi", &vzi, &b_vzi);
  Notify();
  //cout<<"ReadSingleArm()::Init():  all is good! Number of Event is "<<fChain->GetEntries()<<"\n";
}

Bool_t ReadSingleArm::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void ReadSingleArm::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t ReadSingleArm::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

#endif // #ifdef ReadSingleArm_cxx
