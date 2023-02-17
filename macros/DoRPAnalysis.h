#ifndef  DORPANALYSIS_HH
#define  DORPANALYSIS_HH

#include "../flowformat/STFlowInfo.hh"

//@@@@ 
Int_t  seltrackID = 4;
UInt_t selReactionPlanef = 10;

Int_t  seltrack;

// Reading tree
Int_t   iVer[2];

Bool_t  bRedo;
Float_t dMct;


void      SetEnvironment();
Bool_t    DefineVersion();

void      PrintProcess(Int_t ievt);
void      Open();
void      OutputTree();

TString    finname;
TString    foutname;


TTree    *mflw;
TFile    *fout;
TF1      *fv1y;
TF1      *fv2y;

#endif
