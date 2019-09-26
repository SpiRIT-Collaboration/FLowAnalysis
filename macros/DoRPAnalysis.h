#ifndef  DORPANALYSIS_HH
#define  DORPANALYSIS_HH

//@@@@ 
Int_t  seltrackID = 4;
UInt_t selReactionPlanef = 10;

Int_t  seltrack;

// Reading tree
Int_t   iVer[2];
UInt_t  iRun;
TString sRun;
TString sSuf;
TString sVer;
TString dVer;
Bool_t  bRedo;

void      SetEnvironment();
Bool_t    DefineVersion();

void      PrintProcess(Int_t ievt);
void      Open();
void      OutputTree();

TString    finname;
TString    foutname;

TClonesArray     *aBDC;
TClonesArray     *aParticleArray;
TClonesArray     *aFlowArray;
TClonesArray     *anewFlow;

STFlowTask       *aFlowTask;


Long64_t  nEntry;
TTree    *mflw;
TFile    *fout;
TChain   *rChain;

#endif
