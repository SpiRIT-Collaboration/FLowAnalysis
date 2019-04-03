#ifndef  DOFLOWANALYSIS_HH
#define  DOFLOWANALYSIS_HH

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


void      SetEnvironment();
Bool_t    DefineVersion();

void      PrintProcess(Int_t ievt);
void      Open();
void      OutputTree();


TClonesArray     *aBDC;
TClonesArray     *aParticleArray;
TClonesArray     *aFlowArray;
TClonesArray     *anewFlow;

STFlowTask       *aFlowTask;


TString   fname;
Long64_t  nEntry;
TTree    *mflw;
TFile    *fout;
TChain   *rChain;

#endif
