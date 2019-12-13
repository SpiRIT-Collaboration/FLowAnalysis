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
TString oVer;
TString dVer;
Bool_t  bRedo;
Float_t dMct;
UInt_t  isys;


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

STFlowInfo       *aFlowInfo;
UInt_t            beamPID;

Double_t          RPPsi;

Long64_t  nEntry;
TTree    *mflw;
TFile    *fout;
TChain   *rChain;
TF1      *fv1y;
TF1      *fv2y;

#endif
