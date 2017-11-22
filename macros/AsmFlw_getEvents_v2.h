#ifndef  ASMPLW_GETEVENTSHH
#define  ASMPLW_GETEVENTSHH

#include "TVector3.h"
#include "TVector2.h"

// Reading tree
TClonesArray *flowcorrArray    = NULL;
Int_t         ntrack[7];
Int_t         kymult;
Double_t      ProjA;
Double_t      ProjB;
Double_t      aoq;
Double_t      z;
Int_t         snbm;

void      SetEnvironment();
void      PrintProcess(Int_t ievt);
void      SetNumberOfProcess(Int_t nmax);
void      Open();
void      Initialize();
void      OutputTree(Long64_t val);
Long64_t  GetRandomNumTrack();
Long64_t  GetMultiplicityDistribution();
Bool_t    DefineVersion();

STParticle* GetRealTrack(Long64_t ival);
STParticle* GetMixedTrack(Long64_t *ival, Int_t *kval);

void      ResetPID(STParticle *apart);
void      SetPtWeight(STParticle *apart);
Bool_t    CheckPID(STParticle *apart);

void      LoadPIDFile();

Int_t   iVer[3];
TString sRun;
TString sMix;
Bool_t  bMix;  // kTRUE mixing kFALSE real data
TString sVer;
TString sAsm;
Int_t   iAsm;
TString sbRun;
TString sbVer;
TString sBinp;
UInt_t  nBin; 
TString binpara;
Int_t   mxntrk;
Int_t  maxProc;;

//TChain *fChain;
TTree  *fTree;
Long64_t nEntry;

TRandom3 rnd;
vector<UInt_t> trackID;

TFile *fout;
TTree *mflw      = NULL;
TCutG *gProton   = NULL;
TCutG *gDeuteron = NULL;
TCutG *gTriton   = NULL;
TCutG *gPip      = NULL;
TCutG *gPim      = NULL;

vector<TVector2> pt;

TClonesArray *aParticleArray = NULL;
TClonesArray *nParticleArray = NULL;

vector<TString> vfname;

UInt_t  binmapsize = 0;

// Tree out
Int_t   iRun;
Int_t    numGoodTrack;
Int_t    mtrack;
vector<Int_t>    event;

Double_t         aX;
Double_t         bY;

TClonesArray     *p_rot = NULL;
TClonesArray     *p_org = NULL;

// Tree out end


TRandom2 pran;
TH1I *histGT_r = NULL;
TH1I *histGT;
TH1I *histMixEvt;
TH1I *histMixTrack;
TH1D *hRPrapd;
TH1I *hm_t;
TH1I *hm_b;
TH1I *hgtc;
TH1D *hvphi;
TH1D *hvthet;
TH1I *hvmtk;

TFile *mhfile;


// Important parameters.
UInt_t ntr_diff = 5;  //ver. 2.0.14 
//UInt_t ntr_diff = 3;  // ver. 2.0.15


#endif
