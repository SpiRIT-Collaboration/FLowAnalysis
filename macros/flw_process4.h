#ifndef  FLW_PROCESS4_HH
#define  FLW_PROCESS4_HH

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

void      SetEnvironment();
void      PrintProcess(Int_t ievt);
void      SetNumberOfProcess(Int_t nmax);
void      Open();
void      Initialize();
void      OutputTree();
Bool_t    DefineVersion();

void      RotateAsBeamAngle(STParticle *apart, TVector3 *p1, TVector2 *pt);
void      SetPtWeight(STParticle *apart);
void      FlatteningCorrection(STParticle *apart, Int_t ival);


void      SubEventAnalysis();
void      AzmAngleRPTReactionPlane();

void      LoadPIDFile();
Int_t     GetPID(Double_t valx, Double_t valy);
UInt_t    SetDatabaseFiles();
Int_t     GetThetaCorretionIndex(Double_t fval, Int_t ival);
Int_t     GetMultiplicityCorretionIndex(UInt_t ival);
void      CheckPlot(UInt_t ival = 0);

Int_t   iVer[3];
TString sRun;
TString sMix;
TString sRot;
Bool_t  bMix;  // kTRUE mixing kFALSE real data
TString sVer;
TString sAsm;
Bool_t  BeamAngle;
Int_t   iAsm;
TString sbRun;
TString sbVer;
TString scVer;
UInt_t  nBin; 
TString binpara;

Int_t  maxProc;;

//TChain *fChain;
TString foutname;
TTree  *fTree;
Long64_t nEntry;

TRandom3 rnd;
vector<UInt_t> trackID;

TFile *fout;
TTree *mflw = NULL;
TCutG *gProton = NULL;

vector<TVector2> pt;

TClonesArray *aParticleArray = NULL;
TClonesArray *aflowcorrArray = NULL;
STFlowCorrection *flowcorr = NULL;

vector<TString> vfname;
vector< vector<Double_t> >  binmax[2];
vector< vector<Double_t> >  binmin[2];

vector< pair<Double_t, Double_t> > pbinmin;
vector< pair<UInt_t  , Double_t> > binmatrix;
vector< pair<UInt_t  , Double_t> >::iterator itbinmatrix;
vector <vector< pair<UInt_t, Double_t> >::iterator > ppbinmin;

vector<UInt_t> mtkbin;

UInt_t  binmapsize = 0;

// Tree out
Int_t   iRun;
Int_t    numGoodTrack;
Int_t    mtrack;
TVector3 unitP;
TVector2 unitP_lang;
TVector2 unitP_rot;
TVector2 unitP_1;
TVector2 unitP_2;
TVector2 unitP_1r;
TVector2 unitP_2r;
UInt_t   mtrack_1;
UInt_t   mtrack_2;

Double_t bsPhi[2];
Double_t bsPhi_1[2];
Double_t bsPhi_2[2];
Double_t bsPhi_ex[3];

Double_t         aX;
Double_t         bY;

TClonesArray     *p_rot = NULL;
TClonesArray     *p_org = NULL;

// Tree out end

TRandom2 pran;
TFile *mhfile;

TH1D *hvphi;
TH1D *hvthet;
TH1I *hvmtk;

// constant
const Double_t ycm = 0.388568; // 132Sn(278MeV/u) + 124Sn
#endif
