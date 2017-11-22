#ifndef ASSEMBLE_FLWV3_H
#define ASSEMBLE_FLWV3_H

#include <TROOT.h>
#include <TChain.h>
#include <TCutG.h>
#include <TFile.h>
#include <TVector3.h>
#include <TVector2.h>
#include "TSystem.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "KatanaRoot/KatanaRoot_Load.h"
#include "KyotoRoot/STTriggerArray.hh"



//setParameters                                                                                               
Double_t tx_right= -14.7;
Double_t tx_left =  19.4;
Double_t ty_top  =  15.62;
Double_t ty_btm  = -13.74;

//Genie's 132Sn Aug 30 2017
Double_t vrt_Zmin = -12.80121;
Double_t vrt_Zmax =  -9.49569;
Double_t trktgt_right = -15;
Double_t trktgt_left  =  15;
Double_t trktgt_top   = -206.06;
Double_t trktgt_btm   = -246.06 ;

// KATANA offset                                                                     
Double_t KATANA_frame_OffSet     = -1041.+ 25;
Double_t KATANA_paddle_Width     =   100.;
Double_t KATANA_paddle_Height    =   380.;
Double_t KATANA_paddle_Right[12] = {35,135,235,335,435,535,635,855,956,1057,1156,1256};
Double_t KATANA_Max_Dist = 100.; //55.;                            

  //set run2334                                                                          
Double_t beamVx_offset = -17.54;
Double_t beamVx_sigma  =  1.53;
//Double_t beamVy_offset = -2.28e+02;
Double_t beamVy_offset = 2.27e+02;
Double_t beamVy_sigma  =  1.08;
Double_t beamVx_nsig   =  4;
Double_t beamVy_nsig   =  4;

Double_t trackVx_sigma  =  3.292;
Double_t trackVy_sigma  =  2.7;
Double_t trackVx_nsig   =  6;
Double_t trackVy_nsig   =  4;
Double_t trackVy_offset = -224.3;

Double_t beamVx_min = -80.;
Double_t beamVx_max =  40.;

Double_t  BeamBendingDistance = 250;

Int_t  ibeam;

TCutG *g132Sn;
TCutG *g108Sn;
TCutG *g124Sn;
TCutG *g112Sn;
TCutG *gPip;

TChain *fChain;
TChain *ribfChain;
TChain *bdcChain;
TChain *kChain;
TChain *kaChain;

  //ribf data
  //Declaration of leaves types                                                                                         
Double_t         aoq;
Double_t         z;
Double_t         tof;
Double_t         tx;
Double_t         ty;
Double_t         beta;
Double_t         brho;
Double_t         isGood;
Double_t         intZ;
Double_t         intA;

Double_t         bdcax;
Double_t         bdcby;
Double_t         ProjX;
Double_t         ProjY;
Double_t         ProjZ;
Double_t         ProjP;
Double_t         ProjPX;
Double_t         ProjPY;
Double_t         ProjPZ;
Double_t         ProjA;
Double_t         ProjB;



TClonesArray *trackArray  = NULL;
TClonesArray *vertexArray = NULL;
KatanaRoot   *katanaroot  = NULL;
TriggerBox   *triggerbox  = NULL;
TClonesArray *tpcParticle = NULL;

TBranch      *brtrackArray;
TBranch      *brvertexArray;


ULong_t event_number;
Float_t max_veto;
Int_t   katanaM;

vector<Int_t> *bitpat  = 0;
TBranch       *bbitpat = 0;


Int_t  kynHit;
vector<int>   *kyhitch;
vector<float> *kyhitx ;
vector<float> *kyhitz ;

TBranch *bkyhitch;
TBranch *bkyhitx ;
TBranch *bkyhitz ;

Int_t    SnA;
TString  rootDir;
TString  ktnrootDir;
TString  kytrootDir;


void      Setup();
void      BeamPID();
Int_t     GetBeamPID();
void      SetDataDirectory();
void      SetKATANADirectory();
void      SetBeamOnTarget(TVector2 vt);
void      SetBeamOnTarget();

Bool_t    CheckBeamPosition();
Bool_t    CheckVertex(TVector3 vec);
Bool_t    CheckBDCvsVertexCorrelation(TVector2 vxy);
Bool_t    CheckTrackonTarget(TVector3 trackatTarget);
Bool_t    CheckBDCvsTrackCorrelation(TVector3 trackatTarget);

void      SetTPC();
void      SetKATANARoot();
void      SetKATANARoot_bt();
Bool_t    SetKyotoArray();
Bool_t    SetKyotoMultiplicity();
void      SetBigRIPS();

void      Initialize(Int_t ievt);
void      OutputTree(Int_t nmax);
Bool_t    DefineVersion();



// ouput tree
TString foutname;
TFile *fout;
TTree *flw ;

Int_t evtid;
Int_t ntrack[7];
Int_t bmpid;
Int_t kyotomL;
Int_t kyotomR;
Int_t KatanaMult;
Int_t trknKATANA;
Int_t foundKATANA;
TVector2 *BeamonTarget;


std::vector<Double_t> phi;
std::vector<Int_t> gTrack;
std::vector<Int_t>gTgTrack;

TString sRun;
TString sVer;
Int_t   iRun;
Int_t   iVer[2];

Bool_t  BigRIPS;
Bool_t  KyotoArry;
Bool_t  KATANA;
Int_t   KyotoRoot;


#endif
