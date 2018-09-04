#ifndef MC_PROCESS1_H
#define MC_PROCESS1_H

#include <TROOT.h>
#include <TChain.h>
#include <TCutG.h>
#include <TFile.h>
#include <TVector3.h>
#include <TVector2.h>
#include "TSystem.h"
#include "TROOT.h"
#include "TClonesArray.h"


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


TString recoFile;

TChain *fChain;


TClonesarray *trackArray  = NULL;
TClonesArray *vertexArray = NULL;
STMCEventHeader *mceventHeader = NULL;
STMCTriggerResponse *mctriggerResponse = NULL;

ULong_t event_number;
Float_t max_veto;


vector<Int_t> *bitpat  = 0;
TBranch       *bbitpat = 0;


Int_t    SnA;


void      Setup();
void      BeamPID();
Int_t     GetBeamPID();
void      SetBeamOnTarget(TVector2 vt);
void      SetBeamOnTarget();

Bool_t    CheckBeamPosition();
Bool_t    CheckVertex(STParticle *aPart);
Bool_t    CheckBDCvsVertexCorrelation(TVector2 vxy);
Bool_t    CheckTrackonTarget(TVector3 trackatTarget);
Bool_t    CheckBDCvsTrackCorrelation(TVector3 trackatTarget);

void      SetTPC();


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
Int_t nhitnl[9];

std::vector<Double_t> phi;
std::vector<Int_t> gTrack;
std::vector<Int_t>gTgTrack;

TString sRun;
TString sVer;
Int_t   iRun;
Int_t   iVer;

Bool_t  STPC;
Bool_t  BigRIPS;
Bool_t  KyotoArry;
Bool_t  KATANA;
Int_t   KyotoRoot;
Bool_t  NeuLAND;

#endif
