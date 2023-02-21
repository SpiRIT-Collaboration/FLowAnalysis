#ifndef STSpiRITTPCTask_hh
#define STSpiRITTPCTask_hh

// FAIRROOT classes
#include "FairTask.h"
#include "FairRunAna.h"
#include "FairRootManager.h"
#include "STRecoTrack.hh"
#include "STKParticle.hh"
#include "STMassCalculator.hh"
#include "ST_ClusterNum_DB.hh"
#include "STBDC.hh"
#include "STFlowInfo.hh"
#include "STFlowTask.hh"
#include "STLorentzBoostVector.hh"
#include "STBootStrap.hh"
#include "STEventHeader.hh"
#include "TDatime.h"
#include "TRandom3.h"
#include "STFlowCorrection.hh"
#include <TMath.h>
#include "STFlowCorrection.hh"
#include "STRunToBeamA.hh"
#include "TObject.h"
#include "FairLink.h"
#include "STRunToBeamA.hh"
#include "STMassCalSimpleBB.hh"


class STSpiRITTPCTask : public FairTask
{
public:
  STSpiRITTPCTask();

  virtual ~STSpiRITTPCTask();

  InitStatus Init();
  void       Exec(Option_t *opt);
  void       FinishEvent();
  void       Finish(){fActive = kFALSE;}

  void   SetPersistence(Bool_t value = kTRUE);
  void   SetRunInfo(TString vDir, TString tver, TString sver)
  {rootDir = vDir, tVer = tver, sVer=sver;
    LOG(INFO) << " setruninfor " << sVer << FairLogger::endl;}  

  Long64_t GetEventID()  {return fEventIDLast;}
  Long64_t GetEntries() {return nEntry;}
  void     SetProcessingNumberOfEvent(Long64_t val) { val < nEntry ? prcEntry = val : prcEntry = nEntry;
    LOG(INFO) << " process number " << prcEntry << FairLogger::endl;
  }

  Bool_t SetupBeamCut(UInt_t SnA);
  Int_t  GetBeamPID(Int_t SnA, Double_t aoq, Double_t z);

private:
  Bool_t fIsPersistence;         ///< Persistence check variable
  Bool_t fIsFlowAnalysis;        ///< flow analysis flag
  Bool_t fIsFlowCorrection;      ///< reaction plane flattening correction
  Bool_t fIsBootStrap;           ///< bootstrap analysis flag
  Bool_t fIsSubeventAnalysis;    ///< Subevent analysis flag  
  UInt_t selReactionPlanef;      // track quality for reaction plane 


  FairRunAna      *fRunAna;       //!
  FairLogger      *fLogger;       //!
  FairRootManager *fRootManager;  //!
  
  ST_ClusterNum_DB* db;    //!
  
  UInt_t  iRun;
  TString tVer;  //!
  TString sVer;  //!
  UInt_t  bmA;   //!

  TString rootDir; //!       
  TChain *fChain;  //!
  STEventHeader *eventHeader    = NULL;  //!
  STBeamInfo    *beamInfo       = NULL;

  TClonesArray *trackArray     = NULL;  //!
  TClonesArray *trackVAArray   = NULL;  //!
  TClonesArray *vertexArray    = NULL;  //!
  TClonesArray *vertexVAArray  = NULL;  //!
  TClonesArray *vertexBDCArray = NULL;  //!

  TClonesArray   *tpcParticle;  //!

  //--- mass fitter
  STMassCalculator* massCal;    //!
  STMassCalSimpleBB* massCalH;  //!
  STMassCalSimpleBB* massCalHe; //!

  TClonesArray *flowAnalysis; //!
  STFlowInfo   *fflowinfo;    //!
  STFlowTask   *fflowtask;    //!

  TClonesArray *aflowcorrArray[2];  //!

  UInt_t  ntrack[7]; //!
  
  STBDC   *fBeam;  //!
  Int_t    BeamPID;//!
  Double_t ProjA;  //!
  Double_t ProjB;  //!
  Double_t ProjX;  //!
  Double_t ProjY;  //!
  TCutG    *gBeamCut; 

private:
  Long64_t rEventID;      //! read event ID from recotrack data  
  Long64_t fEventID;      //! counted event ID
  Long64_t fEventIDLast;  //!
  Long64_t nEntry;
  Long64_t prcEntry;      //!
  TDatime  dtime;         //! 
  TDatime  beginTime;     //!  

  // for v38(20190804/data/):  
  //STVertex.fPos
  // TVector3 VtxMean[4] = { TVector3( 2.10312, -2.04848e+02, -1.50300e+01), // 132Sn
  // 			  TVector3(-8.25318e-01, -2.06242e+02, -1.44491e+01), //108Sn
  // 			  TVector3(-8.25318e-01, -2.06242e+02, -1.44491e+01), //108Sn (temp)
  // 			  TVector3(-8.25318e-01, -2.06242e+02, -1.44491e+01)}; //108Sn(temp)
  
  // Double_t VtxSigm[4] = { 2.0505,   //! 132Sn 
  // 			  1.38199,  //108Sn
  // 			  1.38199,  //124Sn (temp)
  // 			  1.38199}; //112Sn (temp)

  // Kaneko-kun's after v53
  TVector3 VtxMean[4] = { TVector3( 2.10312,     -2.04848e+02, -1.485675e+01), // 132Sn
			  TVector3(-8.25318e-01, -2.06242e+02, -1.4847947e+01), //108Sn
			  TVector3(-8.25318e-01, -2.06242e+02, -1.4758978e+01), //108Sn (temp)
			  TVector3(-8.25318e-01, -2.06242e+02, -1.4421470e+01)}; //108Sn(temp)
  
  Double_t VtxSigm[4] = { 1.2932625,   //! 132Sn 
			  1.3323066,  //108Sn
			  1.2022664,  //124Sn (temp)
			  0.98318320}; //112Sn (temp)

  Double_t  protonMaxMomentum = 2500.;  //!  
  Double_t  momRange[8][2] = {{0.,1000.},{0.,1000.},{100.,1400.}, {100.,2200.}, {200.,2800.}, {300.,1400.}, {400.,1800.}, {1300.,2500.}};  
  Double_t  momCut[7] = {1000, 1000, 2000., 4000., 4200., 4200.,4500.};

  Double_t MassRegion[7][4] ={{ 127.2,   21.3,      4.,  4.},            //pi  
			      { 911.044, 68.4656,   2.,  2.},            //p  685.3 to 1,165.9
			      { 1874.76, 129.962,   1.5, 1.5},           //d  1,552.3 to 
			      { 2870.62, 212.896,   1.,  1.},            //t 2463
			      { 2760.47, 196.659,   1.,  1.},            //He3 
			      { 3531.15, 278.23,    1.,  1.},            //He4 
			      { 5751.97, 673.339,   0.5, 0.5}};      //! //He6  // read from fitted function
 
  Double_t MassRegionLU_L[7][2] = { {   0.0,  400.0},     // pi
				    { 500.0, 1300.0},     // p
				    {1400.0, 2300.0},     // d
				    {2350.0, 3300.0},     // t
				    {2200.0, 3050.0},     // fBBMassHe  3He
				    {3100.0, 4200.0},     // fBBMassHe  4He
				    {4200.0, 7000.0}};    //!           6He

  Double_t MassRegionLU_N[7][2] = { {   0.0,  350.0},     // pi
				    { 700.0, 1200.0},     // p
				    {1500.0, 2200.0},     // d
				    {2450.0, 3200.0},     // t
				    {2300.0, 2950.0},     // fBBMassHe  3He
				    {3200.0, 4100.0},     // fBBMassHe  4He
				    {4100.0, 7000.0}};    //!           6He


  Double_t MassRegionLU_T[7][2] = { {   0.0,  250.0},     // pi
				    { 800.0, 1050.0},     // p
				    {1600.0, 2100.0},     // d
				    {2550.0, 3100.0},     // t
				    {2400.0, 2850.0},     // fBBMassHe  3He
				    {3200.0, 4000.0},     // fBBMassHe  4He
				    {4000.0, 7000.0}};    //!           6He

  Double_t MassRange_Fit[6][2] = { {700.,1200},      // Kanekokun's pid range
				   {1500.,2300.},
				   {2400.,3400.},
				   {2565.,3200.}, //{2400.,3200.},_flow.v2
				   {3212.,4307.}, //{3250.,4200.},_flow.v2
				   {5200.,6000.} };   //!

  TF1 *f1MassGate[5][4][2][4][2];//!
  TFile *massGateFile; //!


  void   Clear();
  Bool_t SetupParameters();
  Bool_t SetupInputDataFile();
  Bool_t SetupFlowDataBase();

  Bool_t ProceedEvent();
  Bool_t SetupEventInfo();
  void   SetupTrackQualityFlag(STKParticle *apart);
  Int_t  GetPID(Double_t mass[2], Double_t fMom,  Double_t dedx);
  Int_t  GetPIDTight(Double_t mass[2], Double_t fMom, Double_t dedx);
  Int_t  GetPIDLoose(Double_t mass[2], Double_t fMom, Double_t dedx);
  Int_t  GetPIDFit  (Int_t z, Double_t mass, TVector3 vMom, Int_t mult);
  Bool_t SetupPIDFit();
  void   ShowProcessTime();
  Bool_t GetVertexQuality(TVector3 vert);

public:
  void   SetFlowAnalysis(Bool_t val) {fIsFlowAnalysis = val;}
  Bool_t GetFlowAnalysisFlag()       {return fIsFlowAnalysis;}
  void   SetSubEventAnalysis(Bool_t val) {fIsSubeventAnalysis = val;}
  Bool_t GetSubEventAnalysisFlag()   {return fIsSubeventAnalysis;}
  void   SetBootStrap(Bool_t val)    {fIsBootStrap = val;}
  Bool_t GetBootStrapFlag()          {return fIsBootStrap;}


  // return input file chain
  TChain* GetChain()                 {return fChain;}

  ClassDef(STSpiRITTPCTask, 0);

};

#endif


