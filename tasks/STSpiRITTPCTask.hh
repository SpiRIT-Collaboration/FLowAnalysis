#ifndef STSpiRITTPCTask_hh
#define STSpiRITTPCTask_hh

// FAIRROOT classes
#include "FairTask.h"
#include "FairRunAna.h"
#include "FairRootManager.h"
#include "STRecoTrack.hh"
#include "STParticle.hh"
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
  void   SetRunInfo(TString vDir, TString tver, TString sver){rootDir = vDir, tVer = tver, sVer=sver;
    LOG(INFO) << " setruninfor " << sVer << FairLogger::endl;}  

  Long64_t GetEventID()  {return fEventIDLast;}
  Long64_t GetEntries() {return nEntry;}
  void     SetProcessingNumberOfEvent(Long64_t val) { val < nEntry ? prcEntry = val : prcEntry = nEntry;
    LOG(INFO) << " process number " << prcEntry << FairLogger::endl;
  }

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

  TString rootDir; //!       
  TChain *fChain;  //!
  STEventHeader *eventHeader    = NULL;  //!
  TClonesArray *trackArray     = NULL;  //!
  TClonesArray *trackVAArray   = NULL;  //!
  TClonesArray *vertexArray    = NULL;  //!
  TClonesArray *vertexVAArray  = NULL;  //!
  TClonesArray *vertexBDCArray = NULL;  //!

  TClonesArray   *tpcParticle;  //!

  //--- mass fitter
  STMassCalculator* massCal;    //!

  TClonesArray *flowAnalysis; //!
  STFlowInfo   *fflowinfo;    //!
  STFlowTask   *fflowtask;    //!

  TClonesArray *aflowcorrArray[2];  //!

  UInt_t  ntrack[7]; //!
  
  Int_t    BeamPID;//!
  Int_t    BeamIndex; //!
  Double_t ProjA;  //!
  Double_t ProjB;  //!

private:
  Long64_t fEventID;      //! 
  Long64_t fEventIDLast;  //!
  Long64_t nEntry;
  Long64_t prcEntry;      //!
  TDatime  dtime;         //! 
  TDatime  beginTime;     //!  

  // by v37:  
  //vertex mean
  TVector3 VtxMean[4] = { TVector3( 2.25909e+00,  -2.26504e+02,  -1.34982e+01), //! 132Sn
			  TVector3(-9.66045e-01,  -2.29565e+02,  -1.37814e+01), //! 108Sn
			  TVector3( 1.47193e+00,  -2.27352e+02,  -1.34916e+01), //! 124Sn
			  TVector3( 1.55937e+00,  -2.29191e+02,  -1.32485e+01)};//! 112Sn
  // vetex sigma Z);
  Double_t VtxSigm[4] = { 1.91418,  //! 132Sn 
			  1.84580,  //! 108Sn
			  1.56109,  //! 124Sn
			  1.41496}; //! 112Sn


  Double_t MassRegion[7][4] ={{ 127.2,   21.3,      4.,  4.},            //pi  
			      { 911.044, 68.4656,   2.,  2.},            //p  685.3 to 1,165.9
			      { 1874.76, 129.962,   1.5, 1.5},           //d  1,552.3 to 
			      { 2870.62, 212.896,   1.,  1.},            //t 2463
			      { 2760.47, 196.659,   1.,  1.},            //He3 
			      { 3720.77, 255.71,    1.,  1.},            //He4 
			      { 5751.97, 673.339,   0.5, 0.5}};      //! //He6  // read from fitted function
 
  Double_t MassRegionLU_N[7][2] = { {   0.0,  250.0},     // pi
				    { 750.0, 1100.0},     // p
				    {1500.0, 2200.0},     // d
				    {2450.0, 3200.0},     // t
				    {2300.0, 2950.0},     // fBBMassHe  3He
				    {3200.0, 4100.0},     // fBBMassHe  4He
				    {4100.0, 7000.0}};    //!           6He

  Double_t MassRegionLU_L[7][2] = { {   0.0,  250.0},     // pi
				    { 700.0, 1200.0},     // p
				    {1400.0, 2300.0},     // d
				    {2350.0, 3300.0},     // t
				    {2200.0, 3050.0},     // fBBMassHe  3He
				    {3100.0, 4200.0},     // fBBMassHe  4He
				    {4200.0, 7000.0}};    //!           6He

  Double_t MassRegionLU_T[7][2] = { {   0.0,  250.0},     // pi
				    { 800.0, 1050.0},     // p
				    {1600.0, 2100.0},     // d
				    {2550.0, 3100.0},     // t
				    {2400.0, 2850.0},     // fBBMassHe  3He
				    {3300.0, 4000.0},     // fBBMassHe  4He
				    {4000.0, 7000.0}};    //!           6He

  void   Clear();
  Bool_t SetupParameters();
  Bool_t SetupInputDataFile();
  Bool_t SetupFlowDataBase();

  void   ProceedEvent();
  void   SetupEventInfo();
  void   SetupTrackQualityFlag(STParticle *apart);
  Int_t  GetPID(Double_t mass[2], Double_t dedx);
  Int_t  GetPIDTight(Double_t mass[2], Double_t dedx);
  Int_t  GetPIDNorm (Double_t mass[2], Double_t dedx);
  Int_t  GetPIDLoose(Double_t mass[2], Double_t dedx);

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


