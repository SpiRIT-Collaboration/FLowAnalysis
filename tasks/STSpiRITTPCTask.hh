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
  void   SetRunInfo(TString vDir, TString tver, TString sver){rootDir = vDir, tVer = tver, sVer=sver;}  

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
  Double_t ProjA;  //!
  Double_t ProjB;  //!


private:
  Long64_t fEventID;      //! 
  Long64_t fEventIDLast;  //!
  Long64_t nEntry;
  Long64_t prcEntry;      //!
  TDatime  dtime;         //! 
  TDatime  beginTime;     //!  

  Double_t MassRegion[7][4] ={{ 127.2,   21.3,      4.,  4.},            //pi  
			      { 911.044, 68.4656,   2.,  2.},            //p  685.3 to 1,165.9
			      { 1874.76, 129.962,   1.5, 1.5},           //d  1,552.3 to 
			      { 2870.62, 212.896,   1.,  1.},            //t 2463
			      { 2760.47, 196.659,   1.,  1.},            //He3 
			      { 3720.77, 255.71,    1.,  1.},            //He4 
			      { 5751.97, 673.339,   0.5, 0.5}};      //! //He6  // read from fitted function
 
  Double_t MassRegionLU[7][2] = { {   0.0,  250.0},     // pi
				  { 500.0, 1137.0},     // p
				  {1137.0, 2200.0},     // d
				  {2200.0, 3200.0},     // t
				  {2430.0, 2950.0},     // fBBMassHe  3He
				  {2950.0, 4000.0},     // fBBMassHe  4He
				  {4000.0, 7000.0}};    //!           6He

  void   Clear();
  Bool_t SetupParameters();
  Bool_t SetupInputDataFile();
  Bool_t SetupFlowDataBase();

  void   ProceedEvent();
  void   SetupEventInfo();
  void   SetupTrackQualityFlag(STParticle *apart);
  Int_t  GetPID(Double_t mass[2], Double_t dedx);
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


