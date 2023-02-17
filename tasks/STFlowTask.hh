#ifndef STFlowTask_hh
#define STFlowTask_hh

#include <TMath.h>
//#include "STParticle.hh"
#include "STKParticle.hh"
#include "STBDC.hh"
#include "STFlowInfo.hh"
#include "STBootStrap.hh"
#include "TDatime.h"
#include "TRandom3.h"
#include "STFlowCorrection.hh"
#include "STFlowCorrection.hh"
#include "STLorentzBoostVector.hh"
#include "STRunToBeamA.hh"

class STFlowTask
{
public:
  STFlowTask(Bool_t bfltn, Bool_t bsub);
  STFlowTask(Bool_t bfltn, Bool_t bsub, Bool_t bbst);

  virtual ~STFlowTask();

  Bool_t Init(UInt_t irun, TString aver);

private:
  UInt_t iRun;                   //!
  Bool_t fIsFlowCorrection;      //!< reaction plane flattening correction
  Bool_t fIsBootStrap;           //!< bootstrap analysis flag
  Bool_t fIsSubeventAnalysis;    //!< Subevent analysis flag  
  UInt_t selReactionPlanef;      //! track quality for reaction plane 
  UInt_t iSystem;                //!  system ID {"(132Sn + 124Sn)", "(108Sn + 112Sn)", "(124Sn + 112Sn)", "(112Sn + 124Sn)", "(p + p)"}
  Float_t fRPMidCut;             //!  exclude mid-rapidity reation from the RP calculation abs(y_cm)<fRPMidCut

  UInt_t PID_sel;                //!

  TString sVer; //!

  STFlowInfo     *fflowinfo;    //!
  TClonesArray   *tpcParticle;  //!
  TClonesArray   *aflowcorrArray[4];
  
  UInt_t  ntrack[7]; //!

  Bool_t SetupParameters();
  Bool_t SetupInputDataFile();
  Bool_t SetupFlowDataBase();

  void   SetupTrackQualityFlag(STKParticle *apart);
  void   SetupTrackExtraQualityFlag(STKParticle *apart);


  // for flow analysis
  vector<UInt_t> mtkbin[4];		             //!
  vector< TString > vfname[4];	                     //!
  vector< vector<Double_t> >  binmax;	             //!
  vector< vector<Double_t> >  binmin;	             //!
  vector< pair<Double_t, Double_t> > pbinmin[4];     //!
  TRandom3 rnd; //!
  STBootStrap *bs_unitP; //!
  
  Double_t sum_omg2;  //!
  Double_t sum_omg;   //!


  UInt_t *RandomDivide2(const UInt_t npart);
  UInt_t SetFLowDatabaseFiles(UInt_t version);
  Int_t  GetMultiplicityCorretionIndex(UInt_t isel, UInt_t ival);

public:
  void   Clear();
  UInt_t *RandumPickUp(const UInt_t val, const UInt_t npart);
  void   SetNTrack(UInt_t i, UInt_t val);
  void   SetupFlow(STKParticle &apart);
  void   SetSelectionOfReactionPlaneFlag(UInt_t val) {selReactionPlanef = val;}
  void   SetFlowInfo(STFlowInfo *aflowinfo);
  void   DoFlowAnalysis(STKParticle &apart);

  Bool_t DoFlattening();
  Bool_t DoFlatteningSub();
  TVector3 DoFlattening(UInt_t isel, TVector3 mvec, UInt_t ntrk);
  void   DoSubeventAnalysis();
  void   DoSubeventAnalysisFixedMultiplicity(UInt_t val);

  TVector3 Psi_FlatteningCorrection(UInt_t isel, Int_t ival, TVector3 Pvec);
  TVector3 Psi_ReCenteringCorrection(UInt_t isel, Int_t ival, TVector3 Pvec);
  void     DoIndividualReactionPlaneAnalysis();
  void     SetIndividualReactionPlane( STKParticle &apart );
  void     SetIndividualReactionPlane_recal( STKParticle &apart );

  void   SetupEventInfo(Long64_t eventid, UInt_t val);
  void   SetFlowTask( TClonesArray &atpcParticle );
  void   SetFlowTask();
  void   FinishEvent();
  void   SetParticleArray( TClonesArray &atpcParticle ){tpcParticle = &atpcParticle;}

  void   SetFlowCorrection(Bool_t val)   {fIsFlowCorrection = val;}
  Bool_t GetFlowCorrectionFlag()         {return fIsFlowCorrection;}
  void   SetSubEventAnalysis(Bool_t val) {fIsSubeventAnalysis = val;}
  Bool_t GetSubEventAnalysisFlag()       {return fIsSubeventAnalysis;}
  void   SetBootStrap(Bool_t val)        {fIsBootStrap = val;}
  Bool_t GetBootStrapFlag()              {return fIsBootStrap;}



  STFlowInfo  *GetFlowInfo() {return fflowinfo;}
  void    SetGoodEventFlag(UInt_t val){fflowinfo->SetGoodEventFlag(val);}
  UInt_t  GetGoodEventFlag(){return fflowinfo->GetGoodEventFlag();}

  void    SetRPMidRapidityCut(Float_t val) {fRPMidCut = val; fflowinfo->fRPMidCut = fRPMidCut;}

  void   SetPIDSelection(UInt_t ival) {PID_sel = ival;}
  UInt_t GetPIDSelection()            {return PID_sel;}

};

#endif


