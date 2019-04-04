#ifndef STFlowTask_hh
#define STFlowTask_hh

#include <TMath.h>
#include "STParticle.hh"
#include "STBDC.hh"
#include "STFlowInfo.hh"
#include "STBootStrap.hh"
#include "TDatime.h"
#include "TRandom3.h"
#include "STFlowCorrection.hh"
#include "STFlowCorrection.hh"
#include "STLorentzBoostVector.hh"

class STFlowTask
{
public:
  STFlowTask(Bool_t bfltn, Bool_t bsub, Bool_t bbst);

  virtual ~STFlowTask();

  Bool_t Init(UInt_t irun, TString aver);

private:
  Bool_t fIsFlowCorrection;      ///< reaction plane flattening correction
  Bool_t fIsBootStrap;           ///< bootstrap analysis flag
  Bool_t fIsSubeventAnalysis;    ///< Subevent analysis flag  
  UInt_t selReactionPlanef;      // track quality for reaction plane 
  
  TString sVer; //!

  STFlowInfo     *fflowinfo;    //!
  TClonesArray   *tpcParticle;  //!
  TClonesArray   *aflowcorrArray[2];
  
  UInt_t  ntrack[7]; //!

  void   Clear();
  Bool_t SetupParameters();
  Bool_t SetupInputDataFile();
  Bool_t SetupFlowDataBase();

  void   SetupTrackQualityFlag(STParticle *apart);
  void   SetupTrackExtraQualityFlag(STParticle *apart);


  // for flow analysis
  vector<UInt_t> mtkbin[2];		             //!
  vector< TString > vfname[2];	                     //!
  vector< vector<Double_t> >  binmax;	             //!
  vector< vector<Double_t> >  binmin;	             //!
  vector< pair<Double_t, Double_t> > pbinmin[2]; //!
  TRandom3 rnd; //!
  STBootStrap *bs_unitP; //!

  UInt_t *RandomDivide2(const UInt_t npart);
  UInt_t SetFLowDatabaseFiles(UInt_t version);
  Int_t  GetMultiplicityCorretionIndex(UInt_t isel, UInt_t ival);

public:
  void   SetNTrack(UInt_t *nval);
  void   SetupFlow(STParticle &apart);
  void   SetFlowInfo(STFlowInfo *aflowinfo);
  void   DoFlowAnalysis(STParticle &apart);

  Bool_t DoFlattening();
  Bool_t DoFlatteningSub();
  Bool_t DoFlattening(STParticle &part);
  void   DoSubeventAnalysis();
  TVector3 Psi_FlatteningCorrection(UInt_t isel, Int_t ival, TVector3 Pvec);
  TVector3 Psi_ReCenteringCorrection(UInt_t isel, Int_t ival, TVector3 Pvec);
  void     DoIndividualReactionPlaneAnalysis();
  void     SetIndividualReactionPlane( STParticle &apart );

  void   SetupEventInfo(Long64_t eventid, STBDC *aBDC);
  void   SetFlowTask( TClonesArray &atpcParticle );
  void   SetParticleArray( TClonesArray &atpcParticle ){tpcParticle = &atpcParticle;}
  void   FinishEvent();

  void   SetFlowCorrection(Bool_t val)   {fIsFlowCorrection = val;}
  Bool_t GetFlowCorrectionFlag()         {return fIsFlowCorrection;}
  void   SetSubEventAnalysis(Bool_t val) {fIsSubeventAnalysis = val;}
  Bool_t GetSubEventAnalysisFlag()       {return fIsSubeventAnalysis;}
  void   SetBootStrap(Bool_t val)        {fIsBootStrap = val;}
  Bool_t GetBootStrapFlag()              {return fIsBootStrap;}


  STFlowInfo  *GetFlowInfo() {return fflowinfo;}

};

#endif


