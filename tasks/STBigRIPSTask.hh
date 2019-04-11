#ifndef STBigRIPSTask_hh
#define STBigRIPSTask_hh

#include "FairLogger.h"
#include "TSystem.h"
#include "FairTask.h"
#include "FairRunAna.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TCutG.h"
#include "TString.h"
#include "STBDC.hh"

class STBigRIPSTask : virtual public FairTask
{
public:

  STBigRIPSTask();
  ~STBigRIPSTask(){};

  InitStatus Init();
  void Exec(Option_t* opt);

  void     SetPersistence(Bool_t value = kTRUE);
  Long64_t GetEventID() {return fEventIDLast;}
  Bool_t   SetupInputDataFile();

private:
  Bool_t fIsPersistence;  ///< Persistence check variable            
  FairLogger *fLogger;           //!
  FairRootManager* fRootManager; //!

  TClonesArray* fBDCArray;
  TClonesArray* fEvent;
  STBDC*        fBDC;
  
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

  UInt_t iRun;
  UInt_t SnA;
  UInt_t beamPID;

  TString rootDir; //!       

  TChain *ribfChain;
  TChain *bdcChain;

  TCutG *g132Sn;
  TCutG *g108Sn;
  TCutG *g124Sn;
  TCutG *g112Sn;


private:
  Long64_t fEventID;      //!   
  Long64_t fEventIDLast;  //!  
  
  void Clear();
  
  void    ProceedEvent();
  Bool_t  SetupBeamCut();
  void    SetupBeamA(UInt_t vRun);

public:
  Int_t   GetBeamPID();

  ClassDef(STBigRIPSTask, 0);
};


#endif
