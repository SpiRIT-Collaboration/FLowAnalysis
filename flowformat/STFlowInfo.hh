#ifndef STFLOWINFO_H
#define STFLOWINFO_H

#include "vector"
#include <TObject.h>
#include <iostream>
#include <TVector3.h>
#include "FairLogger.h"

class STFlowInfo : public TObject
{

public:
  STFlowInfo();
  STFlowInfo(const STFlowInfo &cp);
  STFlowInfo &operator=(const STFlowInfo &cp);

  virtual ~STFlowInfo(){};
  
public:
  Int_t    run;
  Long64_t evt;
  Int_t    SnA;
  Int_t    beamPID;

  UInt_t   ntrack[7];
  UInt_t   mtrack0;
  UInt_t   mtrack1;
  UInt_t   mtrack2;
  UInt_t   mtrack3;
  UInt_t   mtrack4;
  UInt_t   mtrack5;
  UInt_t   mtrack6;

  
  TVector2 unitP;
  TVector3 unitP_fc;
  TVector3 unitP_rc;

  TVector2 unitP_1;
  TVector2 unitP_2;
  UInt_t   mtrack_1;
  UInt_t   mtrack_2;
  
  Double_t bsPhi[3];
  Double_t bsPhi_1[3];
  Double_t bsPhi_2[3];

public:

  void Clear();
  void SetRun(Int_t val){ run = val; }
  void SetBeamA(Int_t val)  { SnA = val; }
  void SetBeamPID(Int_t val){ beamPID = val; }
  void SetEventID(Long64_t val){ evt = val; }
  void SetNTrack(UInt_t *nval);
  UInt_t *GetNTrack() { return ntrack; }
  UInt_t GetNTrack(UInt_t ival) { 
    if( ival < 7 )
      return ntrack[ival];
    else 
      return 0;
  }

  ClassDef(STFlowInfo,1);
};

#endif
