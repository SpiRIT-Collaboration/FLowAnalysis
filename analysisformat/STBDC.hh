#ifndef STBDC_H
#define STBDC_H

#include "vector"
#include <TObject.h>
#include <iostream>


class STBDC : public TObject
{

public:
  STBDC();
  STBDC(const STBDC &cp);
  STBDC &operator=(const STBDC &cp);

  virtual ~STBDC(){};
  
public:  
  Int_t    run;
  Long64_t evt;
  Int_t    SnA;
  Int_t    beamPID;

  Double_t         aoq;
  Double_t         z;
  Double_t         tof;
  Double_t         beta;      //!
  Double_t         brho;      //!
  Double_t         isGood;
  Double_t         intZ;
  Double_t         intA;

  Double_t         bdcax;     
  Double_t         bdcby;     
  Double_t         ProjX;     
  Double_t         ProjY;     
  Double_t         ProjZ;     
  Double_t         ProjP;     //!
  Double_t         ProjPX;    //!
  Double_t         ProjPY;    //!
  Double_t         ProjPZ;    //!
  Double_t         ProjA;
  Double_t         ProjB;

  void      SetRun(Int_t val){ run = val;}
  void      SetEventID(Long64_t val){ evt = val;}

  Double_t  GetProjA(){return ProjA;}
  Double_t  GetProjB(){return ProjB;}

  void      SetBeamPID(Int_t val){ beamPID = val;}
  Int_t     GetBeamPID() {return beamPID;}
  
  Int_t     GetRunNumber() {return run;}
  

  void Clear();

  ClassDef(STBDC,1);
};

#endif
