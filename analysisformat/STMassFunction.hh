#ifndef STMASSFUNCTION_H
#define STMASSFUNCTION_H

#include "TF1.h"
#include "TVector3.h"
#include "TFile.h"
#include "TSystem.h"


class STMassFunction {

public:
  STMassFunction(){};
  virtual  ~STMassFunction(){};

private:

  const int nTheta = 10;
  const int nPhi   = 18;

  TF1 *fBBMassFunction[10][18];
  TF1 *fLVMassFunction[10][18];


  Double_t fBBMass;
  Double_t fLVMass;

public:
  Double_t GetBBMass(TVector3 mom, Double_t dedx, Double_t chrg);
  Double_t GetLBMass(TVector3 mom, Double_t dedx, Double_t chrg);

  UInt_t GetIndexTheta(Double_t val);
  UInt_t GetIndexPhi(Double_t val);

  Bool_t SetFunction(TString dir, TString fname);
  TF1*   GetBBFunction(Double_t valtheta, Double_t valphi);
  TF1*   GetLVFunction(Double_t valtheta, Double_t valphi);

};

#endif
