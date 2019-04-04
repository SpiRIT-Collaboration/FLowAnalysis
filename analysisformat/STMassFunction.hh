#ifndef STMASSFUNCTION_H
#define STMASSFUNCTION_H

#include "TF1.h"
#include "TVector3.h"
#include "TFile.h"
#include "TSystem.h"
#include "FairLogger.h"
#include "STBetheBlochFittingFunction.hh"
#include "MassEstimator.h"
#include "STPID.hh"

class STMassFunction {

public:
  STMassFunction() 
    : coordinate(1)
  {fLogger = FairLogger::GetLogger();}


  virtual  ~STMassFunction(){}

private:
  FairLogger* fLogger;   //!

  UInt_t coordinate;

  // coordinate 0:
  const UInt_t nTheta = 10;
  const UInt_t nPhi   = 18;

  // cordinate 1:
  const UInt_t nYawBin = 20;
  const UInt_t nPitchBin = 20;

  UInt_t nBin0;
  UInt_t nBin1;

  TF1 *fBBMassFunction[20][20];
  TF1 *fBBMassGaussFit[10];

  Double_t MassRegion[7][4] ={{ 127.2,21.3,  2.,  2.},            //pi  
			      {  1.,   1.,   2.,  2.},            //p  685.3 to 1,165.9
			      {  1.,   1.,   1.5, 1.5},           //d  1,552.3 to 
			      {  1.,   1.,   1.,  1.},            //t 2463
			      {  1.,   1.,   1.,  1.},            //He3 
			      {  1.,   1.,   1.,  1.},            //He4 
			      {  1,    1.,   0.5, 0.5}};      //! //He6  // read from fitted function
 
  Double_t MassRegionLU[7][2] = { {  50.0,  200.0},     // pi
				  { 500.0, 1137.0},     // p
				  {1137.0, 2200.0},     // d
				  {2200.0, 3200.0},     // t
				  {2430.0, 2950.0},     // fBBMassHe  3He
				  {2950.0, 4000.0},     // fBBMassHe  4He
				  {4000.0, 7000.0}};    //!           6He

public:
  Bool_t GetBBMass(TVector3 mom, Double_t dedx, Int_t chrg, Double_t &fmassH, Double_t &fmassHe, Int_t &pid1, Int_t &pid2);
  Int_t  GetPID(Double_t fmassH, Double_t fmassHe, Double_t dedx);
  Int_t  GetPIDLU(Double_t fmassH, Double_t fmassHe, Double_t dedx);

  UInt_t GetIndexTheta(Double_t val);
  UInt_t GetIndexPhi(Double_t val);
  
  UInt_t GetIndexYaw(TVector3 mom);
  UInt_t GetIndexPitch(TVector3 mom);

  Bool_t SetFunction(TString dir, TString fname, TString funcname);

  Bool_t SetMassFitFunction(TString dir, TString fname, TString funcname);
  void   SetMassRegion();
};

#endif
