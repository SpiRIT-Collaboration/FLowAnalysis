#ifndef STBOOTSTRAP_HH
#define STBOOTSTRAP_HH

#include "TH1D.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "TMath.h"
#include <iostream>

class STBootStrap : public TObject {

public:
  STBootStrap(){ clear(); };
  STBootStrap(UInt_t ival);
  STBootStrap(UInt_t ival1, UInt_t ival2, Double_t *sample);
  STBootStrap(UInt_t ival1, UInt_t ival2, TVector2 *sample);
  STBootStrap(UInt_t ival1, std::vector<TVector2> *sample);
  ~STBootStrap(){};

  void clear();

  void SetBootNumber(UInt_t ival) {nboot = ival;}
  void Add(TVector2 sample);
  void Add(Double_t sample);

  Double_t GetMean()           {return cnvMean; }
  Double_t GetMod()            {return cnvMod; }
  Double_t GetCosMean()        {return cnvCosMean;}
  Double_t GetStdDev()         {return cnvStdv;}
  Double_t GetStdDevError();
  Double_t GetStdDev2()        {return cnvStdv2;}
  Double_t GetNElem()          {return (Double_t)numElements;}
  Double_t GetResidualMean(UInt_t ival);
  Double_t GetResidualStdDev(UInt_t ival);
  Double_t GetOrdinaryMean()    {return ordn_Mean;}
  Double_t GetOrdinaryStdDev()  {return ordn_StdDev;}
  TVector2 GetOrdinarySum()     {return org_sum;}

  Double_t GetCLLow()           {return CL_25lw;}
  Double_t GetCLUp()            {return CL_975up;}
  Double_t GetError()           {return Error;}

  std::vector< Double_t > GetReplaceVector()   {return replace; }
  std::vector< Double_t > GetMeanVector()      {return resMean; }
  std::vector< Double_t > GetStdDevVector()    {return resStdv; }

  void StoreConfideneLevel();
  
  std::vector< UInt_t > Resampling(UInt_t ival);
  UInt_t BootStrapping(UInt_t nbt = 0);
  UInt_t BootStrappingDouble(UInt_t nbt);
  UInt_t BootStrappingTVector2(UInt_t nbt);


private:

  std::vector< Double_t > elements;
  std::vector< TVector2 > elementsTV2;

  UInt_t numElements;
  UInt_t numElementsTV2;

  Double_t     ordn_Mean = -99.;
  Double_t     ordn_StdDev = 0.;
  TVector2     org_sum;
  Double_t     r_sum;

  UInt_t nboot = 0;

  std::vector< Double_t >    replace;    // resampling event
  std::vector< Double_t >    resMean;    // <Phi> / resampling event
  std::vector< Double_t >    resStdv;    // std<Phi> / resampling event
  
  Double_t cnvMean;    // <Phi> for bootstrapped events
  Double_t cnvStdv;    // std<Phi> for bootstrapped events
  Double_t cnvCosMean; // cos(std_Phi) for bootstrapped events
  Double_t cnvStdv2;   // std<std<Phi>> for bootstrapped events
  Double_t cnvMod;     // Vector2 Mod for bootstrapped events


  Double_t CL_25lw;
  Double_t CL_975up;
  Double_t Error;

  TRandom3 rnd;

  TH1D *hbsphi; // temporary
public:
  TH1D *GetBSPhiDistribution() {return hbsphi;}
  
  ClassDef(STBootStrap,0);
};

#endif
