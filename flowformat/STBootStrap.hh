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

  Double_t GetPhiOriginal() {return TVector2::Phi_mpi_pi( phi_off );}

  Double_t GetMean()       {return cnvMean; }
  Double_t GetCosMean()    {return cnvCosMean;}
  Double_t GetStdDev()     {return cnvStdv;}
  Double_t GetStdDevError();
  Double_t GetStdDev2()    {return cnvStdv2;}
  Double_t GetNElem(){return (Double_t)numElements;}
  Double_t GetResidualMean(UInt_t ival);
  Double_t GetResidualStdDev(UInt_t ival);


  std::vector< Double_t > GetMeanVector()      {return resMean;   }
  std::vector< Double_t > GetStdDevVector()    {return resStdv;   }
								  
  Double_t  GetMeanCnvVector()   {return cnvMean;}
  Double_t  GetCosMeanVector()   {return cnvCosMean;}
  Double_t  GetStdDevCnvVector() {return cnvStdv;   }
  Double_t  GetStdDev2CnvVector(){return cnvStdv2;  }
  

  void StoreResults();
  
  std::vector< UInt_t > Resampling(UInt_t ival);
  UInt_t BootStrapping(UInt_t nbt = 0);

private:


  std::vector< Double_t > elements;
  std::vector< TVector2 > elementsTV2;

  UInt_t numElements;
  Double_t     phi_off = -99.;

  UInt_t nboot = 0;

  std::vector< Double_t >    replace;    // resampling event
  std::vector< Double_t >    resMean;    // <Phi> / resampling event
  std::vector< Double_t >    resStdv;    // std<Phi> / resampling event
  
  Double_t cnvMean;    // <Phi> for bootstrapped events
  Double_t cnvStdv;    // std<Phi> for bootstrapped events
  Double_t cnvCosMean; // cos(std_Phi) for bootstrapped events
  Double_t cnvStdv2;   // std<std<Phi>> for bootstrapped events

  TRandom3 rnd;

  ClassDef(STBootStrap,0);
};

#endif
