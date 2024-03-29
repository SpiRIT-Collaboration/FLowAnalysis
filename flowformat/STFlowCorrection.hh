#ifndef STFLOWORRECTION_HH
#define STFLOWORRECTION_HH


#include <fstream>
#include <iostream>
#include <TMath.h>
#include "TString.h"
#include "TCollection.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TVector3.h"
#include "FairLogger.h"

class STFlowCorrection : public TObject {

 public:
  STFlowCorrection(){fLogger = FairLogger::GetLogger();}
  STFlowCorrection( TChain *chele, UInt_t ival1, UInt_t ival2) {Initialize(chele, ival1, ival2);}
  STFlowCorrection(UInt_t ival1, UInt_t ival2)                    {Initialize(ival1, ival2);}
  ~STFlowCorrection(){};


 private :
  FairLogger *fLogger;      //!

  UInt_t     irm; // 1 = mix, otherwise real
  UInt_t     iver;
  UInt_t     harm;
  UInt_t     charm;
  UInt_t    *indx;
  Double_t  *An;
  Double_t  *Bn;
  Double_t  *An_rms;
  Double_t  *Bn_rms;
  TString    fname;
  TChain    *ChEle = NULL;
  std::vector<Double_t> vphi;  //! corrected Phi
  std::vector<Double_t> rcphi; //! ReCentering phi
  std::vector<Double_t> bphi;  //! original phi
  std::vector<TVector3> vvec;  //! corrected vector
  std::vector<TVector3> bvec;  //! original vector

  Double_t   constX;
  Double_t   meanX;
  Double_t   sigX;
  Double_t   constY;
  Double_t   meanY;
  Double_t   sigY;


  Double_t   binmax[2];
  Double_t   binmin[2];
  TString    binpara[2];

  std::vector<Double_t> vtheta;
  std::vector<Int_t>    vmtrack;

 public:
  void   Initialize(TChain *chele=NULL, UInt_t ival1=4, UInt_t ival2=0);
  void   Initialize(UInt_t ival1=4, UInt_t ival2=0);
  void   SetHarmonics(UInt_t ival){charm = ival;}
  void   SetRealOrMix(UInt_t ival){irm  = ival;}


  void   Add(Double_t val)                 {vphi.push_back(val);}
  void   Add(Double_t val1, Double_t val2) {vphi.push_back(val1);    vtheta.push_back(val2);}
  void   Add(Int_t ival,    Double_t val1) {vmtrack.push_back(ival); vphi.push_back(val1);}
  void   Add(Int_t ival,    Double_t val1, Double_t val2) 
  { vmtrack.push_back(ival);  vphi.push_back(val1);   vtheta.push_back(val2);  }
  void   Add(Int_t ival, TVector3 vval) 
  { vmtrack.push_back(ival);  bvec.push_back(vval);   bphi.push_back(vval.Phi()); }

  void   clear();


  Int_t                 GetNPhi()         {return bphi.size();}
  std::vector<Double_t> GetOriginalPhi()  {return bphi;}             
  std::vector<Double_t> GetCorrectedPhi() {return vphi;}             
  std::vector<Double_t> GetReCeneringPhi(){return rcphi;}
  std::vector<Double_t> GetTheta();       
  std::vector<Int_t>    GetMTrack()      {return vmtrack;}
  Double_t         GetMTrackMean(){return TMath::Mean(vmtrack.begin(), vmtrack.end());}
  Double_t         GetThetaMean() {return TMath::Mean(vtheta.begin(),  vtheta.end());}
  

  TVector3 ReCentering(TVector3 val);
  UInt_t   ReCenteringFourierCorrection();
  TVector3 ReCenteringFourierCorrection(TVector3 val);
  UInt_t   FourierCorrection();
  void     FourierCorrection(std::vector<Double_t> &val);
  void     FourierCorrection(Float_t hrm, std::vector<Double_t> &val);
  Double_t GetCorrection(Double_t val);
  void     GetCorrection(std::vector<Double_t> &val);


  UInt_t   LoadCorrectionFactor(UInt_t val=0);
  UInt_t   SaveCorrectionFactor(TString comm1=":", TString comm2="");
  void     PrintContents();
  void     PrintRange();

  UInt_t   GetNumberOfParam() {return (UInt_t)vphi.size();};
  UInt_t   SetHarmonics()     {return harm;}

  TVector3 GetReCentering(TVector3 vec);
  Double_t *GetAverageCosin(Int_t ival, std::vector<Double_t> &val);

  
  Int_t    GetNHarmonics() {return harm;}
  void     SetFileName(TString sval);
  TString  GetFileName()   {return fname;}

  void     SetBin_max(UInt_t idx=0, Double_t val=999.) {if( idx < 3) binmax[idx] = val;}
  void     SetBin_min(UInt_t idx=0, Double_t val=0.) {if( idx < 3) binmin[idx] = val;}
  Double_t GetBin_max(UInt_t idx=0)  {if(idx < 3) return binmax[idx]; else return -1.;}
  Double_t GetBin_min(UInt_t idx=0)  {if(idx < 3) return binmin[idx]; else return -1.;}
  TString  GetBinParameter(UInt_t idx=0) {return binpara[idx];}

  void     SetReCenteringParameter(TString cprm, Double_t val[]);

  void   ShowParameters();
  void   ShowBinInformation();

private:

  void   SetFileName();
  void   Init();
  void   SetDirectory();

  ClassDef(STFlowCorrection,0);
};


#endif
