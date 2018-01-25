#ifndef STNEULANDHIT_HH
#define STNEULANDHIT_HH

#include <TObject.h>

class STNeuLANDHit : public TObject
{
 public:
  STNeuLANDHit();
  virtual ~STNeuLANDHit() {}
  void Clear(Option_t * ="");

  void SetBarID (Int_t v){bar_id = v;}

  void SetEdep  (Double_t v){edep = v;}
  void SetTOF   (Double_t v){tof = v;}
  void SetBeta  (Double_t v){beta = v;}

  void SetLocalX(Double_t v){localx = v;}
  void SetLocalY(Double_t v){localy = v;}
  void SetLocalZ(Double_t v){localz = v;}

  Int_t GetBarID (){return bar_id;}

  Double_t GetEdep(){return edep;}
  Double_t GetTOF(){return tof;}
  Double_t GetBeta(){return beta;}

  Double_t GetLocalX(){return localx;}
  Double_t GetLocalY(){return localy;}
  Double_t GetLocalZ(){return localz;}

 private:
  Int_t bar_id;
  Double_t edep;
  Double_t tof;
  Double_t beta;
  Double_t localx;
  Double_t localy;
  Double_t localz;

  ClassDef(STNeuLANDHit, 1);

};

#endif
