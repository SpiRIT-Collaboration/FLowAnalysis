#ifndef STNEULANDCLUSTER_HH
#define STNEULANDCLUSTER_HH

#include <TObject.h>
#include <TVector3.h>
#include <TMath.h>

class STNeuLANDCluster : public TObject
{

 public:
  STNeuLANDCluster();
  virtual ~STNeuLANDCluster() {}
  void Clear(Option_t * ="");

  void SetNHit  (Int_t v){nhit = v;}
  void SetEdep  (Double_t v){edep = v;}
  void SetTOF   (Double_t v){tof = v;}

  void SetLocalX(Double_t v){localx = v;}
  void SetLocalY(Double_t v){localy = v;}
  void SetLocalZ(Double_t v){localz = v;}
  void SetLocalPos();
  void SetLocalPosLast();


  void SetTOFLast   (Double_t v){tof_last = v;}
  void SetLocalXLast(Double_t v){localx_last = v;}
  void SetLocalYLast(Double_t v){localy_last = v;}
  void SetLocalZLast(Double_t v){localz_last = v;}

  void SetVetoHitAll(Int_t v){veto_all = v;}
  void SetVetoHitOne(Int_t v){veto_bar = v;}
  void SetVetoHitMid(Int_t v){veto_mid = v;}
  void SetVetoHitLoose(Int_t v){veto_loose = v;}

  Int_t GetNHit (){return nhit;}
  Double_t GetEdep(){return edep;}
  Double_t GetTOF(){return tof;}
  Double_t GetBeta();
  Double_t GetGamma();

  Double_t GetLocalX(){return localx;}
  Double_t GetLocalY(){return localy;}
  Double_t GetLocalZ(){return localz;}

  // for global coordinate. target is at z=0 
  TVector3 GetGlobalPos()      {return globalPos;}
  TVector3 GetGlobalPos_last() {return globalPos_last;}

  Double_t GetGlobalX();
  Double_t GetGlobalY();
  Double_t GetGlobalZ();
  Double_t GetPx();
  Double_t GetPy();
  Double_t GetPz();

  Double_t GetTOFLast(){return tof_last;}
  Double_t GetLocalXLast(){return localx_last;}
  Double_t GetLocalYLast(){return localy_last;}
  Double_t GetLocalZLast(){return localz_last;}

  Int_t GetVetoHitAll(){return veto_all;}
  Int_t GetVetoHitOne(){return veto_bar;}
  Int_t GetVetoHitMid(){return veto_mid;}
  Int_t GetVetoHitLoose(){return veto_loose;}

  // on the assumption of neutron
  Double_t GetMom(){return nmass * GetBeta() * GetGamma(); }

 private:
  Int_t    nhit;
  Double_t edep;
  Double_t tof;
  Double_t localx;
  Double_t localy;
  Double_t localz;
  Double_t distance;

  Double_t tof_last;
  Double_t localx_last;
  Double_t localy_last;
  Double_t localz_last;

  TVector3 localPos;
  TVector3 localPos_last;
  TVector3 globalPos;
  TVector3 globalPos_last;
  TVector3 global_offset;
  

  Int_t veto_all; // strong cut, if there is one hit on any veto
  Int_t veto_bar; // strong cut, if there is one hit on veto bar in front of the cluster, delta x<30cm
  Int_t veto_mid; // intermediate cut, if there is one hit on veto delta x < 22cm && delta y < 30cm
  Int_t veto_loose; // loose cut, if there is one hit on veto delta x < 18cm && delta y < 20cm

  // constant values
  Double_t distance_to_center; //! distance from target to center of NeuLAND 1st plane
  Double_t angle; //!
  Double_t c;     //!
  Double_t nmass; //! neutron mass

  ClassDef(STNeuLANDCluster, 2);

};

#endif
