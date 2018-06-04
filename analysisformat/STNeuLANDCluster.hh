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

  void SetNHit  (Int_t v)   {ncnhit = v;}
  void SetEdep  (Double_t v){ncedep = v;}
  void SetTOF   (Double_t v){nctof = v;}
  void SetMass  (UInt_t   v);
  void SetMass  (TString v);

  void SetLocalX(Double_t v){nclocalx = v;}
  void SetLocalY(Double_t v){nclocaly = v;}
  void SetLocalZ(Double_t v){nclocalz = v;}
  void SetLocalPos();
  void SetLocalPosLast();
  void SetLocalPos(TVector3 posv);
  void SetLocalPosLast(TVector3 posv);

  void SetTOFLast   (Double_t v){nctof_last = v;}
  void SetLocalXLast(Double_t v){nclocalx_last = v;}
  void SetLocalYLast(Double_t v){nclocaly_last = v;}
  void SetLocalZLast(Double_t v){nclocalz_last = v;}

  void SetVetoHitAll(Int_t v)  {ncveto_all = v;}
  void SetVetoHitOne(Int_t v)  {ncveto_bar = v;}
  void SetVetoHitMid(Int_t v)  {ncveto_mid = v;}
  void SetVetoHitLoose(Int_t v){ncveto_loose = v;}

private:
  void SetMomentum();

public:
  Int_t GetNHit ()  {return ncnhit;}
  UInt_t GetPID()   {return ncPID;}

  Double_t GetEdep(){return ncedep;}
  Double_t GetTOF() {return nctof;}

  Double_t GetBeta() {return ncbeta;}
  Double_t GetGamma(){return ncgamma;}
  Double_t GetMom()  {return ncP.Mag(); }
  TVector3 GetP()    {return ncP;}
  Double_t GetEnergy()   {return ncE; }
  Double_t GetRapidity() {return ncRapidity; }

  Double_t GetLocalX(){return nclocalx;}
  Double_t GetLocalY(){return nclocaly;}
  Double_t GetLocalZ(){return nclocalz;}

  // for global coordinate. target is at z=0 
  TVector3 GetGlobalPos()      {return ncglobalPos;}
  TVector3 GetGlobalPos_last() {return ncglobalPos_last;}

  Double_t GetGlobalX(){return ncglobalPos.X();};
  Double_t GetGlobalY(){return ncglobalPos.Y();};
  Double_t GetGlobalZ(){return ncglobalPos.Z();};
  Double_t GetPx(){return ncP.X();};
  Double_t GetPy(){return ncP.Y();};
  Double_t GetPz(){return ncP.Z();};

  Double_t GetTOFLast()   {return nctof_last;}
  Double_t GetLocalXLast(){return nclocalx_last;}
  Double_t GetLocalYLast(){return nclocaly_last;}
  Double_t GetLocalZLast(){return nclocalz_last;}

  Int_t GetVetoHitAll(){return ncveto_all;}
  Int_t GetVetoHitOne(){return ncveto_bar;}
  Int_t GetVetoHitMid(){return ncveto_mid;}
  Int_t GetVetoHitLoose(){return ncveto_loose;}

  
  // on the assumption of neutron

 private:
  Int_t    ncnhit;
  Double_t ncedep;
  Double_t nctof;
  Double_t nclocalx;
  Double_t nclocaly;
  Double_t nclocalz;
  Double_t ncdistance;

  Double_t nctof_last;
  Double_t nclocalx_last;
  Double_t nclocaly_last;
  Double_t nclocalz_last;

  TVector3 nclocalPos;
  TVector3 nclocalPos_last;
  TVector3 ncglobalPos;
  TVector3 ncglobalPos_last;
  TVector3 ncglobal_offset;
  
  Double_t ncbeta;
  Double_t ncgamma;
  TVector3 ncP;
  Double_t ncmass; 
  UInt_t   ncPID;
  Double_t ncE;
  Double_t ncRapidity;


  Int_t ncveto_all; // strong cut, if there is one hit on any veto
  Int_t ncveto_bar; // strong cut, if there is one hit on veto bar in front of the cluster, delta x<30cm
  Int_t ncveto_mid; // intermediate cut, if there is one hit on veto delta x < 22cm && delta y < 30cm
  Int_t ncveto_loose; // loose cut, if there is one hit on veto delta x < 18cm && delta y < 20cm

  // constant values
  Double_t distance_to_center; //! distance from target to center of NeuLAND 1st plane
  Double_t angle; //!
  TVector3 targetPos; //!
  Double_t c;     //!
  Double_t target_offset; //!
  Double_t tof_offset;  //!
  ClassDef(STNeuLANDCluster, 3);

};

#endif
