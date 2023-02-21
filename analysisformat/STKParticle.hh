/**     
 * @brief STTrack Class         
 *                                   
 * @author Mizuki         
 */

#ifndef STKPARTICLE
#define STKPARTICLE

#include "STBetheBlochFittingFunction.hh"
#include "TCutG.h"
#include "TFile.h"

#include "TVector2.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TLorentzVector.h"

#include "STRecoTrack.hh"
#include "STGenfitVATask.hh"
#include "STPID.hh"
#include "STVertex.hh"

#include <vector>
#include <utility>

class STKParticle : public TObject {
public:
  STKParticle();
  STKParticle(const STKParticle &cp);
  STKParticle &operator=(const STKParticle &cp);

  virtual ~STKParticle(){};

public:

  Int_t    ftrackID;
  UInt_t   fPID;
  UInt_t   fPID_seq;
  UInt_t   fPID_tight;
  UInt_t   fPID_loose;
  Int_t    fpipid;
  Int_t    fChar;
  Int_t    fGFChar;
  Double_t fMass;
  Double_t fBBMass;
  Double_t fBBMassHe;
  Double_t fPIDProbability  ;

  TVector3 fvertex;
  TVector3 forigP3;         // Momemtum vector without any correction.
  TVector3 frmom;           // Momentum vector rotated with respect to the beam angle.
  TLorentzVector fLzvec;    // LorentzVector flzvec(frmom, Etotal) Labo. frame.
  Double_t fRapidity;
  Double_t fRapiditycm;
  Double_t fPt;
  Double_t fPhi;
  Double_t fP;               // Momentum/Q [MeV/c]
  Double_t fdEdx;            // dEdx
  Double_t fNDF;             // STVertex::GetNDF()
  Int_t    fncl;             // number of cluster

  Double_t fwgt;             // Summing up weight
  TVector3 frpv;             // Individual Reaction plane vector (IRPO) for each particle
  Double_t frpphi;           // Individual Reaction plane orientaion (IRPO) for each particle
  Double_t fdeltphi;         // Azimuthal opening angle with respect to IRPO
  TVector3 frpv2;            // Individual v2 Reaction plane vector (IRPO2) for each particle
  Double_t frpphi2;          // Individual v2 Reaction plane orientaion (IRPO2) for each particle
  Double_t fdeltphi2;        // Azimuthal opening angle with respect to IRPO2

  //quality flags
  UInt_t   fVatTargetf      ; //flag for reconstructed vertex XY within the target
  UInt_t   fdistanceatvertexf;
  UInt_t   fNDFf            ;
  UInt_t   fdedxf           ; 
  UInt_t   fnclf            ; // Number Of Cluster Flag
  UInt_t   fmassf           ;
  UInt_t   fmomentumf       ;
  UInt_t   fpidf            ;
  UInt_t   fgoodtrackf      ; // good track flag
  UInt_t   fReactionPlanef  ; // excellent track for flow
  UInt_t   frdEdxPointSizef ;  
  UInt_t   fdoublef         ;

  //STRecoTrack parameters
  STRecoTrack *fRTrack; //!
  STRecoTrack *fVATrack; //!
  Int_t     rVertexID;  
  Int_t     rHelixID;
  Int_t     rNDF;
  Int_t     rncl;
  Double_t  rDist;
  TVector3  rPOCAVertex;  
  Double_t  rChi2;
  Int_t     rdEdxPointSize;
  Int_t     rdEdxPointSize_thr;

private:
  virtual void Clear(Option_t *option = "");
  
  void     Initialize();

  void     SetProperty();

  void     CheckKATANAHit(){};
  void     CheckKYOTOHit(){};


public:
  void     Clean()                       {Clear();}
  void     SetRecoTrack(STRecoTrack *atrack, STRecoTrack *atrackva=NULL);

  Int_t    GetNCL()                      {return rncl;}
  Int_t    GetHelixID()                  {return rHelixID;}

  void     RotateAlongBeamDirection(Double_t valuex, Double_t valuey);
  void     SetP(Double_t value)          {fP = value;}

  void     SetPiPID();
  Int_t    GetPiPID()                    {return fpipid;}

  //  void     SetPID();
  void     SetPID(Int_t value)           { fPID = value; SetMass(fPID);}
  Int_t    GetPID()                      {return fPID;}
  void     SetPIDTight(Int_t value)      {fPID_tight = value;}           
  Int_t    GetPIDTight()                 {return fPID_tight;}           
  void     SetPIDLoose(Int_t value)      {fPID_loose = value;}
  Int_t    GetPIDLoose()                 {return fPID_loose;}           

  Int_t    GetPID_seq()                  {return fPID_seq;}


  void     SetPIDProbability(Double_t value){fPIDProbability = value;}
  Double_t GetPIDProbability()           {return fPIDProbability;}

  void     SetMass(Int_t val);
  //  void     SetMass(Double_t value)       {fMass = value;}
  Double_t GetMass()                     {return fMass;}

  void     SetBBMass(Double_t val);
  Double_t GetBBMass()                   {return fBBMass;}
  void     SetBBMassHe(Double_t val);
  Double_t GetBBMassHe()                 {return fBBMassHe;}

  Double_t GetRapidity()                 {return fRapidity;}

  void     SetRapiditycm(Double_t val)   {fRapiditycm = val;}
  Double_t GetRapiditycm()               {return fRapiditycm;}

  Double_t GetP()                        {return fP;}
  void     SetdEdx(Double_t value)       {fdEdx = value;}
  Double_t GetdEdx()                     {return fdEdx;}


  void     SetIndividualRPVector(TVector3 vec)  {frpv = vec;}
  TVector3 GetIndividualRPVector()              {return frpv;}
  void     SetIndividualRPVector2(TVector3 vec) {frpv2 = vec;}
  TVector3 GetIndividualRPVector2()             {return frpv2;}
  void     SetRotatedMomentum(TVector3 value)   {frmom = value;
    SetLorentzVector();}
  TVector3 GetRotatedMomentum()                 {return frmom;}
  void     SetLorentzVector();
  TLorentzVector GetLorentzVector()             {return fLzvec;}


  Int_t  GetTrackID()                           {return   ftrackID;}
  void   SetTrackID(Int_t ival)                 {ftrackID = ival;}

  // good track flag
private:
  void   SetGoodTrackFlag();

public:
  Int_t  GetGoodTrackFlag()                     {return fgoodtrackf;}

  void   SetVertexAtTargetFlag(Int_t value)     
  { fVatTargetf = value;     SetGoodTrackFlag(); }  
  Int_t  GetVertexAtTargetFlag()                {return fVatTargetf;}

  void   SetDistanceAtVertexFlag(UInt_t value)  
  { fdistanceatvertexf = value;  SetGoodTrackFlag(); }
  UInt_t GetDistanceAtVertexFlag()              {return fdistanceatvertexf;}

  void   SetMomentumFlag(UInt_t value)          {fmomentumf = value; SetGoodTrackFlag(); }
  void   SetPIDFlag(UInt_t value)               {fpidf = value; SetGoodTrackFlag();}

  void   SetNCLFlag(UInt_t value)               {fnclf = value; SetGoodTrackFlag();}
  UInt_t GetNCLFlag()                           {return fnclf;}

  void   SetMassFlag(UInt_t value)              {fmassf = value; SetGoodTrackFlag(); }
  UInt_t GetMassFlag()                          {return fmassf;}

  void   SetdEdxFlag(UInt_t value)              {fdedxf = value; SetGoodTrackFlag(); }
  UInt_t GetdEdxFlag()                          {return fdedxf;}

  void   SetNDFFlag(UInt_t value)               {fNDFf = value; SetGoodTrackFlag();} 
  //  UInt_t GetNDFFlag()                           {return fNDFf;}

  void   SetDoubleFlag(UInt_t value)            {fdoublef = value; SetGoodTrackFlag();} 
  UInt_t GetDoubleFlag()                        {return fdoublef;}

  // --end
  void     SetRPWeight(Double_t value)  {fwgt = value;}
  Double_t GetRPWeight()                {return fwgt;}

  void     SetAzmAngle_wrt_RP(Double_t val)  {fdeltphi = val;}
  Double_t GetAzmAngle_wrt_RP()           {return fdeltphi;}
  void     SetAzmAngle2_wrt_RP(Double_t val) {fdeltphi2 = val;}
  Double_t GetAzmAngle2_wrt_RP()          {return fdeltphi2;}

  void     SetIndividualRPAngle(Double_t val) {frpphi = val;}
  Double_t GetIndividualRPAngle()       {return frpphi;}

  void     SetIndividualRPAngle2(Double_t val) {frpphi2 = val;}
  Double_t GetIndividualRPAngle2()       {return frpphi2;}

  void     SetReactionPlaneFlag(Int_t value)    {fReactionPlanef = value;}
  void     AddReactionPlaneFlag(Int_t value)    {fReactionPlanef += value;}
  Int_t    GetReactionPlaneFlag()               {return fReactionPlanef;}
  

  void     SetMomentumAtTarget(TVector3 value)  {forigP3 = value;}
  TVector3 GetMomentumAtTarget()                {return forigP3;}

  Int_t    GetCharge()                          {return fChar;}
  Int_t    GetGFCharge()                        {return fGFChar;}

  STRecoTrack* GetRecoTrack()                   {return fRTrack;}

  void         SetVertexID(Int_t val)           {rVertexID = val;}
  Int_t        GetVertexID()                    {return rVertexID;}

  Int_t        GetdEdxPointSize()               {return rdEdxPointSize;}
  void         SetdEdxPointSizeCut(Int_t value) {rdEdxPointSize_thr = value;}

  Int_t        GetdEdxPointSizeFlag()            
  {
    if(rdEdxPointSize < rdEdxPointSize_thr)       frdEdxPointSizef = 0;
    return frdEdxPointSizef;
  }

  void         SetNDF(Int_t val)                  {rNDF = val;} 
  Int_t        GetNDF()                           {return rNDF;}
  Int_t        GetNDFvertex()                     {return fNDF;}
  Double_t     GetDistanceAtVertex()              {return rDist;}

  void         SetVertex(TVector3 value);
  TVector3     GetVertex()                        { return fvertex;}
  TVector3     GetPOCAVertex()                    {return rPOCAVertex;}


  ClassDef(STKParticle, 2)

};


#endif
