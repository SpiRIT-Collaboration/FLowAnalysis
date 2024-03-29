/**     
 * @brief STTrack Class         
 *                                   
 * @author Mizuki         
 */

#ifndef STPARTICLE
#define STPARTICLE

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

class STParticle : public TObject {
public:
  STParticle();
  STParticle(const STParticle &cp);
  STParticle &operator=(const STParticle &cp);

  virtual ~STParticle(){};

private:
  Bool_t    bRotated = kFALSE;
  Bool_t    bFlatten = kFALSE;

  Int_t    ftrackID;

  UInt_t   fPID;
  UInt_t   fPID_seq;
  UInt_t   fPID_tight;
  UInt_t   fPID_norm;
  UInt_t   fPID_loose;
  Double_t fRapidity;
  Double_t fRapiditycm;
  Double_t fEtotal;
  Int_t    fChar;
  Int_t    fGFChar;
  Double_t fMass;
  Double_t fBBMass;
  Double_t fBBMassHe;

  Int_t    fpipid;
  TVector3 fvertex;
  Double_t fPIDProbability;

  TVector3 forigP3;;         // Momemtum vector without any correction.
  Double_t fP;               // Momentum/Q [MeV/c]
  Double_t fdEdx;            // dEdx
  Double_t fNDF;             // STVertex::GetNDF()
  Int_t    fnclust;          // number of cluster
  Double_t fclustex;         // expected cluster number
  Double_t fclusterratio;    // rClusterSize/fclustex

  // for flow analysis 
  TVector3 fRotatedP3;       // Momentum vector rotated with respect to the beam angle.
  TVector2 fRotatedPt;       // Transverse momentum vector rotated with respect to the beam angle. 
  TVector2 fPxz;             // TVector2( fRotatedP3.X(), fRotatedP3.Z() )
  TVector2 fPyz;             // TVector2( fRotatedP3.Y(), fRotatedP3.Z() )

  TLorentzVector fLzvec;     // LorentzVector flzvec(fRotatedP3, Etotal);

  Double_t fwgt;             // Summing up weight
  TVector3 frpv;             // Individual Reaction plane vector (IRPO) for each particle
  Double_t frpphi;           // Individual Reaction plane orientaion (IRPO) for each particle
  Double_t fdeltphi;         // Azimuthal opening angle with respect to IRPO
  TVector3 frpv2;            // Individual v2 Reaction plane vector (IRPO2) for each particle
  Double_t frpphi2;          // Individual v2 Reaction plane orientaion (IRPO2) for each particle
  Double_t fdeltphi2;        // Azimuthal opening angle with respect to IRPO2

  // --- mixed partiels
  Int_t    fmxevt  = -1;
  Int_t    fmxntrk = -1;
  Int_t    fmxtrackid = -1;

  //quality flags
  UInt_t   fBeamonTargetf   ; //flag for beam tracked by BDC goes on the target
  UInt_t   fVBDCCorf        ; //flag for reconstructed vertex is correated with BDC at the target
  UInt_t   fBDCCorf         ; //
  UInt_t   fTargetf         ;
  UInt_t   fgotoKatanaf     ;
  UInt_t   fgotoKyotof      ;
  UInt_t   frdEdxPointSizef ;  

  UInt_t   fVatTargetf      ; //flag for reconstructed vertex XY within the target
  UInt_t   fdistanceatvertexf;
  UInt_t   fNDFf            ;
  UInt_t   fmomentumf    ;
  UInt_t   fdedxf           ;
  UInt_t   fnclustf         ; // number of cluster flag
  UInt_t   fclusterratiof   ; // cluster ratio flag
  UInt_t   fmassf           ;

  UInt_t   fgoodtrackf      ; // good track flag
  UInt_t   fReactionPlanef  ; // excellent track for flow



  //STRecoTrack parameters
  STRecoTrack *fRTrack; //!
  STRecoTrack *fVATrack; //!
  Int_t     rVertexID;  
  Int_t     rdEdxPointSize;
  Int_t     rdEdxPointSize_thr = 1;
  Int_t     rNDF;
  Double_t  rDist;
  TVector3  rPOCAVertex;  
  Double_t  rChi2;
  Int_t     rClusterSize;

  Double_t  maxdEdx     = 2000.;  //!
  Double_t  maxMomentum = 4000.;  //!
  Double_t  minMomentum = 20.;    //!
  Double_t  maxbbmass   = 20000.; //!

private:
  virtual void Clear(Option_t *option = "");
  
  void     Initialize();

  void     SetProperty();
  void     SetRotatedPt();                       

  void     CheckKATANAHit(){};
  void     CheckKYOTOHit(){};


public:
  void     Clean()                       {Clear();}
  void     SetRecoTrack(STRecoTrack *atrack, STRecoTrack *atrackva=NULL);

  void     SetNumCluster(STRecoTrack *atrack) 
  {        fnclust = atrack->GetClusterIDArray()->size();}
  Int_t    GetNumCluster()              {return rClusterSize;}

  void     RotateAlongBeamDirection(Double_t valuex, Double_t valuey);
  void     SetP(Double_t value)          {fP = value;}
  void     SetPiPID();
  Int_t    GetPiPID()                    {return fpipid;}

  //  void     SetPID();
  void     SetPID(Int_t value)           { fPID = value; SetMass(fPID);}
  Int_t    GetPID()                      {return fPID;}
  void     SetPIDTight(Int_t value)      {fPID_tight = value;}           
  Int_t    GetPIDTight()                 {return fPID_tight;}           
  void     SetPIDNorm(Int_t value)       {fPID_norm = value; SetMass(fPID_norm);}           
  Int_t    GetPIDNorm()                  {return fPID_norm;}           
  void     SetPIDLoose(Int_t value)      {fPID_loose = value;}
  Int_t    GetPIDLoose()                 {return fPID_loose;}           

  Int_t    GetPID_seq()                  {return fPID_seq;}


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

  Double_t GetEtotal()                   {return fEtotal;}

  Double_t GetP()                        {return fP;}
  void     SetdEdx(Double_t value)       {fdEdx = value;}
  Double_t GetdEdx()                     {return fdEdx;}


  void     SetIndividualRPVector(TVector3 vec)  {frpv = vec;}
  TVector3 GetIndividualRPVector()              {return frpv;}
  void     SetIndividualRPVector2(TVector3 vec) {frpv2 = vec;}
  TVector3 GetIndividualRPVector2()             {return frpv2;}
  void     SetRotatedMomentum(TVector3 value)   {fRotatedP3 = value; SetRotatedPt(); SetLorentzVector();}
  TVector3 GetRotatedMomentum()                 {return fRotatedP3;}
  TVector2 GetRotatedPt()                       {return fRotatedPt;}
  void     SetLorentzVector();
  TLorentzVector GetLorentzVector()             {return fLzvec;}

  Double_t GetYawAngle()                        {return TVector2::Phi_mpi_pi(fPxz.Phi());}
  Double_t GetPitchAngle()                      {return TVector2::Phi_mpi_pi(fPyz.Phi());}

  void     SetExpectedClusterNumber(Double_t val) 
  {
    fclustex = val; 
    if( fclustex > 0 )
      fclusterratio = (Double_t)rClusterSize / fclustex;
  }

  Int_t  GetTrackID()                           {return   ftrackID;}
  void   SetTrackID(Int_t ival)                 {ftrackID = ival;}


  void  SetVertexBDCCorrelationFlag(Int_t value){fVBDCCorf   = value;}
  Int_t GetVertexBDCCorrelationFlag()           {return fVBDCCorf;}

  void  SetBDCCorrelationFlag(Int_t value)      {fBDCCorf     = value;}
  Int_t GetBDCCorrelationFlag()                 {return fBDCCorf;}

  void  SetgotoKATANAFlag(Int_t value)          {fgotoKatanaf   = value;}
  Int_t GetgotoKATANAFlag()                     {return fgotoKatanaf;}

  void  SetgotoKYOTOFlag(Int_t value)           {fgotoKyotof   = value;}
  Int_t GetgotoKYOTOFlag()                      {return fgotoKyotof;}

  // good track flag
private:
  void   SetGoodTrackFlag();

public:
  Int_t  GetGoodTrackFlag()                     {return fgoodtrackf;}
  //  void   SetGoodTrackFlag(Int_t value)          {fgoodtrackf  = fgoodtrackf != 0 ? value : 0;}

  void  SetBeamonTargetFlag(Int_t value)        {fBeamonTargetf = value; SetGoodTrackFlag(); }
  Int_t GetBeamonTargetFlag()                   {return fBeamonTargetf;}

  void   SetVertexAtTargetFlag(Int_t value)     
  { fVatTargetf = value; 
    fTargetf = fVatTargetf*fdistanceatvertexf;
    SetGoodTrackFlag();
  }  
  Int_t  GetVertexAtTargetFlag()                {return fVatTargetf;}

  void   SetDistanceAtVertexFlag(UInt_t value)  
  { fdistanceatvertexf = value; 
    fTargetf = fVatTargetf*fdistanceatvertexf;
    SetGoodTrackFlag();
  }
  UInt_t GetDistanceAtVertexFlag()              {return fdistanceatvertexf;}
  // void   SetFromTargetFlag(Int_t value)         {fTargetf  = value;}
  // Int_t  GetFromTargetFlag()                    {return fTargetf;}

  void   SetdEdxFlag(UInt_t value)              {fdedxf = value; SetGoodTrackFlag(); }
  UInt_t GetdEdxFlag()                          {return fdedxf;}

  void   SetMomentumFlag(UInt_t value)          {fmomentumf = value; SetGoodTrackFlag(); }
  //  UInt_t GetMomentumFlag()                      {return fmomentumf; }

  void   SetNDFFlag(UInt_t value)               {fNDFf = value; SetGoodTrackFlag();} 
  //  UInt_t GetNDFFlag()                           {return fNDFf;}
  void   SetNumClusterFlag(UInt_t value)        {fnclustf = value; SetGoodTrackFlag();}
  UInt_t GetNumClusterFlag()                    {return fnclustf;}
  

  void   SetClusterRatioFlag(UInt_t value)      {fclusterratiof = value;}
  UInt_t GetClusterRatioFlag()                  {return fclusterratiof;}

  void   SetMassFlag(UInt_t value)              {fmassf = value; SetGoodTrackFlag(); }
  //  UInt_t GetMassFlag()                          {return fmassf;}

  // --end

  // for flow analysis
  void     SetMixedEventID(Int_t value) {fmxevt = value;}
  void     SetMixedNtrack(Int_t value)  {fmxntrk= value;}
  void     SetMixTrackID(Int_t ival)    {fmxtrackid = ival;}

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

  TVector3 GetMomentum()                        {return fRotatedP3;}

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

  void         SetClustNum(Int_t value)          {fnclust = value;}
  Int_t        GetClustNum()                     {return fnclust;}

  void         SetNDF(Int_t val)                  {rNDF = val;} 
  Int_t        GetNDF()                           {return rNDF;}
  Int_t        GetNDFvertex()                     {return fNDF;}
  void         SetDistanceAtVertex(Double_t val)  {rDist = val;}
  Double_t     GetDistanceAtVertex()              {return rDist;}

  Double_t     GetClusterRatio() {return fclusterratio;}


  void         SetVertex(TVector3 value);
  void         SetVertex(STVertex *value);
  TVector3     GetVertex()                      { return fvertex;}
  TVector3     GetPOCAVertex()                  {return rPOCAVertex;}


  ClassDef(STParticle, 22)

};


#endif
