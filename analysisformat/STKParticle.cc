/**
 * STKParticle Class
 *
 * @author Mizuki
 */

#include "STKParticle.hh"

#include <iostream>

ClassImp(STKParticle)

STKParticle::STKParticle() : ftrackID(-1)
{
  Clear();
  Initialize();
}


STKParticle::STKParticle(const STKParticle &cp)
{

  ftrackID      = cp.ftrackID;;
  fPID          = cp.fPID;
  fPID_seq      = cp.fPID_seq;
  fPID_tight    = cp.fPID_tight;
  fPID_loose    = cp.fPID_loose;
  fpipid        = cp.fpipid;
  fChar         = cp.fChar;
  fGFChar       = cp.fGFChar;
  fMass         = cp.fMass;
  fBBMass       = cp.fBBMass;
  fBBMassHe     = cp.fBBMassHe;
  fPIDProbability = cp.fPIDProbability;  

  fvertex       = cp.fvertex;
  forigP3       = cp.forigP3;
  frmom         = cp.frmom;
  fLzvec        = cp.fLzvec;
  fP            = cp.fP;
  fdEdx         = cp.fdEdx;
  fNDF          = cp.fNDF;
  fncl          = cp.fncl;
  fRapidity     = cp.fRapidity;
  fRapiditycm   = cp.fRapiditycm;
  fPt           = cp.fPt;
  fPhi          = cp.fPhi;

  fwgt          = cp.fwgt;
  frpv          = cp.frpv;
  frpphi        = cp.frpphi;
  fdeltphi      = cp.fdeltphi;
  frpv2         = cp.frpv2;
  frpphi2       = cp.frpphi2;
  fdeltphi2     = cp.fdeltphi2;


  //flags
  fVatTargetf      = cp.fVatTargetf;
  fdistanceatvertexf = cp.fdistanceatvertexf;
  fNDFf            = cp.fNDFf;
  fnclf            = cp.fnclf;
  fdedxf           = cp.fdedxf;
  fmassf           = cp.fmassf;
  fmomentumf       = cp.fmomentumf;  
  fpidf            = cp.fpidf;
  fgoodtrackf      = cp.fgoodtrackf;
  fReactionPlanef  = cp.fReactionPlanef;
  frdEdxPointSizef = cp.frdEdxPointSizef;
  fdoublef         = cp.fdoublef;

  //STRecoTrack parameters
  fRTrack           = cp.fRTrack;
  fVATrack          = cp.fVATrack;
  rVertexID         = cp.rVertexID;  
  rHelixID          = cp.rHelixID;
  rNDF              = cp.rNDF;
  rncl              = cp.rncl;
  rDist             = cp.rDist;
  rPOCAVertex       = cp.rPOCAVertex;
  rChi2             = cp.rChi2;
  rdEdxPointSize    = cp.rdEdxPointSize;
  rdEdxPointSize_thr= cp.rdEdxPointSize_thr;
  Initialize();
}

STKParticle &STKParticle::operator=(const STKParticle &cp)
{

  if( this != &cp )
    *this = cp;

  return *this;
}


void STKParticle::Initialize()
{
}

void STKParticle::Clear(Option_t *option)
{
  fRTrack  = NULL;
  fVATrack = NULL;

  ftrackID     = -1;
  fPID         = 0;
  fPID_seq     = 99;
  fPID_tight   = 0;
  fPID_loose   = 0;
  fpipid       = 0;
  fChar        = 0;
  fGFChar      = 0;
  fMass        = 0.;
  fBBMass      = 0.;
  fBBMassHe    = 0.;
  fPIDProbability = 0.;

  fvertex      = TVector3(-9999,-9999,-9999);
  forigP3      = TVector3(-9999,-9999,-9999);
  frmom        = TVector3(-9999,-9999,-9999);
  fLzvec       = TLorentzVector(0.,0.,0.,0.);
  fP           = -9999.;
  fdEdx        = -9999.;
  fNDF         = 0;
  fncl         = 0;
  fRapidity    = -9.;
  fRapiditycm  = -9.;
  fPt          = -1.;
  fPhi         = -5.;

  fwgt   = 0.;
  frpv     = TVector3(-9999,-9999,-9999);
  frpphi   = -5.;
  fdeltphi = -5.;
  frpv2    = TVector3(-9999,-9999,-9999);
  frpphi2  = -5.;
  fdeltphi2= -5.;


  // Track quality flag
  fVatTargetf  = 1;   
  fdistanceatvertexf = 1;
  fnclf            = 1;
  fmassf           = 1;
  fmomentumf       = 1;    
  fpidf            = 1;
  fdoublef         = 1;
  fgoodtrackf      = 1;

  fNDFf            = 1;    
  fdedxf           = 1;
  fReactionPlanef  = 1;
  frdEdxPointSizef = 1;

  rVertexID      =  0;
  rHelixID       =  0;
  rdEdxPointSize =  0;
  rNDF           =  0;
  rncl           =  0;
  rDist          =  0.;
  rPOCAVertex    =  TVector3(-9999,-9999,-9999);  
  rChi2          =  0.;
  rdEdxPointSize =  0;
  rdEdxPointSize_thr = 0.;
}

void STKParticle::SetGoodTrackFlag()
{
  fgoodtrackf = fgoodtrackf != 0 ?  
    1     * fdoublef    +
    10    * fmomentumf  + 
    100   * fmassf      +
    1000  * fnclf       +
    10000 * fVatTargetf +
    100000* fdistanceatvertexf
    : 0;
}


void STKParticle::SetRecoTrack(STRecoTrack *atrack, STRecoTrack *atrackva)
{    
  Clear();

  fRTrack  = atrack;
  fVATrack = atrackva;

  fdEdx   = fVATrack -> GetdEdxWithCut(0, 0.7);
  forigP3 = fVATrack -> GetMomentumTargetPlane(); //v56
  fncl    = fVATrack -> GetClusterIDArray() -> size();
  frmom   = forigP3;

  fP    = forigP3.Mag();
  fChar = fVATrack->GetCharge();
  fGFChar = fVATrack->GetGenfitCharge();

  //  fPID  = STPID::GetPDG(fVATrack->GetPID());
  //  fPIDProbability = fVATrack->GetPIDProbability();

  rVertexID      =  fRTrack -> GetVertexID();
  rHelixID       =  fRTrack -> GetHelixID();
  rdEdxPointSize =  fRTrack -> GetdEdxPointArray() -> size();
  rNDF           =  fRTrack -> GetNDF();
  rPOCAVertex    =  fRTrack -> GetPOCAVertex();
  rChi2          =  fRTrack -> GetChi2();
  rncl           =  fRTrack -> GetClusterIDArray() -> size();

  // quality flag
  //  if( fdEdx > maxdEdx  || fdEdx <= 0 ) SetdEdxFlag(0);
  //  if( fGFChar != fChar )               SetMassFlag(0);

  // vertex pos
  Bool_t bgtrack = 
    (abs(fRTrack->GetPOCAVertex().Z() + 14.85 ) <= 3.* 1.33) &&
    (abs(fRTrack->GetPOCAVertex().X()) <= 15.) &&
    (abs(fRTrack->GetPOCAVertex().Y() + 205.) <= 20. );
  if( !bgtrack ) 
    SetVertexAtTargetFlag(0);
  

  SetGoodTrackFlag();
  
}



void STKParticle::SetProperty()
{
  
  //  SetLinearPID();
  //  CheckTrackonTarget();
  CheckKATANAHit();
  CheckKYOTOHit();

  SetPiPID();

  
}


void STKParticle::SetVertex(TVector3 value)  
{
  fvertex = value;
  rDist = (rPOCAVertex - fvertex).Mag();
  
  if( rDist > 20 ) SetDistanceAtVertexFlag(0);
  
}

void STKParticle::SetBBMass(Double_t val)       
{
  if( val > 0 ) 
    fBBMass   = val;
}
void STKParticle::SetBBMassHe(Double_t val)       
{
  if( val > 0 && val < 1.e+8) 
    fBBMassHe   = val;
}


void STKParticle::SetMass(Int_t val)
{
  auto mpi  =    139.57018;
  auto mp   =    938.2720813;
  auto mn   =    939.565346;
  auto md   =   1875.612762;
  auto mt   =   2808.921112;
  auto mhe3 =   2808.39132;
  auto mhe4 =   3727.379378;

  Double_t mass = 0;
  switch (val) {
  case -211:
    mass = mpi;
    fPID_seq = 0;
    break;
  case 211:
    mass = mpi;
    fPID_seq = 1;
    break;
  case 2212:
    mass = mp;
    fPID_seq = 2;
    break;
  case 1000010020:
    mass = md;
    fPID_seq = 3;
    break;
  case 1000010030:
    mass = mt;
    fPID_seq = 4;
    break;
  case 1000020030:
    mass = mhe3;
    frmom.SetMag(fP * 2.);
    fPID_seq = 5;    
    break;
  case 1000020040:
    frmom.SetMag(fP * 2.);
    mass = mhe4;
    fPID_seq = 6;    
    break;

  default:
    mass = 0;
    fPID_seq = 99;
    break;
  }

  fMass = mass;

  if( fMass > 0 ) 
    SetLorentzVector();

  else {
    SetMassFlag(0);
    SetPIDFlag(0);
  }

}


void  STKParticle::SetLorentzVector()
{
  Double_t fEtotal     = sqrt(fMass*fMass + frmom.Mag2());

  fLzvec.SetVect(frmom);
  fLzvec.SetT(fEtotal);

  fRapidity   = fLzvec.Rapidity();

}

void  STKParticle::RotateAlongBeamDirection(Double_t valuex, Double_t valuey)
{
  TVector3 beamDirection(TMath::Tan(valuex), TMath::Tan(valuey), 1.);
  beamDirection = beamDirection.Unit();
  auto rotationAxis  = beamDirection.Cross(TVector3(0,0,1));
  auto rotationAngle = beamDirection.Angle(TVector3(0,0,1));

  frmom.Rotate(rotationAngle, rotationAxis);

  fPhi = frmom.Phi()*TMath::RadToDeg();
  fPt  = frmom.Perp();

}


void STKParticle::SetPiPID()
{
  //pion cut                                                                                                                
  // TFile *gcutPiFile = new TFile("/cache/scr/spirit/mizuki/SpiRITROOT/macros/AnalysisMacro/gcutPip.root");
  // TCutG *gPip = (TCutG*)gcutPiFile->Get("gPi");
  // gcutPiFile->Close();

  // if(gPip->IsInside(fdEdx,fP))
  //   fpipid = 1;
  // else
  //   fpipid = 0;

}


