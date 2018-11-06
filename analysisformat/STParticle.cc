/** 
 * STParticle Class
 *
 * @author Mizuki
 */

#include "STParticle.hh"

#include <iostream>

ClassImp(STParticle)

STParticle::STParticle() : ftrackID(-1)
{
  Clear();
  Initialize();
}


STParticle::STParticle(const STParticle &cp)
{

  ftrackID       = cp.ftrackID;;
  fvertex        = cp.fvertex;

  forigP3       = cp.forigP3;
  fRotatedP3    = cp.fRotatedP3;
  fRotatedPt    = cp.fRotatedPt;
  fP            = cp.fP;
  fdEdx         = cp.fdEdx;
  ftheta        = cp.ftheta;
  frtheta       = cp.frtheta;
  fphi          = cp.fphi;
  frphi         = cp.frphi;
  frpphi        = cp.frpphi;
  fwgt          = cp.fwgt;
  fPID          = cp.fPID;
  fPIDProbability = cp.fPIDProbability;  
  fNDF          = cp.fNDF;

  fRapidity     = cp.fRapidity;
  fRapiditycm   = cp.fRapiditycm;
  fpsudoRapidity= cp.fpsudoRapidity;
  fEtotal       = cp.fEtotal;
  fChar         = cp.fChar;
  fMass         = cp.fMass;
  fBBMass       = cp.fBBMass;

  fpipid        = cp.fpipid;

  //flags
  fBeamonTargetf   = cp.fBeamonTargetf;
  fVBDCCorf        = cp.fVBDCCorf;
  fBDCCorf         = cp.fBDCCorf;
  fTargetXYf       = cp.fTargetXYf;
  fgotoKatanaf     = cp.fgotoKatanaf;
  fgotoKyotof      = cp.fgotoKyotof;
  frdEdxPointSizef = cp.frdEdxPointSizef;

  fVatTargetf      = cp.fVatTargetf;
  fVZatTargetf     = cp.fVZatTargetf;
  fdistanceatvertexf = cp.fdistanceatvertexf;
  fNDFf            = cp.fNDFf;
  fmaxmomentumf    = cp.fmaxmomentumf;
  fmaxthetaf       = cp.fmaxthetaf;
  fmaxdedxf        = cp.fmaxdedxf;
  
  fgoodtrackf      = cp.fgoodtrackf;
  fReactionPlanef  = cp.fReactionPlanef;

  
  //mixed event
  fmxevt           = cp.fmxevt;
  fmxntrk          = cp.fmxntrk;
  fmxtrackid       = cp.fmxtrackid;


  //STRecoTrack parameters
  fRTrack          = cp.fRTrack;
  rVertexID        = cp.rVertexID;  
  rdEdxPointSize   = cp.rdEdxPointSize;
  rdEdxPointSize_thr = cp.rdEdxPointSize_thr;
  rNDF             = cp.rNDF;
  rDist            = cp.rDist;
  rPOCAVertex      = cp.rPOCAVertex;


  Initialize();

}

STParticle &STParticle::operator=(const STParticle &cp)
{

  if( this != &cp )
    *this = cp;

  return *this;
}


void STParticle::Initialize()
{
  gcutHe3BBmass = new TCutG("gHe3",6);
  gcutHe3BBmass->SetVarX("fP");
  gcutHe3BBmass->SetVarY("fBBMass");
  gcutHe3BBmass->SetPoint(0,1657.52,  4754.8);
  gcutHe3BBmass->SetPoint(1,1264.57,  3668.17);
  gcutHe3BBmass->SetPoint(2,1037.6,   3248.5);
  gcutHe3BBmass->SetPoint(3,140,      2550);
  gcutHe3BBmass->SetPoint(4,155,      3130);
  gcutHe3BBmass->SetPoint(5,1657.52,  4754.8);


  gcutHe4BBmass = new TCutG("gHe4",7);
  gcutHe4BBmass->SetVarX("fP");
  gcutHe4BBmass->SetVarY("fBBMass");
  gcutHe4BBmass->SetPoint(0,2165.65, 6965.53);
  gcutHe4BBmass->SetPoint(1,1986.11, 5968.82);
  gcutHe4BBmass->SetPoint(2,1654.13, 4762.29);
  gcutHe4BBmass->SetPoint(3,153.455, 3136.09);
  gcutHe4BBmass->SetPoint(4,248.306, 3825.54);
  gcutHe4BBmass->SetPoint(5,1908.2,  6755.7);
  gcutHe4BBmass->SetPoint(6,2162.26, 6973.02);
}

void STParticle::CheckTrackonTarget()
{
  // Track XY
  Double_t trktgt_right =  -10.2; //!
  Double_t trktgt_left  =   16.2; //!
  Double_t trktgt_top   = -210.; //!
  Double_t trktgt_btm   = -235.; //!


  if(fvertex.X() >= trktgt_right && fvertex.X() <= trktgt_left &&
     fvertex.Y() >= trktgt_btm   && fvertex.Y() <= trktgt_top)
    fTargetXYf = 1;
  else
    fTargetXYf = 0;
}


void STParticle::Clear(Option_t *option)
{

  fRotatedP3 = TVector3(-9999,-9999,-9999);
  forigP3    = TVector3(-9999,-9999,-9999);
  fvertex    = TVector3(-9999,-9999,-9999);

  fP            = -9999.;
  fdEdx         = -9999.;
  ftheta        = -10.;
  frtheta       = -10.;
  fphi          = -10.;
  frphi         = -10.;


  fpipid       = 0;
  fPID         = 0;
  fNDF         = 0.;

  // Track quality flag
  fgoodtrackf  = 0;
  fVatTargetf  = 1;   
  fVZatTargetf = 1;   
  fdistanceatvertexf = 1;
  fNDFf        = 1;    
  fmaxmomentumf= 1;    
  fmaxthetaf   = 1;
  fmaxdedxf    = 1;

  // for flow
  ffltnP3 = TVector3(-9999,-9999,-9999);
  ffltnPt = TVector2(-9999,-9999);

  frpphi   = -10.;
  fdeltphi = -10.;
  fwgt   = 0.;
  fcorrBin[0] = -1;  
  fcorrBin[1] = -1;
  
  fmxevt = -1;
  fmxntrk = -1;

  fcorrPt = ffltnPt;

  fReactionPlanef = 0;

  rChi2   = 0.;
  fBBMass = 0.;
}


void STParticle::SetRecoTrack(STRecoTrack *atrack)
{
  fRTrack = atrack;

  forigP3 = fRTrack->GetMomentumTargetPlane();

  //  forigP3 = fRTrack->GetMomentum();  // modified on 9 Nov. 2017
  //Because of PZ bug for v1.04
  if(forigP3.Z() < 0)
    forigP3 = -forigP3;

  fphi  = forigP3.Phi();
  ftheta= forigP3.Theta();
  fP    = forigP3.Mag();
  fdEdx = fRTrack->GetdEdxWithCut(0, 0.7);
  fChar = fRTrack->GetCharge();

  fPID  = STPID::GetPDG(fRTrack->GetPID());
  fPIDProbability = fRTrack->GetPIDProbability();

  fRotatedP3 = forigP3;

  rVertexID      =  fRTrack -> GetVertexID();
  rdEdxPointSize = (fRTrack -> GetdEdxPointArray()) -> size();
  rNDF           =  fRTrack -> GetNDF();
  rPOCAVertex    =  fRTrack -> GetPOCAVertex();
  rChi2          =  fRTrack -> GetChi2();

  SetMass();
  SetProperty();
}




void STParticle::SetProperty()
{
  
  //  SetLinearPID();
  //  CheckTrackonTarget();
  CheckKATANAHit();
  CheckKYOTOHit();

  SetPiPID();

  
}

void STParticle::SetVertex(STVertex *value)  
{
  fNDF = value->GetNDF();
  SetVertex( value->GetPos() );
}

void STParticle::SetVertex(TVector3 value)  
{
  fvertex = value;
  rDist = (rPOCAVertex - fvertex).Mag();
}

void STParticle::SetPID(Int_t value)
{
  fPID = value;

  SetMass();

  SetRapidity();
}



void STParticle::SetMass()
{
  auto mpi  =    139.57018;
  auto mp   =    938.2720813;
  auto mn   =    939.565346;
  auto md   =   1875.612762;
  auto mt   =   2808.921112;
  auto mhe3 =   2808.39132;
  auto mhe4 =   3727.379378;


  Double_t mass = 0;
  switch (fPID) {
  case 12212:
  case 2212:
    mass = mp;
    break;

  case 211:
  case -211:
    mass = mpi;
    break;

  case 1000010020:
    mass = md;
    break;

  case 1000010030:
    mass = mt;
    break;

  case 1000020030:
    fChar = 2;
    mass = mhe3;
    break;

  case 1000020040:
    fChar = 2;
    mass = mhe4;
    break;

  default:
    mass = 0;
    break;
  }

  fMass = mass;

}


void STParticle::SetpsudoRapidity()
{
  fpsudoRapidity = -log( tan((fRotatedP3.Theta()/2)) );
}

 
void STParticle::SetRapidity()
{
  Double_t P    = fRotatedP3.Mag();
  Double_t Pz   = fRotatedP3.Z();


  if(fMass != 0 ){
    fEtotal   = sqrt(fMass*fMass + P*P);
    fRapidity = 0.5 * log( (fEtotal + Pz) / (fEtotal - Pz) );
  }
  else {
    fEtotal = 0;
    fRapidity = -10.;
  }

  SetpsudoRapidity();
}


void  STParticle::RotateAlongBeamDirection(Double_t valuex, Double_t valuey)
{

  fRotatedP3.RotateY(-valuex);
  fRotatedP3.RotateX(-valuey);

  fRotatedPt = TVector2(fRotatedP3.X(),fRotatedP3.Y());

  fcorrPt = fRotatedPt;

  frphi   = fRotatedP3.Phi();
  frtheta = fRotatedP3.Theta();

  bRotated = kTRUE;
}

void STParticle::Flattening(Double_t value)
{

  ffltnP3  = fRotatedP3;

  ffltnP3.SetPhi(value);

  ffltnPt  = TVector2(ffltnP3.X(), ffltnP3.Y());

  fcorrPt  = ffltnPt;

  frphi    = ffltnP3.Phi();

  frtheta  = ffltnP3.Theta();

  bFlatten = kTRUE;
}


void STParticle::SetBetheBlochMass(Double_t *para)
{
  fitterpara[0] = para[0];
  fitterpara[1] = para[1];
  
  fitterpara[2] = fChar * fP;
  fitterpara[3] = fChar;
  fitterpara[4] = fdEdx;
  
  fBBMass = BetheBlochFitter::CalcMass(fitterpara);

  SetBBPID();
}

void STParticle::SetBBPID()
{


  for(UInt_t i = 0; i < nBBsize; i++){

    Double_t mass_low = BBmassRegion[i][0]-BBmassRegion[i][1]*BBmassRegion[i][2] ;
    Double_t mass_up  = BBmassRegion[i][0]+BBmassRegion[i][1]*BBmassRegion[i][2] ;

    Bool_t triton_cut = 1;
    if( i == 3 && fBBMass > 0.8*fP+2400. ) triton_cut = 0;
      

    if( fBBMass >= mass_low &&	fBBMass <= mass_up && triton_cut) {

      STPID::PID pid = static_cast<STPID::PID>(i);
      fPID = STPID::GetPDG(pid);

      SetMass();
      SetRapidity();
      return;
    }
  }

  // extended p
  if( fBBMass > BBmassRegion[1][0]+BBmassRegion[1][1]*BBmassRegion[1][2] &&
      fBBMass < BBmassRegion[2][0]-BBmassRegion[2][1]*BBmassRegion[2][2] ) {

    STPID::PID pid = static_cast<STPID::PID>(1);
    fPID = STPID::GetPDG(pid);
  }    
  // exptended d
  if( fBBMass > BBmassRegion[2][0]+BBmassRegion[2][1]*BBmassRegion[2][2] &&
      fBBMass < BBmassRegion[3][0]-BBmassRegion[3][1]*BBmassRegion[3][2] ) {

    STPID::PID pid = static_cast<STPID::PID>(2);
    fPID = STPID::GetPDG(pid);
  }    


  if( gcutHe3BBmass == NULL || gcutHe4BBmass == NULL ) return;

  if( gcutHe3BBmass->IsInside(fP, fBBMass) ) {
    STPID::PID pid = static_cast<STPID::PID>(4);
    fPID = STPID::GetPDG(pid);
  }
  
  if( gcutHe4BBmass->IsInside(fP, fBBMass) ) {
    STPID::PID pid = static_cast<STPID::PID>(5);
    fPID = STPID::GetPDG(pid);
  }

  SetMass();

}


void STParticle::SetPiPID()
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
