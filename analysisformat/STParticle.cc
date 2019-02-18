/** 111
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
  bRotated      = cp.bRotated;
  bFlatten      = cp.bFlatten;
  
  ftrackID      = cp.ftrackID;;
  fPID          = cp.fPID;
  fRapidity     = cp.fRapidity;
  fRapiditycm   = cp.fRapiditycm;
  fEtotal       = cp.fEtotal;
  fChar         = cp.fChar;
  fGFChar       = cp.fGFChar;
  fMass         = cp.fMass;
  fBBMass       = cp.fBBMass;
  fpipid        = cp.fpipid;
  fvertex       = cp.fvertex;
  fPIDProbability = cp.fPIDProbability;  
  fNDF          = cp.fNDF;
  fLzvec        = cp.fLzvec;

  forigP3       = cp.forigP3;
  fP            = cp.fP;
  fdEdx         = cp.fdEdx;

  fRotatedP3    = cp.fRotatedP3;
  fRotatedPt    = cp.fRotatedPt;

  frpphi        = cp.frpphi;
  fdeltphi      = cp.fdeltphi;
  fwgt          = cp.fwgt;

  fBLMass       = cp.fBLMass;

  //flags
  fBeamonTargetf   = cp.fBeamonTargetf;
  fVBDCCorf        = cp.fVBDCCorf;
  fBDCCorf         = cp.fBDCCorf;
  fTargetf         = cp.fTargetf;
  fgotoKatanaf     = cp.fgotoKatanaf;
  fgotoKyotof      = cp.fgotoKyotof;
  frdEdxPointSizef = cp.frdEdxPointSizef;
  fclusterratiof   = cp.fclusterratiof;
  fVatTargetf      = cp.fVatTargetf;
  fdistanceatvertexf = cp.fdistanceatvertexf;
  fNDFf            = cp.fNDFf;
  fmomentumf       = cp.fmomentumf;
  fdedxf           = cp.fdedxf;

  fmassf           = cp.fmassf;
  
  fgoodtrackf      = cp.fgoodtrackf;
  fReactionPlanef  = cp.fReactionPlanef;

  
  //mixed event
  fmxevt           = cp.fmxevt;
  fmxntrk          = cp.fmxntrk;
  fmxtrackid       = cp.fmxtrackid;


  //STRecoTrack parameters
  fRTrack          = cp.fRTrack;
  fVATrack         = cp.fVATrack;
  rVertexID        = cp.rVertexID;  
  rdEdxPointSize   = cp.rdEdxPointSize;
  rdEdxPointSize_thr = cp.rdEdxPointSize_thr;
  rNDF             = cp.rNDF;
  rDist            = cp.rDist;
  rPOCAVertex      = cp.rPOCAVertex;
  rChi2            = cp.rChi2;
  rClusterSize     = cp.rClusterSize;    
  fclustex         = cp.fclustex;
  fclusterratio    = cp.fclusterratio;

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

void STParticle::Clear(Option_t *option)
{

  fRotatedP3 = TVector3(-9999,-9999,-9999);
  forigP3    = TVector3(-9999,-9999,-9999);
  fvertex    = TVector3(-9999,-9999,-9999);
  fLzvec     = TLorentzVector(0.,0.,0.,0.);

  fP            = -9999.;
  fdEdx         = -9999.;

  fpipid       = 0;
  fPID         = 0;
  fNDF         = 0.;
  fclustex     = -1.;
  fclusterratio= -1.;

  // Track quality flag
  fgoodtrackf  = 1;
  fVatTargetf  = 1;   
  fdistanceatvertexf = 1;
  fTargetf     = 1;
  fmomentumf   = 1;    
  fdedxf       = 1;

  fNDFf        = 1;    
  fclusterratiof = 1;  
  fmassf       = 1;

  // for flow
  fdeltphi = -10.;
  fwgt   = 0.;
  
  fmxevt = -1;
  fmxntrk = -1;


  fReactionPlanef = 0;

  rChi2   = 0.;
  fBBMass = 0.;
  fBLMass = 0.;
  rClusterSize = 0;

}


void STParticle::SetRecoTrack(STRecoTrack *atrack)
{
  fRTrack = atrack;

  forigP3 = fRTrack->GetMomentumTargetPlane();
  fRotatedP3 = forigP3;

  fP    = forigP3.Mag();
  fdEdx = fRTrack->GetdEdxWithCut(0, 0.7);
  fChar = fRTrack->GetCharge();
  fGFChar = fRTrack->GetGenfitCharge();

  fPID  = STPID::GetPDG(fRTrack->GetPID());
  fPIDProbability = fRTrack->GetPIDProbability();

  rVertexID      =  fRTrack -> GetVertexID();
  rdEdxPointSize = (fRTrack -> GetdEdxPointArray()) -> size();
  rNDF           =  fRTrack -> GetNDF();
  rPOCAVertex    =  fRTrack -> GetPOCAVertex();
  rChi2          =  fRTrack -> GetChi2();
  rClusterSize   = (fRTrack -> GetClusterIDArray()) -> size();

}

void STParticle::SetVATrack(STGenfitVATask *atrack)
{
  fVATrack = atrack;

  SetRecoTrack( (STRecoTrack*)fVATrack );

  std::cout << "SetVATrack " << std::endl;
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
    mass = mhe3;
    break;

  case 1000020040:
    mass = mhe4;
    break;

  default:
    mass = 0;
    break;
  }

  fMass = mass;

  SetLorentzVector();
}

void  STParticle::SetLorentzVector()
{
  fEtotal     = sqrt(fMass*fMass + fRotatedP3.Mag2());

  fLzvec.SetVect(fRotatedP3);
  fLzvec.SetT(fEtotal);

  fRapidity   = fLzvec.Rapidity();
    
  if(fMass == 0 ) {
    fEtotal = 0;
    fRapidity = -10.;
    fLzvec.SetT(0.);
  }
}

void  STParticle::RotateAlongBeamDirection(Double_t valuex, Double_t valuey)
{
 
 fRotatedP3.RotateY(-valuex);
  fRotatedP3.RotateX(-valuey);

  fRotatedPt = TVector2(fRotatedP3.X(),fRotatedP3.Y());

  bRotated = kTRUE;
}

void STParticle::GetMasswithMassFitter(TF1 *afunc)
{
  if( afunc == NULL ) return;

  TString f1name = ((TString)afunc->GetName())(2,2);

  Double_t dx = 0.1;

  if( f1name == "BB") {

    Double_t *bbPar = new Double_t[6];
    bbPar[0] = afunc->GetParameter(0);
    bbPar[1] = afunc->GetParameter(1);
    bbPar[2] = 1.; 
    bbPar[3] = fP;
    bbPar[4] = fdEdx;


    auto funcBB  = [bbPar](double x)    { return MassEstimator::BBMassFinderEq(&x, bbPar);};
    auto dfuncBB = [bbPar, dx](double x){ return MassEstimator::BBMassFinderDeriv(&x, bbPar, dx);};
    ROOT::Math::RootFinder finder;  
    ROOT::Math::GradFunctor1D gradfunc1dBB(funcBB, dfuncBB);
    
    finder.SetMethod(ROOT::Math::RootFinder::kGSL_SECANT);
    finder.SetFunction(gradfunc1dBB, 1000.); 
    finder.Solve();
    fBLMass  = finder.Root();

  }

  else {

    Double_t *lvPar = new Double_t[6];
    lvPar[0] = afunc->GetParameter(0);
    lvPar[1] = afunc->GetParameter(1);
    lvPar[2] = afunc->GetParameter(2);
    lvPar[3] = 1.; 
    lvPar[4] = fP;
    lvPar[5] = fdEdx;

    auto funcLV  = [lvPar](double x)    { return MassEstimator::LVMassFinderEq(&x, lvPar);};
    auto dfuncLV = [lvPar, dx](double x){ return MassEstimator::LVMassFinderDeriv(&x, lvPar, dx);};
    ROOT::Math::RootFinder finder;  
    ROOT::Math::GradFunctor1D gradfunc1dLV(funcLV, dfuncLV);
    
    finder.SetMethod(ROOT::Math::RootFinder::kGSL_SECANT);
    finder.SetFunction(gradfunc1dLV, 1000.); 
    finder.Solve();
    fBLMass  = finder.Root();

  }

  if( fBLMass <= 0 || fBLMass > 6000) 
    fmassf = 0;

}

void STParticle::SetBLPID()
{
  fPID = 0;

  if( fmassf ) {

    for(UInt_t i = 0; i < nBBsize; i++){
      
      Bool_t pcut = 1;
      if( i == 3 && fBLMass > fP*1.1+2282. ) pcut = 0;

      Double_t mass_low = LVmassRegion[i][0]-LVmassRegion[i][1]*LVmassRegion[i][2];
      Double_t mass_up  = LVmassRegion[i][0]+LVmassRegion[i][1]*LVmassRegion[i][3];
      
      if( fBLMass >= mass_low && fBLMass <= mass_up && pcut ) {

	STPID::PID pid = static_cast<STPID::PID>(i);
	fPID = STPID::GetPDG(pid);

	break;

    }
    
      

    } 
  }
}


void STParticle::SetBetheBlochMass(Double_t *para)
{
  fitterpara[0] = para[0];
  fitterpara[1] = para[1];
  
  fitterpara[2] = fGFChar * fP;
  fitterpara[3] = fGFChar;
  fitterpara[4] = fdEdx;
  
  fBBMass = BetheBlochFitter::CalcMass(fitterpara);

  SetBBPID();
}

void STParticle::SetBBPID()
{  
  fPID  = 0;

  if( fdEdx <= maxdEdx && fP <= maxMomentum ) { 

    for(UInt_t i = 0; i < nBBsize; i++){

      Bool_t pcut = 1;
      if( i == 1 && fP > protonMaxMomentum ) pcut = 0;
      else if( i == 3 && fBBMass > 0.8*fP+2400. ) pcut = 0;      

      Double_t mass_low = BBmassRegion[i][0]-BBmassRegion[i][1]*BBmassRegion[i][2] ;
      Double_t mass_up  = BBmassRegion[i][0]+BBmassRegion[i][1]*BBmassRegion[i][3] ;

      if( fBBMass >= mass_low && fBBMass <= mass_up && pcut) {
	
	STPID::PID pid = static_cast<STPID::PID>(i);
	fPID = STPID::GetPDG(pid);

	break;
      }
    }
      
    if( gcutHe3BBmass != NULL && gcutHe4BBmass != NULL ) {

      if( fPID == 0 && gcutHe3BBmass->IsInside(fP, fBBMass) ) {
	STPID::PID pid = static_cast<STPID::PID>(4);
	fPID = STPID::GetPDG(pid);

	fRotatedP3.SetMag(fP * 2.);
      }
  
      if( fPID == 0 && gcutHe4BBmass->IsInside(fP, fBBMass) ) {
	STPID::PID pid = static_cast<STPID::PID>(5);
	fPID = STPID::GetPDG(pid);

	fRotatedP3.SetMag(fP * 2.);
      }
    }
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
