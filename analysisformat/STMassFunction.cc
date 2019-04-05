#include "STMassFunction.hh"

Bool_t STMassFunction::GetBBMass(TVector3 mom, Double_t dedx, Int_t chrg, Double_t &fmassH, Double_t &fmassHe, Int_t &fpid1, Int_t &fpid2)
{
  UInt_t ipitch = GetIndexPitch(mom);
  UInt_t iyaw   = GetIndexYaw(mom);

  if( ipitch >= 20 || iyaw >= 20) {
    LOG(DEBUG) << " Out of range in  mass function. " << ipitch << " : " << iyaw  << FairLogger::endl;
    return 0;
  }

  auto afunc =  fBBMassFunction[ipitch][iyaw];

  if( afunc == NULL ) {
    LOG(ERROR) << " No mass function was found. " << FairLogger::endl;
    return 0;
  }

  Double_t dx = 0.1;
  auto bbPar = new Double_t[7];
  bbPar[0] = afunc->GetParameter(0);
  bbPar[1] = afunc->GetParameter(1);
  bbPar[2] = 1.;
  bbPar[3] = 0.; 
  bbPar[4] = 1.;
  bbPar[5] = mom.Mag();
  bbPar[6] = dedx;

  auto funcBB  = [bbPar](double x)    { return MassEstimator::BBMassFinderEq(&x, bbPar);};
  auto dfuncBB = [bbPar, dx](double x){ return MassEstimator::BBMassFinderDeriv(&x, bbPar, dx);};
  ROOT::Math::RootFinder finder;
  ROOT::Math::GradFunctor1D gradfunc1dBB(funcBB, dfuncBB);

  finder.SetMethod(ROOT::Math::RootFinder::kGSL_SECANT);
  finder.SetFunction(gradfunc1dBB, 1500.);
  finder.Solve();
  fmassH = finder.Root();
  if( std::isnan(fmassH) )
    fmassH = 0;

  bbPar[4] = 2.;
  finder.SetFunction(gradfunc1dBB, 3500.);
  finder.Solve();
  fmassHe = finder.Root();
  if( std::isnan(fmassHe) )
    fmassHe = 0;


  fpid1 = GetPID(fmassH, fmassHe, dedx);
  fpid2 = GetPIDLU(fmassH, fmassHe, dedx);
}

Int_t STMassFunction::GetPID(Double_t fmassH, Double_t fmassHe, Double_t dedx)
{
  // p, d, t
  if( fmassHe < 2500 && fmassH > 0 ) { 

    for(UInt_t i = 1; i < 4; i++) {
      Double_t mass_low = MassRegion[i][0]-MassRegion[i][1]*MassRegion[i][2] ;
      Double_t mass_up  = MassRegion[i][0]+MassRegion[i][1]*MassRegion[i][3] ;
      
      if( fmassH >= mass_low && fmassH <= mass_up ) {
	STPID::PID pid = static_cast<STPID::PID>(i);
	return STPID::GetPDG(pid);
      }
    }
  } 

  // He3, He4, He6
  else if( fmassH >= 3100 && dedx <= 700) {
    for( UInt_t i = 4; i < 7; i++ ){
      Double_t mass_low = MassRegion[i][0]-MassRegion[i][1]*MassRegion[i][2] ;
      Double_t mass_up  = MassRegion[i][0]+MassRegion[i][1]*MassRegion[i][3] ;

      if( fmassHe >= mass_low && fmassHe <= mass_up ) {
	STPID::PID pid = static_cast<STPID::PID>(i);
	return STPID::GetPDG(pid);
      }
    }
  } 
  return 0;
}

Int_t STMassFunction::GetPIDLU(Double_t fmassH, Double_t fmassHe, Double_t dedx)
{
  // p, d, t
  if( fmassHe < MassRegionLU[4][0] && fmassH > 0 ) { 

    for(UInt_t i = 1; i < 4; i++) {

      if( fmassH >= MassRegionLU[i][0] && fmassH < MassRegionLU[i][1] ) {
	STPID::PID pid = static_cast<STPID::PID>(i);
	return STPID::GetPDG(pid);
      }
    }
  } 

  // He3, He4, He6
  else if( fmassH >= 3100 && dedx <= 700) {
    for( UInt_t i = 4; i < 7; i++ ){

      if( fmassHe >= MassRegionLU[i][0] && fmassHe < MassRegionLU[i][1] ) {
	STPID::PID pid = static_cast<STPID::PID>(i);
	return STPID::GetPDG(pid);
      }
    }
  } 
  return 0;
}



Bool_t STMassFunction::SetFunction(TString dir, TString fname, TString funcname)
{
  LOG(INFO) << fname << " is called " << FairLogger::endl;
  
  if( !gSystem->FindFile(dir,fname) ) {
    LOG(INFO) << "STMassFunction:: " << fname << " is not fined." << FairLogger::endl;
    return 0;
  }
  else
    LOG(INFO) << "STMassFunction:: " << fname << " is opened." << FairLogger::endl;
    

  TFile *file = new TFile(fname);

  UInt_t i = 0;
  while( i < 20 ) {

    UInt_t j = 0;
    while( j < 20 ) {

      fBBMassFunction[i][j] = (TF1*)file->Get(funcname+Form("_%d_%d",i,j));

      if( fBBMassFunction[i][j] == NULL ) {
	LOG(ERROR) << funcname+Form("_%d_%d",i,j) << " is not Found. " << FairLogger::endl;
	break;
      }
      else 
	LOG(DEBUG) << funcname+Form("_%d_%d",i,j) << " is loaded. " << FairLogger::endl;

      j++;
    }

    i++;
  }


  file->Close();
  delete file;

  return 1;
}

Bool_t STMassFunction::SetMassFitFunction(TString dir, TString fname, TString funcname)
{
  LOG(INFO) << fname << " is called " << FairLogger::endl;
  
  if( !gSystem->FindFile(dir,fname) ) {
    LOG(INFO) << "STMassFunction:: " << fname << " is not fined." << FairLogger::endl;
    return 0;
  }
  else
    LOG(INFO) << "STMassFunction:: " << fname << " is opened." << FairLogger::endl;
    

  TFile *file = new TFile(fname);

  UInt_t i = 0;
  while( i < 6 ) {
    
    fBBMassGaussFit[i] = (TF1*)file->Get(funcname+Form("_%d",i));

    if( fBBMassGaussFit[i] == NULL ) {
      LOG(ERROR) << funcname+Form("_%d",i) << " is not Found. " << FairLogger::endl;
      return 0;
    }
    
    i++;
  }


  file->Close();
  delete file;

  SetMassRegion();

  return 1;
}


UInt_t STMassFunction::GetIndexTheta(Double_t val)
{
  // from 0 deg to 90 per 9 deg 
  UInt_t itheta = UInt_t(val/TMath::DegToRad() / (90/nTheta) );

  if( itheta == nTheta ) itheta = nTheta - 1;
  
  return itheta;
}

UInt_t STMassFunction::GetIndexPhi(Double_t val)
{
  // from -180 to 180 per 20 deg
  UInt_t iphi = UInt_t((val/TMath::DegToRad () + 180.) / (360/nPhi) );
  if( iphi == nPhi ) iphi = nPhi - 1;

  return iphi;
}


UInt_t STMassFunction::GetIndexYaw(TVector3 mom)
{
  TVector2 mompt(mom.Z(), -mom.X());
  Double_t yaw =  TVector2::Phi_mpi_pi( mompt.Phi() );

  return (UInt_t)(yaw + TMath::PiOver2()*8./9.) * nYawBin/(TMath::Pi()*8./9.);
}

UInt_t STMassFunction::GetIndexPitch(TVector3 mom)
{
  TVector2 mompt(mom.Z(), mom.Y());
  Double_t pitch =  TVector2::Phi_mpi_pi( mompt.Phi() );

  return (UInt_t)(pitch + TMath::PiOver2()*8./9.) * nPitchBin/(TMath::Pi()*8./9.);
}



void STMassFunction::SetMassRegion()
{

  for(UInt_t i = 1; i < 7; i++){
    MassRegion[i][0] =  fBBMassGaussFit[i-1]->GetParameter(1);
    MassRegion[i][1] =  fBBMassGaussFit[i-1]->GetParameter(2);

    // LOG(INFO) << " MassRegion is loaded " << i << " " << MassRegion[i][0] 
    // 	      << " - " <<  MassRegion[i][0]-MassRegion[i][1]*MassRegion[i][2] 
    // 	      << " + " <<  MassRegion[i][0]+MassRegion[i][1]*MassRegion[i][3] 
    // 	      << "("<<MassRegion[i][1] <<")"
    // 	      << FairLogger::endl;


    LOG(INFO) << MassRegion[i][0] << ", " << MassRegion[i][1] << FairLogger::endl;


  }

}
