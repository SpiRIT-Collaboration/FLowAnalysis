#include "STMassFunction.hh"


Bool_t STMassFunction::SetFunction(TString dir, TString fname)
{
  std::cout << fname << " is called " << std::endl;
  
  if( !gSystem->FindFile(dir,fname) ) {
    std::cout << "STMassFunction:: " << fname << " is not fined." << std::endl;
    return 0;
  }
  else
    std::cout << "STMassFunction:: " << fname << " is opened." << std::endl;
    

  TFile *file = new TFile(fname);

  for(UInt_t i = 0; i < nTheta; i++) {

    for(UInt_t j = 0; j < nPhi; j++) {


      fBBMassFunction[i][j] = (TF1*)file->Get(Form("f1BBProton_%d_%d",i,j));
      fLVMassFunction[i][j] = (TF1*)file->Get(Form("f1LVProton_%d_%d",i,j));

      if( fBBMassFunction[i][j] == NULL )
	std::cout << Form("f1BBProton_%d_%d",i,j) << " is not Found. " << std::endl;

      if( fLVMassFunction[i][j] == NULL )
	std::cout << Form("f1LBProton_%d_%d",i,j) << " is not Found. " << std::endl;


    }
  }


  file->Close();
  delete file;

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


TF1* STMassFunction::GetBBFunction(Double_t valtheta, Double_t valphi)
{
  auto itheta = GetIndexTheta(valtheta);
  auto iphi   = GetIndexPhi(valphi);

  return fBBMassFunction[itheta][iphi];
}

TF1* STMassFunction::GetLVFunction(Double_t valtheta, Double_t valphi)
{
  auto itheta = GetIndexTheta(valtheta);
  auto iphi   = GetIndexPhi(valphi);

  return fLVMassFunction[itheta][iphi];
}
