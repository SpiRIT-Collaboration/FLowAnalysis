void BootStrapTest()
{
  Double_t pxy[][2] = { { 300.30820 , -196.1504 },
			{ 151.66651 , -271.1333 },
			{ 177.78397 , -93.96121 },
			{ 48.376354 , 97.354136 },
			{ -93.56772 , -13.82681 },
			{ -101.2725 , -206.6223 },
			{ -330.6198 , -233.9083 },
			{ -134.8822 , -43.27171 },
			{ 31.988070 , -280.8166 },
			{ 1085.5061 , -13.50405 },
			{ -77.49755 , 154.50885 },
			{ 407.17399 , -234.7863 },
			{ 311.75433 , 135.90189 },
			{ 108.63103 , -280.2595 },
			{ -518.5169 , -209.0715 },
			{ 384.42002 , -278.8419 },
			{ -815.0031 , 28.662856 },
			{ -414.8316 , 66.177955 },
			{ 326.15717 , -426.7524 },
			{ 500.45768 , 346.02029 }};



  auto hstrap = new TGraphErrors();
  hstrap->SetTitle(";Number of Trial; Reaction Plane Orientaion [rad]");

  auto bstrap = new STBootStrap(1);

  TVector2 vcsum = TVector2(0.,0.);
  

  for(UInt_t i = 0; i < 22; i++) {
    TVector2 vec = TVector2( pxy[i][0], pxy[i][1] ).Unit(); 

    bstrap->Add(vec);

    vcsum += vec;
  }

  UInt_t itry = 2000000;
  bstrap->BootStrapping(itry);

  for(UInt_t it = 1; it < itry; it++){

    hstrap->SetPoint(it-1, (Double_t)it, bstrap->GetResidualMean(it));
    hstrap->SetPointError(it-1, 0., bstrap->GetResidualStdDev(it));

  }
  
  //  cout << " res mean " << bstrap->GetResidualMean() << endl;
  cout << " BST mean " << bstrap->GetMean() << endl;
  cout << " vsum Phi " << TVector2::Phi_mpi_pi(vcsum.Phi()) << endl;


  hstrap->SetMarkerStyle(20);
  hstrap->SetMarkerColor(2);
  hstrap->Draw("alp");

  auto aLine = new TLine(0., TVector2::Phi_mpi_pi(vcsum.Phi()), 
			 hstrap->GetXaxis()->GetXmax(), TVector2::Phi_mpi_pi(vcsum.Phi()));
  aLine->Draw();

}
