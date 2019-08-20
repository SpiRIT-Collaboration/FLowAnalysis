#include "FlowFunctions.h"
#include "DoFlow_adv.C"


void simFlow(Double_t p1=-0.04)
{
  UInt_t mult = 20;
  STBootStrap *bstrap  = new STBootStrap(1);

  gRandom->SetSeed(0);

  Double_t px;
  Double_t py;

  auto *fcos1 = new TF1("fcos1","1. + 2.*[0]*cos(x)"   ,-1.*TMath::Pi(),TMath::Pi());
  auto *fcos2 = new TF1("fcos2","1. + 2.*[0]*cos(x)+2.*[1]*cos(2.*x)"   ,-1.*TMath::Pi(),TMath::Pi());

  auto hcos0 = new TH1D("hcos0","phi rand "  ,100,-1.*TMath::Pi(),TMath::Pi());
  auto hcos1 = new TH1D("hcos1","2xphi rand ",100,-1.*TMath::Pi(),TMath::Pi());
  auto hcos2 = new TH1D("hcos2","cos(rand) ",100,-1., 1.);
  auto hcos3 = new TH1D("hcos3","cos(2*rand) ",100,-1., 1.);

  auto hcos11 = new TH1D("hcos11","<cos(#Delta #phi)>",100,-1.,1.);
  auto hcos12 = new TH1D("hcos12","<cos(2x#Delta #Phi)>",100,-1.,1.);
  auto hcos13 = new TH1D("hcos13","#phi - #Psi"     ,100,-1.*TMath::Pi(),TMath::Pi());
  auto hcos14 = new TH1D("hcos14","2x(#phi - #Psi)" ,100,-1.*TMath::Pi(),TMath::Pi());
  auto hcos15 = new TH1D("hcos15","abs(2x(#phi - #Psi))" ,100,        0.,TMath::Pi());

  auto hphi12 = new TH1D("hphi12","<cos(phi)>",100,-1.,1.);

  auto hrpphi0  = new TH1D("hrpphi0" ,"RP Phi"     ,100,0.,TMath::Pi());
  auto hrpphi90 = new TH1D("hrpphi90","RP Phi > 90",100,0.,TMath::Pi());

  //--

  std::vector<Double_t> vvcos;
  Double_t *x = new Double_t[2];

  fcos1->SetParameter(0, 0.2);
  fcos2->SetParameter(0, 0.);
  fcos2->SetParameter(1, p1);


  for(UInt_t i = 0; i < 20000; i++) {

    Double_t rphi1 = fcos1->GetRandom();
    Double_t rphi2 = fcos2->GetRandom();

    hcos0->Fill(rphi1);
    hcos1->Fill(rphi2);
    hcos2->Fill(cos( rphi1 ));
    hcos3->Fill(cos( 2.*rphi2 ));


    bstrap->Clear();
    TVector3 pxyA;
    TVector3 pxyB;
    for(UInt_t j = 0; j < mult; j++){
      Double_t rphi = fcos1->GetRandom();
      pxyA += TVector3(cos(rphi), sin(rphi), 0);
      bstrap->Add( TVector2( cos(rphi), sin(rphi) ) );

      rphi = fcos1->GetRandom();
      pxyB += TVector3(cos(rphi), sin(rphi), 0);
      bstrap->Add( TVector2( cos(rphi), sin(rphi) ) );

    }

    TVector3 pxy = pxyA + pxyB;

    auto delta_Phi = TVector2::Phi_mpi_pi( pxyA.Phi() - pxyB.Phi() );
    hrpphi0->Fill(abs( delta_Phi ) ); 
    if( abs( delta_Phi ) > TMath::Pi()/2. ) 
      hrpphi90->Fill(abs( delta_Phi ) );
    
    hcos11->Fill(cos(rphi1 - pxy.Phi()));
    hcos12->Fill(cos(2.*(rphi2 - pxy.Phi())));
    hcos13->Fill(TVector2::Phi_mpi_pi(rphi1 - pxy.Phi()));
    hcos14->Fill(TVector2::Phi_mpi_pi(2.*(rphi2 - pxy.Phi())));
    hcos15->Fill(abs(TVector2::Phi_mpi_pi(2.*(rphi2 - pxy.Phi()))));

    bstrap->BootStrapping(100);

  }

  Double_t *rpres = new Double_t[4];
  rpres = GetRPResolutionwChi(hrpphi0, hrpphi90);
  
  auto cv1  = hcos11 -> GetMean() / rpres[0];
  auto cv1n = hcos11 -> GetStdDev()/sqrt( hcos11->GetEntries() );
  auto cv1e = GetError( cv1, rpres[0], cv1n, rpres[1] );

  auto cv2  = hcos12 -> GetMean() / rpres[2];
  auto cv2n = hcos12 -> GetStdDev()/sqrt( hcos12->GetEntries() );
  auto cv2e = GetError( cv2, rpres[2], cv2n, rpres[3] );
  
  cout << "--- calculated v1 ---" << endl;
  cout << " RP res v1 " << rpres[0] << " and v2 " << rpres[2] << endl;
  cout << " v1 " << cv1 << " +- " << cv1e << " / " << hcos11 -> GetMean() << endl;
  cout << " v2 " << cv2 << " +- " << cv2e << " / " << hcos12 -> GetMean() << endl;

  cout << "--- Original v1 and v2 ---" << endl;
  cout << " <cos(Phi)>  v1 " << hcos2->GetMean() << " +- " << hcos2->GetRMS()/sqrt(hcos2->GetEntries()) << endl;
  cout << " <cos(2Phi)> v2 " << hcos3->GetMean() << " +- " << hcos3->GetRMS()/sqrt(hcos3->GetEntries()) << endl;

  cout << "--- BootStrap v1 ---" << endl;
  cout << " mean = " << bstrap->GetMean() << endl;
  cout << " <cos phi> = " << bstrap->GetCosMean() << endl;
  cout << " StdDev = " << bstrap->GetStdDev() << endl;


  auto cc3 = new TCanvas("cc3","cc3");
  hrpphi0 ->SetLineColor(4);
  hrpphi90->SetLineColor(2);
  hrpphi0 ->Draw();
  hrpphi90->Draw("same");

  auto cc2 = new TCanvas("cc2","cc2");
  cc2->Divide(1,2);
  cc2->cd(1);
  hcos0->SetNormFactor(100);
  hcos0->Draw();
  cc2->cd(2);
  hcos1->SetNormFactor(100);
  hcos1->Draw();

  auto cc4 = new TCanvas("cc4","cc4");
  hcos11->Draw();

  auto cc5 = new TCanvas("cc5","cc5");
  hcos12->Draw();

  auto cc6 = new TCanvas("cc6","cc6");
  fcos2->Draw();

  auto cc7 = new TCanvas("cc7","cc7");
  cc7->Divide(1,3);
  cc7->cd(1);
  hcos13->SetNormFactor(100);
  hcos13->Draw();

  cc7->cd(2);
  hcos14->SetNormFactor(100);
  hcos14->Draw();

  cc7->cd(3);
  hcos15->SetNormFactor(100);
  hcos15->Draw();
}


