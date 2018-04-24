#include "FlowFunctions.h"

void simFlw(Double_t p1=0.2)
{
  UInt_t mult = 15;

  gRandom->SetSeed(0);

  Double_t px;
  Double_t py;

  auto *fcos  = new TF1("fcos" ," 1 +[0]*cos(x)",-1.*TMath::Pi(),TMath::Pi());
  auto *fcos1 = new TF1("fcos1","[0]+[1]*cos(x)",-1.*TMath::Pi(),TMath::Pi());
  fcos->SetParameter(0,p1);
  auto cc0 = new TCanvas("cc0","cc0");
  fcos->Draw();
  
  // auto hcos = new TH1D("hcos","v1 rand",100,-3.2,3.2);
  // hcos->FillRandom("fcos",200000);

  auto hcos1 = new TH1D("hcos1","phi rand ",100,-1.*TMath::Pi(),TMath::Pi());
  auto hcos2 = new TH1D("hcos2","cos rand ",100,-1., 1.);
  auto hcos3 = new TH1D("hcos3","sin rand ",100,-1., 1.);

  auto hcos11 = new TH1D("hcos11","<cos>",100,-1.,1.);

  auto hcos12 = new TH1D("hcos12","<cos(RP)>",100,-1.,1.);
  auto hcos13 = new TH1D("hcos13","RP_Phi"   ,100,-1.*TMath::Pi(),TMath::Pi());

  auto hphi12 = new TH1D("hphi12","<cos(phi)>",100,-1.,1.);
  //--

  std::vector<Double_t> vvcos;
  Double_t *x = new Double_t[2];

  for(UInt_t i = 0; i < 2000000; i++) {
    Double_t rphi = fcos->GetRandom();
    hcos1->Fill(rphi);
    hcos2->Fill(cos( rphi ));
    hcos3->Fill(sin( rphi ));

    std::vector<Double_t> vcos;
    std::vector<Double_t> vphi;
    TVector3 pxy;

    for(UInt_t j = 0; j < mult; j++){
      Double_t rphi = fcos->GetRandom();
      vphi.push_back( rphi );
      vcos.push_back( cos(rphi) );

      pxy += TVector3(cos(rphi), sin(rphi), 0);
    }
    
    x = vMean(vcos);
    
    vvcos.push_back(x[0]);

    rphi = fcos->GetRandom();
    hcos13->Fill(TVector2::Phi_mpi_pi(rphi-pxy.Phi()));
    hcos12->Fill(cos(pxy.Phi()));

    x = vMean(vphi);
    hphi12->Fill(x[0]);

  }

  x = vMean(vvcos);
  cout << " <cos<phi>> " << x[0] << " +- " << x[1] << endl;
  

  auto cc2 = new TCanvas("cc2","cc2");
  hcos2->Draw();

  auto cc3 = new TCanvas("cc3","cc3");
  hphi12->Draw();

  cout << " <cos(Phi)> " << hcos2->GetMean() << " +- " << hcos2->GetRMS();
  cout << " <sin(Phi)> " << hcos3->GetMean() << " +- " << hcos3->GetRMS() << endl;;


  auto cc4 = new TCanvas("cc4","cc4");
  hcos11->Draw();
  cout << " << cos(Phi) >> " << hcos11->GetMean() << endl;

  auto cc5 = new TCanvas("cc5","cc5");
  hcos12->Draw();

  auto cc6 = new TCanvas("cc6","cc6");
  hcos13->SetNormFactor(100);
  hcos13->Draw();
  //  hcos13->Fit("fcos");

  auto cc1 = new TCanvas("cc1","cc1");
  //  hcos1->SetNormFactor(100);
  TH1D *hcosfit = new TH1D( (*hcos1) );
  hcosfit->Draw();
  fcos1->SetParameter(0,1);
  fcos1->SetParameter(1,p1);
  hcos1->Fit("fcos1");

  Double_t fn  = fcos1->GetParameter(0);
  Double_t fv1 = fcos1->GetParameter(1);
  Double_t fne = fcos1->GetParError(0);
  Double_t fv1e= fcos1->GetParError(1);

  cout << " v1 = " << fv1/fn/2. << " +- " << 0.5*sqrt( pow(fv1/fn,2)*( pow(fv1e/fv1,2) + pow(fne/fn,2) ) ) << endl;

}


