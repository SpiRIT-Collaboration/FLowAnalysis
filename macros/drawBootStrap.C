#include "openFlw.C"
#include "SetStyle.C"
TChain* rChain0;
auto *fcos1 = new TF1("fcos1","[0]+2.*[1]*cos(x)"   ,-1.*TMath::Pi(),TMath::Pi());

void sub_Phi1();
void sub_Phi2();
void GetResolution(TH1D* hres, Double_t *mcos);
void GetResolution2(Double_t *mcos);


void sub_Phi5()
{
  auto h2phi_od = new TH2D("h2phi_od",";2#Psi_{od}; 2#Delta(#Psi^{A}_{od}-#Psi^{B}_{od})", 80,-3.14,3.14, 100, -3.14, 3.14);
  auto h2phi_bs = new TH2D("h2phi_bs",";2#Psi_{od}; 2#Delta(#Psi^{A}_{bs}-#Psi^{B}_{bs})", 80,-3.14,3.14, 100, -3.14, 3.14);
  auto h2cos_od = new TH2D("h2cos_od","", 80,-3.14,3.14, 100, -3.14, 3.14);
  auto h2cos_bs = new TH2D("h2cos_bs","", 80,-3.14,3.14, 100, -3.14, 3.14);


  rChain0->Project("h2phi_od","od_fd2Phi:TVector2::Phi_mpi_pi(2.*unitP_fc.Phi())");
  rChain0->Project("h2phi_bs","bs_fd2Phi:TVector2::Phi_mpi_pi(2.*unitP_fc.Phi())");

  rChain0->Project("h2cos_od","od_fd2Phi:TVector2::Phi_mpi_pi(2.*unitP_fc.Phi())");
  rChain0->Project("h2cos_bs","bs_fd2Phi:TVector2::Phi_mpi_pi(2.*unitP_fc.Phi())");


  auto cc0 = new TCanvas("cc0","cc0");
  h2phi_od->Draw("colz");

  auto cc1 = new TCanvas("cc1","cc1");
  h2phi_bs->Draw("colz");

  
}

void sub_Phi4()
{
  TH1D* hh2phi0_od;
  TH1D* hh2phi0_bs;

  Double_t div[4] = {-1.*TMath::Pi(), -1.*TMath::Pi()/2., TMath::Pi()/2., TMath::Pi()};


  auto cc3 = new TCanvas("cc3","cc3",1000,500);
  cc3->Divide(3,1);

  for(UInt_t i = 0; i < 3; i++) {
    
    TString title = Form(" %4.2f < 2#Psi_{od} < %4.2f ; 2#Delta(#Psi^{A}-#Psi^{B})", div[i], div[i+1]);
    hh2phi0_od = new TH1D(Form("hh2phi_od%d",i), title, 40,-3.14,3.14);
    hh2phi0_bs = new TH1D(Form("hh2phi_bs%d",i), title, 40,-3.14,3.14);

    rChain0->Project(Form("hh2phi_od%d",i),"od_fd2Phi", Form("%f <= 2Psi && 2Psi < %f",div[i], div[i+1] ));
    rChain0->Project(Form("hh2phi_bs%d",i),"bs_fd2Phi", Form("%f <= 2Psi && 2Psi < %f",div[i], div[i+1] ));


    cc3->cd(i+1);
    hh2phi0_bs->SetNormFactor(40);
    hh2phi0_bs->SetLineColor(4);
    hh2phi0_bs->Draw("e");

    hh2phi0_od->SetNormFactor(40);
    hh2phi0_od->SetLineColor(2);
    hh2phi0_od->Draw("same e");

  }

  TLegend *aleg0 = new TLegend(0.56,0.8,0.87,0.92,"");
  aleg0->AddEntry("hh2phi_od2" ,"Ordinal");
  aleg0->AddEntry("hh2phi_bs2" ,"BootStrap");
  aleg0->Draw();

}
 
void sub_Phi0()
{
  auto hphio = new TH1D("hphio","Subevent A            ; Azimuthal angle [rad]",80, -3.14, 3.14);
  auto hphib = new TH1D("hphib","Subevent A (BootStrap); Azimuthal angle [rad]",80, -3.14, 3.14);

  auto hphioc = new TH1D("hphioc","Subevent A            ; Azimuthal angle [rad]",80, -3.14, 3.14);
  auto hphibc = new TH1D("hphibc","Subevent A (BootStrap); Azimuthal angle [rad]",80, -3.14, 3.14);

  rChain0->Project("hphio","TVector2::Phi_mpi_pi(unitP_1r.Phi())");
  rChain0->Project("hphib","TVector2::Phi_mpi_pi(bsPhi_1[0])");

  rChain0->Project("hphioc","TVector2::Phi_mpi_pi(unitP_1.Phi())");
  rChain0->Project("hphibc","TVector2::Phi_mpi_pi(bsP_1.Phi())");

  hphio->SetLineColor(2);
  hphib->SetLineColor(4);
  hphioc->SetLineColor(5);
  hphibc->SetLineColor(8);


  auto cc0 = new TCanvas("cc0","cc0");
  hphio->Draw("e");
  hphib->Draw("same e");

  auto aleg = new TLegend(0.7,0.1,0.88 ,0.3,"");
  aleg->AddEntry(hphio, "Ordinal ", "lp");
  aleg->AddEntry(hphib, "BootStrap ", "lp");
  //  aleg->Draw();



  //  auto cc1 = new TCanvas("cc1","cc1");
  //  hphioc->SetMinimum(0);
  hphioc->Draw("same e");
  hphibc->Draw("same e");

  auto alegc = new TLegend(0.7,0.1,0.88 ,0.3,"");
  aleg->AddEntry(hphioc, "Ordinal(corrected)", "lp");
  aleg->AddEntry(hphibc, "BootStrap(corrected)", "lp");
  aleg->Draw();

}


void sub_Phi1()
{
  auto hmult = new TH1D("hmult","; Number of Tracks;",65,0,65);

  auto hstd  = new TH2D("hstd" ,"; #Psi^{A}_{bs}; StdDev",100,-3.14,3.14, 100,0.,0.3);
  auto hcor  = new TH2D("hcor" ,"; #Psi_{od};  #Psi^{A}_{bs}", 100, -3.14, 3.14, 100, -3.14, 3.14);
  auto hdif  = new TH2D("hdif" ,"; #Psi_{od};  ^{*}#Psi^{A}_{bs} - ^{*}#Psi^{A}_{od}", 100, -3.14, 3.14, 100, -0.5, 0.5);
  auto hdphi = new TH2D("hdphi","; #Psi_{od};  #Psi^{A}_{bs} - #Psi^{A}_{od}", 100, -3.14, 3.14, 100, -0.5, 0.5);
  auto hbsmlt = new TH2D("hbsmlt","; ^{*}#Psi; Number of Tracks",80, -3.14, 3.14, 65,0,65);


  rChain0->Project("hmult", "ntrack[4]");
  rChain0->Project("hstd" , "bsPhi_1[1]:TVector2::Phi_mpi_pi(bsP_1.Phi())");
  rChain0->Project("hcor" , "TVector2::Phi_mpi_pi(bsPhi_1[0]):unitP_fc.Phi()");
  rChain0->Project("hdif" , "TVector2::Phi_mpi_pi(bsPhi_1[0]-unitP_1r.Phi()):unitP_fc.Phi()");
  rChain0->Project("hdphi", "TVector2::Phi_mpi_pi(bsP_1.Phi()-unitP_1.Phi()):unitP_fc.Phi()");
  rChain0->Project("hbsmlt","ntrack[4]:unitP_rot.Phi()");

  auto cc1 = new TCanvas("cc1","cc1");
  hstd->Draw("colz");

  auto cc2 = new TCanvas("cc2","cc2");
  hcor->Draw("colz");

  auto cc3 = new TCanvas("cc3","cc3");
  hdif->Draw("colz");

  auto cc4 = new TCanvas("cc4","cc4");
  hdphi->Draw("colz");


  auto cc5 = new TCanvas("cc5","cc5");
  hmult->SetFillColor(4);
  hmult->Draw();


  auto cc6 = new TCanvas("cc6","cc6");
  hbsmlt->Draw("colz");
}

void sub_Phi2()
{
  // resolution
  auto hres_od = new TH2D("hres_od","Resolution Ordinal;   #Psi_{od}; abs(#Psi^{A}_{od} - #Psi^{B}_{od})", 80, -3.14, 3.14, 200, 0., 3.14);
  auto hres_bs = new TH2D("hres_bs","Resolution BootStrap; #Psi_{od}; abs(#Psi^{A}_{bs} - #Psi^{B}_{bs})", 80, -3.14, 3.14, 200, 0., 3.14); 

  rChain0->Project("hres_od","abs(od_fdPhi):unitP_fc.Phi()");
  rChain0->Project("hres_bs","abs(bs_fdPhi):unitP_fc.Phi()");


  auto cc = new TCanvas("cc11","cc11", 1200,1000);
  cc->Divide(2,5);


  TH1D *hres0;
  TH1D *hres1;
  auto grres1_od = new TGraphErrors();
  auto grres1_bs = new TGraphErrors();

  UInt_t div = 80/10;

  for(UInt_t i = 0; i < 10; i++){
    cc->cd(i+1);

    Double_t le = hres_bs->GetXaxis()->GetBinLowEdge(div*i);
    Double_t ue = hres_bs->GetXaxis()->GetBinLowEdge(div*(i+1));
    auto xcenter  = (le + ue)/2.;
    auto xcentere = (ue - le)/2.;

    Double_t mcos0[4];
    Double_t mcos1[4];

    hres0  = (TH1D*)hres_od->ProjectionY((TString)Form("hres0_od%d",i), div*i, div*(i+1), "eo");
    cout << " od divide  " << div*i << " to " << div*(i+1) << endl;
    hres0->SetTitle(Form("%4.2f ~ %4.2f; #Delta( #Psi^{A} - #Psi^{B})",le,ue));
    hres0->SetNormFactor(80);
    hres0->SetLineColor(kRed);
    hres0->Draw("e");

    hres1  = (TH1D*)hres_bs->ProjectionY((TString)Form("hres0_bs%d",i), div*i, div*(i+1), "eo"); 
    cout << " bs divide  " << div*i << " to " << div*(i+1) << endl;
    hres1->SetTitle(Form("%4.2f ~ %4.2f; #Delta( #Psi^{A} - #Psi^{B})",le,ue));
    hres1->SetNormFactor(80);
    hres1->SetLineColor(kBlue);
    hres1->Draw("same e");

    mcos0[0] = hres_od->Integral(div*i, div*(i+1));
    mcos0[1] = hres_od->Integral(div*i, div*(i+1), 101, 201);
    cout << hres0->GetName() << endl;
    GetResolution2(mcos0);
    grres1_od->SetPoint(i, xcenter, mcos0[0]);
    grres1_od->SetPointError(i, xcentere, mcos0[1]);
    grres1_od->SetLineColor(2);
    
    mcos1[0] = hres_bs->Integral(div*i, div*(i+1));
    mcos1[1] = hres_bs->Integral(div*i, div*(i+1), 101, 201);
    cout << hres1->GetName() << endl;
    GetResolution2(mcos1);
    grres1_bs->SetPoint(i, xcenter, mcos1[0]);
    grres1_bs->SetPointError(i, xcentere, mcos1[1]);
    grres1_bs->SetLineColor(4);

  }


  Double_t mcos[4];
  mcos[0] = hres_bs->Integral(-1, 81);
  mcos[1] = hres_bs->Integral(-1, 81, 101, 201);
  GetResolution2(mcos);

  auto averes1_bs = new TLine(-3.53, mcos[0], 3.4, mcos[0]);
  averes1_bs->SetLineColor(4);


  mcos[0] = hres_od->Integral(-1, 81);
  mcos[1] = hres_od->Integral(-1, 81, 101, 201);
  GetResolution2(mcos);

  auto averes1_od = new TLine(-3.53, mcos[0], 3.4, mcos[0]);
  averes1_od->SetLineColor(2);


  auto cc14 = new TCanvas("cc14","cc14");
  auto mgrres1 = new TMultiGraph();
  mgrres1->SetTitle(";#Psi_{od}; <cos #Delta #Psi >");
  mgrres1->Add(grres1_od,"lp");
  mgrres1->Add(grres1_bs,"lp");
  mgrres1->Draw("ALP");
  averes1_bs->Draw();
  averes1_od->Draw();
  TLegend *aleg5= new TLegend(0.5,0.8,0.75,0.9,"");
  aleg5->AddEntry(grres1_od,"Ordinal");
  aleg5->AddEntry(grres1_bs,"BootStrap");
  aleg5->Draw();

}


void sub_Phi3()
{
  auto hfsub_bs = new TH1D("hfsub_bs","; #Delta( #Psi^{A}_{bs} - #Psi^{B}_{bs})",100,-3.14, 3.14);
  auto hfsub_od = new TH1D("hfsub_od","; #Delta( #Psi^{A}_{od} - #Psi^{B}_{od})",100,-3.14, 3.14);
  auto hf_bsod  = new TH2D("hf_bsod" ,"; #Psi; #Psi^{A}_{bs} - #Psi^{A}_{od}", 100, -3.14, 3.14, 100, -0.5, 0.5);
  auto hf_dbsod = new TH2D("hf_dbsod","; #Psi^{A}_{od} - #Psi^{B}_{od}; #Psi^{A}_{bs} - #Psi^{B}_{bs}", 100, -3.14, 3.14, 100, -3.14, 3.14);
  rChain0->Project("hfsub_bs","bs_fdPhi","mtrack_1>1");
  rChain0->Project("hfsub_od","od_fdPhi","mtrack_1>1");
  rChain0->Project("hf_bsod", "TVector2::Phi_mpi_pi(bsP_1.Phi()-unitP_1.Phi()):TVector2::Phi_mpi_pi(unitP_fc.Phi())");
  rChain0->Project("hf_dbsod","od_fdPhi:bs_fdPhi");


  auto hsub_bs = new TH1D("hsub_bs",";#Delta( #Psi^{A}_{bs} - #Psi^{B}_{bs} )",100,-3.14, 3.14);
  auto hsub_od = new TH1D("hsub_od",";#Delta( #Psi^{A}_{od} - #Psi^{B}_{od} )",100,-3.14, 3.14);
  auto h_bsod  = new TH2D("h_bsod" ,";#Psi; #Psi^{A}_{bs} - #Psi^{A}_{od} ", 100, -3.14, 3.14, 100, -0.5, 0.5);
  
  rChain0->Project("hsub_bs","TVector2::Phi_mpi_pi(bsPhi_1[0]-bsPhi_2[0])","mtrack_1>1");
  rChain0->Project("hsub_od","TVector2::Phi_mpi_pi(unitP_1r.Phi()-unitP_2r.Phi())","mtrack_1>1");
  rChain0->Project("h_bsod", "TVector2::Phi_mpi_pi(bsPhi_1[0]-unitP_1r.Phi()):TVector2::Phi_mpi_pi(unitP_fc.Phi())");


  auto hfdif_bs = new TH2D("hfdif_bs","BootStrap; #Psi_{od}; #Delta( #Psi^{A}_{bs} - #Psi^{B}_{bs})",80,-3.14, 3.14, 80,-3.14, 3.14);
  auto hfdif_od = new TH2D("hfdif_od","Ordinal;   #Psi_{od}; #Delta( #Psi^{A}_{od} - #Psi^{B}_{od})",80,-3.14, 3.14, 80,-3.14, 3.14);

  rChain0->Project("hfdif_bs","bs_fdPhi:unitP_fc.Phi()");
  rChain0->Project("hfdif_od","od_fdPhi:unitP_fc.Phi()");

  //---------------------------------------------
  auto cc9 = new TCanvas("cc9","cc9");
  hfdif_bs->Draw("colz");

  auto cc10 = new TCanvas("cc10","cc10");
  hfdif_od->Draw("colz");



  TH1D *hslice;


  UInt_t div = 80/10;

  auto cc = new TCanvas("cc11","cc11", 1200,1000);
  cc->Divide(2,5);


  for(UInt_t i = 0; i < 10; i++){
    cc->cd(i+1);

    Double_t le = hfdif_bs->GetXaxis()->GetBinLowEdge(div*i);
    Double_t ue = hfdif_bs->GetXaxis()->GetBinLowEdge(div*(i+1));
    auto xcenter  = (le + ue)/2.;
    auto xcentere = (ue - le)/2.;
    hslice = (TH1D*)hfdif_bs->ProjectionY((TString)Form("hdiff_bs%d",i),div*i, div*(i+1), "eo");
    hslice->SetTitle(Form("%4.2f ~ %4.2f; #Delta( #Psi^{A} - #Psi^{B})",le,ue));
    hslice->SetNormFactor(80);
    hslice->SetLineColor(kRed);
    //      hslice->SetLineColor(kRed-10+2*i);
    //      if(i == 0)
    hslice->Draw();
    // else
    // 	hslice->Draw("same");
    //    }

    //    cc = new TCanvas("cc12","cc12");
    //    for(UInt_t i = 0; i < 5; i++){
    //      UInt_t il = 16*i;
    hslice = (TH1D*)hfdif_od->ProjectionY((TString)Form("hdiff_od%d",i),div*i, div*(i+1), "eo");
    hslice->SetNormFactor(80);
    hslice->SetLineColor(kBlue);
    //      hslice->SetLineColor(kBlue-10+2*i);
    // if(i == 0)
    // 	hslice->Draw();
    // else
      
    hslice->Draw("same");
  }
    



  auto cc7 = new TCanvas("cc7","cc7");
  hfsub_od->SetNormFactor(100);
  hfsub_od->SetLineColor(2);
  hfsub_od->SetTitle("; #Delta( #Psi^{A} - #Psi^{B})");
  hfsub_od->Draw("e");
  hfsub_bs->SetNormFactor(100);
  hfsub_bs->SetLineColor(4);
  hfsub_bs->Draw("same e");
  
  TLegend *aleg2= new TLegend(0.65,0.75,0.85,0.9,"w Flattening");
  aleg2->AddEntry("hfsub_od","Ordinal");
  aleg2->AddEntry("hfsub_bs","BootStrap");
  aleg2->Draw();


  auto cc5 = new TCanvas("cc5","cc5");
  hsub_od->SetNormFactor(100);
  hsub_od->SetLineColor(2);
  hsub_od->SetTitle("; #Delta( #Psi^{A} - #Psi^{B})");
  hsub_od->Draw("e");
  
  hsub_bs->SetNormFactor(100);
  hsub_bs->SetLineColor(4);
  hsub_bs->Draw("same e");

  TLegend *aleg1 = new TLegend(0.65,0.75,0.85,0.9,"w/o flattening");
  aleg1->AddEntry("hsub_od","Ordinal");
  aleg1->AddEntry("hsub_bs","BootStrap");
  aleg1->Draw();


  auto cc4 = new TCanvas("cc4","cc4");
  hsub_bs->SetLineColor(6);
  hsub_bs->SetTitle("; #Delta( #Psi^{A} - #Psi^{B})");
  hsub_bs->Draw("e");
  hfsub_bs->Draw("same e");
  TLegend *aleg3 = new TLegend(0.65,0.75,0.85,0.9,"BootStrap");
  aleg3->AddEntry("hsub_bs" ,"w/o corr");
  aleg3->AddEntry("hfsub_bs","w   corr");
  aleg3->Draw();

  auto cc6 = new TCanvas("cc6","cc6"); 
  hsub_od->SetLineColor(8);
  hsub_od->SetTitle("; #Delta( #Psi^{A} - #Psi^{B})");
  hsub_od->Draw("e");
  hfsub_od->Draw("same e");
  TLegend *aleg4 = new TLegend(0.65,0.75,0.85,0.9,"Ordinal");
  aleg4->AddEntry("hsub_od" ,"w/o corr");
  aleg4->AddEntry("hfsub_od","w   corr");
  aleg4->Draw();


  auto cc8 = new TCanvas("cc8","cc8");
  h_bsod->Draw("colz");
}

void GetResolution(TH1D* hres, Double_t *mcos)
{
  auto ndiv  = hres->GetNbinsX();
  auto denom = hres->Integral();
  auto numer = hres->Integral(ndiv/2+1,ndiv+1);

  auto chi    = sqrt(-2. * log(2. * numer/denom));
  auto ratioe = numer/denom * sqrt(1./numer + 1./denom);
  auto chie   = chi - sqrt( -2.*log(2. * (numer/denom + ratioe) ) );

  mcos[0] = sqrt(TMath::Pi())/(2.*TMath::Gamma(1))*chi;
  mcos[1] = sqrt(TMath::Pi())/(2.*TMath::Gamma(1))*(chi+chie) - mcos[0];
  mcos[2] = sqrt(TMath::Pi())/(2.*2.*TMath::Gamma(1.5))*pow(chi,2);
  mcos[3] = sqrt(TMath::Pi())/(2.*2.*TMath::Gamma(1.5))*pow(chi+chie,2) - mcos[2];

  std::cout << hres->GetName()  << std::endl;

  std::cout << numer << " / " << denom << " = " << numer/denom << " -> Chi " << chi << " +- " << chie
	    << std::endl;

  std::cout << " <cos(Phi)> = " << mcos[0] << " +- " << mcos[1]
	    << " <cos(2Phi)> = "<< mcos[2] << " +- " << mcos[3]
	    << std::endl;

  return mcos;

}

void GetResolution2(Double_t *mcos)
{
  auto denom = mcos[0];
  auto numer = mcos[1];

  auto chi    = sqrt(-2. * log(2. * numer/denom));
  auto ratioe = numer/denom * sqrt(1./numer + 1./denom);
  auto chie   = chi - sqrt( -2.*log(2. * (numer/denom + ratioe) ) );

  mcos[0] = sqrt(TMath::Pi())/(2.*TMath::Gamma(1))*chi;
  mcos[1] = sqrt(TMath::Pi())/(2.*TMath::Gamma(1))*(chi+chie) - mcos[0];
  mcos[2] = sqrt(TMath::Pi())/(2.*2.*TMath::Gamma(1.5))*pow(chi,2);
  mcos[3] = sqrt(TMath::Pi())/(2.*2.*TMath::Gamma(1.5))*pow(chi+chie,2) - mcos[2];


  std::cout << numer << " / " << denom << " = " << numer/denom << " -> Chi " << chi << " +- " << chie
	    << std::endl;

  std::cout << " <cos(Phi)> = " << mcos[0] << " +- " << mcos[1]
	    << " <cos(2Phi)> = "<< mcos[2] << " +- " << mcos[3]
	    << std::endl;

  return mcos;

}
void drawBootStrap()
{

  SetStyle();

  openFlw();

  rChain0 = (TChain*)gROOT->FindObject("rChain0");

  rChain0->SetAlias("od_fdPhi" ,"TVector2::Phi_mpi_pi(unitP_1.Phi()-unitP_2.Phi())");
  rChain0->SetAlias("bs_fdPhi" ,"TVector2::Phi_mpi_pi(bsP_1.Phi()-bsP_2.Phi())");
  rChain0->SetAlias("od_fd2Phi","TVector2::Phi_mpi_pi(2.*(unitP_1.Phi()-unitP_2.Phi()))");
  rChain0->SetAlias("bs_fd2Phi","TVector2::Phi_mpi_pi(2.*(bsP_1.Phi()-bsP_2.Phi()))");

  rChain0->SetAlias("2Psi"     ,"TVector2::Phi_mpi_pi(2.*(unitP_fc.Phi()))");

  // sub_Phi1();
  //  sub_Phi2();
  sub_Phi5();
  sub_Phi4();
  
}


