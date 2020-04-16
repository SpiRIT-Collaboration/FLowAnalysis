
void PlotTuQMD()
{

  auto gv1_prot = new TGraphErrors();
  gv1_prot -> SetName("gv1_prot");
  auto gv2_prot = new TGraphErrors();
  gv2_prot -> SetName("gv2_prot");

  auto gv1_deut = new TGraphErrors();
  gv1_deut -> SetName("gv1_deut");
  auto gv2_deut = new TGraphErrors();
  gv1_deut -> SetName("gv2_deut");

  auto gv1_trit = new TGraphErrors();
  gv1_trit -> SetName("gv1_trit");
  auto gv2_trit = new TGraphErrors();
  gv1_trit -> SetName("gv2_trit");

  auto gv1_He3 = new TGraphErrors();
  gv1_He3 -> SetName("gv1_He3");
  auto gv2_He3 = new TGraphErrors();
  gv1_He3 -> SetName("gv2_He3");


  for( UInt_t i = 0; i < 30; i++ ) {
    gv1_prot -> SetPoint(i, protTQMD[i][0], -protTQMD[i][3]);
    gv1_prot -> SetPointError(i, 0., protTQMD[i][4]);
    gv2_prot -> SetPoint(i, protTQMD[i][0], protTQMD[i][5]);
    gv2_prot -> SetPointError(i, 0., protTQMD[i][6]);

    gv1_deut -> SetPoint(i, deutTQMD[i][0], -deutTQMD[i][3]);
    gv1_deut -> SetPointError(i, 0., deutTQMD[i][4]);
    gv2_deut -> SetPoint(i, deutTQMD[i][0], deutTQMD[i][5]);
    gv2_deut -> SetPointError(i, 0., deutTQMD[i][6]);


    gv1_trit -> SetPoint(i, tritTQMD[i][0], -tritTQMD[i][3]);
    gv1_trit -> SetPointError(i, 0., tritTQMD[i][4]);
    gv2_trit -> SetPoint(i, tritTQMD[i][0], tritTQMD[i][5]);
    gv2_trit -> SetPointError(i, 0., tritTQMD[i][6]);

    gv1_He3 -> SetPoint(i, He3TQMD[i][0], -He3TQMD[i][3]);
    gv1_He3 -> SetPointError(i, 0., He3TQMD[i][4]);
    gv2_He3 -> SetPoint(i, He3TQMD[i][0], He3TQMD[i][5]);
    gv2_He3 -> SetPointError(i, 0., He3TQMD[i][6]);

  }

  UInt_t ic = 0; UInt_t id = 1;
  TCanvas *cc;

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  auto legv1 = new TLegend(0.03, 0.7, 0.4, 0.98, "");

  gv1_He3  -> SetLineColor(6);
  gv1_He3  -> Draw("ALP");

  gv1_prot -> SetLineColor(2);
  gv1_prot -> Draw("same");

  gv1_deut -> SetLineColor(8);
  gv1_deut -> Draw("same");

  gv1_trit -> SetLineColor(4);
  gv1_trit -> Draw("same");

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  gv2_He3  -> SetLineColor(6);
  gv2_He3  -> Draw("ALP");

  gv2_prot -> SetLineColor(2);
  gv2_prot -> Draw("same");

  gv2_deut -> SetLineColor(8);
  gv2_deut -> Draw("same");

  gv2_trit -> SetLineColor(4);
  gv2_trit -> Draw("same");


  TString fHeader = "tuqmd_100Sn";
  TFile *fopen;;
  
  fopen = TFile::Open(fHeader+"_proton.v0.root","recreate");
  gv1_prot->SetName("gu_v1");
  gv2_prot->SetName("gu_v2");
  gv1_prot->Write();
  gv2_prot->Write();
  fopen->Close();

  fopen = TFile::Open(fHeader+"_deuteron.v0.root","recreate");
  gv1_deut->SetName("gu_v1");
  gv2_deut->SetName("gu_v2");
  gv1_deut->Write();
  gv2_deut->Write();
  fopen->Close();

  fopen = TFile::Open(fHeader+"_triton.v0.root","recreate");
  gv1_trit->SetName("gu_v1");
  gv2_trit->SetName("gu_v2");
  gv1_trit->Write();
  gv2_trit->Write();
  fopen->Close();

  fopen = TFile::Open(fHeader+"_3He.v0.root","recreate");
  gv1_He3->SetName("gu_v1");
  gv2_He3->SetName("gu_v2");
  gv1_He3->Write();
  gv2_He3->Write();
  fopen->Close();
}


