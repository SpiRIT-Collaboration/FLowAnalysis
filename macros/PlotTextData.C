TGraph* GetTextGraph(TString fname);

void PlotTextData()
{
  auto mgrph = new TMultiGraph();

  TString dir = "data/FOPI/NPA876";
  gSystem->cd(dir);

  TGraph* gr = GetTextGraph("Fig7_v1Ut_0.4_proton.txt");
  gr -> SetName("gr_p");
  gr -> SetMarkerStyle(20);
  gr -> SetMarkerColor(2);
  gr -> SetLineColor(2);
  mgrph->Add(gr,"lp");

  gr = GetTextGraph("Fig7_v1Ut_0.4_deuteron.txt");
  gr -> SetName("gr_d");
  gr -> SetMarkerStyle(20);
  gr -> SetMarkerColor(4);
  gr -> SetLineColor(4);
  mgrph->Add(gr,"lp");

  gr = GetTextGraph("Fig7_v1Ut_0.4_A3.txt");
  gr -> SetName("gr_a3");
  gr -> SetMarkerStyle(20);
  gr -> SetMarkerColor(7);
  gr -> SetLineColor(7);
  mgrph->Add(gr,"lp");

  gr = GetTextGraph("Fig7_v1Ut_0.4_4He.txt");
  gr -> SetName("gr_4he");
  gr -> SetMarkerStyle(20);
  gr -> SetMarkerColor(6);
  gr -> SetLineColor(6);
  mgrph->Add(gr,"lp");


  
  auto cc = new TCanvas("cc","cc");
  mgrph->Draw("ALP");

  gSystem->cd("../../../");
}


TGraph* GetTextGraph(TString fname)
{
  std::fstream fread;
  fread.open(fname, std::fstream::in);

  auto vut = new TGraph();

  Double_t x, y;
  TString sget;
  UInt_t  in = 0;
  while( !fread.eof() ) {
    fread >> sget;
    x = (Double_t)atof(sget);
    fread >> sget;
    y = (Double_t)atof(sget);

    if( !std::isnan(x) ) {
      vut -> SetPoint(in, x, y);
      in++;
    }
  }
  vut -> RemovePoint(in-1);

  return vut;
  
}
