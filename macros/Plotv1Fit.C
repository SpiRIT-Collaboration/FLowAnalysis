{
  gROOT->Macro("SetStyle.C");

  Double_t x_d[] = {0.5,1.,2.,3.,4.,4.5};
  Double_t y_d[] = {0.2,0.2,0.2,0.2,0.2,0.8};
  auto gr_dummy = new TGraph(6, x_d, y_d);
  gr_dummy->SetMarkerStyle(0);
  gr_dummy->SetLineColorAlpha(kBlack,0.99);

  TString   x_label[]  = {"1H","2H","3H","3He"};
  Double_t  x_p[]      = {1.,2.,3.,4.};
  Double_t  x_pe[]     = {0.,0.,0.,0.};
  
  //v37.1.4 
  Double_t y_132v1slop[] = {0.336699, 0.498031, 0.677569, 0.757752};
  Double_t y_132v1slope[]= {0.00258905, 0.00316241, 0.00379723, 0.0059261};

  Double_t y_108v1slop[]  = {0.32761, 0.46231, 0.619141, 0.733029};
  Double_t y_108v1slope[] = {0.00268186, 0.00354251, 0.00485736, 0.00663348};

  auto gr_132v1slp = new TGraphErrors(4, x_p, y_132v1slop, x_pe, y_132v1slope);
  auto gr_108v1slp = new TGraphErrors(4, x_p, y_108v1slop, x_pe, y_108v1slope);

  for( UInt_t i = 0; i < 4; i++) {
    gr_132v1slp->GetXaxis()->SetBinLabel( gr_132v1slp->GetXaxis()->FindBin(i+1.), x_label[i]);
    gr_108v1slp->GetXaxis()->SetBinLabel( gr_108v1slp->GetXaxis()->FindBin(i+1.), x_label[i]);
  }


  gr_132v1slp->SetMarkerStyle(20);
  gr_132v1slp->SetMarkerColor(2);
  gr_132v1slp->SetLineColor(2);

  gr_108v1slp->SetMarkerStyle(20);
  gr_108v1slp->SetMarkerColor(4);
  gr_108v1slp->SetLineColor(4);


  auto m_v1slp = new TMultiGraph("",";Fragment; v_{11}");
  auto leg = new TLegend(0.5,0.2,0.8,0.4);

  m_v1slp->Add(gr_132v1slp,"ALP");
  m_v1slp->Add(gr_108v1slp,"LP");
  m_v1slp->Add(gr_dummy,"P");

  leg->AddEntry(gr_132v1slp,"^{132}Sn+^{124}Sn");
  leg->AddEntry(gr_108v1slp,"^{108}Sn+^{112}Sn");

  for( UInt_t i = 0; i < 4; i++) 
    m_v1slp->GetXaxis()->SetBinLabel( m_v1slp->GetXaxis()->FindBin(i+0.9), x_label[i]);
  m_v1slp->GetXaxis()->SetRangeUser(-1.,4.5);


  auto gr_v1slpdiff = new TGraphErrors();
  gr_v1slpdiff->SetTitle(";Fragment; v_{11}(^{132}Sn+^{124}Sn) - v_{11}(^{108}Sn+^{112}Sn)");
  for( UInt_t i = 0; i < 4; i++ ){
    auto y_v1diff  = y_132v1slop[i] - y_108v1slop[i];
    auto y_v1diffe = sqrt( pow(y_132v1slope[i],2) + pow(y_108v1slope[i],2) );
    gr_v1slpdiff->SetPoint(i, x_p[i], y_v1diff);
    gr_v1slpdiff->SetPointError(i, x_pe[i], y_v1diffe);
  }

  gr_v1slpdiff->SetMarkerStyle(20);
  gr_v1slpdiff->SetMarkerColor(2);
  gr_v1slpdiff->SetLineColor(2);
  gr_v1slpdiff->GetYaxis()->SetTitleSize(0.05);
  for( UInt_t i = 0; i < 4; i++) 
    gr_v1slpdiff->GetXaxis()->SetBinLabel( gr_v1slpdiff->GetXaxis()->FindBin(i+1.), x_label[i]);

  UInt_t ic = -1;

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  gr_v1slpdiff->Draw("AP");

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  m_v1slp->Draw("ALP");
  leg->Draw();
  //gr_132v1slp->Draw("ALP");
}
