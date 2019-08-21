{
  gROOT->Macro("SetStyle.C");

  Double_t x_d[] = {0.5,1.,2.,3.,4.,4.5};
  Double_t y_d[] = {-0.14,-0.13,-0.13,-0.13,-0.13,-0.03};
  auto gr_dummy = new TGraph(6, x_d, y_d);
  gr_dummy->SetMarkerStyle(0);
  gr_dummy->SetLineColorAlpha(kBlack,0.99);

  TString   x_label[]  = {"1H","2H","3H","3He"};
  Double_t  x_p[]      = {1.,2.,3.,4.};
  Double_t  x_pe[]     = {0.,0.,0.,0.};
  
  //v37.1.4 
  Double_t y_132v2min[] = {-0.0572906, -0.0902588, -0.122411, -0.110268};
  Double_t y_132v2mine[]= {0.0087086,  0.00781448, 0.00833103, 0.00854178};

  Double_t y_108v2min[]  = {-0.050298, -0.0792591, -0.112002, -0.111895};
  Double_t y_108v2mine[] = {0.00973763, 0.00932302, 0.0102602, 0.0105336};

  auto gr_132v2min = new TGraphErrors(4, x_p, y_132v2min, x_pe, y_132v2mine);
  auto gr_108v2min = new TGraphErrors(4, x_p, y_108v2min, x_pe, y_108v2mine);

  for( UInt_t i = 0; i < 4; i++) {
    gr_132v2min->GetXaxis()->SetBinLabel( gr_132v2min->GetXaxis()->FindBin(i+1.), x_label[i]);
    gr_108v2min->GetXaxis()->SetBinLabel( gr_108v2min->GetXaxis()->FindBin(i+1.), x_label[i]);
  }


  gr_132v2min->SetMarkerStyle(20);
  gr_132v2min->SetMarkerColor(2);
  gr_132v2min->SetLineColor(2);

  gr_108v2min->SetMarkerStyle(20);
  gr_108v2min->SetMarkerColor(4);
  gr_108v2min->SetLineColor(4);


  auto m_v2min = new TMultiGraph("",";Fragment; v2_{min}");
  auto leg = new TLegend(0.5,0.7,0.8,0.9);

  m_v2min->Add(gr_132v2min,"ALP");
  m_v2min->Add(gr_108v2min,"LP");
  m_v2min->Add(gr_dummy,"P");

  leg->AddEntry(gr_132v2min,"^{132}Sn+^{124}Sn");
  leg->AddEntry(gr_108v2min,"^{108}Sn+^{112}Sn");

  for( UInt_t i = 0; i < 4; i++) 
    m_v2min->GetXaxis()->SetBinLabel( m_v2min->GetXaxis()->FindBin(i+0.9), x_label[i]);
  m_v2min->GetXaxis()->SetRangeUser(-1.,4.5);


  auto gr_v2mindiff = new TGraphErrors();
  gr_v2mindiff->SetTitle(";Fragment; v2_{min}(^{132}Sn+^{124}Sn) - v2_{min}(^{108}Sn+^{112}Sn)");
  for( UInt_t i = 0; i < 4; i++ ){
    auto y_v1diff  = y_132v2min[i] - y_108v2min[i];
    auto y_v1diffe = sqrt( pow(y_132v2mine[i],2) + pow(y_108v2mine[i],2) );
    gr_v2mindiff->SetPoint(i, x_p[i], y_v1diff);
    gr_v2mindiff->SetPointError(i, x_pe[i], y_v1diffe);
  }

  gr_v2mindiff->SetMarkerStyle(20);
  gr_v2mindiff->SetMarkerColor(2);
  gr_v2mindiff->SetLineColor(2);
  gr_v2mindiff->GetYaxis()->SetTitleSize(0.05);
  for( UInt_t i = 0; i < 4; i++) 
    gr_v2mindiff->GetXaxis()->SetBinLabel( gr_v2mindiff->GetXaxis()->FindBin(i+1.), x_label[i]);

  UInt_t ic = -1;

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  gr_v2mindiff->Draw("AP");

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  m_v2min->Draw("ALP");
  leg->Draw();
  //gr_132v2min->Draw("ALP");
}
