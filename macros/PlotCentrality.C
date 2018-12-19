void previous()
{
  //132Sn _rf.v4.7.0.cv0.
  Double_t mult[]  = {4., 12., 20., 28., 36.,};
  Double_t multe[] = {4.,  4.,  4.,  4.,  4.,};

  Double_t mcos[]  = {0.06667,0.103186,0.130410,0.13108,0.12151};
  Double_t mcose[] = {0.00219,0.000552,0.000423,0.00075,0.00351};


  auto gv_mcos = new TGraphErrors(5, mult,mcos,multe,mcose);
  gv_mcos->SetName("gv_mcos");
  gv_mcos->SetTitle("; Multiplicity ; <cos #Delta#Phi_sub> ");

  gv_mcos->SetMarkerStyle(20);
  gv_mcos->SetMarkerColor(4);
  gv_mcos->SetLineColor(4);

  auto cc0 = new TCanvas("cc0","cc0");
  gv_mcos->Draw("ALP");

}

void v12cos()
{
  Double_t mult[]  = {39.1,      31.0,    21.2};
  Double_t multe[] = {3.4 ,      2.,       4.3};
  Double_t cos1[]  = {0.818874,0.865476, 0.743198};
  Double_t cos1e[] = {0.002623,0.002272, 0.001989};
  Double_t cos2[]  = {0.426889,0.476859, 0.351632};
  Double_t cos2e[] = {0.002738,0.002507, 0.001885};


  auto gv_mcos1 = new TGraphErrors(3, mult, cos1, multe, cos1e);
  auto gv_mcos2 = new TGraphErrors(3, mult, cos2, multe, cos2e);
  

  auto m_cos = new TMultiGraph("mcos",";Multiplicity; <cos>");
  m_cos->Add(gv_mcos1,"LP");
  m_cos->Add(gv_mcos2,"LP");

  m_cos->Draw("AP");

}

void PlotCentrality()
{
  v12cos();
}

