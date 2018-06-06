void PlotCentrality()
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
