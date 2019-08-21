void PlotSystemDependence()
{
  TString  sysname[] = {"132SN","108Sn","124Sn","112Sn"};
  Double_t sysdlt[]  = {0.22,    0.09,   0.15,   0.15};
  Double_t sysA[]    = {256.,    220.,   236.,   236};

  auto mv2 = new TGraphErrors();

  UInt_t ip = 0;
  mv2->SetPoint(ip, sysdlt[ip], -0.0279831 );
  mv2->SetPointError(ip, 0., 0.000749319); ip++;

  mv2->SetPoint(ip, sysdlt[ip], -0.0265948);
  mv2->SetPointError(ip, 0., 0.000958754); ip++;

  mv2->SetPoint(ip, sysdlt[ip], -0.0293554);
  mv2->SetPointError(ip, 0., 0.00106792); ip++;

  // mv2->SetPoint(ip, sysdlt[ip], -0.0205363);
  // mv2->SetPointError(ip, 0., 0.00602508); ip++;


  TCanvas *cc;
  
  UInt_t ic = 0;
  cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic)); ic++;

  mv2->SetMarkerStyle(20);
  mv2->SetMarkerColor(2);
  mv2->SetLineColor(2);

  //  mv2->SetTitle("; (N-P)/A; v2");
  mv2->SetTitle("; A(Beam+Target); v2");
  mv2->Draw("AP");

}
