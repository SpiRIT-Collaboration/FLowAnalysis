{

  auto c0 = new TCanvas("c0","c0"); 
  c0->Divide(2,2);

  c0->cd(1);
  auto hx = new TH1D("hx","hx",200,-30,30);
  rChain->Project("hx","fvertex.X()","mtrack4>0","");
  hx->Fit("gaus");

  c0->cd(2);
  auto hy = new TH1D("hy","hy",200,-260,-200);
  rChain->Project("hy","fvertex.Y()","","colz");
  hy->Fit("gaus");

  c0->cd(3);
  auto hz = new TH1D("hz","hz",200,-30,30);
  rChain->Project("hz","fvertex.Z()","mtrack4>0","");
  hz->Fit("gaus");
}
