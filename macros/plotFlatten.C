{

  auto hm1 = new TH1D("hm1","",100,-3.15,3.15);
  rChain->Project("hm1","unitP_fc.Phi()","mtrack4>0");
  hm1->SetNormFactor(100);

  auto hm2 = new TH1D("hm2","",100,-3.15,3.15);
  rChain->Project("hm2","unitP.Phi()","mtrack4>0");
  hm2->SetNormFactor(100);
  hm2->SetLineColor(2);

  TCanvas *c1 = new TCanvas("c1","c1");
  hm2->Draw("e");
  hm1->Draw("samee");


  auto hk1 = new TH1D("hk1","",100,-3.15,3.15);
  rChain->Project("hk1","frpv.Phi()","mtrack4>2");
  hk1->SetNormFactor(100);

  auto hk2 = new TH1D("hk2","",100,-3.15,3.15);
  rChain->Project("hk2","frpphi","mtrack4>2");
  hk2->SetNormFactor(100);
  hk2->SetLineColor(2);

  TCanvas *c4 = new TCanvas("c4","c4");
  hk2->Draw("e");
  hk1->Draw("samee");
  

  TCanvas *c2 = new TCanvas("c2","c2");
  c2->Divide(1,2);
  c2->cd(1);
  rChain->Draw("frpphi:mtrack4>>hrpm0(90,0,90,100,-3.15,3.15)","mtrack4>2","colz");
  c2->cd(2);
  rChain->Draw("frpv.Phi();:mtrack4>>hrpm1(90,0,90,100,-3.15,3.15)","mtrack4>2","colz");


  TCanvas *c3 = new TCanvas("c3","c3");
  c3->Divide(1,2);
  c3->cd(1);
  rChain->Draw("unitP.Phi():mtrack4>>hrpm2(90,0,90,100,-3.15,3.15)","mtrack4>2","colz");
  c3->cd(2);
  rChain->Draw("unitP_fc.Phi():mtrack4>>hrpm3(90,0,90,100,-3.15,3.15)","mtrack4>2","colz");

}
