void PlotNLProperty()
{
  Double_t vpt1[50];
  Double_t vpt2[50];
  Double_t vrap1[50];;
  Double_t vrap2[50];;


  Double_t mass = 939.565346;

  Double_t thta[2] = {0.4, 0.65};
  Double_t dist = 8900.;


  Double_t beta;
  Double_t gamm;
  Double_t E;
  Double_t p;
  for(UInt_t i = 0; i < 50; i++){

    beta =  i * 0.012;
    gamm = 1./sqrt(1. - beta*beta);
    
    p  = mass * beta * gamm;
    E  = mass * gamm;
    vpt1[i] = p * sin(thta[0]);
    vpt2[i] = p * sin(thta[1]);

    vrap1[i] = 0.5 * log( (E + p*cos(thta[0]))/(E - p*cos(thta[0])));
    vrap2[i] = 0.5 * log( (E + p*cos(thta[1]))/(E - p*cos(thta[1])));

  }
  

  auto mgr = new TMultiGraph();

  auto gt1 = new TGraph(50, vrap1, vpt1);
  auto gt2 = new TGraph(50, vrap2, vpt2);

  gt1->SetName("gt1");
  gt2->SetName("gt2");


  mgr->Add(gt1,"");
  mgr->Add(gt2,"");

  mgr->Draw("al");

}
