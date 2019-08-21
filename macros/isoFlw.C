{
  UInt_t iccv = 0;
  // auto ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv));
  // rChain0->Draw("eisobm.X():eisobm.Y()>>h132(200,0,4000,200,0,10000)","","colz");

  // iccv++;
  // ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv));
  // rChain1->Draw("eisobm.X():eisobm.Y()>>h108(200,0,4000,200,0,10000)","","colz");




  UInt_t hiv = 0;
  auto h132v = new TH1D(Form("h132v%d",hiv),"Y > Y_{cm} ; iso Trans/Long ", 200, 0., 2.);
  rChain0->Project(Form("h132v%d",hiv),"eisobm.Y()/eisobm.X()");

  auto h108v = new TH1D(Form("h108v%d",hiv),"Y > Y_{cm} ; iso Trans/Long ", 200, 0., 2.);
  rChain1->Project(Form("h108v%d",hiv),"eisobm.Y()/eisobm.X()");

  iccv++;
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv));
  ccv -> SetTitle(" all Y > Y_cm");
  auto f132 = (TH1D*)gROOT->Get(Form("h132v%d",hiv));
  auto f108 = (TH1D*)gROOT->Get(Form("h108v%d",hiv));
  
  f132 -> SetNormFactor(1);
  f108 -> SetNormFactor(1);
  f132 -> SetLineColor(2);
  f108 -> SetLineColor(4);
  f108 -> Draw();
  f132 -> Draw("same");


  hiv++;
  h132v = new TH1D(Form("h132v%d",hiv),"Y < Y_{cm} ; iso Trans/Long ", 200, 0., 5.);
  rChain0->Project(Form("h132v%d",hiv),"eisobm.Y()/eisobm.X()");

  h108v = new TH1D(Form("h108v%d",hiv),"Y < Y_{cm} ; iso Trans/Long ", 200, 0., 5.);
  rChain1->Project(Form("h108v%d",hiv),"eisobm.Y()/eisobm.X()");

  iccv++;
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv));
  ccv -> SetTitle(" all Y<Y_cm");
  
  f132 = (TH1D*)gROOT->Get(Form("h132v%d",hiv));
  f108 = (TH1D*)gROOT->Get(Form("h108v%d",hiv));

  f132 -> SetNormFactor(1);
  f108 -> SetNormFactor(1);
  f132 -> SetLineColor(2);
  f108 -> SetLineColor(4);
  f108 -> Draw();
  f132 -> Draw("same");



  TCut mcut[4];
  mcut[0] = "ntrack[4]>=  0 && ntrack[4] < 20";
  mcut[1] = "ntrack[4]>= 20 && ntrack[4] < 35";
  mcut[2] = "ntrack[4]>= 35 && ntrack[4] < 40";
  mcut[3] = "ntrack[4]>= 40";


  UInt_t shiv = hiv+1;
  for(UInt_t i = 0; i < 4; i++){
    hiv++;

    h132v = new TH1D(Form("h132v%d",hiv),"Y > Y_{cm} ; iso Trans/Long ", 200, 0., 2.);
    rChain0->Project(Form("h132v%d",hiv),"eisobm.Y()/eisobm.X()",mcut[i]);

    h108v = new TH1D(Form("h108v%d",hiv),"Y > Y_{cm} ; iso Trans/Long ", 200, 0., 2.);
    rChain1->Project(Form("h108v%d",hiv),"eisobm.Y()/eisobm.X()",mcut[i]);
  }

  iccv++;
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv));
  ccv->SetTitle("iso ratio Y>Ycm");
  ccv->Divide(2,2);
  UInt_t idd = 1;

  for(UInt_t i = 0; i < 4; i++){
    ccv->cd(idd); idd++;

    f132 = (TH1D*)gROOT->Get(Form("h132v%d",i+shiv));
    f108 = (TH1D*)gROOT->Get(Form("h108v%d",i+shiv));

    f132 -> SetNormFactor(1);
    f108 -> SetNormFactor(1);

    f132 -> SetLineColor(2);
    f108 -> SetLineColor(4);

    f108->SetTitle(mcut[i].GetTitle());

    f108->Draw();
    f132->Draw("same");
  }


  shiv = hiv+1;
  for(UInt_t i = 0; i < 4; i++){
    hiv++;
    
    h132v = new TH1D(Form("h132v%d",hiv),"Y < Y_{cm} ; iso Trans/Long ", 200, 0., 5.);
    rChain0->Project(Form("h132v%d",hiv),"eisotg.Y()/eisotg.X()",mcut[i]);

    h108v = new TH1D(Form("h108v%d",hiv),"Y < Y_{cm} ; iso Trans/Long ", 200, 0., 5.);
    rChain1->Project(Form("h108v%d",hiv),"eisotg.Y()/eisotg.X()",mcut[i]);
  }

  iccv++;
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv));
  ccv->SetTitle("iso ratio Y<Ycm");

  ccv->Divide(2,2);
  idd = 1;

  for(UInt_t i = 0; i < 4; i++){
    ccv->cd(idd); idd++;

    f132 = (TH1D*)gROOT->Get(Form("h132v%d",i+shiv));
    f108 = (TH1D*)gROOT->Get(Form("h108v%d",i+shiv));

    f132 -> SetNormFactor(1);
    f108 -> SetNormFactor(1);

    f132 -> SetLineColor(2);
    f108 -> SetLineColor(4);

    f108->SetTitle(mcut[i].GetTitle());

    f108->Draw();
    f132->Draw("same");
  }


}
