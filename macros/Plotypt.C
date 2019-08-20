{
  gStyle->SetOptStat(0);

  auto _file0 = new TFile("data/cosYPt_132Sn_proton.v23.1.22.root");
  auto _file1 = new TFile("data/cosYPt_132Sn_proton.v23.1.21.root");
  
  TString label[2] = {"v.22 yaw>0","v.21 yaw < 0 "};



  auto leg0 = new TLegend(0.57, 0.64, 0.95, 0.92);
  auto leg1 = new TLegend(0.57, 0.64, 0.95, 0.92);



  TH1D *hevt0 = (TH1D*)_file0->Get("hphi0_180");
  Double_t evt0 = (Double_t)hevt0->GetEntries();

  TH1D *hevt1 = (TH1D*)_file1->Get("hphi0_180");
  Double_t evt1 = (Double_t)hevt1->GetEntries();


  auto cc0 = new TCanvas("cc0","v1",2200,700);
  cc0->Divide(9,2);


  TH1D *hypt1_0[9];
  TH1D *hypt1_1[9];
  TH1D *hrat1[9];
  for(UInt_t i = 0; i < 9; i++ ){
    hypt1_0[i] = (TH1D*)_file0->Get((TString)Form("hypt1_%d",i));
    hypt1_1[i] = (TH1D*)_file1->Get((TString)Form("hypt1_%d",i));
    

    hypt1_0[i]->Scale(1./evt0);
    hypt1_1[i]->Scale(1./evt1);
    hypt1_0[i]->SetLineColor(4);
    hypt1_1[i]->SetLineColor(2);


    cc0->cd(i+1);
    if( i == 8 ) {
      hypt1_1[i]->Draw();
      hypt1_0[i]->Draw("same");
    }
    else {
      hypt1_0[i]->Draw();
      hypt1_1[i]->Draw("same");
    }
    if( i == 4 ) {
      leg0->AddEntry(hypt1_0[i],label[0]);
      leg0->AddEntry(hypt1_1[i],label[1]);
      leg0->Draw();
    }



    hrat1[i] = new TH1D( (*hypt1_0[i]) / (*hypt1_1[i]) );
    cc0->cd(i+10);
    hrat1[i]->Draw();


  }


  TH1D *hypt2_0[6];
  TH1D *hypt2_1[6];
  TH1D *hrat2[6];
  auto cc1 = new TCanvas("cc1","v2",2200,700);
  cc1->Divide(6,2);

  for(UInt_t i = 0; i < 6; i++ ){
    hypt2_0[i] = (TH1D*)_file0->Get((TString)Form("hypt2_%d",i));
    hypt2_1[i] = (TH1D*)_file1->Get((TString)Form("hypt2_%d",i));
    

    hypt2_0[i]->Scale(1./evt0);
    hypt2_1[i]->Scale(1./evt1);
    hypt2_0[i]->SetLineColor(4);
    hypt2_1[i]->SetLineColor(2);


    cc1->cd(i+1);
    hypt2_0[i]->Draw();
    hypt2_1[i]->Draw("same");
    if( i == 3 ) {
      leg1->AddEntry(hypt2_0[i],label[0]);
      leg1->AddEntry(hypt2_1[i],label[1]);
      leg1->Draw();
    }
    
      



    hrat2[i] = new TH1D( (*hypt2_0[i]) / (*hypt2_1[i]) );
    cc1->cd(i+7);
    hrat2[i]->Draw();

  }

}
