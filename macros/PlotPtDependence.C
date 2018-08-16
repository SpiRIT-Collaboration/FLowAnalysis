TCanvas *cc0;
TCanvas *cc1;

// --> configuration
TString fsys[4] = {"132Sn+124Sn","108Sn+112Sn","124Sn+112Sn","112Sn+124Sn"};
TString rsys[4] = {"132",        "108",        "124",        "112"};
TString fpid[3] = {"proton","deuteron","triton"};  
UInt_t  imrk[4] = {20, 21, 22, 23};
Size_t  imsz[4] = {1, 1, 1.3, 1.3};
Color_t icol[3][4]= { {kRed,          kBlue,  kOrange-3,   kGreen+1}, 
		      {kBlue+2,   kOrange+7,  kGreen-3,     kPink+9},
		      {kGreen-3,    kPink+7,  kCyan-1,    kYellow-2} };

void PlotPtDependence()
{
  // --> Plotting selection
  Bool_t bsys[4]  = { 1, 1, 1, 1};
  Bool_t bpid[3]  = { 1, 0, 0};

  TString fname[4][3];
  for(UInt_t is = 0; is < 4; is++){
    for(UInt_t ip = 0; ip < 3; ip++){
      fname[is][ip] = "data/PT"+rsys[is]+"Sn_"+fpid[ip]+".root";
    }
  }

  //  Double_t rrange[3] = {0.27, 0.54, 1.};
  TString  rapRange[] = {"Rapidity < 0.3 ",
			 "0.3 <= Rapidity < 0.54",
			 "0.54 <= Rapidity < 0.74",
			 "0.74 <= Rapidity < 1.2"};

  const UInt_t rbin = 4;

  TFile *fOpen;
  TGraphErrors *gr_v1[rbin*4];
  TGraphErrors *gr_v2[rbin*4];
  TMultiGraph  *mv1[rbin];
  TMultiGraph  *mv2[rbin];
  TLegend      *lg1[rbin];
  TLegend      *lg2[rbin];

  for(UInt_t k = 0; k < rbin; k++){
    mv1[k] = new TMultiGraph((TString)Form("mv1%d",k),";Pt [MeV/c]; v1");
    mv1[k]->SetTitle(rapRange[k]+";Pt [MeV/c]; v1");
    mv2[k] = new TMultiGraph((TString)Form("mv2%d",k),";Pt [MeV/c]; v2");
    mv2[k]->SetTitle(rapRange[k]+";Pt [MeV/c]; v1");
    lg1[k] = new TLegend(0.4 , 0.15, 0.9 ,  0.4,"");
    lg2[k] = new TLegend(0.12, 0.14, 0.58, 0.35,"");
  }
  
  UInt_t igr = 0;

  for(UInt_t is = 0; is < 4; is++){
    
    if( !bsys[is] ) continue;

    for(UInt_t ip = 0; ip < 3; ip++){

      if( !bpid[ip] ) continue;

      fOpen = TFile::Open(fname[is][ip]); 

      if( fOpen == NULL ) continue;
      else
	std::cout << fname[is][ip] << " is opened. " << std::endl;

    
      for(UInt_t k = 0; k < rbin; k++){
	TGraphErrors *gv1 = (TGraphErrors*)fOpen->Get((TString)Form("gPt_v1%d",k));
	TGraphErrors *gv2 = (TGraphErrors*)fOpen->Get((TString)Form("gPt_v2%d",k));
  
	TString gvname = Form("gPt_v1%d",igr);
	gr_v1[igr] = (TGraphErrors*)gv1->Clone(gvname);
	gvname = Form("gPt_v2%d",igr);
	gr_v2[igr] = (TGraphErrors*)gv2->Clone(gvname);

	mv1[k]->Add(gr_v1[igr],"lp");
	mv2[k]->Add(gr_v2[igr],"lp");

	gr_v1[igr]->SetMarkerStyle(imrk[is]);
	gr_v1[igr]->SetMarkerColor(icol[ip][is]);
	gr_v1[igr]->SetMarkerSize(imsz[is]);
	gr_v1[igr]->SetLineColor(icol[ip][is]);
	lg1[k]->AddEntry(gr_v1[igr], rsys[is]+" "+fpid[ip] ,"lp");


	gr_v2[igr]->SetMarkerStyle(imrk[is]);
	gr_v2[igr]->SetMarkerColor(icol[ip][is]);
	gr_v2[igr]->SetMarkerSize(imsz[is]);
	gr_v2[igr]->SetLineColor(icol[ip][is]);
	lg2[k]->AddEntry(gr_v2[igr], rsys[is]+" "+fpid[ip],"lp");

      }
    

      fOpen->Close();
      igr++;
    }
  }


  

  cc1 = new TCanvas("cc1","v2",110,2280,1000,500);
  cc1->Divide(rbin,1);
  for(UInt_t k = 0; k < rbin; k++){
    cc1->cd(k+1);
    mv2[k]->SetMaximum( 0.02);
    mv2[k]->SetMinimum(-0.07);
    mv2[k]->Draw("ALP");
    if( k == rbin-1 ) lg2[k]->Draw();
  }

  cc0 = new TCanvas("cc0","v1",110,1280,1000,500);
  cc0->Divide(rbin,1);
  for(UInt_t k = 0; k < rbin; k++){
    cc0->cd(k+1);
    mv1[k]->SetMaximum(0.4);
    mv1[k]->SetMinimum(-0.25);
    mv1[k]->Draw("ALP");
    if( k == rbin-1 )
      lg1[k]->Draw();
  }

}


void SaveCanvas(TString cnt = "")
{
  cc0->SaveAs("v1PtSn"+cnt+".png");
  cc1->SaveAs("v2PtSn"+cnt+".png");
}

