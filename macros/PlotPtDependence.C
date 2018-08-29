TCanvas *cc0;
TCanvas *cc1;
TCanvas *cc2;

// --> configuration
TString fsys[] = {"132Sn+124Sn","108Sn+112Sn","124Sn+112Sn","112Sn+124Sn"};
TString rsys[] = {"132",        "108",        "124",        "112"};
TString fpid[] = {"proton","deuteron","triton","neutron"};  
UInt_t  imrk[] = {20, 21, 22, 23, 21};
Size_t  imsz[] = {1, 1, 1.3, 1.3, 1.3};
Color_t icol[][4] = { {kRed,          kBlue,  kOrange-3,   kGreen+1}, 
		      {kBlue+2,   kOrange+7,  kGreen-3,     kPink+9},
		      {kGreen-3,    kPink+7,  kCyan-1,    kYellow-2},
		      {kOrange-2, kGreen+3,   kCyan-1,    kBlue} };

Double_t FittingAndIntegral(TGraphErrors *gr)
{
  TF1 *p1 = new TF1("p1","[0]+[1]*x",0.,800.);
  
  gr->Fit("p1","","",200.,700.);
  
  return p1->Integral(0., 800) / 800.;
}

void PlotPtDependence()
{
  // --> Plotting selection
  Bool_t bsys[]  = { 1, 1, 0, 0};
  Bool_t bpid[]  = { 1, 0, 0, 1};


  TString fname[4][4];
  for(UInt_t is = 0; is < 4; is++){
    if( !bsys[is] ) continue;

    for(UInt_t ip = 0; ip < 3; ip++){
      if(bpid[ip])
	fname[is][ip] = "data/YPT"+rsys[is]+"Sn_"+fpid[ip]+".root";
    }

    if(bpid[3])
      fname[is][3] = "data/NL"+rsys[is]+"Sn_neutron.root";
  }

  //  Double_t rrange[3] = {0.27, 0.54, 1.};
  Double_t yrange[] =  {-0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35, 0.45};
  const UInt_t nyrange = sizeof(yrange)/sizeof(Double_t);

  Double_t yrange2[] = {-0.25, -0.1,  0.1, 0.25, 0.45};
  const UInt_t nyrange2 = sizeof(yrange2)/sizeof(Double_t);



  TString rapRange[nyrange];
  for(UInt_t i = 1; i < nyrange - 1; i++)
    rapRange[i] = Form(" %f <= Rapidity < %f", yrange[i-1],yrange[i]);
  rapRange[0] = Form(" Rapidity < %f", yrange[0]);
  rapRange[nyrange-1] = Form(" %f <= Rapidity", yrange[nyrange-2]);

  TString rapRange2[nyrange2];
  for(UInt_t i = 1; i < nyrange2 - 1; i++)
    rapRange2[i] = Form(" %f <= Rapidity < %f", yrange2[i-1],yrange2[i]);
  rapRange2[0] = Form(" Rapidity < %f", yrange2[0]);
  rapRange2[nyrange2-1] = Form(" %f <= Rapidity", yrange2[nyrange2-2]);


  const UInt_t rbin = nyrange;
  const UInt_t rbin2 = nyrange2;

  TFile *fOpen;
  TGraphErrors *gr_v1[rbin*4];
  TGraphErrors *gr_v2[rbin*4];
  TMultiGraph  *mv1[rbin];
  TMultiGraph  *mv2[rbin];
  TLegend      *lg1[rbin];
  TLegend      *lg2[rbin];
  TGraphErrors *gr_v1f;


  for(UInt_t k = 0; k < rbin; k++){
    mv1[k] = new TMultiGraph((TString)Form("mv1%d",k),";Pt [MeV/c]; v1");
    mv1[k]->SetTitle(rapRange[k]+";Pt [MeV/c]; v1");
    mv2[k] = new TMultiGraph((TString)Form("mv2%d",k),";Pt [MeV/c]; v2");
    mv2[k]->SetTitle(rapRange2[k]+";Pt [MeV/c]; v1");
    lg1[k] = new TLegend(0.4 , 0.15, 0.9 ,  0.4,"");
    lg2[k] = new TLegend(0.12, 0.14, 0.58, 0.35,"");

    gr_v1f = new TGraphErrors();
    gr_v1f->SetName("gr_v1f");
  }
  
  UInt_t igr = 0;

  for(UInt_t is = 0; is < 4; is++){
    
    if( !bsys[is] ) continue;

    for(UInt_t ip = 0; ip < 4; ip++){

      if( !bpid[ip] ) continue;

      fOpen = TFile::Open(fname[is][ip]); 

      if( fOpen == NULL ) continue;
      else
	std::cout << fname[is][ip] << " is opened. " << std::endl;

      fOpen->ls();
    
      for(UInt_t k = 0; k < rbin; k++){
	TGraphErrors *gv1 = (TGraphErrors*)fOpen->Get((TString)Form("gPt_v1%d",k));

	if( gv1 == NULL ) {
	  cout << " not found " << k << endl;
	  continue;
	}
	TString gvname = Form("gPt_v1%d",igr);
	gr_v1[igr] = (TGraphErrors*)gv1->Clone(gvname);

	mv1[k]->Add(gr_v1[igr],"lp");

	gr_v1[igr]->SetMarkerStyle(imrk[is]);
	gr_v1[igr]->SetMarkerColor(icol[ip][is]);
	gr_v1[igr]->SetMarkerSize(imsz[is]);
	gr_v1[igr]->SetLineColor(icol[ip][is]);
	lg1[k]->AddEntry(gr_v1[igr], rsys[is]+" "+fpid[ip] ,"lp");

      }

      for(UInt_t k = 0; k < rbin2; k++){

	TGraphErrors *gv2 = (TGraphErrors*)fOpen->Get((TString)Form("gPt_v2%d",k));
	TString gvname = Form("gPt_v2%d",igr);
	gr_v2[igr] = (TGraphErrors*)gv2->Clone(gvname);


	mv2[k]->Add(gr_v2[igr],"lp");

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
  cc1->Divide(rbin2,1);
  for(UInt_t k = 0; k < rbin2; k++){
    cc1->cd(k+1);
    mv2[k]->SetMaximum( 0.2);
    mv2[k]->SetMinimum(-0.2);
    mv2[k]->Draw("ALP");
    if( k == rbin2-1 ) lg2[k]->Draw();
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

  // cc2 = new TCanvas("cc2","v1fit");
  // gr_v1f->SetMarkerStyle(20);
  // gr_v1f->Draw("ALP");
}


void SaveCanvas(TString cnt = "")
{
  cc0->SaveAs("v1PtSn"+cnt+".png");
  cc1->SaveAs("v2PtSn"+cnt+".png");
}

