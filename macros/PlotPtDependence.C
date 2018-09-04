TCanvas *cc0;
TCanvas *cc1;
TCanvas *cc2;
TCanvas *cv[8];

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

//  Double_t rrange[3] = {0.27, 0.54, 1.};
Double_t yrange[] =  {-0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35, 0.45};
const UInt_t rbin = sizeof(yrange)/sizeof(Double_t);

Double_t yrange2[] = {-0.25, -0.1,  0.1, 0.25, 0.45};
const UInt_t rbin2 = sizeof(yrange2)/sizeof(Double_t);

TMultiGraph  *mv1[rbin];
TMultiGraph  *mv2[rbin2];
TLegend      *lg1[rbin];
TLegend      *lg2[rbin2];

void RemovePoints(TString fn, TGraphErrors &gr);

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




  TString rapRange[rbin];
  for(UInt_t i = 1; i < rbin - 1; i++)
    rapRange[i] = Form(" %f <= Y_cm < %f", yrange[i-1],yrange[i]);
  rapRange[0] = Form(" Y_cm < %f", yrange[0]);
  rapRange[rbin-1] = Form(" %f <= Y_cm", yrange[rbin-2]);

  TString rapRange2[rbin2];
  for(UInt_t i = 1; i < rbin2 - 1; i++)
    rapRange2[i] = Form(" %f <= Y_cm < %f", yrange2[i-1],yrange2[i]);
  rapRange2[0] = Form(" Y_cm < %f", yrange2[0]);
  rapRange2[rbin2-1] = Form(" %f <= Y_cm", yrange2[rbin2-2]);



  TFile *fOpen;
  TGraphErrors *gr_v1;
  TGraphErrors *gr_v2;
  TGraphErrors *gr_v1f;


  for(UInt_t k = 0; k < rbin; k++){
    mv1[k] = new TMultiGraph();
    mv1[k]->SetName((TString)Form("mv1%d",k));
    mv1[k]->SetTitle(rapRange[k]+";Pt [MeV/c]; v1");

    //    mv1[k]->GetXaxis()->SetRangeUser(0.,800);
    //    mv1[k]->GetYaxis()->SetLabelSize(0.02);
    
    //    mv1[k]->GetXaxis()->SetLabelSize(0.04);
    // mv1[k]->GetXaxis()->SetLabelOffset(0.0);
    // mv1[k]->GetYaxis()->SetLabelSize(0.05);
    // mv1[k]->GetYaxis()->SetLabelOffset(0.0);
    
    lg1[k] = new TLegend(0.33, 0.14, 0.96, 0.4,"");
    lg1[k]->SetTextSize(0.07);
  }

  for(UInt_t k = 0; k < rbin2; k++){
    mv2[k] = new TMultiGraph();
    mv2[k]->SetName((TString)Form("mv2%d",k));
    mv2[k]->SetTitle(rapRange2[k]+";Pt [MeV/c]; v2");
    lg2[k] = new TLegend(0.33, 0.14, 0.96, 0.4,"");
  }

  gr_v1f = new TGraphErrors();
  gr_v1f->SetName("gr_v1f");
  
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
       	RemovePoints(fname[is][ip], *gv1);

	if( gv1 == NULL ) {
	  cout << " not found " << k << endl;
	  continue;
	}

	TString gvname = Form("gPt_v1%d",igr);
	gr_v1 = (TGraphErrors*)gv1->Clone(gvname);
	gr_v1->SetName(gvname); 

	mv1[k]->Add(gr_v1,"lp");

	gr_v1->SetMarkerStyle(imrk[is]);
	gr_v1->SetMarkerColor(icol[ip][is]);
	gr_v1->SetMarkerSize(imsz[is]);
	gr_v1->SetLineColor(icol[ip][is]);
	gr_v1->GetXaxis()->SetRangeUser(0.,800.);
	lg1[k]->AddEntry(gr_v1, rsys[is]+" "+fpid[ip] ,"lp");

	igr++;
      }

      for(UInt_t k = 0; k < rbin2; k++){

	TGraphErrors *gv2 = (TGraphErrors*)fOpen->Get((TString)Form("gPt_v2%d",k));
	RemovePoints(fname[is][ip], *gv2);
	TString gvname = Form("gPt_v2%d",igr);
	gr_v2 = (TGraphErrors*)gv2->Clone(gvname);

	mv2[k]->Add(gr_v2,"lp");

 	gr_v2->SetMarkerStyle(imrk[is]);
	gr_v2->SetMarkerColor(icol[ip][is]);
	gr_v2->SetMarkerSize(imsz[is]);
	gr_v2->SetLineColor(icol[ip][is]);
	lg2[k]->AddEntry(gr_v2, rsys[is]+" "+fpid[ip],"lp");

	igr++;
      }
      fOpen->Close();
    }
  }


  cc1 = new TCanvas("cc1","v2",1100,500);
  cc1->Divide(rbin2,1);
  for(UInt_t k = 0; k < rbin2; k++){
    cc1->cd(k+1);
    mv2[k]->SetMaximum( 0.2);
    mv2[k]->SetMinimum(-0.2);
    mv2[k]->Draw("ALP");
    if( k == rbin2-1 ) lg2[k]->Draw();
  }

  cc0 = new TCanvas("cc0","v1",1200,700);
  cc0->Divide(rbin/2,2);
  for(UInt_t k = 0; k < rbin; k++){
    cc0->cd(k+1);
    mv1[k]->SetMaximum(0.4);
    mv1[k]->SetMinimum(-0.35);
    mv1[k]->Draw("ALP");
    if( k == rbin-1 )
      lg1[k]->Draw();
  }

  for(UInt_t i = 0; i < 8; i++) {
    cv[i] = new TCanvas(Form("cv%d",i),Form("cv%d",i),400,450);
    cv[i]->SetRightMargin(0.02);
    cv[i]->SetLeftMargin(0.12);
    cv[i]->SetTopMargin(0.05);
    mv1[i]->GetXaxis()->SetLabelSize(0.04);
    mv1[i]->GetYaxis()->SetLabelSize(0.04);
    mv1[i]->GetXaxis()->SetTitleSize(0.05);
    mv1[i]->GetYaxis()->SetTitleSize(0.05);
    mv1[i]->GetYaxis()->SetTitleOffset(1.2);

    mv1[i]->GetYaxis()->SetRangeUser(-0.32,0.38);

    mv1[i]->Draw("ALP");
    if( i == rbin-2 )
      lg1[i]->Draw();

  }
}





void SaveCanvas(TString cnt = "")
{

  Int_t iCanvas = gROOT->GetListOfCanvases()->GetEntries();  
  for(Int_t i = 0; i < iCanvas; i++)
    gROOT->GetListOfCanvases()->At(i)->SaveAs("ypt"+cnt+Form("_%d",i)+".png");

}

void RemovePoints(TString fn, TGraphErrors &gr)
{

  if( fn == "data/NL132Sn_neutron.root" ) {

    cout << gr.GetName() << endl;
    gr.Print();

    if( (TString)gr.GetName() == "gPt_v11") 
      gr.RemovePoint(1);
    else if( (TString)gr.GetName() == "gPt_v12") 
      gr.RemovePoint(2);
    else if( (TString)gr.GetName() == "gPt_v13") 
      gr.RemovePoint(3);
    else if( (TString)gr.GetName() == "gPt_v15") { 
      gr.RemovePoint(4);

    }    
    else if( (TString)gr.GetName() == "gPt_v16" ) {
      for(UInt_t i = 0; i < 2; i++)
	gr.RemovePoint(4);
    }
    else if( (TString)gr.GetName() =="gPt_v17" ) {
      for(UInt_t i = 0; i < 6; i++)
	gr.RemovePoint(1);
      //      gr.RemovePoint(0);
    }


    if( (TString)gr.GetName() == "gPt_v21" ) 
      gr.RemovePoint(2);
    else if( (TString)gr.GetName() == "gPt_v23") {
      for(UInt_t i = 0; i < 3; i++)
	gr.RemovePoint(3);
      gr.RemovePoint(0);
    }
    else if( (TString)gr.GetName() == "gPt_v24")
      gr.RemovePoint(3);

  }

  if( fn == "data/NL108Sn_neutron.root" ) {

    cout << gr.GetName() << endl;
    gr.Print();

    if( (TString)gr.GetName() == "gPt_v10" ) 
      gr.RemovePoint(1);
    if( (TString)gr.GetName() == "gPt_v11" ) 
      gr.RemovePoint(2);
    if( (TString)gr.GetName() == "gPt_v13" ) 
      gr.RemovePoint(3);
    if( (TString)gr.GetName() == "gPt_v14" ) 
      gr.RemovePoint(3);
    if( (TString)gr.GetName() == "gPt_v15" ) { 
      for(UInt_t i = 0; i < 2; i++)
	gr.RemovePoint(3);
    }
    else if( (TString)gr.GetName() == "gPt_v16" ) { 
      for(UInt_t i = 0; i < 3; i++)
	gr.RemovePoint(3);
    }
    else if( (TString)gr.GetName() == "gPt_v17" ) { 
      for(UInt_t i = 0; i < 5; i++)
	gr.RemovePoint(1);
    }
    else if( (TString)gr.GetName() == "gPt_v22" ) 
      gr.RemovePoint(2);
    else if( (TString)gr.GetName() == "gPt_v23" ) {
      for(UInt_t i = 0; i < 2; i++)
	gr.RemovePoint(3);
      gr.RemovePoint(0);
    }
    else if( (TString)gr.GetName() == "gPt_v24" ) {
      for(UInt_t i = 0; i < 2; i++)
	gr.RemovePoint(4);
      gr.RemovePoint(0);
    }
  }
}
