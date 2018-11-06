TCanvas *cc0;
TCanvas *cc1;
TCanvas *cc2;
TCanvas *cv[20];

// --> configuration
TString fsys[] = {"132Sn+124Sn","108Sn+112Sn","124Sn+112Sn","112Sn+124Sn"};
TString rsys[] = {"132",        "108",        "124",        "112"};
TString fpid[] = {"proton","deuteron","triton","neutron"};  
UInt_t  imrk[] = {20, 21, 22, 23, 21};
Size_t  imsz[] = {1, 1, 1.3, 1.3, 1.3};
Color_t icol[][4] = { {kRed,          kBlue,  kOrange-3,   kGreen+3}, 
		      {kBlue+2,   kOrange+7,  kGreen-3,     kPink+9},
		      {kGreen-3,    kPink+7,  kCyan-1,    kYellow-2},
		      {kGreen+3,   kOrange-1, kCyan-1,    kBlue} };

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
  Bool_t bsys[]  = { 0, 1, 0, 0};
  Bool_t bpid[]  = { 1, 1, 1, 1}; //p, d, t, n


  TString fname[4][4];
  for(UInt_t is = 0; is < 4; is++){
    if( !bsys[is] ) continue;

    for(UInt_t ip = 0; ip < 3; ip++){
      if(bpid[ip])
	fname[is][ip] = "data/YPT"+rsys[is]+"Sn_"+fpid[ip]+".root";
    }

    if(bpid[3])
      fname[is][3] = "data/vaNL"+rsys[is]+"Sn_neutron.root";
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

  // Rapidity dependence 
  auto mrv1 = new TMultiGraph("mrv1",";Rapidity; v1");
  auto mrv2 = new TMultiGraph("mrv2",";Rapidity; v2");
  auto lgr1 = new TLegend(0.54, 0.13, 0.87, 0.4,""); 
  auto lgr2 = new TLegend(0.54, 0.14, 0.9, 0.34,"");


  // Pt dependence
  for(UInt_t k = 0; k < rbin; k++){
    mv1[k] = new TMultiGraph((TString)Form("mv1%d",k), rapRange[k]+";Pt [MeV/c]; v1");
    lg1[k] = new TLegend(0.3, 0.13, 0.88, 0.4,"");
    lg1[k]->SetTextSize(0.07);
  }

  for(UInt_t k = 0; k < rbin2; k++) {
    mv2[k] = new TMultiGraph((TString)Form("mv2%d",k), rapRange2[k]+";Pt [MeV/c]; v2");
    lg2[k] = new TLegend(0.55, 0.70, 0.9, 0.9,"");
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

      //fOpen->ls();
    
      // rapidity dependence
      TGraphErrors *yv1 = (TGraphErrors*)fOpen->Get("gv_v1");
      yv1->SetMarkerColor(icol[ip][is]);
      yv1->SetMarkerStyle(imrk[is]);
      yv1->SetMarkerSize(imsz[is]);
      yv1->SetLineColor(icol[ip][is]);
      mrv1->Add(yv1,"lp");
      lgr1->AddEntry(yv1, rsys[is]+"Sn "+fpid[ip] ,"lp");

      TGraphErrors *yv2 = (TGraphErrors*)fOpen->Get("gv_v2");
      yv2->SetMarkerColor(icol[ip][is]);
      yv2->SetMarkerStyle(imrk[is]);
      yv2->SetMarkerSize(imsz[is]);
      yv2->SetLineColor(icol[ip][is]);
      mrv2->Add(yv2,"lp");
      lgr2->AddEntry(yv2, rsys[is]+"Sn "+fpid[ip] ,"lp");
      // --end of rapidity dependence


      // Pt dependence 
      for(UInt_t k = 0; k < rbin; k++){
	TGraphErrors *gv1 = (TGraphErrors*)fOpen->Get((TString)Form("gPt_v1%d",k));
	//       	RemovePoints(fname[is][ip], *gv1);

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
	lg1[k]->AddEntry(gr_v1, rsys[is]+"Sn "+fpid[ip] ,"lp");

	igr++;
      }

      for(UInt_t k = 0; k < rbin2; k++){

	TGraphErrors *gv2 = (TGraphErrors*)fOpen->Get((TString)Form("gPt_v2%d",k));
	//RemovePoints(fname[is][ip], *gv2);
	TString gvname = Form("gPt_v2%d",igr);
	gr_v2 = (TGraphErrors*)gv2->Clone(gvname);

	mv2[k]->Add(gr_v2,"lp");

 	gr_v2->SetMarkerStyle(imrk[is]);
	gr_v2->SetMarkerColor(icol[ip][is]);
	gr_v2->SetMarkerSize(imsz[is]);
	gr_v2->SetLineColor(icol[ip][is]);
	lg2[k]->AddEntry(gr_v2, rsys[is]+"Sn "+fpid[ip],"lp");

	igr++;
      }

      fOpen->Close();
    }
  }


  cc1 = new TCanvas("cc1","v2",1600,500);
  cc1->Divide(rbin2,1);
  for(UInt_t k = 0; k < rbin2; k++){
    cc1->cd(k+1);
    //    mv2[k]->SetMaximum( 0.2);
    //    mv2[k]->SetMinimum(-0.2);
    mv2[k]->Draw("ALP");
    if( k == rbin2-2 ) lg2[k]->Draw();
  }

  cc0 = new TCanvas("cc0","v1",1600,700);
  cc0->Divide(rbin/2,2);
  for(UInt_t k = 0; k < rbin; k++){
    cc0->cd(k+1);
    //    mv1[k]->SetMaximum(0.45);
    //    mv1[k]->SetMinimum(-0.45);
    mv1[k]->Draw("ALP");
    if( k == rbin-2 )
      lg1[k]->Draw();
  }

  if( kFALSE ) {  
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

      if( i < 4 )
	mv1[i]->GetYaxis()->SetRangeUser(-0.4,0.2);
      else
	mv1[i]->GetYaxis()->SetRangeUser(-0.2,0.48);

      mv1[i]->Draw("ALP");
      if( i == rbin-2 )
	lg1[i]->Draw();

    }

    for(UInt_t i = 0; i < 5; i++) {
      cv[i] = new TCanvas(Form("cv%d",i+10),Form("cv%d",i+10),400,450);
      cv[i]->SetRightMargin(0.02);
      cv[i]->SetLeftMargin(0.12);
      cv[i]->SetTopMargin(0.05);
      mv2[i]->GetXaxis()->SetLabelSize(0.04);
      mv2[i]->GetYaxis()->SetLabelSize(0.04);
      mv2[i]->GetXaxis()->SetTitleSize(0.05);
      mv2[i]->GetYaxis()->SetTitleSize(0.05);
      mv2[i]->GetYaxis()->SetTitleOffset(1.2);

      mv2[i]->GetYaxis()->SetRangeUser(-0.2,0.14);

      mv2[i]->Draw("ALP");
      if( i == rbin2-2 )
	lg2[i]->Draw();

    }
  }

  auto cc2 = new TCanvas("cc2","cc2");
  mrv1->Draw("ALP");
  lgr1->Draw();

  auto Ymin = mrv1->GetYaxis()->GetXmin();
  auto Ymax = mrv1->GetYaxis()->GetXmax();
  auto Xmin = mrv1->GetXaxis()->GetXmin();
  auto Xmax = mrv1->GetXaxis()->GetXmax();

  auto aLineX1 = new TLine(Xmin, 0., Xmax, 0.);
  aLineX1->SetLineColor(1);
  aLineX1->SetLineStyle(3);
  aLineX1->Draw();

  auto aLineY1 = new TLine(0., Ymin, 0., Ymax);
  aLineY1->SetLineColor(1);
  aLineY1->SetLineStyle(3);
  aLineY1->Draw();

   

  auto cc3 = new TCanvas("cc3","cc3");
  mrv2->Draw("ALP");
  lgr2->Draw();

  Ymin = mrv2->GetYaxis()->GetXmin();
  Ymax = mrv2->GetYaxis()->GetXmax();
  Xmin = mrv2->GetXaxis()->GetXmin();
  Xmax = mrv2->GetXaxis()->GetXmax();

  auto aLineX2 = new TLine(0., Ymin, 0., Ymax);
  aLineX2->SetLineColor(1);
  aLineX2->SetLineStyle(3);
  aLineX2->Draw();

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
