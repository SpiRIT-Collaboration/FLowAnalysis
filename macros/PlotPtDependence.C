TCanvas *cv[20];
TCanvas *ccv;

TF1 *lslope = new TF1("lslope","[0]+[1]*x",-0.1,0.2);

// --> configuration
TString fsys[] = {"132Sn+124Sn","108Sn+112Sn","124Sn+112Sn","112Sn+124Sn"};
TString rsys[] = {"132",        "108",        "124",        "112"};
TString fpid[] = {"proton","deuteron","triton","neutron","3He","4He"};  
UInt_t  imrk[] = {20, 21, 22, 23, 24, 25, 26};
Size_t  imsz[] = {1, 1, 1.3, 1.3, 1.3, 1.3, 1.3};
Color_t icol[] = { kRed, kBlue, kSpring, kMagenta, kOrange, kViolet};
Color_t icolnn[]={ kPink, kBlue+1, kGreen+2, kViolet-1};

//  Double_t rrange[3] = {0.27, 0.54, 1.};
Double_t yrange[] =  {-0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35, 0.45};
const UInt_t rbin = sizeof(yrange)/sizeof(Double_t);

Double_t yrange2[] = {-0.25, -0.1,  0.1, 0.25, 0.45};
const UInt_t rbin2 = sizeof(yrange2)/sizeof(Double_t);

//UInt_t  cent[] = {70, 40, 35, 20, 0};
UInt_t   cent[]  = {70, 35, 28, 0};
UInt_t  isobin[] = {1000000, 25, 10, 1, 0};

TMultiGraph  *mv1[rbin];
TMultiGraph  *mv2[rbin2];
TLegend      *lg1[rbin];
TLegend      *lg2[rbin2];

void RemovePoints(TString fn, TGraphErrors &gr);
void RemoveYPoint(TGraphErrors *gv);
void RemovePtPoint(UInt_t ip, TGraphErrors *gv);
void RapidityShift(TGraphErrors *gv);

Double_t FittingAndIntegral(TGraphErrors *gr)
{
  TF1 *p1 = new TF1("p1","[0]+[1]*x",0.,800.);
  
  gr->Fit("p1","","",200.,700.);
  
  return p1->Integral(0., 800) / 800.;
}




void PlotPtDependence()
{
  gStyle->SetOptStat(0);

  // --> Plotting selection
  Bool_t bsys[]  = { 1, 1, 0, 0};
  Bool_t bpid[]  = { 1, 0, 1, 0, 1, 0}; //p, d, t, n
  Bool_t bcnt[]  = { 1, 0, 0}; 
  UInt_t cntw = 3;
  UInt_t bazm[]  = { 1, 0, 0}; // 1: pm, 10:m, 11:xsp, //1: <45, 2:>135, 3:45<&<135

  UInt_t ngr = 0;
  std::vector< TString > fname;
  std::vector< std::vector< UInt_t >> index(5);
  TGraphErrors *g_slope[4][4];
  TGraphErrors *g_v2cm[4][4];

  for(UInt_t is = 0; is < 4; is++){

    for(UInt_t ip = 0; ip < (UInt_t)sizeof(bpid)/sizeof(Bool_t); ip++){

      g_slope[is][ip] = new TGraphErrors();
      g_slope[is][ip] ->SetName( Form("g_slope_%d%d",is,ip) );
      g_v2cm[is][ip]  = new TGraphErrors();
      g_v2cm[is][ip]  ->SetName( Form("g_v2cm_%d%d",is,ip) );


      if( !bsys[is] || !bpid[ip] ) continue;

      for(UInt_t it = 0; it < (UInt_t)sizeof(bcnt)/sizeof(Bool_t); it++){
	
	if( !bcnt[it] ) continue;

	for(UInt_t iz = 0; iz < (UInt_t)sizeof(bazm)/sizeof(UInt_t); iz++){
	  
	  if( bazm[iz] == 0 ) continue;

	
	  //	  fname.push_back( Form("data/cosYPT_n%dto%d_",cent[it],cent[it+cntw]) + rsys[is]+"Snyoff_"+fpid[ip]+".root" );
	  //	fname.push_back( Form("data/cosYPT_n%dto%d_",cent[it],cent[it+cntw]) + rsys[is]+"Snyoffp_"+fpid[ip]+".root" );

	  fname.push_back( Form("data/cosYPT_n%dto%d_",cent[it],cent[it+cntw]) + rsys[is]+"Snyoffp_"+Form("az%d",bazm[iz]) + fpid[ip]+".root" );

	  // fname.push_back( Form("data/cosYPT_n%dto%d_",cent[it],cent[it+cntw]) + rsys[is]+"Sn_"+fpid[ip]+".root" );
	  //	  fname.push_back( Form("data/cosYPT_iso%dto%d_",isobin[it],isobin[ik+1]) + rsys[is]+"Sn_"+fpid[ip]+".root" );

	  std::cout << fname.at(ngr) << std::endl;

	  index[0].push_back(is); // system
	  index[1].push_back(ip); // pid
	  index[2].push_back(it); // centrality
	  index[3].push_back(it+cntw);
	  index[4].push_back(iz);
	  ngr++;
	}
      }
    }
  }

  cout << fname.size() << endl;
  

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


  // Rapidity dependence 
  auto mrv1  = new TMultiGraph("mrv1",";Rapidity; v1");
  auto mrv2  = new TMultiGraph("mrv2",";Rapidity; v2");
  auto mv1sl = new TMultiGraph("mv1sl",";Centrality; v1 slope");
  auto mv2cm = new TMultiGraph("mv2cm",";Centrality; v2 at Y_cm");
  auto lgr1 = new TLegend(0.54, 0.13, 0.87, 0.4,""); 
  auto lgr2 = new TLegend(0.54, 0.14, 0.9, 0.34,"");
  auto lgr3 = new TLegend(0.35, 0.13, 0.7, 0.33,"");
  auto lgr4 = new TLegend(0.15, 0.63, 0.5, 0.85,"");
  auto lgr5 = new TLegend(0.15, 0.63, 0.5, 0.85,"");


  // Pt dependence
  for(UInt_t k = 0; k < rbin; k++){
    mv1[k] = new TMultiGraph((TString)Form("mv1%d",k), rapRange[k]+";Pt [MeV/c]; v1");
    lg1[k] = new TLegend(0.5, 0.15, 0.94, 0.43,""); 
    lg1[k]->SetTextSize(0.05);
  }

  for(UInt_t k = 0; k < rbin2; k++) {
    mv2[k] = new TMultiGraph((TString)Form("mv2%d",k), rapRange2[k]+";Pt [MeV/c]; v2");
    lg2[k] = new TLegend(0.17, 0.2, 0.62, 0.4,"");
    lg2[k]->SetTextSize(0.05);
  }



  //  ccv = new TCanvas("cc10","cc10",1000,1000);
  //  ccv->Divide(2, ngr/2);
  UInt_t id = 1;

  Color_t icolor = 100; 
  for(UInt_t igr = 0; igr < ngr; igr++ ) {

    fOpen = TFile::Open(fname.at(igr)); 

    if( fOpen == NULL ) continue;
    else
      std::cout << fname.at(igr) << " is opened. " << std::endl;


    TH1D *fevt = (TH1D*)fOpen->Get("hphi0_180");
    UInt_t tevt = fevt->GetEntries();

    UInt_t is = index[0].at(igr);
    UInt_t ip = index[1].at(igr);
    UInt_t it = index[2].at(igr);
    UInt_t ik = index[3].at(igr);
    UInt_t iz = index[4].at(igr);

    if( igr < 7)
      icolor = icol[igr];
    else
      icolor -= 2;
    //    UInt_t icolor = (icol[is] + 4*ip - 2*it );
    //    if(icolor  >= 912) icolor = 4;

    cout << igr  << " " << ip << " color " << icolor << endl;


    //    icolor = icolnn[igr];

    //acceptance
    TH2D *ff = (TH2D*)fOpen->Get("hyptacp");
    ff->SetName(Form("hyptacp_%d",igr));
    //    TString otitle = rsys[is]+"Sn "+fpid[ip]+": "+ (TString)ff->GetTitle();
    TString otitle = rsys[is]+"Sn "+fpid[ip](0,3) + Form(":n%dto%d_azm%d",cent[it],cent[ik],bazm[iz]);
    ff->SetTitle(otitle);
    ff->SetDirectory(gROOT);

    gROOT->cd();
    ff->ProjectionX("_px");
    TH1D* ff1 = (TH1D*)gROOT->Get((TString)ff->GetName()+"_px");
    
    Double_t xwd = (Double_t)ff1->GetXaxis()->GetNbins()/(ff1->GetXaxis()->GetXmax() - ff1->GetXaxis()->GetXmin());
    
    // cout << " nbin " << ff1->GetXaxis()->GetNbins()
    // 	 << " max " << ff1->GetXaxis()->GetXmax()
    // 	 << " min " << ff1->GetXaxis()->GetXmin()
    // 	 << " xwd " << xwd
    // 	 << endl;

    if( fpid[ip] == "neutron" )
      ff1->Scale( xwd / tevt /0.26);
    else
      ff1->Scale(xwd / tevt);

    ff1->SetTitle("; Y_{cm}; dN/dy");
    ff1->SetLineColor(icolor);
    lgr3->AddEntry(ff1, otitle );
    

    fOpen->cd();
    
    // rapidity dependence
    TGraphErrors *yv1 = (TGraphErrors*)fOpen->Get("gv_v1");
    
    if((TString)fOpen->GetName() == "data/cosYPT_n35to0_108Sn_proton.root")
      RapidityShift(yv1);

    RemoveYPoint(yv1);

    yv1->SetMarkerColor(icolor);

    cout << "imark[" << is << "] " << imrk[is] << endl;
    yv1->SetMarkerStyle(imrk[is]);
    if( ip == 3 )
      yv1->SetMarkerStyle(imrk[is]+4);
    yv1->SetMarkerSize(imsz[is]);
    yv1->SetLineColor(icolor);
	  
    
    mrv1->Add(yv1,"lp");
    lgr1->AddEntry(yv1,  otitle ,"lp");


    // slope
    //    ccv->cd(id); id++;
    //    yv1->Draw("ALP");
    yv1->Fit("lslope","Q0","",-0.11,0.2);
    Double_t constlslope = lslope->GetParameter(0);
    Double_t slope = lslope->GetParameter(1);
    Double_t slopee= lslope->GetParError(1);

    g_slope[is][ip]->SetPoint(it, (Double_t)it+1+0.05*is, slope);
    g_slope[is][ip]->SetPointError(it, 0., slopee);

    //	yv1->SetDirectory(gROOT);

    
    TGraphErrors *yv2 = (TGraphErrors*)fOpen->Get("gv_v2");

    RemoveYPoint(yv2);

    yv2->SetMarkerColor(icolor);
    yv2->SetMarkerStyle(imrk[is]);
    if( ip == 3 )
      yv2->SetMarkerStyle(imrk[is]+4);
    yv2->SetMarkerSize(imsz[is]);
    yv2->SetLineColor(icolor);

    mrv2->Add(yv2,"lp");
    lgr2->AddEntry(yv2,  otitle ,"lp");
    // --end of rapidity dependence

    Double_t v2x, v2y, v2xe;
    yv2->GetPoint(2, v2x, v2y);
    v2xe = yv2->GetErrorY(2);
    g_v2cm[is][ip]->SetPoint(it, (Double_t)it+1+0.05*is, v2y);
    g_v2cm[is][ip]->SetPointError(it, 0., v2xe);


    // Pt dependence 
    for(UInt_t k = 0; k < rbin; k++){

      TGraphErrors *gr_v1 = (TGraphErrors*)fOpen->Get((TString)Form("gPt_v1%d",k));
      

      RemovePtPoint(ip, gr_v1);
      
      if( gr_v1 == NULL ) continue;
      gr_v1->SetMarkerStyle(imrk[is]);


      if( ip == 3 )
	gr_v1->SetMarkerStyle(imrk[is]+4);
      gr_v1->SetMarkerColor(icolor);
      gr_v1->SetMarkerSize(imsz[is]);
      gr_v1->SetLineColor(icolor);
      gr_v1->GetXaxis()->SetRangeUser(0.,800.);


      mv1[k]->Add(gr_v1,"lp");
      lg1[k]->AddEntry(gr_v1,  otitle,"lp");
    }
    
    cout << " ==== " << endl;

    for(UInt_t k = 0; k < rbin2; k++){

      TGraphErrors *gr_v2 = (TGraphErrors*)fOpen->Get((TString)Form("gPt_v2%d",k));
      if( gr_v2 == NULL ) continue;

      RemovePtPoint(ip,gr_v2);

      
      gr_v2->SetMarkerStyle(imrk[is]);
      if( ip == 3 )
	gr_v2->SetMarkerStyle(imrk[is]+4);
      gr_v2->SetMarkerColor(icolor);
      gr_v2->SetMarkerSize(imsz[is]);
      gr_v2->SetLineColor(icolor);

      mv2[k]->Add(gr_v2,"lp");
      lg2[k]->AddEntry(gr_v2, otitle ,"lp");
    }
	
    fOpen->Close();
  }


  auto cc4 = new TCanvas("cc4","aceptance", 700, 1000);
  cc4->Divide(1,ngr);
  for(UInt_t i = 0; i < ngr; i++){
    cc4->cd(i+1);
    TH2D* ff = (TH2D*)gROOT->Get(Form("hyptacp_%d",i));
    ff->Draw("colz");
  }

  auto cc5 = new TCanvas("cc5","cc5");
  Double_t ymax = 0;
  UInt_t   imax = 0;
  for(Int_t i = 0; i < ngr; i++){
    TH1D* ff = (TH1D*)gROOT->Get(Form("hyptacp_%d_px",i));
    if( ymax < ff->GetMaximum() ) {
      ymax = ff->GetMaximum();
      imax = i;
      ff->Draw();
    }
  }
  for(Int_t i = 0; i < ngr; i++) {
    TH1D* ff = (TH1D*)gROOT->Get(Form("hyptacp_%d_px",i));
    if( i != imax ) ff->Draw("same");
  }
  lgr3->Draw();


  auto cc1 = new TCanvas("cc1","v2",1600,500);
  cc1->Divide(rbin2,1);
  for(UInt_t k = 0; k < rbin2; k++){
    cc1->cd(k+1);
    //    mv2[k]->SetMaximum( 0.2);
    //    mv2[k]->SetMinimum(-0.2);
    mv2[k]->Draw("ALP");
    if( k == rbin2-2 ) lg2[k]->Draw();
  }



  auto cc0 = new TCanvas("cc0","v1",1600,700);
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

  for(UInt_t is = 0; is < 4; is++){
    for(UInt_t ip = 0; ip < 4; ip++){
      UInt_t icolor = (icol[is] + 4*ip);
      if( g_slope[is][ip]->GetN() > 0 ) {
	g_slope[is][ip]->SetMarkerStyle(imrk[is]);
	g_slope[is][ip]->SetMarkerColor(icol[is]+2*ip);
	g_slope[is][ip]->SetLineColor(icol[is]+2*ip);
	mv1sl->Add(g_slope[is][ip],"lp");      

	lgr4->AddEntry(g_slope[is][ip],fsys[is]+":"+fpid[ip](0,4),"lp");
      }

      if( g_v2cm[is][ip]->GetN() > 0) {
	g_v2cm[is][ip]->SetMarkerStyle(imrk[is]);
	g_v2cm[is][ip]->SetMarkerColor(icol[is]+2*ip);
	g_v2cm[is][ip]->SetLineColor(icol[is]+2*ip);
	mv2cm->Add(g_v2cm[is][ip],"lp");      

	lgr5->AddEntry(g_v2cm[is][ip], fsys[is]+":"+fpid[ip](0,4), "lp");
      }
    }

  }
  // auto cc6 = new TCanvas("cc6","cc6");
  // mv1sl->Draw("ALP");
  // lgr4->Draw();
  // auto cc7 = new TCanvas("cc7","cc7");
  // mv2cm->Draw("ALP");
  // lgr5->Draw();


  if( kFALSE ){  
    //for NN v1-pt
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadLeftMargin(0.1);
    gStyle->SetTitleYSize(0.04);
    gStyle->SetTitleYOffset(1.1);


    auto cc9 = new TCanvas("cc9","v1_target",500,500);
    mv1[1]->SetMaximum(0.);
    mv1[1]->SetMinimum(-0.27);
    mv1[1]->Draw("ALP");

    auto cc12 = new TCanvas("cc12","v1_mid",500,500);
    mv1[3]->Draw("ALP");

    auto cc11 = new TCanvas("cc11","v1_beam",500,500);
    mv1[5]->Draw("ALP");
    lg1[5]->SetX1NDC(0.5);
    lg1[5]->SetY1NDC(0.16);
    lg1[5]->SetX2NDC(0.95);
    lg1[5]->SetY2NDC(0.40);
    lg1[5]->Draw();

    // v2-pt
    auto cc8 = new TCanvas("cc8","v2 at Y_cm=0",500,500);
    mv2[2]->Draw("ALP");
    lg2[2]->Draw();
  }
  
}


void SaveCanvas(TString cnt = "", Int_t isel=-1)
{
  if(isel > -1)
    gROOT->GetListOfCanvases()->At(isel)->SaveAs("ypt_"+cnt+Form("_%d",isel)+".png");

  else {
    TString mdir = "ypt_"+cnt;
    gSystem->mkdir(mdir);
    gSystem->cd(mdir);
      
    Int_t iCanvas = gROOT->GetListOfCanvases()->GetEntries();  
    for(Int_t i = 0; i < iCanvas; i++)
      gROOT->GetListOfCanvases()->At(i)->SaveAs(mdir+Form("_%d",i)+".png");

    gSystem->cd("..");
    std::cout << "Figures in " + mdir << std::endl;
  }
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


void RapidityShift(TGraphErrors *gv) 
{
  // shifting artificially
  std::cout << "===> SHIFTING " << 0.018 << std::endl;
  Double_t xx, yy;
  for(UInt_t k = 0; k < (UInt_t)gv->GetN(); k++) {
    gv->GetPoint(k, xx, yy);
    xx += 0.018;
    gv->SetPoint(k, xx, yy);
  }
}

void RemoveYPoint(TGraphErrors *gv)
{
  cout << " removey " << gv->GetN() <<endl;
  Double_t xx, yy;
  for(Int_t k = (Int_t)gv->GetN()-1; k > -1; k--) {

    gv->GetPoint(k, xx, yy);
    if( xx < -0.2 ) {
      gv->RemovePoint(k);
      //      std::cout << "===> REMOVEPOINT in Y (" << xx << "," << yy << ")" << std::endl;    
    }
    if( xx > 0.37 ) {
      gv->RemovePoint(k);
      //      std::cout << "===> REMOVEPOINT in Y (" << xx << "," << yy << ")" << std::endl;    
    }
  }

  for(Int_t k = (Int_t)gv->GetN()-1; k > -1; k--) {
    yy = gv->GetErrorY(k);
    if( yy > 0.15 )
      gv->RemovePoint(k);

  }  
}

void RemovePtPoint(UInt_t ip, TGraphErrors *gv)
{
  if( gv == NULL ) return;

  //  cout << "removept " << gv->GetN() <<  endl;
  Double_t ptcut = 700.;
  if( ip == 3 ) ptcut = 400.;

  Double_t xx, yy;
  for(Int_t k = (Int_t)gv->GetN()-1; k > -1; k--) {
    gv->GetPoint(k, xx, yy);
    if( xx > ptcut ) {
      gv->RemovePoint(k);
      //      std::cout << "===> REMOVEPOINT in Pt (" << xx << "," << yy << ")" <<  std::endl;
    }
  }  
  for(Int_t k = (Int_t)gv->GetN()-1; k > -1; k--) {
    yy = gv->GetErrorY(k);
    if( yy > 0.12 )
      gv->RemovePoint(k);

  }  
}



void integ()
{
  UInt_t idx = 0;
  while(1) {
    TH1D* ff = (TH1D*)gROOT->Get(Form("hyptacp_%d_px",idx));
    if( ff != NULL ) {
      ff->Print();
      cout << " -> " <<  ff->Integral()/200. << endl;;
      //      cout << " -> " <<  ff->Integral()/(Double_t)ff->GetXaxis()->GetNbins() << endl;;
      idx++;
    }
    else break;

    //    hyptacp_0_px->Integral()/200+hyptacp_1_px->Integral()/200+hyptacp_2_px->Integral()/200;
    //    hyptacp_3_px->Integral()/200+hyptacp_4_px->Integral()/200+hyptacp_5_px->Integral()/200;
  }
}
  
