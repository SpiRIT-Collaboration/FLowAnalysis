#include "DoFlowAnalysis.h"

  // --> Plotting selection
Bool_t bsys[]  = { 1, 0, 0, 0};
Bool_t bpid[]  = { 1, 0, 0, 0, 0, 0}; //0:p, 1:d, 2:t, 3:3He, 4:4He, 5:n
Bool_t bcnt[]  = { 1, 0, 0}; 
UInt_t cntw = 3;

UInt_t bver[]  = { 1, 1, 1, 1, 1, 1};
TString sVer[] = {".v23.1.7",".v23.1.10", ".v23.1.11", ".v23.1.12", ".v23.1.13"};
TString sName = "cosYPt_132Sn_";



UInt_t   ccvid = 0;
TF1 *lslope = new TF1("lslope","[0]+[1]*x",-0.1,0.2);

// --> configuration
TString fsys[] = {"132Sn+124Sn","108Sn+112Sn","124Sn+112Sn","112Sn+124Sn"};
TString rsys[] = {"132",        "108",        "124",        "112"};
TString fpid[] = {"proton","deuteron","triton","3He","4He","neutron"};  

Size_t  imsz[] = {1, 1, 1.3, 1.3, 1.3, 1.3, 1.3};
Color_t icol[] = { kRed, kBlue, kSpring, kMagenta, kOrange, kViolet};
Color_t icolnn[]={ kPink, kBlue+1, kGreen+2, kViolet-1};

TMultiGraph  *mv1[ybin1];
TMultiGraph  *mv2[ybin2];
TLegend      *lg1[ybin1];
TLegend      *lg2[ybin2];

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


  Bool_t bplot = 0;
  TString ltitle;


  UInt_t ngr = 0;
  std::vector< TString > fname;
  std::vector< std::vector< UInt_t >> index(5);


  TGraphErrors *g_slope[4][6];
  TGraphErrors *g_v2cm[4][6];
  TGraphErrors *g_v1slp[4];
  TGraphErrors *g_v2mid[4];
  

  for(UInt_t is = 0; is < 4; is++){
    g_v1slp[is] = new TGraphErrors();
    g_v1slp[is]->SetName( Form("g_v1slp%d",is) );

    g_v2mid[is] = new TGraphErrors();
    g_v2mid[is]->SetName( Form("g_v2mid%d",is) );

    for(UInt_t ip = 0; ip < (UInt_t)sizeof(bpid)/sizeof(Bool_t); ip++){

      g_slope[is][ip] = new TGraphErrors();
      g_slope[is][ip] ->SetName( Form("g_slope_%d%d",is,ip) );
      g_v2cm[is][ip]  = new TGraphErrors();
      g_v2cm[is][ip]  ->SetName( Form("g_v2cm_%d%d",is,ip) );
      
      
      if( !bsys[is] || !bpid[ip] ) continue;

      for(UInt_t it = 0; it < (UInt_t)sizeof(bcnt)/sizeof(Bool_t); it++){
	
	if( !bcnt[it] ) continue;

	for(UInt_t iz = 0; iz < (UInt_t)sizeof(bver)/sizeof(UInt_t); iz++){
	  
	  if( bver[iz] == 0 ) continue;

	  fname.push_back( "data/"+ sName + fpid[ip] + sVer[iz] + ".root" );

	  std::cout << fname.at(ngr) << std::endl;
	  ltitle = Form("mult %d ~ %d; azm%d",cent[it],cent[it+cntw],bver[iz]);
	  std::cout << " label title " << ltitle << std::endl;

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
  

  TString rapRange[ybin1];
  for(UInt_t i = 1; i < ybin1 - 1; i++)
    rapRange[i] = Form(" %f <= Y_cm < %f", yrange1[i-1],yrange1[i]);
  rapRange[0] = Form(" Y_cm < %f", yrange1[0]);
  rapRange[ybin1-1] = Form(" %f <= Y_cm", yrange1[ybin1-2]);

  TString rapRange2[ybin2];
  for(UInt_t i = 1; i < ybin2 - 1; i++)
    rapRange2[i] = Form(" %f <= Y_cm < %f", yrange2[i-1],yrange2[i]);
  rapRange2[0] = Form(" Y_cm < %f", yrange2[0]);
  rapRange2[ybin2-1] = Form(" %f <= Y_cm", yrange2[ybin2-2]);


  TFile *fOpen;

  // Rapidity dependence 
  auto mrv1  = new TMultiGraph("mrv1"  ,";Rapidity; v1");
  auto mrv2  = new TMultiGraph("mrv2"  ,";Rapidity; v2");
  auto mv1sl = new TMultiGraph("mv1sl" ,";Centrality; v1 slope");
  auto mv1slp= new TMultiGraph("mv1slp","; Particle ; v1 slope");
  //  mv1slp->GetXaxis()->SetAlphanumeric(kTRUE);
  auto mv2cm = new TMultiGraph("mv2cm" ,";Centrality; v2 at Y_cm");
  auto mv2mid= new TMultiGraph("mv2mid","; Particle ; v2 at Y_cm");


  auto lgr1 = new TLegend(0.54, 0.13, 0.87, 0.4, ltitle); 
  auto lgr2 = new TLegend(0.21, 0.68, 0.57, 0.9, ltitle);
  auto lgr3 = new TLegend(0.35, 0.13, 0.7, 0.33, ltitle);
  auto lgr4 = new TLegend(0.15, 0.63, 0.5, 0.85, ltitle);
  auto lgr5 = new TLegend(0.15, 0.63, 0.5, 0.85, ltitle);
  auto lgr6 = new TLegend(0.72, 0.72, 0.94, 0.88,ltitle);
  auto lgr7 = new TLegend(0.16, 0.70, 0.46, 0.85,ltitle);


  // Pt dependence
  for(UInt_t k = 0; k < ybin1; k++){
    mv1[k] = new TMultiGraph((TString)Form("mv1%d",k), rapRange[k]+";Pt [MeV/c]; v1");
    lg1[k] = new TLegend(0.57, 0.64, 0.95, 0.92); 
    lg1[k]->SetTextSize(0.05);
  }

  for(UInt_t k = 0; k < ybin2; k++) {
    mv2[k] = new TMultiGraph((TString)Form("mv2%d",k), rapRange2[k]+";Pt [MeV/c]; v2");
    lg2[k] = new TLegend(0.5, 0.64, 0.95, 0.92); 
    lg2[k]->SetTextSize(0.05);
  }



  //  UInt_t id = 1;

  UInt_t iss = 9;
  UInt_t ipp = 0;
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
    
    if( is != iss ) {
      ipp = 0;
      iss = is;
    }

    if( igr < 7)
      icolor = icol[igr];
    else
      icolor -= 2;

    cout << igr  << " " << ip << " color " << icolor << endl;


    //acceptance
    TH2D *ff = (TH2D*)fOpen->Get("hyptacp");
    ff->SetName(Form("hyptacp_%d",igr));

    TString otitle = rsys[is]+"Sn "+fpid[ip]+sVer[iz];
    ff->SetTitle(otitle);
    ff->SetDirectory(gROOT);

    gROOT->cd();
    ff->ProjectionX("_px");
    TH1D* ff1 = (TH1D*)gROOT->Get((TString)ff->GetName()+"_px");
    
    Double_t xwd = (Double_t)ff1->GetXaxis()->GetNbins()/(ff1->GetXaxis()->GetXmax() - ff1->GetXaxis()->GetXmin());
    
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
    
    //    RemoveYPoint(yv1);

    yv1->SetMarkerColor(icolor);

    cout << "imark[" << is << "] " << imark[is] << endl;
    yv1->SetMarkerStyle(imark[is]);
    if( ip == 5 )
      yv1->SetMarkerStyle(imark[is]+4);
    yv1->SetMarkerSize(imsz[is]);
    yv1->SetLineColor(icolor);
	  
    
    mrv1->Add(yv1,"lp");
    lgr1->AddEntry(yv1,  otitle ,"lp");


    //slope
    yv1->Fit("lslope","Q0","",-0.11,0.2);
    Double_t constlslope = lslope->GetParameter(0);
    Double_t slope = lslope->GetParameter(1);
    Double_t slopee= lslope->GetParError(1);

    g_slope[is][ip]->SetPoint(it, (Double_t)it+1+0.05*is, slope);
    g_slope[is][ip]->SetPointError(it, 0., slopee);

    g_v1slp[is]->SetPoint(ipp, ip+1, slope);
    g_v1slp[is]->SetPointError(ipp,  0, slopee);
    
    TGraphErrors *yv2 = (TGraphErrors*)fOpen->Get("gv_v2");
    //    RemoveYPoint(yv2);

    yv2->SetMarkerColor(icolor);
    yv2->SetMarkerStyle(imark[is]);
    if( ip == 5 )
      yv2->SetMarkerStyle(imark[is]+4);
    yv2->SetMarkerSize(imsz[is]);
    yv2->SetLineColor(icolor);

    mrv2->Add(yv2,"lp");
    lgr2->AddEntry(yv2,  otitle ,"lp");
    // --end of rapidity dependence

    Double_t v2x, v2y, v2xe;
    yv2->GetPoint(1, v2x, v2y);
    v2xe = yv2->GetErrorY(1);
    g_v2cm[is][ip]->SetPoint(it, (Double_t)it+1+0.05*is, v2y);
    g_v2cm[is][ip]->SetPointError(it, 0., v2xe);

    g_v2mid[is]->SetPoint(ipp, ip+1, v2y);
    g_v2mid[is]->SetPointError(ipp, 0., v2xe); ipp++;

    cout << fsys[is] << " : " << fpid[ip] << " -> " << v2y << " +- " << v2xe << " @ " << ipp <<  endl;

    // Pt dependence 
    for(UInt_t k = 0; k < ybin1; k++){

      TGraphErrors *gr_v1 = (TGraphErrors*)fOpen->Get((TString)Form("gPt_v1%d",k));
      
      //      RemovePtPoint(ip, gr_v1);
      
      if( gr_v1 == NULL ) continue;
      gr_v1->SetMarkerStyle(imark[is]);

      if( ip == 5 )
	gr_v1->SetMarkerStyle(imark[is]+4);
      gr_v1->SetMarkerColor(icolor);
      gr_v1->SetMarkerSize(imsz[is]);
      gr_v1->SetLineColor(icolor);
      gr_v1->GetXaxis()->SetRangeUser(0.,800.);


      mv1[k]->Add(gr_v1,"lp");
      lg1[k]->AddEntry(gr_v1,  otitle,"lp");
    }
    
    cout << " ==== " << endl;

    for(UInt_t k = 0; k < ybin2; k++){

      TGraphErrors *gr_v2 = (TGraphErrors*)fOpen->Get((TString)Form("gPt_v2%d",k));
      if( gr_v2 == NULL ) continue;

      //      RemovePtPoint(ip,gr_v2);

      
      gr_v2->SetMarkerStyle(imark[is]);
      if( ip == 5 )
	gr_v2->SetMarkerStyle(imark[is]+4);
      gr_v2->SetMarkerColor(icolor);
      gr_v2->SetMarkerSize(imsz[is]);
      gr_v2->SetLineColor(icolor);

      mv2[k]->Add(gr_v2,"lp");
      lg2[k]->AddEntry(gr_v2, otitle ,"lp");
    }
	
    fOpen->Close();
  }

  if( bplot || kFALSE ) {
    ic++; cc = new TCanvas(Form("cc%d",ic),"cc6");
    cc->Divide(1,ngr);
    for(UInt_t i = 0; i < ngr; i++){
      cc->cd(i+1);
      TH2D* ff = (TH2D*)gROOT->Get(Form("hyptacp_%d",i));
      ff->Draw("colz");
    }
  
 
    ic++; cc = new TCanvas(Form("cc%d",ic),"cc6");
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
  }

  if( bplot || kTRUE ) { // gothered plots
    ic++; cc = new TCanvas(Form("cc%d",ic),"cc6");
    cc->Divide(ybin1/2,2);
    for(UInt_t k = 0; k < ybin1; k++){
      cc->cd(k+1);
      mv1[k]->Draw("ALP");
      if( k == ybin1-2 )
	lg1[k]->Draw();
    }

    ic++; cc = new TCanvas(Form("cc%d",ic),"cc6");
    cc->Divide(ybin2,1);
    for(UInt_t k = 0; k < ybin2; k++){
      cc->cd(k+1);
      //    mv2[k]->SetMaximum( 0.2);
      //    mv2[k]->SetMinimum(-0.2);
      mv2[k]->Draw("ALP");
      if( k == ybin2-2 ) lg2[k]->Draw();
    }
  }


  if( bplot || kFALSE ) {
    ic++; cc = new TCanvas(Form("cc%d",ic),"cc6");
    mv1[2]->Draw("ALP");
    lg1[2]->Draw();
    
    ic++; cc = new TCanvas(Form("cc%d",ic),"cc6");
    mv1[3]->Draw("ALP");
    lg1[3]->Draw();

    ic++; cc = new TCanvas(Form("cc%d",ic),"cc6");
    mv1[4]->Draw("ALP");
    lg1[4]->Draw();
  }




  if( bplot || kTRUE ) {
    ic++; cc = new TCanvas(Form("cc%d",ic),"cc6");
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
  
    ic++; cc = new TCanvas(Form("cc%d",ic),"cc6");
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


  for(UInt_t is = 0; is < 4; is++){

    for(UInt_t ip = 0; ip < 4; ip++){
      TString otitle = rsys[is]+"Sn "+fpid[ip];

      UInt_t icolor = (icol[is] + 4*ip);
      if( g_slope[is][ip]->GetN() > 0 ) {
	g_slope[is][ip]->SetMarkerStyle(imark[is]);
	g_slope[is][ip]->SetMarkerColor(icol[is]+2*ip);
	g_slope[is][ip]->SetLineColor(icol[is]+2*ip);
	mv1sl->Add(g_slope[is][ip],"lp");      
	
	lgr4->AddEntry(g_slope[is][ip], otitle, "lp");
      }

      if( g_v2cm[is][ip]->GetN() > 0) {
	g_v2cm[is][ip]->SetMarkerStyle(imark[is]);
	g_v2cm[is][ip]->SetMarkerColor(icol[is]+2*ip);
	g_v2cm[is][ip]->SetLineColor(icol[is]+2*ip);
	mv2cm->Add(g_v2cm[is][ip],"lp");      

	lgr5->AddEntry(g_v2cm[is][ip], otitle, "lp");
      }
    }
    
    if( g_v1slp[is]->GetN() > 0 ){
      g_v1slp[is]->SetMarkerStyle(imark[is]);
      g_v1slp[is]->SetMarkerColor(icol[is]);
      g_v1slp[is]->SetLineColor(icol[is]);
      mv1slp->Add(g_v1slp[is], "lp");
      lgr7->AddEntry(g_v1slp[is], rsys[is]+"Sn","lp");
    }

    if( g_v2mid[is]->GetN() > 0) {
      g_v2mid[is]->SetMarkerStyle(imark[is]);
      g_v2mid[is]->SetMarkerColor(icol[is]);
      g_v2mid[is]->SetLineColor(icol[is]);
      mv2mid->Add(g_v2mid[is],"lp");      
      lgr6->AddEntry(g_v2mid[is], rsys[is]+"Sn", "lp");
    }


    mv1slp->GetXaxis()->SetLimits(0.5, 5.5);
    mv2mid->GetXaxis()->SetLimits(0.5, 5.5);
    for( UInt_t ippp = 0; ippp < 5; ippp++ ) {
      mv1slp->GetXaxis()->SetBinLabel(mv1slp->GetXaxis()->FindBin(ippp+1.), fpid[ippp]);
      mv2mid->GetXaxis()->SetBinLabel(mv2mid->GetXaxis()->FindBin(ippp+1.), fpid[ippp]);
    }
    mv1slp->GetXaxis()->LabelsOption("h");
    mv2mid->GetXaxis()->LabelsOption("h");

  }


  if( bplot || kFALSE) {
    //for NN v1-pt
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadLeftMargin(2.8);
    gStyle->SetPadLeftMargin(0.1);
    gStyle->SetTitleYSize(0.04);
    gStyle->SetTitleYOffset(1.1);

    ic++; cc = new TCanvas(Form("cc%d",ic),"cc6");
    mv2mid->GetYaxis()->SetRangeUser(-0.06, -0.01);
    mv2mid->Draw("ALP");
    lgr6->Draw();

    ic++; cc = new TCanvas(Form("cc%d",ic),"cc6");
    mv1slp->Draw("ALP");

    lgr7->SetX1NDC(0.15);
    lgr7->SetY1NDC(0.64);
    lgr7->SetX2NDC(0.41);
    lgr7->SetY2NDC(0.88);
    lgr7->Draw();
  }

  if( bplot || kFALSE ) {  // individual pt plots  
    for(UInt_t i = 0; i < 8; i++) {
      cc = new TCanvas(Form("cv%d",i),Form("cv%d",i),400,450);
      cc->SetRightMargin(0.02);
      cc->SetLeftMargin(0.12);
      cc->SetTopMargin(0.05);
      mv1[i]->GetXaxis()->SetLabelSize(0.04);
      mv1[i]->GetYaxis()->SetLabelSize(0.04);
      mv1[i]->GetXaxis()->SetTitleSize(0.05);
      mv1[i]->GetYaxis()->SetTitleSize(0.05);
      mv1[i]->GetYaxis()->SetTitleOffset(1.2);

      if( i == 0 )
	mv1[0]->GetYaxis()->SetRangeUser(-0.4,0.);
      else if( i == 1)
	mv1[1]->GetYaxis()->SetRangeUser(-0.4,0.);
      else if( i == 2 )
	mv1[2]->GetYaxis()->SetRangeUser(-0.2,0.15);
      else if( i == 3 )
	mv1[3]->GetYaxis()->SetRangeUser(-0.1,0.15);
      else
	mv1[i]->GetYaxis()->SetRangeUser(-0.1,0.5);

      mv1[i]->Draw("ALP");
      lg1[i]->Draw();

    }

    for(UInt_t i = 0; i < 5; i++) {
      cc = new TCanvas(Form("cv%d",i+10),Form("cv%d",i+10),400,450);
      cc->SetRightMargin(0.02);
      cc->SetLeftMargin(0.12);
      cc->SetTopMargin(0.05);
      mv2[i]->GetXaxis()->SetLabelSize(0.04);
      mv2[i]->GetYaxis()->SetLabelSize(0.04);
      mv2[i]->GetXaxis()->SetTitleSize(0.05);
      mv2[i]->GetYaxis()->SetTitleSize(0.05);
      mv2[i]->GetYaxis()->SetTitleOffset(1.2);

      mv2[i]->GetYaxis()->SetRangeUser(-0.07,0.04);

      mv2[i]->Draw("ALP");
      lg2[i]->Draw();

    }
  }

  if( bplot || kFALSE ){  
    ic++; cc = new TCanvas(Form("cc%d",ic),"v1_target",500,500);
    mv1[1]->SetMaximum(0.);
    mv1[1]->SetMinimum(-0.27);
    mv1[1]->Draw("ALP");

    ic++; cc = new TCanvas(Form("cc%d",ic),"v1_mid",500,500);
    mv1[3]->Draw("ALP");

    ic++; cc = new TCanvas(Form("cc%d",ic),"v1_beam",500,500);
    mv1[5]->Draw("ALP");
    lg1[5]->SetX1NDC(0.5);
    lg1[5]->SetY1NDC(0.16);
    lg1[5]->SetX2NDC(0.95);
    lg1[5]->SetY2NDC(0.40);
    lg1[5]->Draw();
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
  if( ip == 5 ) ptcut = 400.;

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
  
