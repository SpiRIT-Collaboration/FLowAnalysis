#include "DoFlow.h"
#include "SetStyle.C"
#include "FlowFunction.C"
#include "SetColor.C"

struct gplot{
  TString Version;
  TString fileHeader;
  TString yv1;
  TString config1;
  TString config2;
  TString config3;
};

//Size_t   imsz[]    = {1, 1, 1.3, 1.3, 1.3, 1.3, 1.3};

std::vector< TString > ssPart = {"proton","deuteron","triton","3He","4He"};

gplot gnames[] = {
  {".v52.15.122","finYPt_","gu_v1","5fm","|#phi|<45" ,"acr"},
  {".v52.15.123","finYPt_","gu_v1","5fm","|#phi|>135","acr"},
    // {".v52.15.82","finYPt_","gu_v1","m42to56","|#phi|<45" ,"acr"},
    // {".v52.15.83","finYPt_","gu_v1","m42to56","|#phi|>135","acr"},
    // {".v52.15.70","finYPt_","gu_v1","m55to80","|#phi|<45" ,"acr"},
    // {".v52.15.71","finYPt_","gu_v1","m55to80","|#phi|>135","acr"},
   // {".v52.15.72","finYPt_","gu_v1","m50to65","|#phi|<45" ,"acr"},
   // {".v52.15.73","finYPt_","gu_v1","m50to65","|#phi|>135","acr"},
  //  {".v52.15.74","finYPt_","gu_v1","m40to55","|#phi|<45" ,"acr"},
  //  {".v52.15.75","finYPt_","gu_v1","m40to55","|#phi|>135","acr"},
   // {".v52.15.76","finYPt_","gu_v1","m35to35","|#phi|<45" ,"acr"},
   // {".v52.15.77","finYPt_","gu_v1","m35to35","|#phi|>135","acr"},
   // {".v52.15.78","finYPt_","gu_v1","m20to40","|#phi|<45" ,"acr"},
   // {".v52.15.79","finYPt_","gu_v1","m20to40","|#phi|>135","acr"},
  //   {".v52.15.80","finYPt_","gu_v1","m0to35" ,"|#phi|<45" ,"acr"},
  //   {".v52.15.81","finYPt_","gu_v1","m0to35" ,"|#phi|>135","acr"},
  //
  //{".v52.15.51.v52.15.52" ,"finYPt_" ,"gu_v1","m55to80","ave","y_nn"},
  //{".v52.15.53.v52.15.54" ,"finYPt_" ,"gu_v1","m50to65","ave","y_nn"},
  //{".v52.15.55.v52.15.56" ,"finYPt_" ,"gu_v1","m40to55","ave","y_nn"},
  //{".v52.15.57.v52.15.58" ,"finYPt_" ,"gu_v1","m35to45","ave","y_nn"},
  //{".v52.15.59.v52.15.60" ,"finYPt_" ,"gu_v1","m20to40","ave","y_nn"},
  //{".v52.15.61.v52.15.62" ,"finYPt_" ,"gu_v1","m0to35" ,"ave","y_nn"},
  //
};


const UInt_t ncp = sizeof(gnames)/sizeof(gplot);

std::vector<UInt_t> sqSys = {0,1,3};

TCanvas *ccv; UInt_t iccv = 0;
std::vector< TString > g_v1label;
std::vector< TString > gev1Label;
std::vector< TString > gev2Label;
std::vector<Double_t> v1_slp ;
std::vector<Double_t> v1_slpe;
std::vector<Double_t> v2_max ;
std::vector<Double_t> v2_maxe;
Double_t v1_sigm = 0.;

TLine *meanV1 = new TLine();
TBox  *stdV1  = new TBox();

TMultiGraph* mv1off; 
TLegend*    lv1off; 

UInt_t  sysID = 0; // 0:(132Sm) 

FairLogger *logger = FairLogger::GetLogger();
void GetMinimumv2(TGraphErrors *gr, Double_t &min, Double_t &mer);
void Plot_MultiGraph(TMultiGraph *mgr);
void Plot_v1offset(TGraphErrors *g_v1off);
void Plot_v1(TMultiGraph *mrv2, TLegend* lrv2, TGraphErrors *gev1);
void Plot_v2(TMultiGraph *mrv2, TLegend* lrv2, TGraphErrors *gev2);
void Plot(UInt_t partID);
void TakeAverage(UInt_t partID);
Double_t* GetMean(TGraphErrors *gr);
TGraphErrors* LoadData(TString fname, TString gname);
TH1D*         LoadHistgram(TString fname, TString gname);

void PlotSystematicError()
{
  gStyle->SetOptStat(0);
  SetStyle();
  SetColor();

  for(auto isys: sqSys) {
    sysID = isys;

    mv1off = new TMultiGraph();
    lv1off = new TLegend(0.2, 0.7, 0.5, 1.,"");    
    mv1off -> SetName( Form("mv1off_%d",sysID) );
    mv1off -> SetTitle( lsys[sysID] );
    auto gtmp = new TGraph();
    gtmp -> SetPoint(0,-1.,0.);
    gtmp -> SetPoint(1,(Double_t)ncp,0.);
    mv1off -> Add(gtmp,"");
    

    std::vector<UInt_t> sqPart  = {0,1,2,3,4};//{0,1,2,3,4};
    for(auto ipart: sqPart) {

      TakeAverage(ipart);

        Plot(ipart);

    }

    //    Plot_MultiGraph(mv1off);
  }

}

void Plot(UInt_t partID)
{
  TString partName = ssPart.at(partID);

  Color_t fcolor = 2;

  TString gtitle = bName[sysID]+partName;

  TGraphErrors *yv1;
  TGraphErrors *yv2;


  auto mrv1 = new TMultiGraph("mrv1_"+partName  ,gtitle+";y_{nrm}; v1");
  auto mrv2 = new TMultiGraph("mrv2_"+partName  ,gtitle+";y_{nrm}; v2");
  auto lrv1 = new TLegend(0.6 , 0.20, 0.85, 0.65, "");
  auto lrv2 = new TLegend(0.25 , 0.55, 0.44, 0.9, "");

  auto gev1    = new TGraphErrors();
  auto g_v1off = new TGraphErrors();
  g_v1off -> SetTitle(gtitle+";;v_{11}off");


  gev1 -> SetName("gev1_"+partName);
  gev1 -> SetTitle(gtitle+";;v_{11}");
  UInt_t iv1off = 0;

  auto gev2    = new TGraphErrors();
  gev2 -> SetName("gev2_"+partName);
  gev2 -> SetTitle(";;v2_{max}");

  TFile *fOpen;

  UInt_t iv1 = 0, iv2 = 0;

  for( UInt_t icp = 0; icp < ncp; icp++ ) {

    TString fname = gnames[icp].fileHeader + bName[sysID] + partName + gnames[icp].Version + ".root";

    if( !gSystem->FindFile("data", fname) ) {
      LOG(ERROR) <<gnames[icp].fileHeader + bName[sysID] + partName + gnames[icp].Version + ".root" 
		 << " file is not found " << FairLogger::endl;
      continue;
    }
    else {
      fOpen = TFile::Open( fname );
      LOG(INFO) << fname << " is opened. " << FairLogger::endl;
    }

    if( gnames[icp].yv1 == "gy_v1" ) {
      yv1 = (TGraphErrors*)fOpen->Get("gy_v1");
      yv2 = (TGraphErrors*)fOpen->Get("gy_v2"); 
    }
    else {
      yv1 = (TGraphErrors*)fOpen->Get("gu_v1");
      yv2 = (TGraphErrors*)fOpen->Get("gu_v2");
    }


    if( yv1 ) {

      yv1 -> SetMarkerColor(fcolor);
      yv1 -> SetLineColor(fcolor);
      yv1 -> SetMarkerStyle(20);

      mrv1 -> Add( yv1,"P");
      lrv1 -> AddEntry( yv1, gnames[icp].config1 +"&"+ gnames[icp].config2 );


      TGraphErrors *yv1_rev;
      if( sysID==1 || sysID==3 ) {	
	yv1_rev = new TGraphErrors();
	for( UInt_t i = 0; i < yv1->GetN(); i++ ) {
	  Double_t x = 0, y = 0;
	  yv1->GetPoint( yv1->GetN()-1-i, x, y );
	  yv1_rev->SetPoint( i, -x, -y );
	}
	
	yv1_rev->SetMarkerColor(fcolor);
	yv1_rev->SetMarkerStyle(24);

	mrv1 -> Add( yv1_rev,"P");
	lrv1 -> AddEntry( yv1, gnames[icp].config1+ gnames[icp].config2+"_rev" );
      }

      //central M>55
      fv1fit->SetParameter(1,0.2);
      fv1fit->SetParameter(2,-2.4e-02);
      fv1fit->SetParameter(3,-5.8e-02);
      //mid-central M<50
      fv1fit->SetParameter(1,0.3);
      fv1fit->SetParameter(2,-8.4e-09);
      fv1fit->SetParameter(3,-6.5e-05);

      fv1fit->SetLineColor(fcolor);


      yv1->Fit("fv1fit","","", v1fit[0], v1fit[1]); //"Q0","");    
      

      gev1 -> SetPoint(iv1, iv1, fv1fit->GetParameter(1) );
      gev1 -> SetPointError(iv1, 0., fv1fit->GetParError(1) );
      gev1Label.push_back( gnames[icp].config1 + gnames[icp].config2 );
      iv1++;
      
      g_v1label.push_back( gnames[icp].config1 + gnames[icp].config2 ) ;
      g_v1off -> SetPoint(iv1off, iv1off, fv1fit->GetParameter(3) );
      g_v1off -> SetPointError(iv1off, 0, fv1fit->GetParError(3) );
      iv1off++;

	
      if( sysID==1 || sysID==3 ) {	
	yv1_rev->Fit("fv1fit","","",v1fit[0], v1fit[1]); //"Q0","");    
	
	gev1 -> SetPoint(iv1, iv1, fv1fit->GetParameter(1) );
	gev1 -> SetPointError(iv1, 0., fv1fit->GetParError(1) );
	gev1Label.push_back(gnames[icp].config1 + gnames[icp].config2 + "_rev" );
	iv1++;
      }
    }

    
  
      
    if( yv2 ) {
      for( Int_t iip = (Int_t)yv2->GetN()-1; iip >= 0; iip-- ){

	Double_t xpnt, ypnt;
	ypnt = yv2->GetErrorY(iip);
	if( ypnt > 0.05 )
	  yv2->RemovePoint(iip);
      }

      Double_t v2x, v2y, v2ye;
      fv2fit->SetLineColor( fcolor );

      fv2fit->SetParameter(0,v2para0[partID][0]);
      fv2fit->SetParameter(1,v2para0[partID][1]);
      fv2fit->SetParameter(2,v2para0[partID][2]);

      yv2 -> Fit("fv2fit","","",v2para0[partID][3], v2para0[partID][4]);
      v2y  = fv2fit->GetParameter(0);
      v2ye = fv2fit->GetParError(0);
 
      yv2 -> SetMarkerColor( fcolor );
      yv2 -> SetMarkerStyle(20);
      yv2 -> SetLineColor( fcolor );

      mrv2 -> Add( yv2, "P" );
      lrv2 -> AddEntry( yv2, gnames[icp].config1 );

      gev2 -> SetPoint( iv2, iv2, v2y );
      gev2 -> SetPointError( iv2, 0., v2ye );
      gev2Label.push_back( gnames[icp].config1 + gnames[icp].config2);
      iv2++;
    }

    fcolor++;
    if( fcolor == 10 ) fcolor++;
      
    fOpen->Close();
  }

  if( g_v1off->GetN() > 0 ) {
    g_v1off -> SetMarkerStyle(20);
    g_v1off -> SetMarkerColor(icol[partID]);
    g_v1off -> SetLineColor(icol[partID]);
    mv1off -> Add( g_v1off, "lp");
    lv1off -> AddEntry( g_v1off, partName);
  }

  if( gev1 && gev2 ) {
    for( UInt_t i = 0; i < gev1->GetN(); i++ ){
      auto fbin = gev1 -> GetXaxis() -> FindBin( i );
      if( fbin != 101 ) 
	gev1 -> GetXaxis() -> SetBinLabel( fbin, gev1Label.at(i) );
    }

    for( UInt_t i = 0; i < gev2->GetN(); i++ ){
      auto fbin = gev2 -> GetXaxis() -> FindBin( i );
      if( fbin != 101 ) 
	gev2 -> GetXaxis() -> SetBinLabel( fbin, gev2Label.at(i) );
    }
    

    Plot_v1(mrv1,lrv1, gev1);

    //	Plot_v1offset(g_v1off);
    Plot_v2(mrv2,lrv2, gev2);

  }
}


void Plot_v1(TMultiGraph* mrv1, TLegend* lrv1, TGraphErrors *gev1)
{
  // Draw -------------
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 600, 1000); iccv++;
  ccv -> Divide(1,2);
  //  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;

  
  ccv -> cd(1);
  mrv1 -> Draw("ALP");
  lrv1 -> Draw();

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


  //------
  ccv -> cd(2);
  //  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;

  // set style
  ccv  -> GetPad(2)->SetBottomMargin(0.3);
  //  gev1 -> SetLineColor(2);
  gev1 -> SetLineWidth(1504);
  gev1 -> SetFillStyle(3004);
  gev1 -> SetMarkerColor(2);
  gev1 -> SetMarkerSize(1.5);
  gev1 -> SetMarkerStyle(21);

  Double_t* ave = GetMean(gev1);
  auto aLineV1  = new TLine( gev1->GetXaxis()->GetXmin(), ave[0],
			     gev1->GetXaxis()->GetXmax(), ave[0]);
  auto aBoxV1   = new TBox(  gev1->GetXaxis()->GetXmin(), ave[0]+ave[1],
			     gev1->GetXaxis()->GetXmax(), ave[0]-ave[1] );
  aBoxV1->SetFillStyle(3008);
  aBoxV1->SetFillColor(4);

  TString fave= Form("mean %4.2f+-%5.2f(%f)",ave[0],ave[1],ave[1]/ave[0]);
  auto latex = new TLatex(0.5, ave[0], fave );


  LOG(INFO) << "v1 slope = " << ave[0] << " +- " << ave[1] << FairLogger::endl;

  // auto aLineY1 = new TLine( 0., mrv1->GetYaxis()->GetXmin(), 
  // 			    0., mrv1->GetYaxis()->GetXmax());

  aLineV1->SetLineColor(1);
  aLineV1->SetLineStyle(3);

  meanV1 -> SetX1( gev1->GetXaxis()->GetXmin() );
  meanV1 -> SetX2( gev1->GetXaxis()->GetXmax() );
  stdV1  -> SetX1( gev1->GetXaxis()->GetXmin() );
  stdV1  -> SetX2( gev1->GetXaxis()->GetXmax() );
  meanV1 -> SetLineColor(kGreen);
  stdV1  -> SetFillStyle(4061);
  stdV1  -> SetFillColor( kGreen -6 );

  Double_t ymax =  ave[0]+ave[1] > gev1->GetYaxis()->GetXmax() ?  ave[0]+ave[1]*1.2  : gev1->GetYaxis()->GetXmax();
  Double_t ymin =  ave[0]-ave[1] < gev1->GetYaxis()->GetXmin() ?  ave[0]-ave[1]*0.97 : gev1->GetYaxis()->GetXmin();

  ccv -> cd(2);
  gev1->GetYaxis()->SetRangeUser(ymin, ymax);
  gev1 -> Draw("ACFP");
  aLineV1->Draw();
  aBoxV1 ->Draw();
  meanV1 -> Draw();
  latex -> Draw();
  //stdV1  -> Draw();

}


void Plot_v2(TMultiGraph* mrv2, TLegend *lrv2,  TGraphErrors *gev2)
{
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 600, 1000); iccv++;
  ccv -> Divide(1,2);
  //ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;

  ccv -> cd(1);
  mrv2 -> Draw("AP");
  lrv2 -> Draw();

  //ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
  ccv  -> cd(2);
  ccv  -> GetPad(2)->SetBottomMargin(0.25);
  gev2 -> SetLineWidth(1504);
  gev2 -> SetFillStyle(3004);
  gev2 -> SetMarkerColor(2);
  gev2 -> SetMarkerSize(1.5);
  gev2 -> SetMarkerStyle(21);

  gev2 -> Draw("AP");

  gev2 -> Print();
  cout << " gname " << gev2->GetName() << endl;

  Double_t* ave = GetMean(gev2);
  cout << gev2->GetName() << " mean "  << ave[0] << " std " << ave[1] << endl;

  auto aLineV  = new TLine( gev2->GetXaxis()->GetXmin(), ave[0],
			     gev2->GetXaxis()->GetXmax(), ave[0]);
  auto aBoxV   = new TBox(  gev2->GetXaxis()->GetXmin(), ave[0]+ave[1],
			     gev2->GetXaxis()->GetXmax(), ave[0]-ave[1] );
  aLineV->SetLineColor(8);
  aBoxV->SetFillStyle(3008);
  aBoxV->SetFillColor(4);

  aLineV->Draw();
  aBoxV->Draw();

}

Double_t* GetMean(TGraphErrors *gr)
{
  Double_t *ave = new Double_t[2];

  Double_t weight = 0.;
  Double_t vals = 0.;
  Double_t vale = 0.;
  for(auto i : ROOT::TSeqI(gr->GetN()) ) { //@@@ move the last one
    Double_t x,y;
    gr->GetPoint(i, x, y);
    Double_t ye = gr->GetErrorY(i);
    if( ye == 0 ) ye = 1.;
    Double_t weight = 1./pow(ye, 2);
    vals += y * weight;
    vale += weight; 
  }

  ave[0] = vals / vale; 

  vals = 0;
  for(auto i : ROOT::TSeqI(gr->GetN()) ) {
    Double_t x,y;
    gr->GetPoint(i, x, y);
    Double_t ye = gr->GetErrorY(i);
    if( ye == 0 ) ye = 1.;
    Double_t weight = 1./pow(ye, 2);
    Double_t delt = y - ave[0];
    vals += pow( delt, 2) * weight;
  }

  ave[1] = sqrt( vals /( (Double_t)(gr->GetN() -1)* vale ) );  

  return ave;
}


void GetMinimumv2(TGraphErrors *gr, Double_t &min, Double_t &mer)
{
  Double_t x, y;
  min = 9.;
  UInt_t mid = 0;
  for(UInt_t i = 0; i < gr->GetN(); i++) {
    gr->GetPoint(i, x, y);
    if( abs(y) > 0.2 ) {
      gr->RemovePoint(i); 
      i = 0;
      continue;
    }

    if( y < min && y < 0) {
      min = y;
      mid = i;
    }
  }


  mer = gr->GetErrorY(mid);
}

void Plot_MultiGraph(TMultiGraph *mgr)
{

  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
  for(UInt_t i = 0; i < g_v1label.size(); i++ ){
    auto fbin = mgr -> GetXaxis() -> FindBin(i);
    if( fbin != 101 )
      mgr -> GetXaxis() -> SetBinLabel( fbin, g_v1label.at(i) );
  }

  ccv -> GetPad(0)->SetBottomMargin(0.3);
  mgr -> Draw("AP");
  lv1off -> Draw();
}

void Plot_v1offset(TGraphErrors *g_v1off)
{
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
  for(UInt_t i = 0; i < g_v1label.size(); i++ ){
    auto fbin = g_v1off -> GetXaxis() -> FindBin(i);
    if( fbin != 101 )
      g_v1off -> GetXaxis() -> SetBinLabel( fbin, g_v1label.at(i) );
  }

  ccv  -> GetPad(0)->SetBottomMargin(0.3);
  g_v1off -> Draw("AP");
}


void TakeAverage(UInt_t partID)
{
  LOG(INFO) << " TakeAverage " << partID << FairLogger::endl;

  TString partName = ssPart.at(partID);

  fv1fit->SetParameter(1,0.2);
  fv1fit->SetParameter(2,-2.4e-02);
  fv1fit->SetParameter(3,-5.8e-02);
  
  fv2fit->SetParameter(0,v2para0[partID][0]);
  fv2fit->SetParameter(1,v2para0[partID][1]);
  fv2fit->SetParameter(2,v2para0[partID][2]);


  auto mrv1 = new TMultiGraph();
  auto mrv2 = new TMultiGraph();

  TGraphErrors* yv1a[2];
  TGraphErrors* yv2a[2];
  TGraphErrors* utv1_0[2];
  TGraphErrors* utv1_1[2];
  TGraphErrors* utv2[2];
  TH1D*         hmult[2];
  Color_t fcolor = 3;

  UInt_t iv1 = 0, iv2 = 0;
  for( UInt_t icp = 0; icp < 2; icp++ ) {

    TString fname = gnames[icp].fileHeader + bName[sysID] + partName + gnames[icp].Version + ".root";

    yv1a[icp] = LoadData(fname, "gu_v1");
    if( yv1a[icp] ) { 
      yv1a[icp] -> SetName(Form("gu_v1%d",icp));

      yv1a[icp] -> SetMarkerColor(fcolor);
      yv1a[icp] -> SetLineColor(fcolor);
      yv1a[icp] -> SetMarkerStyle(20);

      fv1fit -> SetLineColor(fcolor);
      yv1a[icp] -> Fit("fv1fit","","",-0.5,0.6); //"Q0","");    
      mrv1 -> Add(yv1a[icp],"p");
    }

    yv2a[icp] = LoadData(fname, "gu_v2");
    if( yv2a[icp] ) {
      yv2a[icp] -> SetName(Form("gu_v2%d",icp));

      yv2a[icp] -> SetMarkerColor(fcolor);
      yv2a[icp] -> SetLineColor(fcolor);
      yv2a[icp] -> SetMarkerStyle(20);
      
      fv2fit -> SetLineColor(fcolor);
      yv2a[icp] -> Fit("fv2fit","","",v2para0[partID][3], v2para0[partID][4]);
      mrv2 -> Add(yv2a[icp],"p");
    }

    hmult[icp] =  LoadHistgram(fname,"hmult");

    utv1_0[icp] = LoadData(fname, "g_utv1_0");
    utv1_0[icp] -> SetName(Form("g_utv1_0%d",icp));
    utv1_1[icp] = LoadData(fname, "g_utv1_1");
    utv1_1[icp] -> SetName(Form("g_utv1_1%d",icp));
    utv2[icp]   = LoadData(fname, "g_utv2");
    utv2[icp]   -> SetName(Form("g_utv2%d",icp));
    fcolor++;
  }

  TString foutname = "data/ut_" + gnames[0].fileHeader + bName[sysID] + partName + gnames[0].Version + gnames[1].Version + ".root";
  TFile *fFile = TFile::Open(foutname,"recreate");

  TGraphErrors *yv1 = new TGraphErrors();
  yv1 -> SetName("gu_v1");
  yv1 -> SetTitle( yv1a[0]->GetTitle());
  TGraphErrors *yv2 = new TGraphErrors();
  yv2 -> SetName("gu_v2");
  yv2 -> SetTitle(yv2a[0]->GetTitle());
  TGraphErrors *g_utv1_0 = new TGraphErrors();
  g_utv1_0 -> SetName("g_utv1_0");
  g_utv1_0 -> SetTitle(utv1_0[0]->GetTitle());
  TGraphErrors *g_utv1_1 = new TGraphErrors();
  g_utv1_1 -> SetName("g_utv1_1") ;
  g_utv1_1 -> SetTitle(utv1_1[0]->GetTitle());
  TGraphErrors *g_utv2 = new TGraphErrors();
  g_utv2   -> SetName("g_utv2");
  g_utv2   -> SetTitle(utv2[0]->GetTitle());

  if( yv1a[0] && yv1a[1]) {
    for(auto ive : ROOT::TSeqI(2) )for( auto ipo : ROOT::TSeqI(yv1a[0]->GetN()) ) {
	Double_t x0,y0,y0e;
	Double_t x1,y1,y1e;
	yv1a[0]->GetPoint(ipo, x0, y0);
	y0e = yv1a[0]->GetErrorY(ipo);
	yv1a[1]->GetPoint(ipo, x1, y1);
	y1e = yv1a[1]->GetErrorY(ipo);
	
	if( y0 != 0 && y1 != 0) {
	  Double_t wave = ( y0/pow(y0e,2) + y1/pow(y1e,2) )/ ( 1./pow(y0e,2) + 1./pow(y1e,2) );
	  Double_t werr = sqrt(1./ ( 1./pow(y0e,2) + 1./pow(y1e,2) ) );
	  
	  yv1->SetPoint(ipo, x0, wave);
	  yv1->SetPointError(ipo, 0, werr);
	}
      }


    fcolor = 2;
    yv1 -> SetMarkerColor(fcolor);
    yv1 -> SetLineColor(fcolor);
    yv1 -> SetMarkerStyle(20);

    fv1fit -> SetLineColor(fcolor);
    yv1->Fit("fv1fit","","",-0.5,0.6); //"Q0","");    

    mrv1 -> Add(yv1,"p");

  }

  if( yv2a[0] && yv2a[1] ) {
    for(auto ive : ROOT::TSeqI(2) )for( auto ipo : ROOT::TSeqI(yv2a[0]->GetN()) ) {
	Double_t x0,y0,y0e;
	Double_t x1,y1,y1e;
	yv2a[0]->GetPoint(ipo, x0, y0);
	y0e = yv2a[0]->GetErrorY(ipo);
	yv2a[1]->GetPoint(ipo, x1, y1);
	y1e = yv2a[1]->GetErrorY(ipo);

	if( y0 != 0 && y1 != 0) {
	  Double_t wave = ( y0/pow(y0e,2) + y1/pow(y1e,2) )/ ( 1./pow(y0e,2) + 1./pow(y1e,2) );
	  Double_t werr = sqrt(1./ ( 1./pow(y0e,2) + 1./pow(y1e,2) ) );
	  
	  yv2->SetPoint(ipo, x0, wave);
	  yv2->SetPointError(ipo, 0, werr);
	}
      }

    yv2 -> SetMarkerColor(fcolor);
    yv2 -> SetLineColor(fcolor);
    yv2 -> SetMarkerStyle(20);

    fv2fit -> SetLineColor(fcolor);
    yv2 -> Fit("fv2fit","","",v2para0[partID][3], v2para0[partID][4]);
    mrv2 -> Add(yv2,"p");
  }

  if( utv1_0[0] &&utv1_0[1] ) {
    for(auto ive : ROOT::TSeqI(2) )for( auto ipo : ROOT::TSeqI(utv1_0[0]->GetN()) ) {
	Double_t x0,y0,y0e;
	Double_t x1,y1,y1e;
	utv1_0[0]->GetPoint(ipo, x0, y0);
	y0e = utv1_0[0]->GetErrorY(ipo);
	utv1_0[1]->GetPoint(ipo, x1, y1);
	y1e = utv1_0[1]->GetErrorY(ipo);

	if( y0 != 0 && y1 != 0) {
	  Double_t wave = ( y0/pow(y0e,2) + y1/pow(y1e,2) )/ ( 1./pow(y0e,2) + 1./pow(y1e,2) );
	  Double_t werr = sqrt(1./ ( 1./pow(y0e,2) + 1./pow(y1e,2) ) );
	  
	  g_utv1_0->SetPoint(ipo, x0, wave);
	  g_utv1_0->SetPointError(ipo, 0, werr);
	}
      }	
  }
  if( utv1_1[0] &&utv1_1[1] ) {
    for(auto ive : ROOT::TSeqI(2) )for( auto ipo : ROOT::TSeqI(utv1_1[0]->GetN()) ) {
	Double_t x0,y0,y0e;
	Double_t x1,y1,y1e;
	utv1_1[0]->GetPoint(ipo, x0, y0);
	y0e = utv1_1[0]->GetErrorY(ipo);
	utv1_1[1]->GetPoint(ipo, x1, y1);
	y1e = utv1_1[1]->GetErrorY(ipo);

	if( y0 != 0 && y1 != 0) {
	  Double_t wave = ( y0/pow(y0e,2) + y1/pow(y1e,2) )/ ( 1./pow(y0e,2) + 1./pow(y1e,2) );
	  Double_t werr = sqrt(1./ ( 1./pow(y0e,2) + 1./pow(y1e,2) ) );
	  
	  g_utv1_1->SetPoint(ipo, x0, wave);
	  g_utv1_1->SetPointError(ipo, 0, werr);
	}
      }	
  }
  if( utv2[0] &&utv2[1] ) {
    for(auto ive : ROOT::TSeqI(2) )for( auto ipo : ROOT::TSeqI(utv2[0]->GetN()) ) {
	Double_t x0,y0,y0e;
	Double_t x1,y1,y1e;
	utv2[0]->GetPoint(ipo, x0, y0);
	y0e = utv2[0]->GetErrorY(ipo);
	utv2[1]->GetPoint(ipo, x1, y1);
	y1e = utv2[1]->GetErrorY(ipo);

	if( y0 != 0 && y1 != 0) {
	  Double_t wave = ( y0/pow(y0e,2) + y1/pow(y1e,2) )/ ( 1./pow(y0e,2) + 1./pow(y1e,2) );
	  Double_t werr = sqrt(1./ ( 1./pow(y0e,2) + 1./pow(y1e,2) ) );
	  
	  g_utv2->SetPoint(ipo, x0, wave);
	  g_utv2->SetPointError(ipo, 0, werr);
	}
      }	
  }

  // ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 600, 1000); iccv++;
  // ccv -> Divide(1,2);

  // ccv -> cd(1);
  // mrv1 -> Draw("ALP");

  // ccv ->cd(2);
  // mrv2 -> Draw("ALP");

  hmult[0] -> SetName("hmult");
  hmult[0] -> Write();
  hmult[1] -> SetName("hmult_1");
  hmult[1] -> Write();
  yv1 -> Write();
  yv2 -> Write();
  g_utv1_0 -> Write();
  g_utv1_1 -> Write();
  g_utv2   -> Write();

  for( auto i : ROOT::TSeqI(2) ) {
    yv1a[i] -> Write();
    yv2a[i] -> Write();
    utv1_0[i]->Write();
    utv1_1[i]->Write();
    utv2[i]->Write();
  }

  LOG(INFO) << fFile->GetName() << " is created. " << FairLogger::endl;
  fFile -> Close();

}

TGraphErrors* LoadData(TString fname, TString gname)
{
  TGraphErrors* grv = NULL;

  TFile *fOpen;

  TString ffname = fname;
  if( !gSystem->FindFile("data", ffname) ) {
    LOG(ERROR) << fname << " is not found " << FairLogger::endl;
    return NULL;
  }
  else {
    fOpen = TFile::Open( ffname );
    LOG(INFO) << fname << " is opened. " << FairLogger::endl;
  }

  grv =  (TGraphErrors*)fOpen->Get(gname);
  if( grv == NULL ) return NULL;

  fOpen->Close();

  return grv;
}

TH1D* LoadHistgram(TString fname, TString gname)
{
  TH1D* mhst = NULL;

  TFile *fOpen;

  TString ffname = fname;
  if( !gSystem->FindFile("data", ffname) ) {
    LOG(ERROR) << fname << " is not found " << FairLogger::endl;
    return NULL;
  }
  else {
    fOpen = TFile::Open( ffname );
    LOG(INFO) << fname << " is opened. " << FairLogger::endl;
  }

  mhst =  (TH1D*)fOpen->Get(gname);
  if( mhst == NULL ) return NULL;

  mhst -> SetDirectory(gROOT);

  fOpen->Close();

  return mhst;
}
