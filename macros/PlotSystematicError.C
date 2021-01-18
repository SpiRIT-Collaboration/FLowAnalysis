#include "DoFlow.h"
#include "SetStyle.C"
#include "FlowFunction.C"
#include "SetColor.C"

struct gplot{
  TString Version;
  TString fileHeader;
  TString comment;
  TString yv1;
};

//sysID                 0        1        2        3         4
TString  bName[]   = {"132Sn_","108Sn_","124Sn_","112Sn_","100Sn_"};
Double_t sysD[]    = {0.22,    0.09,      0.15,   0.15   , 0.22};
Double_t sysA[]    = {256.,    220.,      236.,   236.   ,  256.};
Double_t sysBN[]   = {82.,      58.,       74.,    62.   ,   50.};
Double_t sysN[]    = {156.,    110.,      136.,   136.   ,  156.};
Size_t   imsz[]    = {1, 1, 1.3, 1.3, 1.3, 1.3, 1.3};

UInt_t  sysID = 1;

gplot gnames[] = {
  //  {".v52.15.36" ,"finYPt_","w/oC&ndf50&#phi<45","gy_v1"},
  // {".v52.15.39" ,"finYPt_","w/oC&ndf20&#phi<45","gy_v1"},
  // {".v52.15.40" ,"finYPt_","w/oC&ndf20&#phi<135","gy_v1"},
  // {".v52.15.35" ,"finYPt_","w/oC&ndf50&#phi<45/135","gy_v1"},
  // {".v52.15.38" ,"finYPt_","w/oC&ndf20&","gy_v1"},
  //  {".v52.15.39" ,"finYPt_","w/oC&ndf20&#phi<45","gy_v1"},
  //  {".v52.15.38" ,"finYPt_","ndf20","gu_v1"},
  {".v52.15.39" ,"finYPt_","ndf20&#phi<45","gu_v1"},
  {".v52.15.40" ,"finYPt_","ndf20&#phi>135","gu_v1"},
  {".v52.15.41" ,"finYPt_","ndf20&45<#phi<135","gu_v1"},
  // {".v52.15.42" ,"finYPt_","ndf50","gu_v1"},
  // {".v52.15.36" ,"finYPt_","ndf50&#phi<45","gu_v1"},
  // {".v52.15.35" ,"finYPt_","ndf50&#phi<45or>135","gu_v1"},
};

const UInt_t ncp = sizeof(gnames)/sizeof(gplot);

TCanvas *ccv; UInt_t iccv = 0;
std::vector< TString > g_v1label;
std::vector< TString > gev1Label;
std::vector<Double_t> v1_slp ;
std::vector<Double_t> v1_slpe;
std::vector<Double_t> v2_max ;
std::vector<Double_t> v2_maxe;
Double_t v1_sigm = 0.;

TLine *meanV1 = new TLine();
TBox  *stdV1  = new TBox();


FairLogger *logger = FairLogger::GetLogger();
void GetMinimumv2(TGraphErrors *gr, Double_t &min, Double_t &mer);
void Plot_v1offset(TGraphErrors *g_v1off);
void Plot_v1(TMultiGraph *mrv2, TLegend* lrv2, TGraphErrors *gev1);
void Plot_v2(TMultiGraph *mrv2, TLegend* lrv2, TGraphErrors *gev2);
void Plot(UInt_t partID);
Double_t* GetMean(TGraphErrors *gr);

void PlotSystematicError()
{

  std::vector<UInt_t> sqPart  = {0,1,2,3};
  for(auto ipart: sqPart) 
    Plot(ipart);
}

void Plot(UInt_t partID)
{
  std::vector< TString > ssPart = {"proton","deuteron","triton","3He"};
  TString partName = ssPart.at(partID);

  gStyle->SetOptStat(0);
  SetStyle();
  SetColor();


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

  gev1 -> SetName("gev1_"+partName);
  gev1 -> SetTitle(gtitle+";;v_{11}");
  UInt_t iv1off = 0;

  auto gev2    = new TGraphErrors();
  gev2 -> SetName("gev2_"+partName);
  gev2 -> SetTitle(";;v2_{max}");

  TFile *fOpen;

  UInt_t iv1 = 0, iv2 = 0;
  for( UInt_t icp = 0; icp < ncp; icp++ ) {
    TString otitle = lsys[sysID] + gnames[icp].comment;

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

      //central M>55
      fv1fit->SetParameter(1,0.2);
      fv1fit->SetParameter(2,-2.4e-02);
      fv1fit->SetParameter(3,-5.8e-02);
      //mid-central M<50
      fv1fit->SetParameter(1,0.3);
      fv1fit->SetParameter(2,-8.4e-09);
      fv1fit->SetParameter(3,-6.5e-05);

      fv1fit->SetLineColor(fcolor);
      yv1->Fit("fv1fit","","",-0.5,0.6); //"Q0","");    
      cout << " fcolor " << fcolor << endl; 
      
      //      yv1->Fit("fv1fitorg","","",-0.7,1.1); // w/o offset

      yv1 -> SetMarkerColor(fcolor);
      yv1 -> SetLineColor(fcolor);
      yv1 -> SetMarkerStyle(20);


      mrv1 -> Add( yv1,"P");
      lrv1 -> AddEntry( yv1, gnames[icp].comment );

      // if( icp == 0 ) {
      // 	meanV1 -> SetY1(fv1fit->GetParameter(1));
      // 	meanV1 -> SetY2(fv1fit->GetParameter(1));
      // 	stdV1  -> SetY1(fv1fit->GetParameter(1) +  fv1fit->GetParError(1) );
      // 	stdV1  -> SetY2(fv1fit->GetParameter(1) -  fv1fit->GetParError(1) );
      // }
      // else {
      // 	v1_slp.push_back ( fv1fit->GetParameter(1) / pow( fv1fit->GetParError(1), 2) );
      // 	v1_slpe.push_back( 1./ pow( fv1fit->GetParError(1), 2 ) ); 
      // 	v1_sigm += 1./ pow( fv1fit->GetParError(1), 2 );
      // }

      gev1 -> SetPoint(iv1, iv1, fv1fit->GetParameter(1) );
      gev1 -> SetPointError(iv1, 0., fv1fit->GetParError(1) );
      gev1Label.push_back( gnames[icp].comment );
      iv1++;
      
      // gev1 -> SetPoint(iv1, iv1, fv1fitorg->GetParameter(1) );
      // gev1 -> SetPointError(iv1, 0., fv1fitorg->GetParError(1) );
      // gev1Label.push_back( "org" + gnames[icp].comment );
      // iv1++;

      g_v1label.push_back( otitle ) ;
      g_v1off -> SetPoint(iv1off, iv1off, fv1fit->GetParameter(3) );
      g_v1off -> SetPointError(iv1off, 0, fv1fit->GetParError(3) );
      iv1off++;
    }
      
    if( yv2 ) {
      for( Int_t iip = (Int_t)yv2->GetN()-1; iip >= 0; iip-- ){

	Double_t xpnt, ypnt;
	ypnt = yv2->GetErrorY(iip);
	if( ypnt > 0.05 )
	  yv2->RemovePoint(iip);
      }

      Double_t v2x, v2y, v2ye;
      //    GetMinimumv2(yv2, v2y, v2ye);
      //      v2_max.push_back( v2y );
      //      v2_maxe.push_back( v2ye );

      fv2fit->SetLineColor( fcolor );

      Double_t v2para0[5][5]={{-0.03, 0.05, 0.02, -0.5, 0.5},
			      { -0.04,0.05, 0.02, -0.5, 0.5},
			      { -0.05,0.05, 0.02, -0.5, 0.5},
			      { -0.07, 0.05, 0.1, -0.4, 0.8},
			      { -0.08, 0.05, 0.1, -0.5, 0.5}};
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
      lrv2 -> AddEntry( yv2, gnames[icp].comment );

      gev2 -> SetPoint( iv2, iv2, v2y );
      gev2 -> SetPointError( iv2, 0., v2ye );

      iv2++;
    }

    fcolor++;
    if( fcolor == 10 ) fcolor++;
      
    fOpen->Close();
  }


  if( gev1 && gev2 ) {
    for( UInt_t i = 0; i < ncp; i++ ){
      auto fbin = gev1 -> GetXaxis() -> FindBin( i );
      if( fbin != 101 ) {
	gev1 -> GetXaxis() -> SetBinLabel( fbin, gev1Label.at(i) );
	gev2 -> GetXaxis() -> SetBinLabel( fbin, gev1Label.at(i) );
      }
    }
    
    Plot_v1(mrv1,lrv1, gev1);
    //  Plot_v1offset(g_v1off);

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
  ccv  -> GetPad(2)->SetBottomMargin(0.30);
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
  ccv  -> GetPad(2)->SetBottomMargin(0.3);
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
  for(auto i : ROOT::TSeqI(gr->GetN()) ) {
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

void Plot_v1offset(TGraphErrors *g_v1off)
{
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
  for(UInt_t i = 0; i < g_v1label.size(); i++ ){
    auto fbin = g_v1off -> GetXaxis() -> FindBin(i);
    if( fbin != 101 )
      g_v1off -> GetXaxis() -> SetBinLabel( fbin, g_v1label.at(i) );
  }

  g_v1off -> SetMarkerStyle(20);
  g_v1off -> SetMarkerColor(3);
  g_v1off -> Draw("AP");
}
