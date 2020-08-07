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

UInt_t  sysID = 0;
TString partName = "proton";
UInt_t  ip = 2;
gplot gnames[] = {
  // {".v52.10.19" ,"finYPt_","M55to80 no corr","gy_v1"},
  // {".v52.10.19" ,"finYPt_","w/o corr","gu_v1"},
  {".v52.10.21" ,"finYPt_","corr","gu_v1"},
  {".v52.10.20" ,"finYPt_","corr yaw<0","gu_v1"},
  {".v52.10.22" ,"finYPt_","corr yaw>0","gu_v1"},
  // {".v52.10.8" ,"finYPt_","PID_Tight","gu_v1"},
  // {".v52.10.9" ,"finYPt_","PID_Norm","gu_v1"},
  // {".v52.10.14" ,"finYPt_","yaw < 0","gu_v1"},
  // {".v52.10.15" ,"finYPt_","yaw > 0","gu_v1"},
  // {".v52.10.10" ,"finYPt_","|#Phi| < 0","gu_v1"},
  // {".v52.10.11" ,"finYPt_","|#Phi| > 0","gu_v1"},
  // {".v52.10.17" ,"finYPt_","aft yaw > 0","gu_v1"},
  // {".v52.10.5" ,"finYPt_","cor",""}, 
 // {".v52.10.1" ,"finYPt_","w/o cor"},
  // {".v52.10.2" ,"finYPt_","cor"},
  // {".v52.10.3" ,"finYPt_","Ut>0.4 cor"},
  // {".v52.10.4" ,"finYPt_","Ut>0.4 w/o cor"},
  // {".v52.9.5" ,"advYPt_","M15-20"},
  // {".v52.9.6" ,"advYPt_","M20-40"},
  // {".v52.9.7" ,"advYPt_","M20-40"},
  // {".v52.9.8" ,"advYPt_","M30-50"},
  // {".v52.9.9" ,"advYPt_","M40-50"},
  // {".v52.9.10","advYPt_","M60-64"},
  // {".v52.9.11","advYPt_","M60-80"},
  // {".v52.9.12","advYPt_","M65-50"},
  // {".v52.9.13","advYPt_","M50-55"},
  //--
  // {".v52.9.0" ,"advYPt_","y_nn*/y_nn*"},
  // {".v52.8.4" ,"advYPt_","y_nn/y_nn"},
  // {".v52.7.2" ,"advYPt_","#theta <>40-50"},
  // {".v52.8.5" ,"advYPt_","#phi>0"},
  // {".v52.8.6" ,"advYPt_","#phi<0"},
  // {".v52.8.0" ,"advYPt_","#theta,#phi<0"},
  // {".v52.8.1" ,"advYPt_","#theta,#phi>0"},
  //  {".v52.8.2" ,"advYPt_","#theta,|#phi|<30"},
  //  {".v52.8.3" ,"advYPt_","#theta,|#phi|>150"}
  //  {".v52.3.0" ,"advYPt_","y_AA/y_AA"},
  //  {".v52.3.1" ,"advYPt_","AA:|#phi|<30"},
  //  {".v52.3.2" ,"advYPt_","AA:|#phi|>150"}
};

const UInt_t ncp = sizeof(gnames)/sizeof(gplot);

TCanvas *ccv; UInt_t iccv = 0;
std::vector< TString > g_v1label;
auto gev1 = new TGraphErrors();
auto g_v1off = new TGraphErrors();
auto gev2 = new TGraphErrors();
TMultiGraph *mrv1;
TMultiGraph *mrv2;
TLegend *lrv1;
TLegend *lrv2;

FairLogger *logger = FairLogger::GetLogger();
void GetMinimumv2(TGraphErrors *gr, Double_t &min, Double_t &mer);
void Plot_v1offset();
void Plot_v2();


void PlotSystematicError()
{
  gStyle->SetOptStat(0);
  SetStyle();
  SetColor();

  Color_t fcolor = 2;

  TString gtitle = bName[sysID]+partName;

  TGraphErrors *yv1;
  TGraphErrors *yv2;

  std::vector<Double_t> v1_slp ;
  std::vector<Double_t> v1_slpe;
  std::vector<Double_t> v2_max ;
  std::vector<Double_t> v2_maxe;

  mrv1 = new TMultiGraph("mrv1"  ,gtitle+";y_{nrm}; v1");
  mrv2 = new TMultiGraph("mrv2"  ,gtitle+";y_{nrm}; v2");
  lrv1 = new TLegend(0.6 , 0.20, 0.85, 0.65, "");
  lrv2 = new TLegend(0.25 , 0.55, 0.44, 0.9, "");

  TLine *meanV1 = new TLine();
  TBox  *stdV1  = new TBox();

  gev1 -> SetTitle(gtitle+";;v_{11}");
  std::vector< TString > gev1Label;
  UInt_t iv1off = 0;

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

    if( gnames[icp].yv1 == "gy_v1" )
      yv1 = (TGraphErrors*)fOpen->Get("gy_v1");
    else
      yv1 = (TGraphErrors*)fOpen->Get("gu_v1");

    yv2 = (TGraphErrors*)fOpen->Get("gy_v2");

    if( yv1 ) {

      //mid-central M<50
      fv1fit->SetParameter(1,0.3);
      fv1fit->SetParameter(2,-8.4e-09);
      fv1fit->SetParameter(3,-6.5e-05);
      //central M>55
      fv1fit->SetParameter(1,0.2);
      fv1fit->SetParameter(2,-2.4e-02);
      fv1fit->SetParameter(3,-5.8e-02);
      fv1fit->SetLineColor(fcolor);
      yv1->Fit("fv1fit","","",-0.5,0.6); //"Q0","");    
      cout << " fcolor " << fcolor << endl; 
      
      //      yv1->Fit("fv1fitorg","","",-0.7,1.1); // w/o offset

      yv1 -> SetMarkerColor(fcolor);
      yv1 -> SetLineColor(fcolor);
      yv1 -> SetMarkerStyle(20);


      mrv1 -> Add( yv1,"P");
      lrv1 -> AddEntry( yv1, gnames[icp].comment );

      if( icp == 0 ) {
	meanV1 -> SetY1(fv1fit->GetParameter(1));
	meanV1 -> SetY2(fv1fit->GetParameter(1));
	stdV1  -> SetY1(fv1fit->GetParameter(1) +  fv1fit->GetParError(1) );
	stdV1  -> SetY2(fv1fit->GetParameter(1) -  fv1fit->GetParError(1) );
      }
      else {
	v1_slp.push_back ( fv1fit->GetParameter(1) / pow( fv1fit->GetParError(1), 2) );
	v1_slpe.push_back( 1./ pow( fv1fit->GetParError(1), 2 ) ); 
      }

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
      GetMinimumv2(yv2, v2y, v2ye);
      v2_max.push_back( v2y );
      v2_maxe.push_back( v2ye );

      yv2 -> SetMarkerColor( fcolor );
      yv2 -> SetMarkerStyle(20);
      yv2 -> SetLineColor( fcolor );

      mrv2 -> Add( yv2, "LP" );
      lrv2 -> AddEntry( yv2, gnames[icp].comment );

      gev2 -> SetPoint( iv2, iv2, v2y );
      gev2 -> SetPointError( iv2, 0., v2ye );

      iv2++;
    }

    fcolor++;
    if( fcolor == 10 ) fcolor++;
      
    fOpen->Close();
  }
  

  // Draw -------------
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
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
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;

  for( UInt_t i = 0; i < ncp; i++ ){
    auto fbin = gev1 -> GetXaxis() -> FindBin( i*2 );
    if( fbin != 101 ) 
      gev1 -> GetXaxis() -> SetBinLabel( fbin, gev1Label.at(i*2) );

    fbin = gev1 -> GetXaxis() -> FindBin( i*2+1 );
    if( fbin != 101 ) 
      gev1 -> GetXaxis() -> SetBinLabel( fbin, gev1Label.at(i*2+1));
    
  }
  // set style
  ccv  -> GetPad(0)->SetBottomMargin(0.30);
  //  gev1 -> SetLineColor(2);
  gev1 -> SetLineWidth(1504);
  gev1 -> SetFillStyle(3004);
  gev1 -> SetMarkerColor(2);
  gev1 -> SetMarkerSize(1.5);
  gev1 -> SetMarkerStyle(21);
  gev1 -> Draw("ACFP");

  Double_t m_v1slp = TMath::Mean( v1_slp.begin(),  v1_slp.end() );
  Double_t s_v1slp = TMath::Mean( v1_slpe.begin(), v1_slpe.end());
  m_v1slp /= s_v1slp;

  s_v1slp = sqrt(1./s_v1slp * v1_slpe.size());
  


  auto aLineV1  = new TLine( gev1->GetXaxis()->GetXmin(), m_v1slp,
			     gev1->GetXaxis()->GetXmax(), m_v1slp);
  auto aBoxV1   = new TBox(  gev1->GetXaxis()->GetXmin(), m_v1slp+s_v1slp,
			     gev1->GetXaxis()->GetXmax(), m_v1slp-s_v1slp );
  aBoxV1->SetFillStyle(3008);
  aBoxV1->SetFillColor(4);


  LOG(INFO) << "v1 slope = " << m_v1slp << " +- " << s_v1slp << FairLogger::endl;

  // auto aLineY1 = new TLine( 0., mrv1->GetYaxis()->GetXmin(), 
  // 			    0., mrv1->GetYaxis()->GetXmax());

  aLineV1->SetLineColor(1);
  aLineV1->SetLineStyle(3);
  aLineV1->Draw();
  aBoxV1 ->Draw();


  meanV1 -> SetX1( gev1->GetXaxis()->GetXmin() );
  meanV1 -> SetX2( gev1->GetXaxis()->GetXmax() );
  stdV1  -> SetX1( gev1->GetXaxis()->GetXmin() );
  stdV1  -> SetX2( gev1->GetXaxis()->GetXmax() );
  meanV1 -> SetLineColor(kGreen);
  stdV1  -> SetFillStyle(4061);
  stdV1  -> SetFillColor( kGreen -6 );
  meanV1 -> Draw();
  stdV1  -> Draw();


  //  Plot_v2();
  //  Plot_v1offset();
}

void Plot_v1offset()
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

void Plot_v2()
{
  //------
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
  mrv2 -> Draw("ALP");
  lrv2 -> Draw();
}
  //----->>>  v1_0 off set




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
