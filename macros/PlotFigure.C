#include "DoFlow.h"
#include "SetStyle.C"
#include "FlowFunction.C"
#include "SetColor.C"

struct gplot{
  TString Version;
  TString fileHeader;
  TString comment;
};

TString  bName[]   = {"132Sn_","108Sn_","124Sn_","112Sn_","pp",    "100Sn_"};
Double_t sysD[]    = {0.22,    0.09,      0.15,   0.15   ,   0,      0.22};
Double_t sysA[]    = {256.,    220.,      236.,   236.   ,   2,    256.};
Double_t sysBN[]   = {82./50., 58./50.,   74./50.,62./50.,   0,     50./50.};
Double_t sysN[]    = {156.,    110.,      136.,   136.   ,   0,    156.};
Double_t sysNA[]   = {156./100.,110./100.,136./100.,136./100., 0,  156./100.};
Size_t   imsz[]    = {1, 1, 1.3, 1.3, 1.3, 1.3, 1.3};
Color_t  pcolor[4] = {2, 3, 4, 5}; // p,d,t,3He


gplot gnames[] = {
  //{".v52.10.16" ,"finYPt_", "M0-50"},
  {".v52.10.28" ,"finYPt_", "b3fm"},
  // {".v52.9.15" ,"advYPt_","M55-80"},
  // {".v52.9.14" ,"advYPt_","M0-50"},
  //  {".v52.8.4" ,"advYPt_","y_sys/y_sys"},
};

TCanvas *ccv; UInt_t iccv = 0;
TString xlabel;
TLatex  plabel;

Double_t *syslabel;
TGraphErrors* g_v1slp[4];
TGraphErrors* g_v2max[4];
TGraphErrors* g_v2n[4];
TGraphErrors* g_v1pslp[4];
TMultiGraph* mrv1[4];
TMultiGraph* mrv2[4];
TLegend*     lrv1[4];
TLegend*     lrv2[4];
TMultiGraph* g_v1sysD;
TLegend*     l_v1sysD;
TMultiGraph* g_v2sysD;
TLegend*     l_v2sysD;
TMultiGraph* g_v1psysD;
TLegend*     l_v1psysD;
TMultiGraph* g_v2nsysD;
TLegend*     l_v2nsysD;

FairLogger *logger = FairLogger::GetLogger();
void GetMinimumv2(TGraphErrors *gr, Double_t &min, Double_t &mer);
void Draw_v1Ratio();
void Draw_v1v2SystemD();
void Draw_v1v2SystemD_fit();
void Draw_Indiv_v1SystemD();
void Draw_Indiv_v2SystemD();
void Draw_Indiv_v2SystemD_Three();
void Draw_Indiv_v1SystemD_Three();
void Draw_v2SystemD();
void Draw_v2nSystemD();
void Draw_v2SystemD_one();

void PlotFigure()
{
  gStyle->SetOptStat(0);
  SetStyle();
  SetColor();
  Color_t fcolor = 2;


  UInt_t ix = 0;
  switch(ix) {
  case 0:
    xlabel = "(N-P)/A";
    syslabel = sysD;
    break;
  case 1:
    xlabel = "N";
    syslabel = sysN;
    break;
  case 2:
    xlabel = "N/Z(Beam)";
    syslabel = sysBN;
    break;
  case 3:
    xlabel = "N/Z";
    syslabel = sysNA;
    break;
  }

  TGraphErrors *yv1;
  TGraphErrors *yv2;

  Double_t v1_slp [4][4];
  Double_t v1_slpe[4][4];
  Double_t v2_max [4][4];
  Double_t v2_maxe[4][4];
  Double_t v20[4][4];
  Double_t v20e[4][4];
  Double_t v21[4][4];
  Double_t v21e[4][4];
  Double_t v2n[4][4];
  Double_t v2ne[4][4];

  for(UInt_t i = 0; i < 4; i++ ){
    mrv1[i] = new TMultiGraph(Form("mrv1_%d",i)  ,";y/y_{nn}-1; v1");
    mrv2[i] = new TMultiGraph(Form("mrv2_%d",i)  ,";y/y_{nn}-1; v2");
    lrv1[i] = new TLegend(0.6 , 0.20, 0.90, 0.40, lsys[i]);
    lrv2[i] = new TLegend(0.4 , 0.70, 0.55, 0.9,  lsys[i]);
  }

  auto g_v1off = new TGraphErrors();
  std::vector< TString > g_v1label;
  UInt_t iv1off = 0;

  TFile *fOpen;

  for( UInt_t is = 0; is < 4; is++ ) {

    if( is == 2 ) continue;

    for( UInt_t ip = 0; ip < 4; ip++ ) {

      TString fname = gnames[0].fileHeader + bName[is] + fpid[ip] + gnames[0].Version + ".root";

      if( !gSystem->FindFile("data", fname) ) {
	LOG(ERROR) << fname << " is not found " << FairLogger::endl;
	continue;
      }
      else {
	fOpen = TFile::Open( fname );
	LOG(INFO) << fname << " is opened. " << FairLogger::endl;
      }


      yv1 = (TGraphErrors*)fOpen->Get("gy_v1");
      yv2 = (TGraphErrors*)fOpen->Get("gy_v2");
      //      yv1 = (TGraphErrors*)fOpen->Get("gu_v1");
      //      yv2 = (TGraphErrors*)fOpen->Get("gu_v2");

      if( yv1 && yv2 ) {
	fcolor = pcolor[ip];
	yv1->SetMarkerColor(fcolor);
	yv1->SetMarkerStyle(imark[is]);
	yv1->SetLineColor(fcolor);

	yv2->SetMarkerColor(fcolor);
        yv2->SetMarkerStyle(imark[is]);
        yv2->SetMarkerSize(imsz[is]);
        yv2->SetLineColor(fcolor);

	for( Int_t j = (Int_t)yv2->GetN()-1; j >= 0; j-- ) {
	  Double_t xx, yy;
	  yv2->GetPoint(j, xx, yy);
	  if( xx > 1. )
	    yv2->RemovePoint(j);
	}
      }	


      mrv1[is] -> Add(yv1,"p");
      lrv1[is] -> AddEntry(yv1,  lpid[ip] ,"lp");
      // mrv1[ip] -> Add(yv1,"p");
      // lrv1[ip] -> AddEntry(yv1,  lsys[is] ,"lp");

      fv1fit->SetLineColor(fcolor);
      fv1fit->SetParameter(1,0.2);
      yv1->Fit("fv1fit","","",-0.6,1.2); //"Q0","");    

      v1_slp [is][ip]   = fv1fit->GetParameter(1);
      v1_slpe[is][ip]   = fv1fit->GetParError(1);

      if( ip != 3 ) {
	g_v1label.push_back( bName[is] + lpid[ip] ) ;
	g_v1off -> SetPoint(iv1off, iv1off, fv1fit->GetParameter(3) );
	g_v1off -> SetPointError(iv1off, 0, fv1fit->GetParError(3) );
	iv1off++;
      }

      if( yv2 != NULL ) {
        for( Int_t iip = (Int_t)yv2->GetN()-1; iip >= 0; iip-- ){

          Double_t xpnt, ypnt;
          ypnt = yv2->GetErrorY(iip);
          if( ypnt > 0.05 )
            yv2->RemovePoint(iip);
        }

	Double_t v2x, v2y, v2ye;
        GetMinimumv2(yv2, v2_max[is][ip], v2_maxe[is][ip]);	
	yv2   ->SetLineColor(fcolor);
	fv2fit->SetLineColor(fcolor);

	fv2fit->SetParameter(0, -0.01);
	fv2fit->SetParameter(1,  0.15);
	fv2fit->SetParameter(2,  0.);

	if( ip == 3 )
	  yv2->Fit("fv2fit","","",-0.48,1.);
	else
	  yv2->Fit("fv2fit","","",-0.5,0.5);

	mrv2[is]->Add(yv2,"p");
	lrv2[is]->AddEntry(yv2, lpid[ip], "p");
	// mrv2[ip]->Add(yv2,"p");
	// lrv2[ip]->AddEntry(yv2, lsys[is], "p");

	v20[is][ip]  = fv2fit->GetParameter(0);
	v20e[is][ip] = fv2fit->GetParError(0);

	v21[is][ip]  = fv2fit->GetParameter(1);
	v21e[is][ip] = fv2fit->GetParError(1);
	
	v2n[is][ip]  = abs(v21[is][ip]) + abs(v20[is][ip]);
	v2ne[is][ip] = 0.; sqrt( pow(v21[is][ip],2) + pow(v20[is][ip],2) );

      }
      
      fOpen->Close();
    }
    fcolor++;
    if( fcolor == 10 ) fcolor++;
    
  }

  //----------
  for( UInt_t ip = 0; ip < 4; ip++ ){
    g_v1slp[ip]  = new TGraphErrors();
    g_v1slp[ip]  -> SetTitle(";"+xlabel+";");
    g_v2max[ip]  = new TGraphErrors();
    g_v2max[ip]  -> SetTitle(";"+xlabel+"; -v20");
    g_v2n[ip]    = new TGraphErrors();
    g_v2n[ip]    -> SetTitle(";"+xlabel+"; -v2n");

    g_v1pslp[ip] = new TGraphErrors();
  }

  
  g_v1sysD = new TMultiGraph("g_v1slpD", ";"+xlabel+"; v_{11}");
  l_v1sysD = new TLegend(0.2, 0.7, 0.5, 0.9, "");
  g_v2sysD = new TMultiGraph("g_v2sysD", ";"+xlabel+"; -v20");
  l_v2sysD = new TLegend(0.35, 0.13, 0.7, 0.33, "");

  g_v1psysD = new TMultiGraph("g_v1pslpD", ";"+xlabel+"; v_{11}/ v_{11}(proton)");
  l_v1psysD = new TLegend(0.2, 0.7, 0.5, 0.9, "");

  g_v2nsysD = new TMultiGraph("g_v2bsysD", ";"+xlabel+"; -v2n");
  l_v2nsysD = new TLegend(0.35, 0.13, 0.7, 0.33, "");

  for( UInt_t ip = 0; ip < 4; ip++ ) {

    UInt_t iss = 0;
    for( UInt_t is = 0; is < 4; is++ ) {
      
      if( is == 2 ) continue;

      if( kTRUE ) {
	LOG(INFO) << " v1slp " << v1_slp[is][ip]
		  << " v2max " << v2_max[is][ip]
		  << FairLogger::endl;

	g_v1slp[ip] -> SetPoint     ( iss,  *(syslabel+is), v1_slp[is][ip] ); 
	g_v1slp[ip] -> SetPointError( iss, 0,   v1_slpe[is][ip] );

	// g_v2max[ip] -> SetPoint     ( iss,  *(syslabel+is), -v2_max[is][ip] ); 
	// g_v2max[ip] -> SetPointError( iss,  0., v2_maxe[is][ip] ); 
	g_v2max[ip] -> SetPoint     ( iss,  *(syslabel+is), -v20[is][ip] ); 
	g_v2max[ip] -> SetPointError( iss,  0., v20e[is][ip] ); 

	g_v2n[ip]   -> SetPoint     ( iss,  *(syslabel+is), v2n[is][ip] );
	g_v2n[ip]   -> SetPointError( iss,  0., v2ne[is][ip] );
	

	if( ip > 0 ) {
	  auto ratio = v1_slp[is][ip]/v1_slp[is][0];
	  auto ratie = ratio * sqrt( pow(v1_slpe[is][ip]/v1_slp[is][ip],2) + pow(v1_slpe[is][0]/v1_slp[is][0],2) );
	  g_v1pslp[ip] -> SetPoint     ( iss, *(syslabel+is), ratio );
	  g_v1pslp[ip] -> SetPointError( iss,        0., ratie );
	}

	iss++;
      }
    }

    g_v1slp[ip]->SetMarkerStyle(20);
    g_v1slp[ip]->SetMarkerSize(1.5);
    g_v1slp[ip]->SetMarkerColor(icol[ip]);
    g_v1slp[ip]->SetLineColor(icol[ip]);

    g_v2max[ip]->SetMarkerStyle(20);
    g_v2max[ip]->SetMarkerSize(1.5);
    g_v2max[ip]->SetMarkerColor(icol[ip]);
    g_v2max[ip]->SetLineColor(icol[ip]);

    g_v2n[ip]->SetMarkerStyle(20);
    g_v2n[ip]->SetMarkerSize(1.5);
    g_v2n[ip]->SetMarkerColor(icol[ip]);
    g_v2n[ip]->SetLineColor(icol[ip]);

    
    g_v1sysD -> Add( g_v1slp[ip], "p" );
    l_v1sysD -> AddEntry( g_v1slp[ip], lpid[ip], "P" );

    g_v2sysD -> Add( g_v2max[ip],"p");
    l_v2sysD -> AddEntry( g_v2max[ip], lpid[ip] ,"P");

    g_v2nsysD -> Add( g_v2n[ip],"p");
    l_v2nsysD -> AddEntry( g_v2n[ip], lpid[ip] ,"P");

    g_v1pslp[ip]->SetMarkerStyle(20);
    g_v1pslp[ip]->SetMarkerSize(1.5);
    g_v1pslp[ip]->SetMarkerColor(icol[ip]);
    g_v1pslp[ip]->SetLineColor(icol[ip]);

    
    if( ip > 1 ) {
      g_v1psysD -> Add( g_v1pslp[ip], "p" );
      l_v1psysD -> AddEntry( g_v1slp[ip], lpid[ip], "P" );
    }

  }

  //
  // Draw function

  plabel.SetTextAlign(13);
  plabel.SetTextSize(0.08);


  // Draw_Indiv_v1SystemD();

  //  Draw_v2SystemD_one();
  //Draw_v1Ratio();

  // fitting results for v1 and v2 vs y
  Draw_v1v2SystemD_fit();

  //  Draw_v1v2SystemD();
  //  Draw_v2SystemD();

  Draw_Indiv_v1SystemD_Three(); 
  Draw_Indiv_v2SystemD_Three();

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


//------- Draw ---------------
void Draw_v2SystemD_one()
{
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 700, 500); iccv++;
  auto gmv2n = new TMultiGraph("gmv2n", ";"+xlabel+"; v2n");
  auto lmv2n = new TLegend(0.35,0.13,0.7,0.33,"");

  gmv2n -> Add( g_v2n[0],  "p");
  gmv2n -> Add( g_v2n[1],  "p");
  gmv2n -> Add( g_v2n[2],  "p");
  gmv2n -> Add( g_v2n[3],  "p");

  lmv2n -> AddEntry( g_v2n[0], lpid[0], "p");
  lmv2n -> AddEntry( g_v2n[1], lpid[1], "p");
  lmv2n -> AddEntry( g_v2n[2], lpid[2], "p");
  lmv2n -> AddEntry( g_v2n[3], lpid[3], "p");

  gmv2n -> Draw("ALP");
  lmv2n -> Draw();
  
}
void Draw_v2SystemD()
{
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 500, 800); iccv++;

  Double_t BYs = 0.1;            // Bottom Y space
  Double_t BYm = 0.1;            // Bottom Y mergin
  Double_t TYm = 0.01;            // Top Y mergin

  Double_t Ny  = 4;               // Number of pads along Y  
  Double_t H   = (1.0 - (TYm+BYm+BYs))/Ny; // pad height
  cout << " H= " << H << endl;

  Double_t Yl = BYs;  Double_t Yu = Yl + H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *p1 = new TPad("p1", "p1", 0., Yl, 1.,  Yu, 0, 0, 0);
  p1->SetTopMargin(0);
  p1->SetBottomMargin(BYm);
  p1->Draw();


  Yl = Yu;  Yu = Yl + H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *p2 = new TPad("p2", "p2", 0., Yl, 1., Yu, 0, 0, 0);
  p2->SetTopMargin(0);
  p2->SetBottomMargin(0);
  p2->Draw();

  Yl = Yu;  Yu = Yl + H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *p3 = new TPad("p3", "p3", 0., Yl, 1., Yu, 0, 0, 0);  
  p3->SetTopMargin(0);
  p3->SetBottomMargin(0);
  p3->Draw();

  Yl = Yu;  Yu = Yl +  H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *p4 = new TPad("p4", "p4", 0., Yl, 1., Yu, 0, 0, 0);  
  p4->SetTopMargin(TYm);
  p4->SetBottomMargin(0);
  p4->Draw();


  //  gStyle->SetLabelSize(0.08);
  g_v2max[0]->SetTitle(";"+xlabel+";"); 
  g_v2max[3]->SetTitle(";; v_{20} ");

  p1->cd(); 
  g_v2max[0] -> GetYaxis() -> SetLabelOffset(0.01);
  g_v2max[0] -> GetYaxis() -> SetLabelSize(0.08);
  g_v2max[0] -> GetXaxis() -> SetLabelSize(0.1);
  g_v2max[0] -> GetXaxis() -> SetTitleOffset(1.);
  g_v2max[0] -> GetXaxis() -> SetTitleSize(1);
  //  g_v2max[0] -> GetXaxis() -> SetTitleOffSet(0.5);
  g_v2max[0] -> Draw("AP");  
  plabel.DrawLatexNDC(0.2, 0.9, fpid[0]);


  p2->cd();
  g_v2max[1] -> GetYaxis() -> SetLabelOffset(0.01);
  g_v2max[1] -> GetYaxis() -> SetLabelSize(0.1);
  g_v2max[1] -> Draw("AP");  
  plabel.DrawLatexNDC(0.2,0.9,fpid[1]);

  p3->cd(); 
  g_v2max[2] -> GetYaxis() -> SetLabelOffset(0.01);
  g_v2max[2] -> GetYaxis() -> SetLabelSize(0.1);
  g_v2max[2] -> Draw("AP");  
  plabel.DrawLatexNDC(0.2,0.9,fpid[2]);

  p4->cd(); 
  g_v2max[3] -> GetYaxis() -> SetTitleOffset(1.);
  g_v2max[3] -> GetYaxis() -> SetTitleSize(0.1);
  g_v2max[3] -> GetYaxis() -> SetLabelOffset(0.01);
  g_v2max[3] -> GetYaxis() -> SetLabelSize(0.1);
  g_v2max[3] -> Draw("AP");  
  plabel.DrawLatexNDC(0.2,0.9,fpid[3]);

}

void Draw_v2nSystemD()
{
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 500, 800); iccv++;

  Double_t BYs = 0.1;            // Bottom Y space
  Double_t BYm = 0.1;            // Bottom Y mergin
  Double_t TYm = 0.01;            // Top Y mergin

  Double_t Ny  = 4;               // Number of pads along Y  
  Double_t H   = (1.0 - (TYm+BYm+BYs))/Ny; // pad height
  cout << " H= " << H << endl;

  Double_t Yl = BYs;  Double_t Yu = Yl + H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *p1 = new TPad("p1", "p1", 0., Yl, 1.,  Yu, 0, 0, 0);
  p1->SetTopMargin(0);
  p1->SetBottomMargin(BYm);
  p1->Draw();


  Yl = Yu;  Yu = Yl + H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *p2 = new TPad("p2", "p2", 0., Yl, 1., Yu, 0, 0, 0);
  p2->SetTopMargin(0);
  p2->SetBottomMargin(0);
  p2->Draw();

  Yl = Yu;  Yu = Yl + H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *p3 = new TPad("p3", "p3", 0., Yl, 1., Yu, 0, 0, 0);  
  p3->SetTopMargin(0);
  p3->SetBottomMargin(0);
  p3->Draw();

  Yl = Yu;  Yu = Yl +  H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *p4 = new TPad("p4", "p4", 0., Yl, 1., Yu, 0, 0, 0);  
  p4->SetTopMargin(TYm);
  p4->SetBottomMargin(0);
  p4->Draw();


  //  gStyle->SetLabelSize(0.08);
  g_v2n[0]->SetTitle(";"+xlabel+";"); 
  g_v2n[3]->SetTitle(";; v_{21} ");

  p1->cd(); 
  g_v2n[0] -> GetYaxis() -> SetLabelOffset(0.01);
  g_v2n[0] -> GetYaxis() -> SetLabelSize(0.08);
  g_v2n[0] -> GetXaxis() -> SetLabelSize(0.1);
  g_v2n[0] -> GetXaxis() -> SetTitleOffset(1.);
  g_v2n[0] -> GetXaxis() -> SetTitleSize(1);
  //  g_v2n[0] -> GetXaxis() -> SetTitleOffSet(0.5);
  g_v2n[0] -> Draw("AP");  
  plabel.DrawLatexNDC(0.2, 0.9, fpid[0]);


  p2->cd();
  g_v2n[1] -> GetYaxis() -> SetLabelOffset(0.01);
  g_v2n[1] -> GetYaxis() -> SetLabelSize(0.1);
  g_v2n[1] -> Draw("AP");  
  plabel.DrawLatexNDC(0.2,0.9,fpid[1]);

  p3->cd(); 
  g_v2n[2] -> GetYaxis() -> SetLabelOffset(0.01);
  g_v2n[2] -> GetYaxis() -> SetLabelSize(0.1);
  g_v2n[2] -> Draw("AP");  
  plabel.DrawLatexNDC(0.2,0.9,fpid[2]);

  p4->cd(); 
  g_v2n[3] -> GetYaxis() -> SetTitleOffset(1.);
  g_v2n[3] -> GetYaxis() -> SetTitleSize(0.1);
  g_v2n[3] -> GetYaxis() -> SetLabelOffset(0.01);
  g_v2n[3] -> GetYaxis() -> SetLabelSize(0.1);
  g_v2n[3] -> Draw("AP");  
  plabel.DrawLatexNDC(0.2,0.9,fpid[3]);

}


void Draw_v1Ratio()
{
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 500, 800); iccv++;

  Double_t BYs = 0.1;            // Bottom Y space
  Double_t TYm = 0.05;            // Top Y mergin
  Double_t BYm = 0.1;            // Bottom Y mergin
  Double_t Ny  = 3;               // Number of pads along Y  
  Double_t H   = (1.0 - (TYm+BYm+BYs))/Ny; // pad height

  Double_t Yl = BYs;   Double_t Yu = Yl + H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *pp1 = new TPad("pp1", "pp1", 0., Yl, 1.,  Yu, 0, 0, 0);
  pp1->SetTopMargin(0);
  pp1->Draw();

  Yl = Yu;  Yu = Yl + H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *pp2 = new TPad("pp2", "pp2", 0., Yl, 1., Yu, 0, 0, 0);
  pp2->SetTopMargin(0);
  pp2->SetBottomMargin(0);
  pp2->Draw();

  Yl = Yu;  Yu = Yl + H ;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *pp3 = new TPad("pp3", "pp3", 0., Yl, 1., Yu, 0, 0, 0);  
  pp3->SetBottomMargin(0);
  pp3->Draw();

  plabel.SetTextSize(0.1);

  g_v1pslp[1]->SetTitle(";"+xlabel+";"); 
  g_v1pslp[3]->SetTitle(";; v11/ v11_{p} ");

  pp1->cd(); 
  plabel.DrawLatexNDC(0.2,0.9,fpid[1]);
  g_v1pslp[1] -> GetYaxis() -> SetLabelSize(0.05);
  g_v1pslp[1] -> GetXaxis() -> SetLabelSize(0.05);
  g_v1pslp[1] -> Draw("AP");  

  pp2->cd();
  plabel.DrawLatexNDC(0.2,0.9,fpid[3]);
  g_v1pslp[2] -> GetYaxis() -> SetLabelSize(0.05);
  g_v1pslp[2] -> Draw("AP");  

  pp3->cd(); 
  plabel.DrawLatexNDC(0.2,0.9,fpid[4]);
  g_v1pslp[3] -> GetYaxis() -> SetLabelSize(0.05);
  g_v1pslp[3] -> Draw("AP");  
  //---
} 


void Draw_v1v2SystemD_fit()
{

  plabel.SetTextAlign(13);
  plabel.SetTextSize(0.08);

  for(UInt_t i = 0; i < 4; i++ ){

    if( mrv1[i]->GetListOfGraphs() != NULL ) {
      ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;

      mrv1[i]->SetTitle();

      mrv1[i]->GetXaxis()->SetRangeUser(-0.8, 1.2);
      mrv1[i]->GetYaxis()->SetRangeUser(-0.5, 0.6);
      mrv1[i]->Draw("ALP");
      lrv1[i]->Draw();
    }

    auto aLineX1 = new TLine( mrv1[i]->GetXaxis()->GetXmin(), 0.,
			      mrv1[i]->GetXaxis()->GetXmax(), 0.);
    auto aLineY1 = new TLine( 0., mrv1[i]->GetYaxis()->GetXmin(), 
			      0., mrv1[i]->GetYaxis()->GetXmax());

    aLineX1->SetLineColor(1);
    aLineX1->SetLineStyle(3);
    aLineX1->Draw();
    
    aLineY1->SetLineColor(1);
    aLineY1->SetLineStyle(3);
    aLineY1->Draw();
  }



  for(UInt_t i = 0; i <4; i++ ){

    if( mrv2[i]->GetListOfGraphs() != NULL ) {
      ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;

      mrv2[i]->GetYaxis()->SetRangeUser(-0.13, 0.02);
      mrv2[i]->Draw("ALP");
      lrv2[i]->Draw();
    }
  }
}

void Draw_v1v2SystemD()
{
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
  g_v2sysD -> Draw("ALP");
  l_v2sysD -> Draw();


  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
  g_v1psysD -> Draw("ALP");
  l_v1psysD -> Draw();


  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv)); iccv++;
  g_v1sysD -> Draw("ALP");
  l_v1sysD -> Draw();
  //<<<---
}

void Draw_Indiv_v1SystemD_Three() 
{
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 500, 800); iccv++;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  plabel.SetTextAlign(13);
  plabel.SetTextSize(0.08);

  Double_t BYs = 0.1;            // Bottom Y space
  Double_t TYm = 0.01;            // Top Y mergin
  Double_t BYm = 0.2;            // Bottom Y mergin
  Double_t Ny  = 3;               // Number of pads along Y  
  Double_t H   = (1.0 - (TYm+BYs))/Ny; // pad height
  cout << " H= " << H << endl;

  Double_t Yl = 0;  Double_t Yu = Yl + H ;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *p1 = new TPad("p1", "p1", 0., Yl, 1.,  Yu, 0, 0, 0);
  p1->SetTopMargin(0);
  p1->SetBottomMargin(BYm);
  p1->Draw();


  Yl = Yu;  Yu = Yl + H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *p2 = new TPad("p2", "p2", 0., Yl, 1., Yu, 0, 0, 0);
  p2->SetTopMargin(0);
  p2->SetBottomMargin(0);
  p2->Draw();

  Yl = Yu;  Yu = Yl + H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *p3 = new TPad("p3", "p3", 0., Yl, 1., Yu, 0, 0, 0);  
  p3->SetTopMargin(0);
  p3->SetBottomMargin(0);
  p3->Draw();


  //  gStyle->SetLabelSize(0.08);
  p1->cd(); 
  g_v1slp[0] -> GetYaxis() -> SetLabelOffset(0.01);
  g_v1slp[0] -> GetYaxis() -> SetLabelSize(0.08);
  g_v1slp[0] -> GetXaxis() -> SetLabelSize(0.1);
  g_v1slp[0] -> GetXaxis() -> SetTitleOffset(0.8);
  g_v1slp[0] -> GetXaxis() -> SetTitleSize(0.1);
  g_v1slp[0] -> GetXaxis() -> SetTitle(xlabel);
  g_v1slp[0] -> Draw("AP");  
  plabel.DrawLatexNDC(0.2, 0.4, fpid[0]);

  p2->cd();
  g_v1slp[1] -> GetYaxis() -> SetLabelOffset(0.01);
  g_v1slp[1] -> GetYaxis() -> SetLabelSize(0.1);
  g_v1slp[1] -> Draw("AP");  
  plabel.DrawLatexNDC(0.2,0.9,fpid[1]);


  p3->cd(); 
  auto mg_T3He = new TMultiGraph("mg_T3He",";;v_{11}");
  auto l_T3He  = new TLegend(0.2,0.7,0.4,0.9,"");
  mg_T3He->Add(g_v1slp[2], "AP");
  mg_T3He->Add(g_v1slp[3], "AP");
  l_T3He->AddEntry(g_v1slp[2],"Triton");
  l_T3He->AddEntry(g_v1slp[3],"^{3}He");

  mg_T3He -> GetXaxis() -> SetLimits(g_v1slp[0]->GetXaxis()->GetXmin(), g_v1slp[0]->GetXaxis()->GetXmax());
  mg_T3He -> SetMaximum(0.794);
  mg_T3He -> GetYaxis() -> SetTitleOffset(0.75);
  mg_T3He -> GetYaxis() -> SetTitleSize(0.1);
  mg_T3He -> GetYaxis() -> SetLabelOffset(0.01);
  mg_T3He -> GetYaxis() -> SetLabelSize(0.1);
  mg_T3He -> Draw("AP");  
  l_T3He  -> Draw();
}

void Draw_Indiv_v1SystemD() 
{
  UInt_t pixH = 800;
  UInt_t pixW = 500;
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), pixW, pixH); iccv++;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  plabel.SetTextAlign(13);
  plabel.SetTextSize(0.08);

  Double_t BYs = 0.1;            // Bottom Y space
  Double_t BYm = 0.1;            // Bottom Y mergin
  Double_t TYm = 0.01;            // Top Y mergin

  Double_t Ny  = 4;               // Number of pads along Y  
  Double_t H   = (1.0 - (TYm+BYm+BYs))/Ny; // pad height
  cout << " H= " << H << endl;

  Double_t Yl = BYs;  Double_t Yu = Yl + H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *p1 = new TPad("p1", "p1", 0., Yl, 1.,  Yu, 0, 0, 0);
  p1->SetTopMargin(0);
  p1->SetBottomMargin(BYm);
  p1->Draw();


  Yl = Yu;  Yu = Yl + H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *p2 = new TPad("p2", "p2", 0., Yl, 1., Yu, 0, 0, 0);
  p2->SetTopMargin(0);
  p2->SetBottomMargin(0);
  p2->Draw();

  Yl = Yu;  Yu = Yl + H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *p3 = new TPad("p3", "p3", 0., Yl, 1., Yu, 0, 0, 0);  
  p3->SetTopMargin(0);
  p3->SetBottomMargin(0);
  p3->Draw();

  Yl = Yu;  Yu = Yl +  H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *p4 = new TPad("p4", "p4", 0., Yl, 1., Yu, 0, 0, 0);  
  p4->SetTopMargin(TYm);
  p4->SetBottomMargin(0);
  p4->Draw();


  //  gStyle->SetLabelSize(0.08);
  g_v1slp[0]->SetTitle(";"+xlabel+";"); 
  g_v1slp[3]->SetTitle(";; v_{11} ");

  p1->cd(); 
  g_v1slp[0] -> GetYaxis() -> SetLabelOffset(0.01);
  g_v1slp[0] -> GetYaxis() -> SetLabelSize(0.08);
  g_v1slp[0] -> GetXaxis() -> SetLabelSize(0.1);
  g_v1slp[0] -> GetXaxis() -> SetTitleOffset(1.);
  g_v1slp[0] -> GetXaxis() -> SetTitleSize(1);
  //  g_v1slp[0] -> GetXaxis() -> SetTitleOffSet(0.5);
  g_v1slp[0] -> Draw("AP");  
  plabel.DrawLatexNDC(0.2, 0.9, fpid[0]);


  p2->cd();
  g_v1slp[1] -> GetYaxis() -> SetLabelOffset(0.01);
  g_v1slp[1] -> GetYaxis() -> SetLabelSize(0.1);
  g_v1slp[1] -> Draw("AP");  
  plabel.DrawLatexNDC(0.2,0.9,fpid[1]);

  p3->cd(); 
  g_v1slp[2] -> GetYaxis() -> SetLabelOffset(0.01);
  g_v1slp[2] -> GetYaxis() -> SetLabelSize(0.1);
  g_v1slp[2] -> Draw("AP");  
  plabel.DrawLatexNDC(0.2,0.9,fpid[2]);

  p4->cd(); 
  g_v1slp[3] -> GetYaxis() -> SetTitleOffset(1.);
  g_v1slp[3] -> GetYaxis() -> SetTitleSize(0.1);
  g_v1slp[3] -> GetYaxis() -> SetLabelOffset(0.01);
  g_v1slp[3] -> GetYaxis() -> SetLabelSize(0.1);
  g_v1slp[3] -> Draw("AP");  
  plabel.DrawLatexNDC(0.2,0.9,fpid[3]);
}

void Draw_Indiv_v2SystemD_Three() 
{
  ccv = new TCanvas(Form("ccv%d",iccv),Form("ccv%d",iccv), 500, 800); iccv++;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  plabel.SetTextAlign(13);
  plabel.SetTextSize(0.08);

  Double_t BYs = 0.1;            // Bottom Y space
  Double_t TYm = 0.01;            // Top Y mergin
  Double_t BYm = 0.2;            // Bottom Y mergin
  Double_t Ny  = 3;               // Number of pads along Y  
  Double_t H   = (1.0 - (TYm+BYs))/Ny; // pad height
  cout << " H= " << H << endl;

  Double_t Yl = 0;  Double_t Yu = Yl + H ;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *p1 = new TPad("p1", "p1", 0., Yl, 1.,  Yu, 0, 0, 0);
  p1->SetTopMargin(0);
  p1->SetBottomMargin(BYm);
  p1->Draw();


  Yl = Yu;  Yu = Yl + H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *p2 = new TPad("p2", "p2", 0., Yl, 1., Yu, 0, 0, 0);
  p2->SetTopMargin(0);
  p2->SetBottomMargin(0);
  p2->Draw();

  Yl = Yu;  Yu = Yl + H;
  cout << " low " << Yl << " up " << Yu << endl;
  TPad *p3 = new TPad("p3", "p3", 0., Yl, 1., Yu, 0, 0, 0);  
  p3->SetTopMargin(0);
  p3->SetBottomMargin(0);
  p3->Draw();


  //  gStyle->SetLabelSize(0.08);
  p1->cd(); 
  g_v2max[0] -> GetYaxis() -> SetLabelOffset(0.01);
  g_v2max[0] -> GetYaxis() -> SetLabelSize(0.08);
  g_v2max[0] -> GetXaxis() -> SetLabelSize(0.1);
  g_v2max[0] -> GetXaxis() -> SetTitleOffset(0.8);
  g_v2max[0] -> GetXaxis() -> SetTitleSize(0.1);
  g_v2max[0] -> GetXaxis() -> SetTitle(xlabel);
  g_v2max[0] -> Draw("AP");  
  plabel.DrawLatexNDC(0.2, 0.4, fpid[0]);

  p2->cd();
  g_v2max[1] -> GetYaxis() -> SetLabelOffset(0.01);
  g_v2max[1] -> GetYaxis() -> SetLabelSize(0.1);
  g_v2max[1] -> Draw("AP");  
  plabel.DrawLatexNDC(0.2,0.9,fpid[1]);


  p3->cd(); 
  auto mg_T3He = new TMultiGraph("mg_T3He",";;-v_{20}");
  auto l_T3He  = new TLegend(0.2,0.7,0.4,0.9,"");
  mg_T3He->Add(g_v2max[2], "AP");
  mg_T3He->Add(g_v2max[3], "AP");
  l_T3He->AddEntry(g_v2max[2],"Triton");
  l_T3He->AddEntry(g_v2max[3],"^{3}He");

  mg_T3He -> GetXaxis() -> SetLimits(g_v1slp[0]->GetXaxis()->GetXmin(), g_v1slp[0]->GetXaxis()->GetXmax());
  //  mg_T3He -> SetMaximum(0.794);
  mg_T3He -> GetYaxis() -> SetTitleOffset(0.75);
  mg_T3He -> GetYaxis() -> SetTitleSize(0.1);
  mg_T3He -> GetYaxis() -> SetLabelOffset(0.01);
  mg_T3He -> GetYaxis() -> SetLabelSize(0.1);
  mg_T3He -> Draw("AP");  
  l_T3He  -> Draw();
  
}
