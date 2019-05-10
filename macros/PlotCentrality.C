#include "DoFlow.h"

//==================================================
//-- plot configuration
//--------------------------------------------------
  // --> Plotting selection
Bool_t bsys[]  = { 1, 1, 1, 1};
Bool_t bpid[]  = { 0, 0, 0, 0, 0, 0, 1}; //0:p, 1:d, 2:t, 3:3He, 4:4He, 5:n 6:H
Bool_t bcnt[]  = { 1, 0, 0}; 
UInt_t cntw = 3;

UInt_t  bver[]  = { 1, 0, 0, 0};
TString sVer[]  = {".v25.0.0",".v25.0.1", ".v25.0.3", ".v25.0.4", ".v23.1.13",".AMD:"};
TString sName   = "mlt_"; //"cosYPt_132Sn_";
TString bName[] = {"132Sn","108Sn","124Sn","112Sn"};

Bool_t amdEOS[]= {0, 0};
TString amdName[] = {"SLy4",
		     "SLy4-L108"};

TString amdHeader[] = {"amd_132Sn124Sn270AMeV_cluster_",
		       "amd_108Sn112Sn270AMeV_cluster_"};
//==================================================

UInt_t   ccvid = 0;
TF1 *lslope = new TF1("lslope","[0]+[1]*x",-0.1,0.2);

// --> configuration
TString fsys[] = {"132Sn+124Sn","108Sn+112Sn","124Sn+112Sn","112Sn+124Sn"};
TString rsys[] = {"132",        "108",        "124",        "112"};
TString fpid[] = {"proton","deuteron","triton","3He","4He","neutron","H"};  

Size_t  imsz[]   = {1, 1, 1.3, 1.3, 1.3, 1.3, 1.3};
Color_t icol[]   = { kRed, kBlue, kSpring, kMagenta, kOrange, kViolet};
Color_t icolnn[] = { kPink, kBlue+1, kGreen+2, kViolet-1};

TH1I *hhmult;

void PlotCentrality()
{
  gStyle->SetOptStat(0);

  TFile *fOpen;
  TMultiGraph *mcos1 = new TMultiGraph("mcos1","");
  TMultiGraph *mcos2 = new TMultiGraph("mcos2","");

  auto  lgr0 = new TLegend(0.65, 0.7 , 0.90, 0.9);
  auto  lgr1 = new TLegend(0.54, 0.13, 0.87, 0.4);
  auto  lgr2 = new TLegend(0.54, 0.13, 0.87, 0.4);



  UInt_t icc = 0;
  TCanvas *cc = new TCanvas(Form("cc%d",icc),Form("cc%d",icc)); icc++;

  UInt_t ip = 0;
  for( UInt_t is = 0; is < 4; is++ ){
    if( !bsys[is] ) continue;

    for(UInt_t iz = 0; iz < (UInt_t)sizeof(bver)/sizeof(UInt_t); iz++){
      if( !bver[iz] ) continue;

      TString label = bName[is] + sVer[iz];
      TString fname = "data/"+ sName + label + ".root";

      fOpen = TFile::Open(fname);
    
      if(fOpen == NULL ) continue;

      std::cout << fname << " is opened. " << std::endl;

      hhmult    = (TH1I*)fOpen->Get("hmult");
      hhmult -> SetName(Form("hmult_%d",is));
      hhmult -> SetNormFactor(1);
      hhmult -> SetLineColor(icol[ip]);
      hhmult -> SetMaximum(100000);
      hhmult -> Draw( iopt[ip] );
      lgr0   -> AddEntry( hhmult, label, "lp" );

      auto hgv_mcos1 = (TGraphErrors*)fOpen->Get("gv_mcos1");
      auto hgv_mcos2 = (TGraphErrors*)fOpen->Get("gv_mcos2");

      hgv_mcos1 -> SetName(Form("hgv_mcos1_%d",is));
      hgv_mcos1 -> SetMarkerStyle(20);
      hgv_mcos1 -> SetMarkerColor(icol[ip]);
      hgv_mcos1 -> SetLineColor(icol[ip]);

      hgv_mcos2 -> SetName(Form("hgv_mcos2_%d",is));
      hgv_mcos2 -> SetMarkerStyle(20);
      hgv_mcos2 -> SetMarkerColor(icol[ip]);
      hgv_mcos2 -> SetLineColor(icol[ip]);

      mcos1->Add( hgv_mcos1 );
      lgr1 ->AddEntry( hgv_mcos1, label, "lp");

      mcos2->Add( hgv_mcos2 );      
      lgr2 ->AddEntry( hgv_mcos2, label, "lp");
      
      //      fOpen->Close();
      ip++;
    }
  }
  lgr0->Draw();

  cc = new TCanvas(Form("cc%d",icc),Form("cc%d",icc)); icc++;
  mcos1->Draw("ALP");
  lgr1 ->Draw();

  // cc = new TCanvas(Form("cc%d",icc),Form("cc%d",icc)); icc++;
  // mcos2->Draw("ALP");

  
}
