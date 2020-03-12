#include "DoFlow.h"
#include "SetStyle.C"

struct gplot{
  TString Version;
  TString fileHeader;
  TString comment;
  TString label;
};

TString  bName[]   = {"132Sn", "108Sn", "124Sn", "112Sn", "100Sn"};
Double_t sysdlt[]  = {0.22,    0.09,      0.15,   0.15   , 0.22};
Double_t sysA[]    = {256.,    220.,      236.,   236.   ,  256.};


// --> Plotting selection
//--- Data
Bool_t bsys[]  = { 0, 0, 0, 0, 1};    //132Sn, 108Sn, 124Sn, 112Sn, 100Sn
//-----------
UInt_t  bver[]  = {1, 1, 1};
const UInt_t nmax = (UInt_t)sizeof(bver)/sizeof(UInt_t);
gplot gnames[] = { 
  {".v50.7"    ,"cpsi_", ".m00to80",""},
  {".v50.7"    ,"epsi_", ".m00to80",""},
  {".v50.7"    ,"dpsi_",".m00to80",""},
  {".v50.7"    ,"ccpsi_",".m00to80",""},
  {".v41.2"    ,"bpsi_",".m00to40",""},
  {".v41.2"    ,"bpsi_",".m00to60",""},
  {".v41.0"    ,"bpsi_","",""},
};

TString sVer[nmax];
TString sName[nmax];
TString cmnt[nmax];

//==================================================

// --> configuration
Size_t  imsz[] = {1, 1, 1.3, 1.3, 1.3, 1.3, 1.3};
Color_t icol[] = {  kRed, kBlue, kSpring, kMagenta, kOrange, kViolet};
Color_t icolnn[]={ kPink, kBlue+1, kGreen+2, kViolet-1};


void PlotPsiCorrection(UInt_t bmp = 0) // 0: phi, 1:mlt
{
  
  gStyle->SetOptStat(0);
  SetStyle();

  for(UInt_t i = 0; i < nmax; i++){
    if( bver[i] ) {
      sVer[i]  = gnames[i].Version;
      sName[i] = gnames[i].fileHeader;
      cmnt[i]  = gnames[i].comment;
    }
  }


  TGraphErrors *gv_psi1[5];
  TGraphErrors *gv_psi2[5];
  auto mrpsi1  = new TMultiGraph("mrpsi1"  ,"");
  auto mrpsi2  = new TMultiGraph("mrpsi2"  ,"");
  auto lgr1 = new TLegend(0.20, 0.75, 0.65, 0.9, ""); 
  auto lgr2 = new TLegend(0.20, 0.75, 0.65, 0.9, "");

  auto lgm = new TLegend(0.2,0.5,0.5,0.9,"");
  

  
  UInt_t itt = 0;

  for(UInt_t is = 0; is < 5; is++){

    if( !bsys[is] ) continue;

    for(UInt_t it = 0; it < nmax; it++){

      if( !bver[it] ) continue;

      TString fname = "data/"+ sName[it] + bName[is] + sVer[it] + cmnt[it] + ".root";
      if( bmp == 1 )
	fname = "data/mlt_" + bName[is] + sVer[it] + ".root";
      
      
      TString ltitle = sVer[it] + ":" + sName[it];

      auto fOpen = TFile::Open(fname); 
	
      if( fOpen == NULL ) continue;
      else
	std::cout << fname << " is opened. " << std::endl;

      if( bmp == 0 ) {
	gv_psi1[itt] = (TGraphErrors*)fOpen->Get("gv_psi1");
	gv_psi2[itt] = (TGraphErrors*)fOpen->Get("gv_psi2");
	gv_psi2[itt]->Print();

	auto mlt = (TH1I*)fOpen->Get("hmult");

	cout << " itt " << itt << " --> # of events = " << mlt->GetEntries() << endl;

      }
      else {
	gv_psi1[itt] = (TGraphErrors*)fOpen->Get("gv_mcos1m1");
	gv_psi2[itt] = (TGraphErrors*)fOpen->Get("gv_mcos2m1");
      }

      gv_psi1[itt]->SetMarkerStyle(20);
      gv_psi1[itt]->SetMarkerColor(icol[itt]);
      gv_psi1[itt]->SetLineColor(icol[itt]);
      gv_psi2[itt]->SetMarkerStyle(20);
      gv_psi2[itt]->SetMarkerColor(icol[itt]);
      gv_psi2[itt]->SetLineColor(icol[itt]);

      mrpsi1->Add( gv_psi1[itt] );
      mrpsi2->Add( gv_psi2[itt] );

      lgr1->AddEntry( gv_psi1[itt], ltitle);
      lgr2->AddEntry( gv_psi2[itt], ltitle);
	  
      itt++;

      fOpen->Close();
    }
  }
  
  if( bmp < 2) {
    mrpsi1->SetTitle( gv_psi1[0]->GetTitle() );
    mrpsi2->SetTitle( gv_psi2[0]->GetTitle() );
    if( bmp == 0 )
      mrpsi2->SetTitle("; #Psi;<cos(2#Delta #Psi)>");
  
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    if( bmp == 0 )
      //      mrpsi1->GetYaxis()->SetRangeUser(0.62,0.76);
    mrpsi1->Draw("ALP");
    lgr1->Draw();

    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    if( bmp == 0)
      //      mrpsi2->GetYaxis()->SetRangeUser(0.26,0.45);
    mrpsi2->Draw("ALP");
    lgr2->Draw();
  }


}

