#include "DoFlow.h"
#include "SetStyle.C"

struct gplot{
  TString Version;
  TString fileHeader;
  TString comment;
};


TString  bName[]   = {"132Sn","108Sn","124Sn","112Sn"};
Double_t sysdlt[]  = {0.22,    0.09,      0.15,   0.15};
Double_t sysA[]    = {256.,    220.,      236.,   236};

// --> Plotting selection
//--- Data
Bool_t bsys[]  = { 1, 0, 0, 0};    //132Sn, 108Sn, 124Sn, 112Sn
//-----------
UInt_t  bver[]  = {1, 1, 0, 0};
const UInt_t nmax = (UInt_t)sizeof(bver)/sizeof(UInt_t);
gplot gnames[] = { 
  {".v41.0"  ,"advYPt_",""},//"ExB&S.C."},
  {".v41.0"  ,"advYPt_",""},//"ExB&S.C."},
  {".v41.1"  ,"advYPt_",""},//"ExB&S.C."},
  {".v42.0"  ,"advYPt_","ExB&S.C. NDF>20"},
  {".v40.0"  ,"advYPt_","ExB"},
  {".v38.0"  ,"advYPt_","no corr*"},
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

  TString ltitle;


  TGraphErrors *gv_psi1[5];
  TGraphErrors *gv_psi2[5];
  auto mrpsi1  = new TMultiGraph("mrpsi1"  ,"");
  auto mrpsi2  = new TMultiGraph("mrpsi1"  ,"");
  auto lgr1 = new TLegend(0.20, 0.75, 0.65, 0.9, ""); 
  auto lgr2 = new TLegend(0.20, 0.75, 0.65, 0.9, "");

  TH1I *hmult1[5];
  auto lgm = new TLegend(0.2,0.5,0.5,0.9,"");
  

  
  UInt_t itt = 0;

  for(UInt_t is = 0; is < 4; is++){

      if( !bsys[is] ) continue;

      for(UInt_t it = 0; it < nmax; it++){

	if( !bver[it] ) continue;

	TString fname = "data/psi_" + bName[is] + sVer[it] + ".root";
	if( bmp >= 1 )
	  fname = "data/mlt_" + bName[is] + sVer[it] + ".root";
	// if( it == 0 )
	//   fname = "data/bpsi_"+ bName[is] + sVer[it] + ".root";


	auto fOpen = TFile::Open(fname); 
	
	if( fOpen == NULL ) continue;
	else
	  std::cout << fname << " is opened. " << std::endl;

	if( bmp == 0 ) {
	  gv_psi1[itt] = (TGraphErrors*)fOpen->Get("gv_psi1");
	  gv_psi2[itt] = (TGraphErrors*)fOpen->Get("gv_psi2");
	  
	  cout << " itt " << itt << endl;
	  gv_psi1[itt]->Print();

	}
	else if( bmp == 1 ) {
	  gv_psi1[itt] = (TGraphErrors*)fOpen->Get("gv_mcos1m1");
	  gv_psi2[itt] = (TGraphErrors*)fOpen->Get("gv_mcos2m1");
	}
	else if( bmp == 2 ) {
	  hmult1[itt] = (TH1I*)fOpen->Get("hmult1");
	  hmult1[itt]->SetName(Form("hmult1_%d",itt));
	  if( hmult1[itt] != NULL )
	    lgm->AddEntry( hmult1[itt], bName[is] );
	}

	if( bmp < 2 ) {
	  gv_psi1[itt]->SetMarkerStyle(20);
	  gv_psi1[itt]->SetMarkerColor(icol[itt]);
	  gv_psi1[itt]->SetLineColor(icol[itt]);
	  gv_psi2[itt]->SetMarkerStyle(20);
	  gv_psi2[itt]->SetMarkerColor(icol[itt]);
	  gv_psi2[itt]->SetLineColor(icol[itt]);

	  mrpsi1->Add( gv_psi1[itt] );
	  mrpsi2->Add( gv_psi2[itt] );

	  lgr1->AddEntry( gv_psi1[itt], bName[is] +cmnt[it] );
	  lgr2->AddEntry( gv_psi2[itt], bName[is] +cmnt[it] );
	  
	}
	else if( bmp == 2 ) {
	  hmult1[itt]->SetLineColor(icol[itt]);
	  hmult1[itt]->SetNormFactor(1);
	}

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
    mrpsi1->Draw("ALP");
    lgr1->Draw();

    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    mrpsi2->Draw("ALP");
    lgr2->Draw();
  }
  else if( bmp == 2 ) {
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    hmult1[0]->Draw();
    // for( UInt_t i = 1; i < 5; i++ ) {
    //   if( hmult1[i] != NULL )
    // 	hmult1[i]->Draw("same");
    // }
    // lgm->Draw();
  }

}

