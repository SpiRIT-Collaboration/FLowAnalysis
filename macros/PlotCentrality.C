#include "DoFlow.h"
#include "SetStyle.C"

//==================================================
//-- plot configuration
//--------------------------------------------------
  // --> Plotting selection
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
UInt_t  bver[]  = {1, 0, 0, 0};
const UInt_t nmax = (UInt_t)sizeof(bver)/sizeof(UInt_t);
gplot gnames[] = { 
  {".v52.10" ,"mlt_",""},
  {".v40.0"  ,"mlt_",""},//"ExB&S.C."},
  {".v29.1"  ,"mlt_",""},//"ExB&S.C."},
  {".v41.1"  ,"mlt_",""},//"ExB&S.C."},
  {".v42.0"  ,"mlt_","ExB&S.C. NDF>20"},
  {".v40.0"  ,"mlt_","ExB"},
  {".v38.0"  ,"mlt_","no corr*"},
};

TString sVer[nmax];
TString sName[nmax];
TString cmnt[nmax];

//==================================================
// --> configuration

Size_t  imsz[]   = {1, 1, 1.3, 1.3, 1.3, 1.3, 1.3};
//Color_t icol[]   = { kRed, kBlue, kSpring, kMagenta, kOrange, kViolet};
Color_t icolnn[] = { kPink, kBlue+1, kGreen+2, kViolet-1};

TH1I *hhmult;
TH1I *hhmult1;
TH1I *hhmult2;

void PlotCentrality()
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

  TFile *fOpen;
  auto  lgr1 = new TLegend(0.7, 0.7, 0.95, 0.9);
  auto  lgr2 = new TLegend(0.7, 0.7, 0.95, 0.9);
  auto  lgr4 = new TLegend(0.6, 0.7, 0.90, 0.9);
 
  TCanvas *cc1 = new TCanvas("cc1","cc1"); 
  TCanvas *cc2 = new TCanvas("cc2","cc2"); 
  TCanvas *cc4 = new TCanvas("cc4","cc4"); 

  UInt_t ip = 0;
  for( UInt_t is = 0; is < 4; is++ ){
    if( !bsys[is] ) continue;

    for(UInt_t iz = 0; iz < (UInt_t)sizeof(bver)/sizeof(UInt_t); iz++){
      if( !bver[iz] ) continue;

      TString label = fsys[is];
      TString fname = "data/"+ sName[iz] + bName[is] + sVer[iz]+ ".root";

      fOpen = TFile::Open(fname);
    
      if(fOpen == NULL ) continue;

      std::cout << fname << " is opened. " << std::endl;

      hhmult    = (TH1I*)fOpen->Get("hmult");
      hhmult -> SetName(Form("hmult_%d",is));
      hhmult -> SetTitle("mtrack4;Multiplicity; Normalized");
      hhmult -> GetYaxis()->SetNdivisions(505);
      hhmult -> SetNormFactor(1);
      hhmult -> SetLineColor(icol[ip]);
      hhmult -> SetMaximum(85000);
      lgr4   -> AddEntry( hhmult, label, "lp" );
      cc4->cd();
      hhmult -> Draw( iopt[ip] );

      hhmult1 = (TH1I*)fOpen->Get("hmult1");
      hhmult1 -> SetName(Form("hmult1_%d",is));
      hhmult1 -> SetTitle("mtrack1;Multiplicity; Normalized");
      hhmult1 -> GetYaxis()->SetNdivisions(505);
      hhmult1 -> SetNormFactor(1);
      hhmult1 -> SetLineColor(icol[ip]);
      hhmult1 -> SetMaximum(65000);
      lgr1    -> AddEntry( hhmult1, label, "lp" );
      cc1->cd();
      hhmult1 -> Draw( iopt[ip] );

      hhmult2 =  (TH1I*)fOpen->Get("hmult2");
      hhmult2 -> SetName(Form("hmult2_%d",is));
      hhmult2 -> SetTitle("mtrack2;Multiplicity; Normalized");
      hhmult2 -> GetYaxis()->SetNdivisions(505);
      hhmult2 -> SetNormFactor(1);
      hhmult2 -> SetLineColor(icol[ip]);
      hhmult2 -> SetMaximum(70000);
      lgr2    -> AddEntry( hhmult, label, "lp" );
      cc2->cd();
      hhmult2 -> Draw( iopt[ip] );

      ip++;
    }
  }
  cc1->cd();
  lgr1->Draw();
  cc2->cd();
  lgr2->Draw();
  cc4->cd();
  lgr4->Draw();
}
