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
Bool_t bsys[]  = { 1, 1, 0, 1};    //132Sn, 108Sn, 124Sn, 112Sn
//-----------
UInt_t  bver[]  = {1, 0, 0, 0};
const UInt_t nmax = (UInt_t)sizeof(bver)/sizeof(UInt_t);
gplot gnames[] = { 
  {".v41.2"  ,"mlt_",""},//"ExB&S.C."},
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
Color_t icol[]   = { kRed, kBlue, kSpring, kMagenta, kOrange, kViolet};
Color_t icolnn[] = { kPink, kBlue+1, kGreen+2, kViolet-1};

TH1I *hhmult;

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
  auto  lgr0 = new TLegend(0.2, 0.7, 0.50, 0.9);
  UInt_t icc = 0;
  TCanvas *cc = new TCanvas(Form("cc%d",icc),Form("cc%d",icc)); icc++;

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

      hhmult    = (TH1I*)fOpen->Get("hmult1");
      hhmult -> SetName(Form("hmult1_%d",is));
      hhmult -> SetTitle(";Multiplicity; Normalized");
      hhmult -> GetYaxis()->SetNdivisions(505);
      hhmult -> SetNormFactor(1);
      hhmult -> SetLineColor(icol[ip]);
      hhmult -> SetMaximum(65000);
      hhmult -> Draw( iopt[ip] );
      lgr0   -> AddEntry( hhmult, label, "lp" );

      ip++;
    }
  }
  lgr0->Draw();
}
