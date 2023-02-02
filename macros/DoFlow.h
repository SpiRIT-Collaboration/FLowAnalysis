TString fsys[] = {"^{132}Sn+^{124}Sn","^{108}Sn+^{112}Sn","^{124}Sn+^{112}Sn(Rev)","^{112}Sn+^{124}Sn", "pp", "RPSim"};
TString lsys[] = {"^{132}Sn",        "^{108}Sn",            "^{124}Sn",           "^{112}Sn"     , "pp","^{100}Sn"  };
TString rsys[] = {"132",        "108",        "124",        "112"};
TString tsys[] = {"124",        "112",        "112",        "124"};
TString asys[] = {"Sn132Sn124", "Sn108Sn112", "Sn124Sn112", "Sn112Sn124"};
TString fpid[] = {"proton","deuteron","triton","3He","4He","HHe","neutron","H"};  
TString Fpid[] = {"Proton","Deuteron","Triton","3He","4He","HHe","neutron","H"};  
TString lpid[] = {"1H",    "2H",      "3H"    ,"3He","4He","HHe",  "H"    ,"N"};  

UInt_t  Asys[] = { 132,   108,  124, 112}; 
std::map< Int_t, Int_t > beamAToSys = { {132, 0}, {108, 1}, {124, 2}, {112,3}};

TString  bName[]   = {"132Sn_","108Sn_","124Sn_","112Sn_","100Sn_"};
Double_t sysdlt[]  = {0.22,    0.09,      0.15,   0.15   , 0.22};
Double_t sysA[]    = {256.,    220.,      236.,   236.   ,  256.};
Double_t sysN[]    = {156.,    110.,      136.,   136.   ,  156.};

Double_t u_p       = 0.355151 * 1.06974; 
Double_t mass[]    = {938.2720813, 1875.612762, 2808.921112, 2808.39132, 3727.379378};
Double_t charge[]  = {1.,          1.,          1.,          2.,         2.};
Double_t atmnum[]  = {1.,          1.,          1.,          2.,         2.};

const UInt_t nsys = 4;
const UInt_t nprt = 5;
TString  iopt[]     = {"","same","same","same","same", "same", "same"};
UInt_t   imark[]    = {20, 21, 22, 23, 32, 24, 21, 22, 23};  
Color_t  icol[] = {  kRed, kBlue, kGreen+2, kYellow+1, kMagenta,  kViolet};
Size_t   imsz[]    = {1, 1, 1.3, 1.3, 1.3, 1.3, 1.3};
Color_t  pcolor[]  = {kRed, kBlue, kGreen+2, kYellow+1, kMagenta}; // p,d,t,3He

Color_t  mdlcol[]  = {kRed,  kOrange+1,kCyan-7,  kMagenta-9,kBlue-9, kOrange,kGreen-6 }; //data, AMD(Soft,Stiff), pBUU(Soft,Stiff), ImQMD(Soft,Stiff)
UInt_t   mstyle[]  = {kFullCircle,  kOpenCircle,  kOpenSquare, kOpenTriangleUp}; //data, AMD, pBUU, ImQMD

const UInt_t dTmy = 0;
const UInt_t dKnk = 1;
struct PlotStyle {
  UInt_t  index;
  Color_t fColor;
  UInt_t  mStyle;
  Float_t  mSize;
  TString comment;
};
PlotStyle DStyle[] = {{0,kRed,     kFullCircle,      1, "^{132}Sn"},
		      {0,kBlue,    kFullSquare,      1, "^{108}Sn"},
		      {0,kGreen+3, kOpenTriangleUp,1.5, "^{124}Sn"},
		      {0,kGreen+3, kFullTriangleUp,1.5, "^{112}Sn"}};

PlotStyle CStyle[] = { {0,kRed,      kFullCircle,    1, "DATA"},
		       {1,kBlue,     kFullCircle,    1, "DATA"},
		       {2,kOrange+1 ,kFullCircle,    1, "DATA"},
		       {3,kGreen-3,  kFullCircle,    1, "DATA"},
		       {4,kMagenta-9,kOpenSquare,    1, "pBUU"},
		       {5,kBlue-9,   kOpenSquare,    1, "pBUU"},
		       {6,kViolet-1, kOpenTriangleUp,1, "ImQMD"},
		       {7,kViolet-1, kOpenTriangleUp,1, "ImQMD"},
		       {8,kRed,      kOpenCircle,    1, "DATA_TmyRev"}};

PlotStyle PStyle[] = { {0,kRed,      kFullCircle,    1,    "Proton"},
		       {1,kBlue,     kFullSquare,    1,    "Deuteron"},
		       {2,kOrange+1, kFullTriangleUp,1.5,  "Triton"},
		       {3,kViolet-1, kFullDiamond,   1.5,  "3He"},
		       {4,kGreen-3,  kFullCross,     1.5,  "4He"},
		       {5,kRed,      24,             1,    "Proton"},
		       {6,kBlue,     25,             1,    "Deuteron"},
		       {7,kOrange+1, 26,             1.5,  "Triton"},
		       {8,kViolet-1, 27,             1.5,  "3He"},
		       {9,kGreen-3,  28,             1.5,  "4He"}};

UInt_t ic = 0;
TCanvas *cc;
const Int_t nbinx = 30;
UInt_t id = 0;

//                        0        1       2        3      4       5        6       7       8       9      
UInt_t trkcut[][2] = {{55, 80},{46, 55},{0, 46}, {46,80}};

Double_t yrange1[]    = {-0.5,  -0.25,  -0.20, -0.15,  -0.1,  -0.05, 0.05,  0.1,  0.15,  0.2,  0.25,  0.3,  0.35,  0.45, 0.5}; //by v52.10
Double_t yrange1nrm[] = {-1.5, -0.75, -0.6, -0.45, -0.3, -0.15, -0.07, 0.07, 0.15, 0.3, 0.45, 0.6, 0.75, 1.0,  3.5};
const UInt_t ybin1 = sizeof(yrange1nrm)/sizeof(Double_t);

Double_t yrange2[]    = { -0.5,  -0.25,   -0.1,   0.,  0.1,   0.25 , 0.35,  0.5};
Double_t yrange2nrm[] = { -1.5, -1.0, -0.5, -0.3,  -0.15, 0.15, 0.3, 0.5, 1.0, 3.5};
const UInt_t ybin2 = sizeof(yrange2nrm)/sizeof(Double_t);


Double_t y_bm[]  = { 0.360199, 0.377779, 0.354065, 0.390301, 0.371326, 0.360199};
Double_t y_cm[] =  { 0.382006, 0.36444,  0.389983, 0.353648, 0.371326, 0.382006,   0.37093, 0.37093, 0.37192, 0.371788, 0.371326, 0.37093};
//                    132AA    108AA     124AA     112AA     pp         100(sim)    132nn    108nn    124nn    112nn    pp         100nn
//Kaneko's version -> 132,108,124,112
const double yAA[]   = {0.3822, 0.3647, 0.3902, 0.3538};
const double yNN[]   = {0.3696, 0.3697, 0.3706, 0.3705};
const double yBeam[] = {0.7421, 0.7423, 0.7441, 0.7439};
const double uNN[]   = {0.3709, 0.3709, 0.3719, 0.3718};

// pt bin
const UInt_t   ptbin1  = 8; //16
const UInt_t   ptbin2  = 5; //10
Double_t pt_max = 800.;

const Double_t ut_max[] = {2.2, 1.6, 1.2, 1.2, 1.0, 2.2} ;
Double_t dpt1 = 0.275;
Double_t dpt2 = 0.44;
// Psi bin
const UInt_t  psibin = 12;
				   
TString  Partname[] = {"pi-","pi+","Proton","Deuteron"  ,"Triton",      "3He",      "4He", "Neutron",        "H", "HHe"};
TString  partname[] = {"pi-","pi+","proton","deuteron"  ,"triton",      "3He",      "4He", "neutron",        "H", "HHe"};
UInt_t   partid[]   = {211,  211,    2212, 1000010020, 1000010030, 1000020030, 1000020040,      2112, 1000010040,  2212};
Double_t partA[]    = {0.15, 0.15,     1.,         2.,         3.,         3.,         4.,        1.,         1.,    1.};
Double_t partZ[]    = { -1.,   1.,     1.,         1.,         1.,         2.,         2.,        0.,         1.,    1.};


struct nucleus{
  UInt_t index;
  UInt_t pdg;
  Double_t A;
  Double_t Z;
  Double_t mass;
  TString name;
  TString Name;
  TString sName;
};

std::vector< nucleus > ncls = {
  {0,      2212, 1., 1.,   938.272081, "proton"  ,"Proton",  "1H"},
  {1,1000010020, 2., 1.,  1875.612762, "deuteron","Deuteron","2H"},
  {2,1000010030, 3., 1.,  2808.921112, "triton",  "Triton",  "3H"},
  {3,1000020030, 3., 2.,  2808.39132,  "3He",     "3He",     "3He"}, 
  {4,1000020020, 4., 2.,  3727.379378, "4He",     "4He",     "4He"}, 
  {6,      2212, 1., 1.,   938.272081, "HHe",     "HHe",     "HHe"}, 
  {5,      2212, 1., 1.,   938.272081, "H",       "H",       "H"}, 
  {7,      2112, 1., 0.,   939.5654,   "neutron", "Neutron", "n"}
};

std::map< TString, UInt_t > pToIDX = {{"1H",0},{"2H",1},{"3H",2},{"3He",3},{"4He",4},{"H",5},{"HHe",6},{"n",7}};
std::map< TString, UInt_t > ptoPDG  = {{"1H",2212},{"1H",1000010020},{"3H",1000010030},{"3He",1000020030},
				       {"4He",1000020040},{"H",2212},{"HHe",2212},{"n",2112}};


//  UInt_t mrange[] = {70, 60, 55, 50, 45, 40, 35, 30, 25, 20, 0}; //mtrack4
UInt_t mrange[] = {800, 100, 90, 85, 80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5, 0};
const UInt_t mbin = sizeof(mrange)/sizeof(UInt_t);

//  Int_t  Seltrk = 2;
Int_t  Lcent  = 0;  // {75, 45, 30, 20, 0};
Int_t  Ucent  = 80;
Int_t  RPBase;

const UInt_t npsi = 12; // Number of bin in Psi
Double_t dphi = 2.*TMath::Pi()/ (Double_t)npsi ;

TString amdpartname[] = {"prt","deut","trit","3He","4He","neut","H"};

Bool_t bCorr = kFALSE;

Double_t y_norm = 1.;

Float_t  v1fit[2] = {-0.5, 0.6};
Float_t  v1para[5] = {0.35, 0.4, 0.55, 0.55, 0.8};


Double_t v2para0[5][5]={{ -0.001,0.05, 0.0, -0.6, 0.6},
			{ -0.05, 0.05, 0.0, -0.6, 0.6},
			{ -0.07, 0.05, 0.0, -0.6, 0.6},
			{ -0.07, 0.05, 0.0, -0.6, 0.8},
			{ -0.09, 0.04, 0.0, -0.6, 0.6}};

TString mcphicut = "womc";


TFile* GetAcceptanceCorrectionFile(UInt_t mtconfig=0, UInt_t phicutID=5, TString comment="" )
{
  //  TString phiname[] = { "45or135", "45", "135", "45to135","all","-20to30"};
  TString phiname[] = { "45or135", "45", "135", "45to135","all","-30to20"};

  //  TString trklabel = Form("_%dto%02d_vapri",trkcut[seltrk][0],trkcut[seltrk][1]);  
  TString trklabel = Form("_%dto%02d_reco",Lcent, Ucent);
  
  //  TString fname =  "data/rootfiles/UnfoldedLCPSpectra_"+ phiname[phicutID]+ "_" + comment + trklabel+".slim.root";
  //  TString fname =  "data/rootfiles/UnfoldedLCPSpectra"+ trklabel+".slim.root";
  //  TString fname =  "data/rootfiles/UnfoldedLCPSpectra"+ trklabel+".slim.hpc.root";

  //  TString fname =  "data/rootfiles/UnfoldedLCPSpectra_"+ phiname[phicutID]+ "_" + comment + trklabel+".slim.root";

  //TString fname = "data/rootfiles/UnfoldedLCPSpectra_-30to20_wmc_55to80_reco.slim.root";   // Kaneko's high stat.
  //  TString fname = "data/rootfiles/UnfoldedLCPSpectra_-30to20_wmc_55to80_reco.slim.m.root"; // low stat.   
  //  TString fname = "data/rootfiles/UnfoldedLCPSpectra_-30to20_20220101_55to80_reco.slim.root"; // 20220101 version

  TString fname;
  if( mtconfig == 1 )
    fname = "data/rootfiles/UnfoldedLCPSpectra_-30to20_46to55_left1rad_vapri.slim.root";
  else
    fname = "data/rootfiles/UnfoldedLCPSpectra_-30to20_55to80_left1rad_vapri.slim.root";

  TFile* fopen = TFile::Open(fname);
  if( fopen ) 
    LOG(INFO) << fname << " is opened. " << FairLogger::endl;
  else 
    LOG(ERROR) << fname << " is not found " << FairLogger::endl;
  
  return fopen;
}

TH2D *GetAcpCorrection(TH2D &h2, TH2D *hacpcorr)
{
  TH2D * h2c = new TH2D(h2);
  //  h2c->Sumw2();

  if( hacpcorr != NULL )
    h2c->Divide(hacpcorr);
  else {
    LOG(ERROR) << " Acceptance correction is failed." << FairLogger::endl;
    return NULL;
  }
  
  return h2c;
}


void SetXTextLabel(TMultiGraph& obj)
{
  std::vector<TString> slabel={"1H","2H","3H","3He","4He"};
  obj.GetXaxis()->SetLimits(0.5,5.5);

  for(auto i : ROOT::TSeqI(slabel.size()) )
    obj.GetXaxis()->SetBinLabel( obj.GetXaxis()->FindBin( i+1 ), slabel.at(i));

}

void SetXTextLabel(THStack& obj)
{
  std::vector<TString> slabel={"1H","2H","3H","3He","4He"};
  //  obj.GetXaxis()->SetRangeUser(0.5,5.5);

  for(auto i : ROOT::TSeqI(slabel.size()) )
    obj.GetXaxis()->SetBinLabel( obj.GetXaxis()->FindBin( i+1 ), slabel.at(i));
}

