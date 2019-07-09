
TString fsys[] = {"^{132}Sn+^{124}Sn","^{108}Sn+^{112}Sn","^{124}Sn+^{112}Sn","^{112}Sn+^{124}Sn"};
TString lsys[] = {"^{132}Sn",        "^{108}Sn",            "^{124}Sn",           "^{112}Sn"};
TString rsys[] = {"132",        "108",        "124",        "112"};
TString tsys[] = {"124",        "112",        "112",        "124"};
TString fpid[] = {"proton","deuteron","triton","3He","4He","neutron","H"};  

const UInt_t nsys = 4;
const UInt_t nprt = 5;
TString  iopt[]     = {"","same","same","same","same", "same"};
UInt_t   imark[]    = {20, 21, 22, 23, 25, 26, 24};  

UInt_t ic = -1;
TCanvas *cc;
const Int_t nbinx = 30;
UInt_t id = 0;


Double_t yrange1[] = { -0.8, -0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.5, 0.8};
const UInt_t ybin1 = sizeof(yrange1)/sizeof(Double_t);
Double_t yrange2[] = { -0.8, -0.2, -0.05,  0.05, 0.2, 0.35, 0.5, 0.8};
const UInt_t ybin2 = sizeof(yrange2)/sizeof(Double_t);

Double_t y_cm[]  = { 0.382453, 0.364873, 0.390302, 0.354066};
Double_t y_bm[]  = { 0.360199, 0.377779, 0.354065, 0.390301};

// pt bin
const UInt_t   ptbin1  = 10; //16
const UInt_t   ptbin2  = 8; //10
Double_t pt_max = 800.;
Double_t dpt1;
Double_t dpt2;

				   
TString  partname[] = {"pi-","pi+","proton","deuteron" ,"triton", "3He", "4He", "neutron", "H",   "triton"};
UInt_t   partid[]   = {211,    211,    2212, 1000010020, 1000010030, 1000020030, 1000020040, 2112, 1000010040, 1000010030};
Double_t partA[]    = {1.,      1.,      1.,         2.,         3.,         3.,         4.,   1.,        1.,   3.};

//  UInt_t mrange[] = {70, 60, 55, 50, 45, 40, 35, 30, 25, 20, 0}; //mtrack4
UInt_t mrange[] = {95, 90, 85, 80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5, 0};
const UInt_t mbin = sizeof(mrange)/sizeof(UInt_t);

Int_t  Lcent = 0;  // {75, 45, 30, 20, 0};
Int_t  Ucent = 80;


TString amdpartname[] = {"prt","deut","trit","3He","4He","neut","H"};

