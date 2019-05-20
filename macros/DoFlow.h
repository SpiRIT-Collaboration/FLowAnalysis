
const UInt_t nsys = 4;
const UInt_t nprt = 5;
TString  iopt[]     = {"","same","same","same","same", "same"};
UInt_t   imark[]    = {20, 21, 22, 23, 25, 26, 24};  

UInt_t ic = -1;
TCanvas *cc;
const Int_t nbinx = 30;
UInt_t id = 0;

Double_t yrange1[] = { -0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.5};
const UInt_t ybin1 = sizeof(yrange1)/sizeof(Double_t);

Double_t yrange2[] = {-0.2, -0.05,  0.05, 0.2, 0.35, 0.5};
const UInt_t ybin2 = sizeof(yrange2)/sizeof(Double_t);

TString  partname[] = {"pi-","pi+","proton","deuteron" ,"triton", "3He", "4He", "neutron", "H"};
UInt_t   partid[]   = {211,    211,    2212, 1000010020, 1000010030, 1000020030, 1000020040, 2112, 1000010040};

//  UInt_t mrange[] = {70, 60, 55, 50, 45, 40, 35, 30, 25, 20, 0}; //mtrack4
UInt_t mrange[] = {95, 90, 85, 80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5, 0};

UInt_t   cent[]  = {75, 35, 28, 0};

TString amdpartname[] = {"prt","deut","trit","3He","4He","neut","H"};
