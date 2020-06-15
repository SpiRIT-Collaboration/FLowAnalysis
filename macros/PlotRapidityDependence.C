TCanvas *cc[5];

const UInt_t nsys = 4;
const UInt_t nprt = 5;
 
// --> configuration
Double_t y_cm[nsys]= {0.382453,  0.364873, 0.390302, 0.354066};
Double_t y_bc[nsys]= {0.360199,  0.377779, 0.354065, 0.390301};
Double_t y_bl[nsys]= {0.742652,  0.742652, 0.744367, 0.744367};

TString fsys[nsys] = {"132Sn+124Sn","108Sn+112Sn","124Sn+112Sn","112Sn+124Sn"};
TString rsys[nsys] = {"132",        "108",        "124",        "112"};
TString fpid[nprt] = {"pi-","pi+","proton","deuteron","triton",};  
UInt_t  imrk[nsys] = {20, 21, 22, 23};
Size_t  imsz[nsys] = {1, 1, 1.3, 1.3};
Color_t icol[nprt][nsys]= { {kRed,          kBlue,  kOrange-3,   kGreen+1}, 
			    {kBlue+2,   kOrange+7,  kGreen-3,     kPink+9},
			    {kGreen-3,    kPink+7,  kCyan-1,    kYellow-2}, 
			    {kRed-2,      kBlue-2,  kOrange-2,   kGreen-1},
			    {kBlue-3,   kOrange+5,  kGreen+3,     kPink-9} };
TString iopt[]     = {"","same","same","same","same"};

Double_t ycmoff[nprt][nsys] = {{0.,0.,0.,0}, {0.,0.,0.,0},{0.,0.,0.,0},
			    {0.,0.,0.,0}, {0.,0.,0.,0}};

//Double_t ycmoff[3][nsys] = { {-0.0810542, -0.0686041, -0.0942895, -0.0773481},
//			  {0.0930758,  0.108757,   0.0723492,  0.0813742},
//			  {0.,0.,0.,0}};



void PlotRapidityDependence()
{
  // --> Plotting selection                                                                                                                           
  Bool_t bsys[nsys]  = { 1, 1, 1, 1};  //{"132","108","124","112"};
  Bool_t bpid[nprt]  = { 1, 0, 0, 0, 0};     //{"pi-","pi+","proton","deuteron","triton"};

  Bool_t bCM = kTRUE; // cm frame
  //------------------------------

  //----- Booking
  TFile* fOpen;
  TH1D* hrap[nsys][nprt];
  TH1D* hnpart[nsys][nprt];

  TLegend* lg[nprt];
  for(UInt_t ip = 0; ip < nprt; ip++){
    if( !bpid[ip] ) continue;
    lg[ip] = new TLegend(0.1,0.7,0.35,0.9,"");
    cc[ip] = new TCanvas(Form("cc%d",ip),fpid[ip]);
  }
  //  auto lg2 = new TLegend(0.1,0.7,0.35,0.9,"");


  for(UInt_t is = 0; is < nsys; is++){

    if( !bsys[is] ) continue;

    fOpen = TFile::Open("data/dNdySn"+rsys[is]+".root");
    
    if( fOpen == NULL ) continue;
    else
      std::cout << fOpen->GetName() << " is opened. " << std::endl;
     

    for(UInt_t ip = 0; ip < nprt; ip++){

      if( !bpid[ip] ) continue;

      hrap[is][ip]   = (TH1D*)fOpen->Get(Form("hrap%d_%d",is,ip));
      hnpart[is][ip] = (TH1D*)fOpen->Get(Form("hnpart%d_%d",is,ip));

      cout <<  hrap[is][ip] -> GetName() << endl;

      if(hrap[is][ip] != NULL && hnpart[is][ip] != NULL ){
	cc[ip]->cd();
	hrap[is][ip]->SetLineColor( icol[0][is] );
	//	Int_t    nybin = hnpart[is][ip]->GetNBin();
	//	Double_t dybin =  
	hrap[is][ip]->Scale(1./hnpart[is][ip]->GetEntries());
	hrap[is][ip]->Draw(iopt[is]);

	lg[ip]->AddEntry( hrap[is][ip], fsys[is]+" : "+fpid[ip],"lp");
      }
    }
  }
  
  for(UInt_t ip = 0; ip < nprt; ip++) {
    if( !cc[ip] ) continue;
    cc[ip]->cd();
    lg[ip]->Draw();
  }

}



