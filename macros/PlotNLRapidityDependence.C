//----------------------------------------------------------------------
// ROOT macro: PlotRapidityDependence()
// (c) Mizuki Nishimura (RIKEN)
//
// This macro plot v1 and v2 as a function of rapidity for 4 systems.
// In advance, each root file which is data/NL###Sn.root must be 
// produced using calcFlw.C.
//----------------------------------------------------------------------

TCanvas *cc0;
TCanvas *cc1;

const UInt_t nsys = 4;  // number of system 132Sn, 108Sn, 124Sn, 112Sn
const UInt_t nprt = 2;  // number of particle species
 
// --> configuration
Double_t y_cm[nsys]= {0.382453,  0.364873, 0.390302, 0.354066};
Double_t y_bc[nsys]= {0.360199,  0.377779, 0.354065, 0.390301};
Double_t y_bl[nsys]= {0.742652,  0.742652, 0.744367, 0.744367};

TString fsys[nsys] = {"132Sn+124Sn","108Sn+112Sn","124Sn+112Sn","112Sn+124Sn"};
TString rsys[nsys] = {"132",        "108",        "124",        "112"};
TString fpid[nprt] = {"neutron","proton"};  
UInt_t  imrk[nsys] = {20, 21, 22, 23};
Size_t  imsz[nsys] = {1, 1};
Color_t icol[nprt][nsys]= { {kRed,          kBlue,  kOrange-3,   kGreen+1}, 
			    {kBlue+2,   kOrange+7,  kGreen-3,     kPink+9} };

Double_t ycmoff[nprt][nsys] = {{0.,0.,0.,0}, {0.,0.,0.,0}};

//Double_t ycmoff[3][nsys] = { {-0.0810542, -0.0686041, -0.0942895, -0.0773481},
//			  {0.0930758,  0.108757,   0.0723492,  0.0813742},
//			  {0.,0.,0.,0}};


void Save(TString cnt = "")
{
  cc0->SaveAs("NLv1RpSn"+cnt+".png");
  cc1->SaveAs("NLv2RpSn"+cnt+".png");
}


void PlotNLRapidityDependence()
{

  // --> Plotting selection                                                                                                                           
  Bool_t bsys[nsys]  = { 1, 1, 1, 1};  //{"132","108","124","112"};
  Bool_t bpid[nprt]  = { 1, 0};     //{"neutron", "proton"}

  Bool_t bCM = kTRUE; // cm frame
  //Bool_t bCM = kFALSE; // lab frame

  //------------------------------

  TString fname[nsys];
  for(UInt_t is = 0; is < nsys; is++){
    fname[is] = "data/NLv1v2"+rsys[is]+"Sn.root";
  }
  
  TFile *fOpen;
  TGraphErrors *grn_v1[12];
  TGraphErrors *grn_v2[12];

  auto mv1 = new TMultiGraph("mv1",";Rapidity_lab; v1");
  auto mv2 = new TMultiGraph("mv2",";Rapidity_lab; v2");
  auto lg1 = new TLegend(0.125,0.7,0.35,0.9,"v1 ");
  auto lg2 = new TLegend(0.125,0.7,0.35,0.9,"v2 ");


  UInt_t igr = 0;

  for(UInt_t is = 0; is < nsys; is++){

    if( !bsys[is] ) continue;
    
    fOpen = TFile::Open(fname[is]);

    if( fOpen == NULL ) continue;
    else
      std::cout << fname[is] << " is opened. " << std::endl;



    for(UInt_t ip = 0; ip < nprt; ip++){

      if( !bpid[ip] ) continue;


      // Getting plots from the file.
      TGraphErrors *gvn1 = (TGraphErrors*)fOpen->Get("gvn_v1");
      TGraphErrors *gvn2 = (TGraphErrors*)fOpen->Get("gvn_v2");
      TGraphErrors *gvp1 = (TGraphErrors*)fOpen->Get("gvp_v1");
      TGraphErrors *gvp2 = (TGraphErrors*)fOpen->Get("gvp_v2");


      // **
      // Plot in Rapidity_cm
      if( bCM ) {

	// -----> v1
	mv1->SetTitle("; y_cm/y_beam; v1");

	const UInt_t npoint1 = gvn1->GetN(); 
	Double_t *xp1 = new Double_t[npoint1];
	Double_t *xe1 = new Double_t[npoint1];
	Double_t *yp1 = new Double_t[npoint1];
	Double_t *ye1 = new Double_t[npoint1];

	xp1 = gvn1->GetX();
	xe1 = gvn1->GetEX();
	yp1 = gvn1->GetY();
	ye1 = gvn1->GetEY();

	for(UInt_t ik = 0; ik < npoint1; ik++)
	  xp1[ik] = (xp1[ik] - y_cm[is])/y_bc[is] - ycmoff[ip][is];

	grn_v1[igr] = new TGraphErrors(npoint1, xp1, yp1, xe1, ye1);

	// -----> v2
	mv2->SetTitle("; y_cm/y_beam; v2");

	const UInt_t npoint2 = gvn2->GetN(); 
	Double_t *xp2 = new Double_t[npoint2];
	Double_t *xe2 = new Double_t[npoint2];
	Double_t *yp2 = new Double_t[npoint2];
	Double_t *ye2 = new Double_t[npoint2];

	xp2 = gvn2->GetX();
	xe2 = gvn2->GetEX();
	yp2 = gvn2->GetY();
	ye2 = gvn2->GetEY();

	for(UInt_t ik = 0; ik < npoint2; ik++)
	  xp2[ik] = (xp2[ik] - y_cm[is])/y_bc[is] - ycmoff[ip][is];

	grn_v2[igr] = new TGraphErrors(npoint2, xp2, yp2, xe2, ye2);

	delete[] xp1, xe1, yp1, ye1;
	delete[] xp2, xe2, yp2, ye2;
      }

      else {	      
	TString gvname = Form("grn_v1%d",igr);
	grn_v1[igr] = (TGraphErrors*)gvn1->Clone(gvname);
	gvname = Form("grn_v2%d",igr);
	grn_v2[igr] = (TGraphErrors*)gvn2->Clone(gvname);
      }

      mv1->Add(grn_v1[igr],"lp");
      mv2->Add(grn_v2[igr],"lp");

      grn_v1[igr]->SetMarkerStyle(imrk[is]);
      grn_v1[igr]->SetMarkerColor(icol[ip][is]);
      grn_v1[igr]->SetMarkerSize(imsz[is]);
      grn_v1[igr]->SetLineColor(icol[ip][is]);
      lg1->AddEntry(grn_v1[igr], rsys[is]+" "+fpid[ip] ,"lp");


      grn_v2[igr]->SetMarkerStyle(imrk[is]);
      grn_v2[igr]->SetMarkerColor(icol[ip][is]);
      grn_v2[igr]->SetMarkerSize(imsz[is]);
      grn_v2[igr]->SetLineColor(icol[ip][is]);
      lg2->AddEntry(grn_v2[igr], rsys[is]+" "+fpid[ip],"lp");

      
      fOpen->Close();
      igr++;
    }
  }


  

  cc1 = new TCanvas("cc1","v2");
  //  mv2->SetMaximum(0.1);
  mv2->Draw("ALP");
  auto aLineX2 = new TLine(mv2->GetXaxis()->GetXmin(), 0., mv2->GetXaxis()->GetXmax(), 0.);
  auto aLineY2 = new TLine(0., mv2->GetYaxis()->GetXmin(), 0., mv2->GetYaxis()->GetXmax());
  aLineX2->SetLineStyle(3);
  aLineY2->SetLineStyle(3);
  aLineX2->Draw();
  aLineY2->Draw();
  lg2->Draw();

  cc0 = new TCanvas("cc0","v1");
  mv1->Draw("ALP");
  auto aLineX1 = new TLine(mv1->GetXaxis()->GetXmin(), 0., mv1->GetXaxis()->GetXmax(), 0.);
  auto aLineY1 = new TLine(0., mv1->GetYaxis()->GetXmin(), 0., mv1->GetYaxis()->GetXmax());
  //auto aLineY1 = new TLine(y_cm[0], mv1->GetYaxis()->GetXmin(), y_cm[0], mv1->GetYaxis()->GetXmax());
  aLineX1->SetLineStyle(3);
  aLineY1->SetLineStyle(3);
  aLineX1->Draw();
  aLineY1->Draw();
  lg1->Draw();
}



