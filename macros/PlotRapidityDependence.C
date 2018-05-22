TCanvas *cc0;
TCanvas *cc1;

// --> configuration
Double_t y_cm[4]= {0.382453,  0.364873, 0.390302, 0.354066};
Double_t y_bc[4]= {0.360199,  0.377779, 0.354065, 0.390301};
Double_t y_bl[4]= {0.742652,  0.742652, 0.744367, 0.744367};

TString fsys[4] = {"132Sn+124Sn","108Sn+112Sn","124Sn+112Sn","112Sn+124Sn"};
TString rsys[4] = {"132",        "108",        "124",        "112"};
TString fpid[3] = {"proton","deuteron","triton"};  
UInt_t  imrk[4] = {20, 21, 22, 23};
Size_t  imsz[4] = {1, 1, 1.3, 1.3};
Color_t icol[3][4]= { {kRed,          kBlue,  kOrange-3,   kGreen+1}, 
		      {kBlue+2,   kOrange+7,  kGreen-3,     kPink+9},
		      {kGreen-3,    kPink+7,  kCyan-1,    kYellow-2} };


Double_t ycmoff[3][4] = {{0.,0.,0.,0}, {0.,0.,0.,0},{0.,0.,0.,0}};

//Double_t ycmoff[3][4] = { {-0.0810542, -0.0686041, -0.0942895, -0.0773481},
//			  {0.0930758,  0.108757,   0.0723492,  0.0813742},
//			  {0.,0.,0.,0}};



void Save(TString cnt = "")
{
  cc0->SaveAs("v1RpSn"+cnt+".png");
  cc1->SaveAs("v2RpSn"+cnt+".png");
}


void PlotRapidityDependence()
{

  // --> Plotting selection                                                                                                                           
  Bool_t bsys[4]  = { 0, 1, 0, 0};  //{"132","108","124","112"};
  Bool_t bpid[3]  = { 1, 1, 1};     //{"proton","deuteron","triton"};

  Bool_t bCM = kTRUE; // cm frame
  //------------------------------

  TString fname[4][3];
  for(UInt_t is = 0; is < 4; is++){
    for(UInt_t ip = 0; ip < 3; ip++){
      fname[is][ip] = "data/VN"+rsys[is]+"Sn_"+fpid[ip]+".root";
    }
  }
  
  TFile *fOpen;
  TGraphErrors *gr_v1[12];
  TGraphErrors *gr_v2[12];

  auto mv1 = new TMultiGraph("mv1",";Rapidity_lab; v1");
  auto mv2 = new TMultiGraph("mv2",";Rapidity_lab; v2");
  auto lg1 = new TLegend(0.1,0.7,0.35,0.9,"v1 ");
  auto lg2 = new TLegend(0.1,0.7,0.35,0.9,"v2 ");


  UInt_t igr = 0;

  for(UInt_t is = 0; is < 4; is++){

    if( !bsys[is] ) continue;

    for(UInt_t ip = 0; ip < 3; ip++){

      if( !bpid[ip] ) continue;

      fOpen = TFile::Open(fname[is][ip]);

      if( fOpen == NULL ) continue;
      else
	std::cout << fname[is][ip] << " is opened. " << std::endl;


      TGraphErrors *gv1 = (TGraphErrors*)fOpen->Get("gv_v1");
      TGraphErrors *gv2 = (TGraphErrors*)fOpen->Get("gv_v2");

      gr_v1[igr] = new TGraphErrors();
      gr_v2[igr] = new TGraphErrors();


      if( bCM ) {
	// mv1->SetTitle("; y_cm/y_beam; v1");
	// mv2->SetTitle("; y_cm/y_beam; v2");
	mv1->SetTitle("; y_lab; v1");
	mv2->SetTitle("; y_lab; v2");
	
	const UInt_t npoint = gv1->GetN(); 
	Double_t xp;
	Double_t xe;
	Double_t yp;
	Double_t ye;

	std::vector<Double_t> xpp;
	std::vector<Double_t> ypp;

      
	for(UInt_t ik = 0; ik < npoint; ik++){
	  gv1->GetPoint(ik, xp, yp);
	  xe = gv1->GetErrorX(ik);
	  ye = gv1->GetErrorY(ik);
	
	  //xp = (xp - y_cm[is])/y_bc[is] - ycmoff[ip][is];
	  
	  //xp = (xp - y_bc[is]); ///y_bl[is] ;

	  gr_v1[igr]->SetPoint(ik, xp, yp);
	  gr_v1[igr]->SetPointError(ik, xe, ye);
	  xpp.push_back(xp);
	  ypp.push_back(yp);


	  gv2->GetPoint(ik, xp, yp);
	  xe = gv2->GetErrorX(ik);
          ye = gv2->GetErrorY(ik);

	  xp = (xp - y_cm[is])/y_bc[is];
          gr_v2[igr]->SetPoint(ik, xp, yp);
          gr_v2[igr]->SetPointError(ik, xe, ye);
	}

	//	auto inter = new ROOT::Math::Interpolator(ypp,xpp,ROOT::Math::Interpolation::kCSPLINE);
	//	std::cout << rsys[is] << " : " << fpid[ip] << " " << inter->Eval(0) << std::endl;
      }
      else {
	      
	TString gvname = Form("gPt_v1%d",igr);
	gr_v1[igr] = (TGraphErrors*)gv1->Clone(gvname);
	gvname = Form("gPt_v2%d",igr);
	gr_v2[igr] = (TGraphErrors*)gv2->Clone(gvname);
      }

      mv1->Add(gr_v1[igr],"lp");
      mv2->Add(gr_v2[igr],"lp");

      gr_v1[igr]->SetMarkerStyle(imrk[is]);
      gr_v1[igr]->SetMarkerColor(icol[ip][is]);
      gr_v1[igr]->SetMarkerSize(imsz[is]);
      gr_v1[igr]->SetLineColor(icol[ip][is]);
      lg1->AddEntry(gr_v1[igr], rsys[is]+" "+fpid[ip] ,"lp");


      gr_v2[igr]->SetMarkerStyle(imrk[is]);
      gr_v2[igr]->SetMarkerColor(icol[ip][is]);
      gr_v2[igr]->SetMarkerSize(imsz[is]);
      gr_v2[igr]->SetLineColor(icol[ip][is]);
      lg2->AddEntry(gr_v2[igr], rsys[is]+" "+fpid[ip],"lp");

      
      fOpen->Close();
      igr++;
    }
  }


  

  cc1 = new TCanvas("cc1","v2");
  mv2->SetMaximum(0);
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
  //  auto aLineY1 = new TLine(0., mv1->GetYaxis()->GetXmin(), 0., mv1->GetYaxis()->GetXmax());
  auto aLineY1 = new TLine(y_cm[0], mv1->GetYaxis()->GetXmin(), y_cm[0], mv1->GetYaxis()->GetXmax());
  aLineX1->SetLineStyle(3);
  aLineY1->SetLineStyle(3);
  aLineX1->Draw();
  aLineY1->Draw();
  lg1->Draw();
}



