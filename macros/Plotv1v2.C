void Plotv1v2()
{
  Bool_t TPC     = 0; //kTRUE; //kTRUE;
  Bool_t NeuLAND = 1; // kTRUE;
  Bool_t AMD     = 0;
  Bool_t bsys[]  = { 1, 0, 0, 0};   //{"132","108","124","112"};
  Bool_t bpid[]  = { 1, 0, 0, 1};   //{"proton","deuteron","triton","neutron"};                           
  Bool_t bamd[]  = { 1, 1, 1, 1};   // amd configulation


  Bool_t bYCM = 0; //kFALSE; //kTRUE;

  Double_t y_cm[]= {0.382453,  0.364873, 0.390302, 0.354066, 0.00, 0.371326};
  Double_t y_bc[]= {0.360199,  0.377779, 0.354065, 0.390301, 0.01, 0.371326};
  Double_t y_bl[]= {0.742652,  0.742652, 0.744367, 0.744367, 0.01, 0.371326};


  TCanvas *cc[2];

  std::vector< TString > label;
  std::vector< UInt_t > isys;

  TString expbeamconfig[] = { "132Sn",
			      "108Sn",
			      "124Sn",
			      "112Sn"};

  TString amdbeamconfig[] = {"132Sn124Sn270AMeV",
			     "108Sn112Sn270AMeV",
			     "124Sn112Sn270AMeV",
			     "112Sn124Sn270AMeV"};

  TString amdconfig[] = { "cluster_SLy4",
                          "cluster_SLy4-L108",
                          "nocluster_SLy4",
			  "nocluster_SLy4-L108"};

  TString dfile[] = {"VN132Sn_proton.root",
		     "NLv1v2132Sn.root"};


  TString pamdname[]= {"prt","neut"};
  TString pname[]   = {"proton","deutron","trit","neut"};
  TString pfname[]  = {"proton","deuteron","triton","neutron"};

  UInt_t  imrk[] = {20, 21, 22, 23};
  Size_t  imsz[] = {1, 1, 1.3, 1.3};
  Color_t icol[] = { kRed,          kBlue,  kOrange-3,   kGreen+1,
		     kBlue+2,   kOrange+7,  kGreen-3,     kPink+9,
		     kGreen-3,    kPink+7,  kCyan-1,    kYellow-2,
		     kRed-2,      kBlue-2,  kOrange-2,   kGreen-1,
		     kBlue-3,   kOrange+5,  kGreen+3,     kPink-9};




  TString foname[20];
  UInt_t ico = 0;

  if( TPC ) {
    // Data:
    // i0 : particle
    // i1 : system
    cout << ico << endl;
  
    for(UInt_t i1 = 0; i1 < 4; i1++){

      if( !bsys[i1] ) continue;
 
      for(UInt_t i0 = 0; i0 < 3; i0++) {
	
	if( !bpid[i0] ) continue;

	//foname[ico] = "VN" + expbeamconfig[i1]+"_"+pfname[i0]+".root"; 
       	foname[ico] = "YPT" + expbeamconfig[i1]+"_"+pfname[i0]+".root"; 
	label.push_back( "DATA:" + expbeamconfig[i1] + " " + pname[i0] );
	isys.push_back(i1);
	ico++;
      }
    }
  }

  if( NeuLAND ) {
    // NeuLAND data
    // i1 : system
    for(UInt_t i1 = 0; i1 < 4; i1++ ){
	
      if( !bsys[i1] ) continue;

      // //      foname[ico] = "NLv1v2cm" + expbeamconfig[i1]+".root";
      // foname[ico] = "NL" + expbeamconfig[i1]+"_neutron.root";
      // label.push_back( "DATA:" + expbeamconfig[i1] + "neut" );
      // isys.push_back(i1);
      // ico++;

      // without beam angle rotation
      // foname[ico] = "noRotationNL" + expbeamconfig[i1]+"_neutron.root";
      // label.push_back( "DATA: no Rot." + expbeamconfig[i1] + "neut" );
      // isys.push_back(i1);
      // ico++;

      // veto_all == 0
      foname[ico] = "vaNL" + expbeamconfig[i1]+"_neutron.root";
      label.push_back( "DATA: no veto_all" + expbeamconfig[i1] + "neut" );
      isys.push_back(i1);
      ico++;

    }
  }
  

  if ( AMD ) {
    // AMD
    // i0 : particle  
    // i1 : system  
    // i2 : amdconfig 

    UInt_t i0 = 0;
    for(UInt_t i1 = 0; i1 < 4; i1++){

      if( !bsys[i1] ) continue;

      for(UInt_t i2 = 0; i2 < 4; i2++) {
	
	if( !bamd[i2] ) continue;
	
	//	foname[ico] = amdbeamconfig[i1] +"_"+amdconfig[i2]+"_ts"+pamdname[i0]+".root"; 
	foname[ico] = amdbeamconfig[i1] +"_"+amdconfig[i2]+"_"+pamdname[i0]+".root"; 
	label.push_back( amdbeamconfig[i1](0,5) + amdconfig[i2] + pamdname[i0] );
	isys.push_back(i1);
	ico++;
      }
    }
  }

  cout << " total file " << ico << endl;

  for(UInt_t lm = 0; lm < ico; lm++){
    std::cout << ico << " - " << lm << " : " << foname[lm] << std::endl;
  }

  gStyle->SetGridColor(7);
  gStyle->SetGridStyle(1);



  UInt_t   iicol = 0;
  TGraphErrors *gv_v1[20];
  TGraphErrors *gv_v2[20];
  TGraphErrors *ginv_v1[20];

  auto aleg1 = new TLegend(0.6,0.13, 0.89 ,0.4,"");
  auto aleg2 = new TLegend(0.6,0.65 ,0.9  ,0.9,"");
  auto mv1 = new TMultiGraph();
  auto mv2 = new TMultiGraph();

  mv1->SetTitle("; Rapidity_cm; v1");
  mv2->SetTitle("; Rapidity_cm; v2");

 
  for(UInt_t lm = 0; lm < ico; lm++){

    if( !gSystem->FindFile("data/",foname[lm]) ){
      cout << "File doesn't exist " << foname[lm] <<  endl;
      continue;
    }

    TFile *fopen = TFile::Open(foname[lm]);
    if( fopen->IsZombie() ) {
      cout << " File not opened " << foname[lm] << endl;
      continue;
    }


    cout << lm << " " << foname[lm] << endl;

    gv_v1[lm] = (TGraphErrors*)fopen->Get("gv_v1");

    if( gv_v1[lm] == NULL) {
      cout << "Error in getting graph " << foname[lm] << " " << lm << endl;
      continue;
    }

    gv_v1[lm] -> SetName(Form("gv_v1_%d",lm));

    gv_v1[lm] -> SetMarkerColor(icol[iicol]);
    gv_v1[lm] -> SetLineColor(icol[iicol]);
    gv_v1[lm] -> SetMarkerStyle(imrk[iicol]);

    // if( bYCM ) {
    //   ginv_v1[lm] = new TGraphErrors(*gv_v1[lm]);
    //   UInt_t npoint = gv_v1[lm]->GetN();
    //   for(UInt_t i = 0; i < npoint; i++) {
    // 	Double_t v1;
    // 	Double_t rapid;
    // 	gv_v1[lm]->GetPoint(i, rapid, v1);
	
    // 	ginv_v1[lm]->SetPoint(i, -rapid, -v1);

    // 	mv1->Add(ginv_v1[lm], "lp");
    // 	aleg1->AddEntry(ginv_v1[lm], "Inv."+label[lm]);
	
    //   }
    // }


    mv1->Add(gv_v1[lm],"lp");
    aleg1->AddEntry(gv_v1[lm], label[lm] );


    gv_v2[lm] = (TGraphErrors*)fopen->Get("gv_v2");

    gv_v1[lm]->Print();
    gv_v2[lm] -> SetName(Form("gv_v2_%d",lm));

    gv_v2[lm] -> SetMarkerColor(icol[iicol]);
    gv_v2[lm] -> SetLineColor(icol[iicol]);
    gv_v2[lm] -> SetMarkerStyle(imrk[iicol]);

    // if( !bYCM ) {
    //   UInt_t npoint = gv_v1[lm]->GetN();
    //   for(UInt_t i = 0; i < npoint; i++) {
    // 	Double_t v2;
    // 	Double_t rapid;
    // 	gv_v2[lm]->GetPoint(i, rapid, v2);
	
    // 	rapid = (rapid + 0.05);
    // 	gv_v2[lm]->SetPoint(i, rapid, v2);

    //   }
    // }

    iicol++;

    // gv_v2[lm]->RemovePoint(5);
    // gv_v2[lm]->RemovePoint(5);
    // gv_v2[lm]->RemovePoint(5);
    
    mv2->Add(gv_v2[lm],"lp");
    aleg2->AddEntry(gv_v2[lm], label[lm] );

    fopen->Close();
  }


  UInt_t ic = 0;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  mv1->Draw("alp");
  aleg1->Draw();


  auto LineV = new TLine(0.,mv1->GetYaxis()->GetXmin(), 0., mv1->GetYaxis()->GetXmax());
  auto LineH = new TLine(mv1->GetXaxis()->GetXmin(),    0., mv1->GetXaxis()->GetXmax(), 0.);
  LineV->SetLineStyle(3);
  LineH->SetLineStyle(3);
  LineV->Draw();
  LineH->Draw();


  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  mv2->Draw("alp");

  if( TPC ) {
    mv2->SetMaximum(0.08);
    mv2->SetMinimum(-0.08);
  }

  aleg2->Draw();

  LineV = new TLine(0.,mv2->GetYaxis()->GetXmin(), 0., mv2->GetYaxis()->GetXmax());
  LineH = new TLine(mv2->GetXaxis()->GetXmin(),    0., mv2->GetXaxis()->GetXmax(), 0.);
  LineV->SetLineStyle(3);
  LineH->SetLineStyle(3);
  LineV->Draw();
  LineH->Draw();

}



void plotAMDFlow()
{
  TCanvas *cc[2];

  TFile *fopen;

  TString beamconfig[] = {"132Sn124Sn270AMeV",
                          "108Sn112Sn270AMeV"};

  TString amdconfig[4] = { "cluster_SLy4",
                           "cluster_SLy4-L108",
                           "nocluster_SLy4",
                           "nocluster_SLy4-L108"};

  TString pname[2] = {"prt","neut"};
  TString pfname[2] = {"proton","neutron"};

  gStyle->SetGridColor(7);
  gStyle->SetGridStyle(1);


  UInt_t   icol[]     = { 2, 58, 78, 68, 53, 96, 156, 52, 100, 226, 229, 108};
  UInt_t iicol = 0;
  TGraphErrors *gv_v1[4];
  TGraphErrors *gv_v2[4];


  auto aleg1 = new TLegend(0.1,0.5, 0.4 ,0.9,"");
  auto aleg2 = new TLegend(0.3,0.65,0.62 ,0.9,"");
  auto mv1 = new TMultiGraph();
  auto mv2 = new TMultiGraph();

  // j : particle                                                                                                                                                           
  // k : beam                                                                                                                                                               
  // i : amdconfig                                                                                                                                      

  gSystem->cd("data/");                    
  UInt_t i = 0;
  for(UInt_t k = 0; k < 2; k++){
    for(UInt_t j = 0; j < 2; j++) {

      std::cout << " file " << beamconfig[k]+"_"+amdconfig[i]+"_"+pname[j]+".root" << " i " << i << " j " << j << " k " << k <<std::endl;
      fopen = new TFile(beamconfig[k]+"_"+amdconfig[i]+"_"+pname[j]+".root");
      if (fopen == NULL) return;

      mv1->SetTitle(pfname[j]+"; Rapidity_lab; v1");
      mv2->SetTitle(pfname[j]+"; Rapidity_lab; v2");

      gv_v1[i] = (TGraphErrors*)fopen->Get("gv_v1");
      gv_v1[i] -> SetName("gv_v1" + amdconfig[i]);

      gv_v1[i] -> SetMarkerColor(icol[iicol]);
      gv_v1[i] -> SetLineColor(icol[iicol]);

      mv1->Add(gv_v1[i],"lp");
      aleg1->AddEntry(gv_v1[i],beamconfig[k](0,5)+amdconfig[i]+"."+pname[j] );

      gv_v2[i] = (TGraphErrors*)fopen->Get("gv_v2");
      gv_v2[i] -> SetName("gv_v2" + amdconfig[i]);

      gv_v2[i] -> SetMarkerColor(icol[iicol]);
      gv_v2[i] -> SetLineColor(icol[iicol]);
      iicol++;

      mv2->Add(gv_v2[i],"lp");
      aleg2->AddEntry(gv_v2[i],beamconfig[k](0,5)+amdconfig[i]+"."+pname[j] );

      fopen->Close();
      delete fopen;
    }
  }
  UInt_t ic = 0;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  mv1->Draw("alp");
  aleg1->Draw();

  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  mv2->Draw("alp");
  aleg2->Draw();

  gSystem->cd("..");
}



void SaveCanvas(TString fopt = "", Int_t isel=-1)
{
  TString printHeader = "dNdy";

  if(isel > -1)
    gROOT->GetListOfCanvases()->At(isel)->SaveAs("v1v2y"+fopt+Form("_%d",isel)+".png");

  else {
    Int_t iCanvas = gROOT->GetListOfCanvases()->GetEntries();
    for(Int_t i = 0; i < iCanvas; i++)
      gROOT->GetListOfCanvases()->At(i)->SaveAs("v1v2y"+fopt+Form("_%d",i)+".png");
  }
}
