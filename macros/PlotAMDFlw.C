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


  Bool_t bbeam[] = {1, 0};
  Bool_t beos[]  = {1, 1, 1, 1};
  Bool_t bpart[] = {1, 0};
  Bool_t rc = 1;

  gStyle->SetGridColor(7);
  gStyle->SetGridStyle(1);


  UInt_t   icol[]     = { 2, 58, 78, 68, 53, 96, 156, 52, 100, 226, 229, 108};
  UInt_t iicol = 0;
  TGraphErrors *gv_v1[4];
  TGraphErrors *gv_v2[4];


  auto aleg1 = new TLegend(0.1,0.5,0.4 ,0.9,"");
  auto aleg2 = new TLegend(0.1,0.5,0.4 ,0.9,"");
  auto mv1 = new TMultiGraph();
  auto mv2 = new TMultiGraph();

  // j : particle
  // k : beam
  // i : amdconfig

  UInt_t il = 0;
  for(UInt_t i = 0; i < 4; i++){
    if( !beos[i] ) continue;
    
    for(UInt_t k = 0; k < 2; k++){
      if( !bbeam[k] ) continue;
      
      for(UInt_t j = 0; j < 2; j++) {
	if( !bpart[j] ) continue;
	mv1->SetTitle(pfname[j]+"; Rapidity_cm; v1");
	mv2->SetTitle(pfname[j]+"; Rapidity_cm; v2");

	std::cout << " file " << beamconfig[k]+"_"+amdconfig[i]+"_"+pname[j]+".root" << " i " << i << " j " << j << " k " << k <<std::endl;

	fopen = new TFile("data/"+beamconfig[k]+"_"+amdconfig[i]+"_"+pname[j]+".root");
	if (fopen == NULL) return;


	gv_v1[il] = (TGraphErrors*)fopen->Get("gv_v1");
	gv_v1[il] -> SetName("gv_v1" + amdconfig[i]);
    
	gv_v1[il] -> SetMarkerColor(icol[iicol]);
	gv_v1[il] -> SetLineColor(icol[iicol]); 

	mv1->Add(gv_v1[il],"lp");
	aleg1->AddEntry(gv_v1[il],beamconfig[k](0,5)+amdconfig[i]+"."+pname[j] );

	gv_v2[il] = (TGraphErrors*)fopen->Get("gv_v2");
	gv_v2[il] -> SetName("gv_v2" + amdconfig[i]);
    
	gv_v2[il] -> SetMarkerColor(icol[iicol]);
	gv_v2[il] -> SetLineColor(icol[iicol]); 
	iicol++;

	mv2->Add(gv_v2[il],"lp");
	aleg2->AddEntry(gv_v2[il],beamconfig[k](0,5)+amdconfig[i]+"."+pname[j] );
    
	il++;
	fopen->Close();
    
      }
    }
  }

  if( rc ){
    fopen = new TFile("data/rcamd_132Sn124Sn270AMeV_cluster_SLy4_proton.root");
    if (fopen != NULL) {
      
      gv_v1[il] = (TGraphErrors*)fopen->Get("gv_v1");
      gv_v1[il] -> SetName("gv_v1rc");
      
      gv_v1[il] -> SetMarkerStyle(20);
      gv_v1[il] -> SetMarkerColor(2);
      gv_v1[il] -> SetLineColor(2);
      
      mv1->Add(gv_v1[il],"lp");
      aleg1->AddEntry(gv_v1[il],"RC AMD 132Sn Proton" );
      il++;
    }

    fopen = new TFile("data/rcamd_132Sn124Sn270AMeV_cluster_SLy4_R_proton.root");
    if (fopen != NULL) {
      
      gv_v1[il] = (TGraphErrors*)fopen->Get("gv_v1");
      gv_v1[il] -> SetName("gv_v1rc");
      
      gv_v1[il] -> SetMarkerStyle(20);
      gv_v1[il] -> SetMarkerColor(4);
      gv_v1[il] -> SetLineColor(4);
      
      mv1->Add(gv_v1[il],"lp");
      aleg1->AddEntry(gv_v1[il],"RC AMD 132Sn Proton" );
      il++;
    }

    fopen = new TFile("data/rcamd_132Sn124Sn270AMeV_cluster_SLy4_pR_proton.root");
    if (fopen != NULL) {
      
      gv_v1[il] = (TGraphErrors*)fopen->Get("gv_v1");
      gv_v1[il] -> SetName("gv_v1rc");
      
      gv_v1[il] -> SetMarkerStyle(20);
      gv_v1[il] -> SetMarkerColor(8);
      gv_v1[il] -> SetLineColor(8);
      
      mv1->Add(gv_v1[il],"lp");
      aleg1->AddEntry(gv_v1[il],"RC AMD 132Sn Proton" );
      il++;
    }


  }


  UInt_t ic = 0;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  mv1->Draw("alp");
  aleg1->Draw();
  
  auto aLineX1 = new TLine(mv1->GetXaxis()->GetXmin(), 0., mv1->GetXaxis()->GetXmax(), 0.);
  aLineX1->SetLineColor(1);
  aLineX1->SetLineStyle(3);
  aLineX1->Draw();
  auto aLineY1 = new TLine(0., mv1->GetYaxis()->GetXmin(), 0., mv1->GetYaxis()->GetXmax());
  aLineY1->SetLineColor(1);
  aLineY1->SetLineStyle(3);
  aLineY1->Draw();


  return;
  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  mv2->Draw("alp");
  aleg2->Draw();

  auto aLineX2 = new TLine(mv2->GetXaxis()->GetXmin(), 0., mv2->GetXaxis()->GetXmax(), 0.);
  aLineX2->SetLineColor(1);
  aLineX2->SetLineStyle(3);
  aLineX2->Draw();
  auto aLineY2 = new TLine(0., mv2->GetYaxis()->GetXmin(), 0., mv2->GetYaxis()->GetXmax());
  aLineY2->SetLineColor(1);
  aLineY2->SetLineStyle(3);
  aLineY2->Draw();



}
