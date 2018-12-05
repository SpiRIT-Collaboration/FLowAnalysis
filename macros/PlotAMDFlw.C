void PlotAMDFlw()
{

  TFile *fOpen;

  TString beamconfig[] = {"132Sn",
			  "108Sn"};

  TString amdconfig[4] = { "cluster_SLy4",
			   "cluster_SLy4-L108",
			   "nocluster_SLy4",
			   "nocluster_SLy4-L108"};

  TString pfname[] = {"proton","deuteron" ,"triton","neutron"};
  TString pname[] = {"prot","deut" ,"trit","neut","neut"};


  Bool_t bbeam[] = {1, 1};
  Bool_t beos[]  = {0, 0, 1, 1};
  Bool_t bpart[] = {0, 0, 0, 1}; //p, d, t, n
  Bool_t rc = 1;

  gStyle->SetGridColor(7);
  gStyle->SetGridStyle(1);


  UInt_t icol[] = { 2, 58, 78, 68, 53, 96, 156, 52, 100, 226, 229, 108};
  UInt_t iicol  = 0;


  auto aleg1 = new TLegend(0.12,0.55 ,0.43 ,0.9,"");
  auto aleg2 = new TLegend(0.1 ,0.13 ,0.4  ,0.3,"");
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
      
      for(UInt_t j = 0; j < 4; j++) {
	if( !bpart[j] ) continue;
	mv1->SetTitle(pfname[j]+"; Rapidity_cm; v1");
	mv2->SetTitle(pfname[j]+"; Rapidity_cm; v2");


	TString fname =  beamconfig[k]+"-"+amdconfig[i]+"-"+pname[j]+".root" ;
	std::cout << " file " << fname ;
	if( !gSystem->FindFile("data/", fname) ) continue; 

	fOpen = TFile::Open(fname);

	gROOT->ls();

	if (fOpen == NULL) return;
	std::cout << " opened " << std::endl;

	TGraphErrors *yv1 = (TGraphErrors*)fOpen->Get("gv_v1");
	yv1 -> SetMarkerColor(icol[iicol]);
	yv1 -> SetLineColor(icol[iicol]); 

	mv1->Add(yv1,"lp");
	aleg1->AddEntry(yv1,beamconfig[k](0,5)+amdconfig[i]+"."+pname[j] );

	TGraphErrors *yv2 = (TGraphErrors*)fOpen->Get("gv_v2");
	yv2 -> SetMarkerColor(icol[iicol]);
	yv2 -> SetLineColor(icol[iicol]); 
	iicol++;

	mv2->Add(yv2,"lp");
	aleg2->AddEntry(yv2,beamconfig[k](0,5)+amdconfig[i]+"."+pname[j] );
    
	fOpen->Close();
    
      }
    }
  }



  UInt_t ic = 0;
  auto cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
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



  ic++;
  cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
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

}
