Double_t      RPPsi;
TClonesArray *aFlowArray;
TClonesArray *aArray;




void mc4()
{

  gROOT->Macro("SetStyle.C");

  auto hcos1 = new TH1D("hcos1",";cos(  #delta #phi )",100, -1., 1.);
  auto hcos2 = new TH1D("hcos2",";cos( 2#delta #phi )",100, -1., 1.);

  auto rChain = new TChain("cbmsim");
  
  //  rChain->Add("data/run0022_rpsim.v22.root");
  rChain->Add("data/run0400_BTt.v50.root");

  if( rChain == NULL ) {
    std::cout << " No file " << std::endl;
    exit(0);
  }

  //  rChain->Print();

  rChain->SetBranchAddress("STParticle",&aArray);  
  //rChain->SetBranchAddress("STFlow"    ,&aFlowArray);
  rChain->SetBranchAddress("RPPsi"     ,&RPPsi);

  Long64_t nEntry = rChain->GetEntries();

  cout << " Entry : " << nEntry << endl;

  for(Long64_t i = 0; i < nEntry; i++ ) {

    rChain->GetEntry(i);

    
    //    STFlowInfo *aflow = (STFlowInfo*)aFlowArray->At(0);

    TIter next(aArray);
    STParticle *aPart = NULL;

    UInt_t mtk = 0;
    while( (aPart = (STParticle*)next()) ) {

      auto ycm = aPart->GetRapiditycm()/0.36;

      if( ycm > 0.8 ) 
	hcos1->Fill( cos( aPart->GetRotatedMomentum().Phi() - RPPsi ) );

      if( abs(ycm) < 0.1 )
	hcos2->Fill( cos(2.*(aPart->GetRotatedMomentum().Phi() - RPPsi) ) );

    } 
  }

  TCanvas *cc1 = new TCanvas("cc1","cc1");
  cc1->Divide(1,2);
  cc1->cd(1);
  hcos1->Draw();

  cc1->cd(2);
  hcos2->Draw();
  
  std::cout << " <cos>" << hcos1->GetMean()
	    << " <cos2> "<< hcos2->GetMean() << std::endl;
}


void mc3()
{
  gROOT->Macro("SetStyle.C");

  TFile *file0 = TFile::Open("data/run0400_BTt.v50.root");
  TTree *cbmsim = (TTree*)file0->Get("cbmsim");


  auto hpsit    = new TH1D("hpsit" ,"; #Psi; 2#pi/NdN/d#Psi",100,-TMath::Pi(), TMath::Pi());
  auto hpsic    = new TH1D("hpsic"  ,"; #Psi; 2#pi/NdN/d#Psi",100,-TMath::Pi(), TMath::Pi());
  auto hpsif    = new TH1D("hpsif"  ,"; #Psi; 2#pi/NdN/d#Psi",100,-TMath::Pi(), TMath::Pi());

  cbmsim->Project("hpsit"  ,"unitP.Phi()");
  cbmsim->Project("hpsic" ,"unitP_fc.Phi()");


  TFile *file1 = TFile::Open("data/run0022_rpsim.v22.root");
  TTree *cbmsimf = (TTree*)file1->Get("cbmsim");
  cbmsimf->Project("hpsif"  ,"unitP.Phi()");


  hpsit ->SetNormFactor(100);
  hpsif ->SetNormFactor(100);
  hpsic->SetNormFactor(100);

  TCanvas *cc1 = new TCanvas("cc1","cc1");
  hpsit -> SetLineColor(2);
  hpsic -> SetLineColor(4);
  hpsif -> SetLineColor(7);

  hpsit -> Draw("e");
  //  hpsic -> Draw("same e");
  hpsif -> Draw("same e");

}

void mc2()
{

  gROOT->Macro("SetStyle.C");

  TFile *file0 = TFile::Open("data/run0022_rpsim.v22.root");

  TTree *cbmsim = (TTree*)file0->Get("cbmsim");


  TH2D *hacp = new TH2D("hacp",";y_{cm}; Pt[MeV/c]",200,-0.42,0.42, 200,0.,800);
  cbmsim->Project("hacp","fRotatedP3.Pt():fRapiditycm");

  TCanvas *cc1 = new TCanvas("cc1","cc1");
  hacp->Draw("colz");
}




void mc1()
{
  gROOT->Macro("SetStyle.C");

  TFile *file0 = TFile::Open("data/run0400_BTt.v50.root");
  
  TF1 *fv2y = (TF1*)file0->Get("fv2y");
  fv2y->SetTitle("; y_{cm}/y_{beam}; v2");

  TLatex *Lv2 = new TLatex(-0.6,-0.01,"-0.08 + 0.1y^2 -0.02y^4");

  TCanvas *cc2 = new TCanvas("cc2","cc2");
  fv2y->SetLineColor(2);
  fv2y->Draw();
  Lv2->Draw();


  TF1 *fv1y = (TF1*)file0->Get("fv1y");
  fv1y->SetTitle("; y_{cm}/y_{beam}; v1");

  TLatex *Lv1 = new TLatex(-0.6,0.2,"0.52*y - 0.18*y^3");


  TCanvas *cc1 = new TCanvas("cc1","cc1");
  fv1y->SetLineColor(2);
  fv1y->Draw();
  Lv1->Draw();


}

void PlotRPSim()
{
  mc3();

}
