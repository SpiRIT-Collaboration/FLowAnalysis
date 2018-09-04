#include "GetLorentzBoost.C"

TCanvas *cc[10];
UInt_t ic = -1;
TChain *rChain0 = new TChain("cbmsim");
TClonesArray *aTrackArray;
STVertex     *avertex;
STMCEventHeader *mceventHeader;
TClonesArray *mctriggerResponse;

const UInt_t nybin = 10;
std::vector< std::vector< Double_t > >vrapd(nybin);
TH1D* hphi1[nybin];
TH1D* hphi2[nybin];


TCutG *gProton;
TCutG *gDeuteron;
TCutG *gTriton;
TH2D  *effProtonI;
TH2D  *effProtonR;
TAxis *xaxis;
TAxis *yaxis;

Double_t amu  = 931.4940954; 
Double_t mass_p = 1.00727646688*amu;

TString printHeader = "";

void mc_open()
{
  TString  MCDIR     = "/cache/scr/spirit/kaneko/rootfiles/recoDataMC/20180711dev_20180713dev.RC/";
  TString  STVERSION = "1660.2c1f3c3";
  TString  MCSYSTEM  = "amd_132Sn124Sn270AMeV_cluster_SLy4_";
  //  TString  MCSYSTEM  = "amd_132Sn124Sn270AMeV_nocluster_SLy4_";
  printHeader = "data/rc"+MCSYSTEM;

  UInt_t inum = 0;
  while(1) {
    
    TString fname = Form(MCSYSTEM+"%d.reco.develop."+STVERSION+".root", inum);
    if(gSystem->FindFile(MCDIR, fname)) {
      cout << fname << endl;
      rChain0->Add(fname);
    }
    else
      break;

    inum++;
    
       if( inum > 2 ) break;
  }

  cout << inum << " files opend. " << endl;


  TFile *gcutFile = new TFile("data/gcutPID132Sn.ROOT");
  gProton   = (TCutG*)gcutFile->Get("gcutProton132Sn2");
  gDeuteron = (TCutG*)gcutFile->Get("gcutDeutron132Sn");
  gTriton   = (TCutG*)gcutFile->Get("gcutTriton132Sn");
  gcutFile->Close();
  delete gcutFile;

  TFile *effFile = new TFile("singProton_Eff.root");
  effFile->ls();
  effProtonI = (TH2D*)effFile->Get("hanglaccI");
  effProtonR = (TH2D*)effFile->Get("hanglaccR");
  if( effProtonI != NULL ) {
    xaxis = effProtonI->GetXaxis();
    yaxis = effProtonI->GetYaxis();
    gcutFile->Close();
  }
  else
    cout << " efficiency file not found. " << endl;
  effFile->Close();
  delete effFile;
}

Double_t GetEfficiency(Double_t theta, Double_t phi)
{
  Int_t xbin = xaxis->FindBin(theta);
  Int_t ybin = yaxis->FindBin(phi);

  cout << theta << " " << phi 
       << " eff " << xbin << " : " << ybin << " -- " ;

  if( xbin < 0 || ybin < 0 ) return 0.;

  Double_t initial = effProtonI->GetBinContent(xbin, ybin);
  Double_t reco    = effProtonR->GetBinContent(xbin, ybin);
  Double_t eff = (Double_t)reco/(Double_t)initial;

  cout << reco << " / " << initial << " = "
       << eff;
  
  if( eff > 1.)
    return 1.;

  else 
    return eff;
}

UInt_t SetBranch()
{

  if(rChain0 == NULL) {
    std::cout << " no file is loaded " << std::endl;
    return 0;
  }


  rChain0 -> SetBranchAddress("STRecoTrack",&aTrackArray);
  rChain0 -> SetBranchAddress("STVertex"   , &avertex);
  rChain0 -> SetBranchAddress("STMCEventHeader", &mceventHeader);
  rChain0 -> SetBranchAddress("STMCTriggerResponse",&mctriggerResponse);


  return rChain0->GetEntries();
}

void GetFittingParameters(TH1D &h1, Double_t pp[6])
{
  Double_t nbin = h1.GetNbinsX();
  Double_t scf  = h1.GetEntries();
  h1.Scale(nbin/scf);
  auto *fcos1 = new TF1("fcos1","[0]+2.*[1]*cos(x)"   ,-1.*TMath::Pi(),TMath::Pi());
  h1.Fit("fcos1","Q","",-3.1,3.1);

  pp[0] = fcos1->GetParameter(0);
  pp[1] = fcos1->GetParameter(1);
  pp[2] = fcos1->GetParError(0);
  pp[3] = fcos1->GetParError(1);

}


void mc_calcFlw()
{
  mc_open();
}

void getv1v2()
{
  //booking
  auto hdedxp    = new TH2D("hdedxp","; Momentum; dEdx", 200, 0., 1500., 200, 0., 800.);
  auto hdedxpprt = new TH2D("hdedxpprt","; Momentum; dEdx", 200, 0., 1500., 200, 0., 800.);
  auto hyphi     = new TH2D("hyphi",";",100, -3.14, 3.14, nybin, -0.4, 0.45);

  auto hangle    = new TH2D("hangle","; #theta; #phi", 100, 0., 1.6, 100, -3.2, 3.2);

  for(UInt_t i = 0 ; i < nybin; i++ ){
    hphi1[i] = new TH1D(Form("hphi1_%d",i),Form("hphi1_%d",i), 100, -3.2, 3.2);
    hphi2[i] = new TH1D(Form("hphi2_%d",i),Form("hphi2_%d",i), 100,  0. , 3.2);
  }


  UInt_t nevt = SetBranch();

  TVector3 boostVec = LorentzBoost(4);

  UInt_t pPID = 0;
  for(UInt_t i = 0; i < nevt; i++){
    
    rChain0->GetEntry(i);

    //    STMCEventHeader *aEvt = (STMCEventHeader*)(mceventHeader->At(0));
    Double_t rp   = mceventHeader->GetReactionPlane();
    Double_t angX = mceventHeader->GetBeamAngle().X();

    //    angX = 0.1;

    TIter next( aTrackArray );
    STRecoTrack *aTrack = NULL;

    while( (aTrack = (STRecoTrack*)next() ) ) {
      Double_t p    = aTrack->GetMomentum().Mag();
      Double_t dedx = aTrack->GetdEdxWithCut(0, 0.7);

      //      if( avertex->GetPos().Mag() < 10

      hdedxp->Fill(p, dedx);

      if(gProton->IsInside(p, dedx)) {
	//      if(gDeuteron->IsInside(p, dedx)) {
	//      if(gTriton->IsInside(p, dedx)) {
	hdedxpprt->Fill(p, dedx);
	pPID = 2212;

	TVector3 rotateP = aTrack->GetMomentum();
	rotateP.RotateY(-angX);
	
	//	if( rotateP.Theta() < 0.6) continue;
	
	cout << " Eff " << GetEfficiency( rotateP.Theta(), rotateP.Phi()) << endl;

	auto rapidity = GetRapidity_cm( rotateP, mass_p, -boostVec);

	hyphi->Fill(TVector2::Phi_mpi_pi( rotateP.Phi() - rp), rapidity);

	hangle->Fill( rotateP.Theta(), rotateP.Phi() );

        for(UInt_t i = 0; i < nybin; i++) {
          Double_t upbin = hyphi->GetYaxis()->GetBinUpEdge(i);
          if( rapidity < upbin ) {
            vrapd[i].push_back(rapidity);
            hphi1[i]->Fill( TVector2::Phi_mpi_pi( rotateP.Phi() - rp) );
            hphi2[i]->Fill( TVector2::Phi_mpi_pi(2. * rotateP.Phi() - rp) );
            break;
          }
        }
      }
    }
  }
  

  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  hdedxpprt->Draw("colz");

  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);  
  //  cc[ic]->SetLogz();
  //  hyphi->Draw("colz");
  hangle->Draw("colz");

  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);  

  TFile *froot = new TFile(printHeader+"test_proton.root","recreate");
  auto gv_v1 = new TGraphErrors();
  gv_v1->SetName("gv_v1");
  gv_v1->SetTitle("AMD proton; Rapidity ; v1 ");


  UInt_t ip1 = 0;
  for(UInt_t jn = 1; jn < nybin; jn++){

    if( vrapd[jn].size() > 0 ){
      Double_t rpd  = TMath::Mean(vrapd[jn].begin(), vrapd[jn].end());
      Double_t rpde = TMath::StdDev(vrapd[jn].begin(), vrapd[jn].end());

      //    auto ybin = (TH1D*)hyphi->ProjectionX((TString)Form("hydphi1_%d",jn+1),jn+1, jn+1,"eo");
    
      Double_t para[6];
      GetFittingParameters(*hphi1[jn], para);

      Double_t v1  = para[1];
      Double_t v1e = para[3];

      gv_v1->SetPoint(ip1, rpd, v1);
      gv_v1->SetPointError(ip1, rpde, v1e);

      ip1++;
    }
  }


  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  gv_v1->Draw("ALP");
  gv_v1->Write();

  auto aLineX1 = new TLine(gv_v1->GetXaxis()->GetXmin(), 0., gv_v1->GetXaxis()->GetXmax(), 0.);
  aLineX1->SetLineColor(1);
  aLineX1->SetLineStyle(3);
  aLineX1->Draw();


  auto aLineY1 = new TLine(0., gv_v1->GetYaxis()->GetXmin(), 0., gv_v1->GetYaxis()->GetXmax());
  //    auto aLineY1 = new TLine(ybeam_cm[isys[m]], Ymin, ybeam_cm[isys[m]], Ymax);                                                       
  aLineY1->SetLineColor(1);
  aLineY1->SetLineStyle(3);
  aLineY1->Draw();


}


void compInitP()
{
  auto hinpz  = new TH1D("hinpz","initial Pz",200,0.,1500.);
  auto hrcpz  = new TH1D("hrcpz","Reconst. Pz",200,0.,1500.);
  auto hirpz  = new TH2D("hirpz","Initial and Reconst. Pz",200,0.,1500.,200,0.,1500.);


  UInt_t pPID = 0;
  for(UInt_t i = 0; i < nevt; i++){

    rChain0->GetEntry(i);

    //    STMCEventHeader *aEvt = (STMCEventHeader*)(mceventHeader->At(0));                                                                
    Double_t rp   = mceventHeader->GetReactionPlane();
    Double_t angX = mceventHeader->GetBeamAngle().X();
    TIter next( aTrackArray );
    STRecoTrack *aTrack = NULL;

    while( (aTrack = (STRecoTrack*)next() ) ) {
      Double_t p    = aTrack->GetMomentum().Mag();
      Double_t dedx = aTrack->GetdEdxWithCut(0, 0.7);

      if(gProton->IsInside(p, dedx)) {

	hinpz->Fill(



	ic++;
	cc[ic] = new TCanvas(Form("cc%d",ic),"cc");
  cc[ic]->Divide(2,2);
  cc[ic]->cd(1);
  rChain0->Draw("fPz*1000.>>hinpz","gProton","");
  cc[ic]->cd(2);
  rChain0->Draw("fMomentum.Z()>>hrcpz","gProton","");
  cc[ic]->cd(3);
  rChain0->Draw("fMomentum.Z():fPz*1000.>>hirpz","gProton","colz");
  


}
