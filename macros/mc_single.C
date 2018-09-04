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
Double_t amu  = 931.4940954; 
Double_t mass_p = 1.00727646688*amu;

TString printHeader = "";

void mc_open()
{
  TString  MCDIR     = "/cache/scr/spirit/kaneko/rootfiles/recoDataMC/20180711dev_20180713dev.RC/";
  TString  STVERSION = "1660.2c1f3c3";
  TString  MCSYSTEM  = "singleProton_";
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


void mc_single()
{
  mc_open();

  //booking
  auto hangle    = new TH2D("hangle","; #theta_init; #theta_rc", 100, 0., 1.6, 100, 0., 1.6);
  auto hangle2   = new TH2D("hangle2","; #theta_init; #theta_diff", 100, 0., 1.6, 100, -1.6, 1.6);
  auto hangy     = new TH2D("hangy","; Rapidity; #theta_diff", 100, -0.4, 0.4, 100, -1.6, 1.6);
  auto hyy       = new TH2D("hyy"  ,"; inti Rapidity; RC Rapidity", 100, -0.4, 0.4, 100, -0.4, 0.4);
  auto hyydiff   = new TH2D("hyydiff","; inti Rapidity; diff Rapidity", 100, -0.4, 0.4, 100, -0.4, 0.4);

  auto hangleaccI = new TH2D("hanglaccI","Generated;#Theta, #Phi", 100,0., 1.6, 100, -3.2, 3.2);
  auto hangleaccR = new TH2D("hanglaccR","Reconstructed ;#Theta, #Phi", 100,0., 1.6, 100, -3.2, 3.2);

  auto hirp       = new TH2D("hirp",    "proton; Initial Momentum; Reconstructed P",200,0.,1500.,200,0.,1500.);

  auto hirpd      = new TH2D("hirpd",   "proton; Initial Momentum; Rec - Init P",200,0.,2000.,200,-500., 500.);

  auto hirpz      = new TH2D("hirpz",   "proton; Initial Momentum Z; Rec - Init Pz",200,0.,1500.,200,0., 1500.);
  auto hirpzd     = new TH2D("hirpzd",  "proton; #Theta ; Pz_Rec - Pz_Init ",100,0., 1.6, 200,-50.,50.);


  UInt_t nevt = SetBranch();

  TVector3 boostVec = LorentzBoost(4);

  UInt_t pPID = 0;
  for(UInt_t i = 0; i < nevt; i++){
    
    rChain0->GetEntry(i);

    STMCTrack* primTrack = (STMCTrack*)mceventHeader->GetPrimaryTrack(0); 
    TVector3 PrimTrackVect;
    primTrack->GetMomentum(PrimTrackVect);
    PrimTrackVect.SetMag( PrimTrackVect.Mag()*1000.);

    auto pRapidity =  GetRapidity_cm( PrimTrackVect, mass_p, -boostVec);

    hangleaccI->Fill(PrimTrackVect.Theta(), PrimTrackVect.Phi() );

    //    cout << " prapid " << pRapidity << " P " << PrimTrackVect.Mag()<< endl;

    //    angX = 0.1;

    TIter next( aTrackArray );
    STRecoTrack *aTrack = NULL;

    while( (aTrack = (STRecoTrack*)next() ) ) {

      hangle->Fill(PrimTrackVect.Theta(), aTrack->GetMomentum().Theta());

      hangle2->Fill(PrimTrackVect.Theta(), (aTrack->GetMomentum().Theta() - PrimTrackVect.Theta()));

      hirp->Fill(PrimTrackVect.Mag(), aTrack->GetMomentum().Mag());
      hirpd->Fill(PrimTrackVect.Mag(), aTrack->GetMomentum().Mag()-PrimTrackVect.Mag());

      hirpz->Fill(PrimTrackVect.Pz(), aTrack->GetMomentum().Pz());
      hirpzd->Fill(aTrack->GetMomentumTargetPlane().Theta(), aTrack->GetMomentumTargetPlane().Pz()-PrimTrackVect.Pz());


      auto rapidity = GetRapidity_cm( aTrack->GetMomentum(), mass_p, -boostVec);

      hangy->Fill(rapidity, (aTrack->GetMomentum().Theta() - PrimTrackVect.Theta()));

      hyy->Fill( pRapidity, rapidity);
      hyydiff->Fill( pRapidity, rapidity-pRapidity);

      hangleaccR->Fill(PrimTrackVect.Theta(), PrimTrackVect.Phi() );
    }
  }
  
  TFile *fout = new TFile("singProton_Eff.root","recreate");
  hangleaccI->Write();
  hangleaccR->Write();
  fout->Close();
  delete fout;



  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  hirp->Draw("colz");

  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  hirpd->Draw("colz");

  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  hirpz->Draw("colz");

  ic++;
  cc[ic] = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  hirpzd->FitSlicesY();
  TH1D* hirpzd_1 = (TH1D*)gROOT->Get("hirpzd_1");
  hirpzd->Draw("colz");
  hirpzd_1->SetLineColor(2);
  hirpzd_1->Draw("same");
}


