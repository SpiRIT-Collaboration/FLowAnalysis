#include  "../analysisformat/STRunToBeamA.hh"

std::map< UInt_t, UInt_t > beamAToSys = { {132, 0}, {108, 1}, {124, 2}, {112,3}};

TString lpid[] = {"1H",    "2H",      "3H"    ,"3He","4He","N"      ,"H"};  
TString pidName[]={"Proton","Deuteron","Triton","Helium3","Alpha"}; 
TString fpid[] = {"proton","deuteron","triton","3He","4He","neutron","H"};  

TString  bName[]   = {"132Sn_","108Sn_","124Sn_","112Sn_","100Sn_"};
TString rsys[] = {"132",        "108",        "124",        "112"};
const double yAA[]   = {0.3822, 0.3647, 0.3538, 0.3902};
const double yNN[]   = {0.3696, 0.3697, 0.3705, 0.3706};
const double yBeam[] = {0.7421, 0.7423, 0.7439, 0.7441};
const double uNN[]   = {0.3709, 0.3709, 0.3718, 0.3719};

void makeFlowBranch()
{
  gROOT->Reset();

  TDatime dtime;
  TRandom3 grnd(dtime.GetSecond());
  gRandom->SetSeed(dtime.GetSecond());

  TString sRun = gSystem -> Getenv("RUN");
  TString sVer = gSystem -> Getenv("VER");
  UInt_t  iRun = atoi(sRun);
  UInt_t  bmA =  STRunToBeamA::GetBeamA(iRun);
  //  UInt_t systemID = beamAToSys[bmA];
  UInt_t systemID = 1;

  cout << " beam is " << bmA << " and systemID is " << systemID <<  endl;

  //@@@@@Input++++++++++++++++++++++++++++++++++++++++++++++++++
  auto fChain = new TChain("cbmsim");
  
  TString filename = Form("run%04d_ph.v%s.root",iRun,sVer.Data());
  if( gSystem->FindFile("data",filename) ){
    fChain->Add(filename);
    cout << " LOAD " << filename << endl;
  }
  else {
    cout << " No data is found with RUN= " << iRun << endl;
    return;
  }

  Int_t           beamPID;
  Int_t           run;
  Int_t           event;
  Int_t           mtrack0;
  Int_t           mtrack1;
  Int_t           mtrack2;
  auto ParticleArray = new TClonesArray("STParticle",80);
  auto aFlowArray    = new TClonesArray("STFlowInfo",1);

  // Set branch addresses.
  fChain->SetBranchAddress("beamPID",&beamPID);
  fChain->SetBranchAddress("run",&run);
  fChain->SetBranchAddress("event",&event);
  fChain->SetBranchAddress("mtrack0",&mtrack0);
  fChain->SetBranchAddress("mtrack1",&mtrack1);
  fChain->SetBranchAddress("mtrack2",&mtrack2);
  fChain->SetBranchAddress("Particle",&ParticleArray);

  //@@@@@Output++++++++++++++++++++++++++++++++++++++++++++++++++
  TFile* outfile = new TFile(Form("data/run%04d_flow.v%s.root",iRun,sVer.Data()),"recreate");
  
  auto outtree = new TTree("cbmsim","Flow");
  auto aflow = new STFlowTask(kFALSE, kTRUE, kFALSE); //flattening, subevent
  STFlowInfo* aflowInfo = new STFlowInfo();

  // outtree->Branch("run",&run); 
  // outtree->Branch("event",&event);
  // outtree->Branch("mtrack0",&mtrack0);
  // outtree->Branch("mtrack1",&mtrack1);
  // outtree->Branch("mtrack2",&mtrack2);
  outtree->Branch("STParticle",&ParticleArray);
  outtree->Branch("STFlow"    ,&aFlowArray);

  aflow -> Init(iRun, sVer);

  Int_t evt_prev = -1; 
  Long64_t nentries = fChain->GetEntries();
  //  nentries = 10;
  cout << nentries << " will be analized. " << endl;

  for ( auto ievt : ROOT::TSeqL(nentries)) {
    fChain->GetEntry(ievt);

    aflow->SetFlowTask( *ParticleArray );

    aflow->SetupEventInfo(ievt, bmA);
    aflow->FinishEvent();

    aflowInfo = aflow->GetFlowInfo();
    aflowInfo -> SetNTrack(0,mtrack0);
    aflowInfo -> SetNTrack(1,mtrack1);
    aflowInfo -> SetNTrack(2,mtrack2);

    new( (*aFlowArray)[0] ) STFlowInfo(*aflowInfo);
    
    outtree->Fill();
  }


  outfile->Write();
  cout << outfile->GetName() << " is saved. " << endl;
  
}
