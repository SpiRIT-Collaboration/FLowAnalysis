#include "flw_process1.h"
//----------------------------------------------------------------------
// Macro: flw_process1
//  (c) Mizuki Kurata-Nishimura
//   2017 Aug. 30 
//----------------------------------------
// Assemble TPC reconstructed data, Kyoto,  Katana, and BigRIPS root
//----------------------------------------------------------------------
void Setup()
{
  sRun = gSystem -> Getenv("RUN");
  sVer = gSystem -> Getenv("VER");

  if( sRun=="" || sVer=="" || !DefineVersion() ) {
    cout << "Plase type " << endl;
    cout << "$ RUN=#### VER=# root flw_process1.C" << endl;
    exit(0);
  }

  iRun = atoi(sRun);

  SnA = 0;
  if(iRun >= 2174 && iRun <= 2509)
    SnA = 108;
  else if( iRun >= 2836 && iRun <= 3039)
    SnA = 132;
  else if( iRun >= 3058 && iRun <= 3184)
    SnA = 124;
  else if( iRun >= 2520 && iRun <= 2653)
    SnA = 112;

  BigRIPS  = kTRUE;  //kFALSE;
  KyotoArry= kFALSE; //kTRUE;
  KATANA   = kFALSE; //kTRUE;
}

void flw_process1(Int_t nevt = -1)
{
  //////////////////////////////////////////////////////////
  // The calibrated beam_run#.ridf.root can be loaded.
  //////////////////////////////////////////////////////////

  Setup();

  //Reset ROOT and connect tree file
  gROOT->Reset();

  //----- Beam PID ------------------------
  BeamPID();

  //----- TPC data ------------------------
  SetTPC();
  Int_t nEvtTPC = 0;
  if(fChain) nEvtTPC = fChain -> GetEntries();
  cout << "Number of events in TPC: " << nEvtTPC << endl;
  if(nEvtTPC == 0) return;
  
  //----- BigRIPS data --------------------
  Long64_t nEvents = 0;
  if(BigRIPS) {
    SetBigRIPS();
    if( ribfChain) {
      nEvents   = ribfChain->GetEntries();
      Int_t nEventsBDC = bdcChain->GetEntries();

      if(nEvents != nEventsBDC) {
	cout << "Inconsistent event number in bigRIS RIDF (quit)" << endl;
	exit(0);
      }
    }
    cout << "Number of events in RIDF: " << nEvents << endl;
  }
  //----- KATANA Array --------------------  
  Long64_t nEvtKTN = 0;
  if(KATANA)   {
    SetKATANARoot_bt();
    nEvtKTN = kChain -> GetEntries();
    cout << "Number of events in KATANA: " << nEvtKTN << endl;
  }
  //----- Kyoto Array ---------------------
  Long64_t nEvtKyt = 0;
  if(KyotoArry) {

    KyotoRoot = 0;
    if( !SetKyotoMultiplicity()) {
      KyotoRoot = 1;
      if( !SetKyotoArray() )
	KyotoRoot = 2;
    }

    nEvtKyt = kaChain -> GetEntries();
    cout << "Number of events in KyotoArray: " << nEvtKyt << endl;
  }

  //-------------------- event number check --------------------
  Int_t nEntry = nEvtTPC;
  if(nEvtTPC > nEvents && nEvents != 0)
    nEntry = nEvents;

  if(nEvtKTN != 0 && nEntry > nEvtKTN)
    nEntry = nEvtKTN;
  
  if(nEvtKyt != 0 && nEntry > nEvtKyt)
    nEntry = nEvtKyt;

  //--------------------output root--------------------
  OutputTree(nevt);

  if(nevt > 0)  nEntry = nevt;
  // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  //-------------------- event loop --------------------
  // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  TDatime dtime;
  TDatime btime(dtime);
 
  for (Int_t i=0; i<nEntry;i++) {
    Initialize(i);

    if(i%1000 == 0) {
      dtime.Set();
      Int_t ptime = dtime.Get() - btime.Get();
      std::cout << "Process " << i << "/"<< nEntry << " = "
		<< ((Double_t)(i)/(Double_t)nEntry)*100. << " % --->"
                << dtime.AsString() << " ---- "
		<< (Int_t)ptime/60 << " [min] "
                << std::endl;
    }


    // --------------- RIDF --------------
    if(BigRIPS) {
      ribfChain->GetEntry(i);
      bdcChain ->GetEntry(i);
    }

    // --------------- TPC ---------------
    fChain -> GetEntry(i);


    // --------------- KATANA ------------
    if(KATANA && i < nEvtKTN) {
      kChain -> GetEntry(i);

      max_veto = katanaroot->max_veto;
    }
    // --------------- KytoArray ----------
    if(KyotoArry && i < nEvtKyt){
      kaChain -> GetEntry(i);

      if(KyotoRoot == 1){
	bkyhitx -> GetEntry(i);
	bkyhitz -> GetEntry(i);

	for(Int_t ik = 0; ik < (Int_t)kyhitz->size(); ik++){
	  if(kyhitch->at(ik) <= 31) kyotomL++;
	  else kyotomR++;
	}
      }
    }
    

    // --- BDC offset --------------------
    SetBeamOnTarget();
    
    // --- Beam on the target  --------------------
    bmpid = GetBeamPID();
    

    //-------------------- User Analysis --------------------

    Int_t numTracksFromArray = trackArray -> GetEntries();
    ntrack[0] = numTracksFromArray;

    UInt_t mtrack = 0;
    tpcParticle->Clear();

    TIter next(trackArray);
    STRecoTrack *trackFromArray = NULL;

    //    cout << " mtrack " << mtrack << " ?? " << ntrack[0] << " at " << i << endl;
    
    while( (trackFromArray = (STRecoTrack*)next()) ) {
      
      TClonesArray &ptpcParticle = *tpcParticle;

      auto parentvid = trackFromArray->GetVertexID();
      STVertex* vertex;

      if (parentvid > -1) {
	vertex = (STVertex *) vertexArray -> At(parentvid);
	
	STParticle *aParticle = new STParticle();
	aParticle->SetRecoTrack(trackFromArray);

	if( CheckVertex(vertex->GetPos()) )   {

	  //--- Rotate tracks along beam direction ---;
	  if(ProjA > -1000 && ProjB > -1000)
	    aParticle->RotateAlongBeamDirection(ProjA/1000., ProjB/1000.);


	  //--- check origin of the track ---;
	  aParticle->SetTrackAtTarget(vertex->GetPos()); 
	  if( aParticle->GetNDF() > 30 && aParticle->GetDistanceAtVergtex() < 5) {
	    aParticle->SetBestTrackFlag(1);
	    ntrack[2]++;
	  }

	  aParticle->SetTrackID(mtrack);      
	  new(ptpcParticle[mtrack]) STParticle(*aParticle);      
    
	  mtrack++;
	}
      }
    
      //-------------------- end of track LOOP User Analysis --------------------
    }
    
    ntrack[1] = mtrack;
    if(ntrack[1] > 0)
      flw->Fill();
  }

  cout << " Writing " << endl;

  flw->Write();
  
  cout << " Output root file is : " << fout->GetName() << endl;

  if(gROOT->IsBatch()) {
    fout->Close();
    
    cout << " is closed " << endl;
    exit(0);
  }
}

  //##################################################//
void OutputTree(Int_t nmax)
{
  TString sdeb = ".s";
  if(nmax < 0)  sdeb = "";

  TString foutname = "data/run"+sRun+"_f0.v"+sVer+sdeb+".root";

  fout = new TFile(foutname,"recreate");
  flw  = new TTree("flw","Beam and TPC track");
  
  BeamonTarget = new TVector2();

  tpcParticle = new TClonesArray("STParticle",150);

  cout << "Output file is " << foutname << endl;

  //-- output
  flw->Branch("irun",&iRun,"irun/I");
  flw->Branch("aoq",&aoq,"aoq/D");
  flw->Branch("z",&z,"z/D");
  flw->Branch("snbm",&bmpid,"snbm/I");
  flw->Branch("intZ",&intZ,"intZ/I");
  flw->Branch("intA",&intA,"intA/I");
  flw->Branch("BDCtc",&BeamonTarget);
  flw->Branch("bdcax",&bdcax,"bdcax/D");
  flw->Branch("bdcby",&bdcby,"bdcby/D");
  flw->Branch("ProjA",&ProjA,"ProjA/D");
  flw->Branch("ProjB",&ProjB,"ProjB/D");
  flw->Branch("STVertex",&vertexArray);  
  flw->Branch("STParticle",&tpcParticle);
  flw->Branch("ntrack",ntrack,"ntrack[7]/I");
 
  if(KyotoArry){  //Kyoto Array
    flw->Branch("kymult",&kynHit,"kymult/I");
    
    if(KyotoRoot == 1) {
      flw->Branch("kyhit",&kyhitch);           
      flw->Branch("kyxpos",&kyhitx);
      flw->Branch("kyzpos",&kyhitz);
    }
  }
  
  if(KATANA){// KATANA
    flw->Branch("kamult",&katanaM,"kamult/I");
    flw->Branch("max_veto",&max_veto,"max_veto/F");
  }

}

void Initialize(Int_t ievt)
{
  evtid = ievt;

  for (Int_t m = 0; m < 7; m++) ntrack[m] = 0;

  phi.clear();
  gTrack.clear();
  gTgTrack.clear();

  BeamonTarget->SetX(-999.);
  BeamonTarget->SetY(-999.);

}


Int_t GetBeamPID(){
  Int_t pid = 0;
  if(g132Sn){
    if(g132Sn->IsInside(aoq,z) && z < 50.536)
      pid = 132;
  }
  else if(g108Sn->IsInside(aoq,z)){
    if(g108Sn->IsInside(aoq,z))
      pid = 108;
  }

  return pid;
}

void BeamPID()
{

  if(SnA == 132){
    auto gcutFile = new TFile("data/gcut132Sn.ROOT");
    g132Sn=(TCutG*)gcutFile->Get("sigma20");
    gcutFile->Close();
  }

  else if(SnA == 108){
    auto gcutFile = new TFile("data/gcut108Sn.ROOT");
    g108Sn=(TCutG*)gcutFile->Get("sigma20");
    gcutFile->Close();
  }
  // else if(SnA == 124){
  //   auto gcutFile = new TFile("data/gcut124Sn.ROOT");
  //   g124Sn=(TCutG*)gcutFile->Get("g124Sn");
  //   gcutFile->Close();
  // }    
}


void SetDataDirectory()
{
  // upto 9 layer
  //  rootDir = "/cache/scr/spirit/DataAskedByMizuki/Sn132-All-LayerCut90-GC-DS/";
  if(SnA == 132)
    rootDir = "/cache/scr/spirit/DataAskedByMizuki/Sn132-All-LayerCut90-GC-DS-GiordanoCommentOut/";
  else if(SnA == 108)
    rootDir = "/cache/scr/spirit/DataAskedByMizuki/Sn108-All-LayerCut90-GC-DS-GiordanoCommentOut/";
  else if(SnA == 124)
    rootDir = "/cache/scr/spirit/DataAskedByMizuki/Sn124-All-LayerCut90-GC-DS-GiordanoCommentOut/";

}

void SetKATANADirectory()
{
  ktnrootDir = "/data/spdaq01/katana/root/katana/";
}


void SetBeamOnTarget(TVector2 vt)
{
  Double_t Xa = 1.;
  Double_t Xb = 0.;
  Double_t Ya = 1.;
  Double_t Yb = 0.;;

  if(SnA == 108){
    Xa = 1.002;
    Xb = -20.39;
    Ya = 1.002;
    Yb = -226.24;
  }
  else if(SnA == 112){
    Xa = 0.975;
    Xb = -18.90;
    Ya = 0.986;
    Yb = 225.23;
  }
  else if(SnA == 124){
    Xa = 0.935;
    Xb = -15.96;
    Ya = 0.987;
    Yb = -226.07;
  }
  else if(SnA == 132){
    Xa = 0.976;
    Xb = -15.88;
    Ya = 1.002;
    Yb = -225.23;
  }
  
  BeamonTarget->SetX(vt.X()*Xa + Xb);
  BeamonTarget->SetY(vt.Y()*Ya + Yb);
    
   
}

void SetBeamOnTarget()
{
  if(ProjX > -990)
    BeamonTarget->SetX(ProjX);

  if(ProjY > -990) 
    BeamonTarget->SetY(ProjY);
  
}

Bool_t CheckBeamPosition()
{
  TVector2 txy = *BeamonTarget;

  if((txy.X() >= tx_right &&txy.X() <= tx_left) &&  
     (txy.Y() >= ty_btm  && txy.Y() <= ty_top ))   
    return kTRUE;
  else
    return kFALSE;
}

Bool_t CheckVertex(TVector3 vec)
{    if( (vec.Z() > vrt_Zmin     && vec.Z() < vrt_Zmax) &&
	 (vec.X() > trktgt_right && vec.X() < trktgt_left) &&
	 (vec.Y() > trktgt_btm   && vec.Y() < trktgt_top) )	 
    return kTRUE;
  else 
    return kFALSE;
}
				   

Bool_t CheckBDCvsVertexCorrelation(TVector2 vxy)
{

  Double_t diffx = BeamonTarget->X() - vxy.X();
  Double_t diffy = BeamonTarget->Y() - vxy.Y() - beamVy_offset;


  if((abs(diffx) < beamVx_nsig*beamVx_sigma) &&
     (abs(diffy) < beamVy_nsig*beamVy_sigma)) {
    return kTRUE;
  }
  else
    return kFALSE;
}

Bool_t CheckBDCvsTrackCorrelation(TVector3 trackatTarget)
{
  Double_t diffx = BeamonTarget->X() - trackatTarget.X();
  Double_t diffy = BeamonTarget->Y() - trackatTarget.Y() + trackVy_offset;


  if(abs(diffx) < trackVx_nsig*trackVx_sigma &&
     abs(diffy) < trackVy_nsig*trackVy_sigma)
    return kTRUE;
  else
    return kFALSE;
}


void SetKATANARoot()
{
  //----- KATANA data --------------------
  kChain = new TChain("kat");

  SetKATANADirectory();

  kChain->Add(ktnrootDir+"run"+sRun+".katana.root");

  kChain -> SetBranchAddress("nhit",&katanaM);
  kChain -> SetBranchAddress("evt",&event_number);
  kChain -> SetBranchAddress("bitpat",&bitpat,&bbitpat);

}

void SetKATANARoot_bt()
{
  //----- KATANA data --------------------
  //  gSystem->Load("KatanaRoot/KatanaRoot_Load_cpp.so");

  kChain = new TChain("tree");

  SetKATANADirectory();

  kChain->Add(ktnrootDir+"run"+sRun+".root");

  kChain -> SetBranchAddress("Katana",&katanaroot);
  kChain -> SetBranchAddress("STTriggerBox",&triggerbox);

}

Bool_t SetKyotoArray()
{
  kaChain = new TChain("kyotoM");

  TString kytDir = "/cache/scr/spirit/kaneko/rootfile/kyoto/";
  TString kytFile = "run"+sRun+".kyotopos.root";
  if( !gSystem->FindFile(kytDir, kytFile)) 
    return kFALSE;

  kaChain-> Add(kytFile);

  kaChain-> SetBranchAddress("multiplicity",&kynHit);
  kaChain-> SetBranchAddress("hitch",&kyhitch,&bkyhitch);
  kaChain-> SetBranchAddress("hitxpos",&kyhitx,&bkyhitx);
  kaChain-> SetBranchAddress("hitzpos",&kyhitz,&bkyhitz);
  cout << "Set Kyoto Array " << kytFile << endl;

  return kTRUE;
}

Bool_t SetKyotoMultiplicity()
{
  kaChain = new TChain("tree");

  TString kytDir = "/cache/scr/spirit/kaneko/rootfile/kyoto_re/mult/";
  TString kytFile = "run"+sRun+".mult.root";
  if( !gSystem->FindFile(kytDir, kytFile) ) 
    return kFALSE;

  kaChain-> Add(kytFile);

  kaChain-> SetBranchAddress("mult_ctcut",&kynHit);
  cout << "Set Kyoto Multiplicity " << kytFile << endl;
  
  return kTRUE;
}

void SetBigRIPS()
{
  //----- BigRIPS data --------------------

  TString beamDir;
  if(SnA == 132) 
    beamDir  = "/cache/scr/spirit/DataAskedByMizuki/beam.Sn132_all/";
  
  else if(SnA == 108) 
    beamDir  = "/cache/scr/spirit/DataAskedByMizuki/beam.Sn108/";

  auto beamFile = "beam_run"+sRun+".ridf.root";


  if( gSystem->FindFile(beamDir, beamFile)) 
    std::cout << " BigRIPS data : " << beamFile << " is opened." << std::endl;
  
  else {
    std::cout << "ERROR:  RIDF data was not found. " << std::endl;
    exit(0);
  }


  ribfChain = new TChain("TBeam");
  bdcChain = new TChain("TBDC");

      
  ribfChain ->Add(beamFile);

  //----- Set branch addresses.
  ribfChain->SetBranchAddress("aoq",&aoq);
  ribfChain->SetBranchAddress("z",&z);
  ribfChain->SetBranchAddress("beta",&beta);
  ribfChain->SetBranchAddress("brho",&brho);
  ribfChain->SetBranchAddress("isGood",&isGood);
  ribfChain->SetBranchAddress("intZ",&intZ);
  ribfChain->SetBranchAddress("intA",&intA);

  bdcChain -> Add(beamFile);  
  bdcChain->SetBranchAddress("bdcax",&bdcax);
  bdcChain->SetBranchAddress("bdcby",&bdcby);
  bdcChain->SetBranchAddress("ProjA",&ProjA);
  bdcChain->SetBranchAddress("ProjB",&ProjB);
  bdcChain->SetBranchAddress("ProjX",&ProjX);
  bdcChain->SetBranchAddress("ProjY",&ProjY);
  bdcChain->SetBranchAddress("ProjZ",&ProjZ);
  bdcChain->SetBranchAddress("ProjP",&ProjP);
  bdcChain->SetBranchAddress("ProjPX",&ProjPX);
  bdcChain->SetBranchAddress("ProjPY",&ProjPY);
  bdcChain->SetBranchAddress("ProjPZ",&ProjPZ);

}

void SetTPC()
{
  fChain = new TChain("cbmsim");
  SetDataDirectory();
   
  Int_t i = 0;
  while(kTRUE){
    
    TString recoFile = Form("run"+sRun+"_s%d.reco.v1.04.root",i);
    cout << " recoFile " << rootDir+recoFile << endl;
    
    if(gSystem->FindFile(rootDir,recoFile)){
      fChain -> Add(recoFile);
    }
    else
      break;
    i++;

  }

  fChain -> SetBranchAddress("STRecoTrack", &trackArray, &brtrackArray);
  fChain -> SetBranchAddress("STVertex"   , &vertexArray,&brvertexArray);
}

Bool_t DefineVersion()
{
  Bool_t bfound = kFALSE;

  Ssiz_t end = sVer.First(".");

  cout << " end " << end << endl;
  if( end == -1) {

    iVer = atoi(sVer);

    bfound = kTRUE;
  }
  //  TString ver = sVer + ".";
  
  // for ( Int_t i = 0; i < 2; i++) {
  //   if( ver.First(".") < 0 ) break;

  //   Ssiz_t end = ver.First(".")  ;
  //   TString ver1 = ver(0, end);

  //   ver = ver(end+1, ver.Length());

  //   iVer[i] = atoi(ver1);

  //   if(i==1) bfound = kTRUE;
  // }

  
  
  if(!bfound)
    cout << " missing version number v# " << endl;

  return bfound;

}

