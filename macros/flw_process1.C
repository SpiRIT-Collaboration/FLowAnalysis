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
    std::cout << "Plase type " << std::endl;
    std::cout << "$ RUN=#### VER=# root flw_process1.C" << std::endl;
    exit(0);
  }

  iRun = atoi(sRun);

  SnA = 0;
  if(iRun >= 2174 && iRun <= 2509)
    SnA = 108;
  else if( iRun >= 2520 && iRun <= 2653)
    SnA = 112;
  else if( iRun >= 2836 && iRun <= 3039)
    SnA = 132;
  else if( iRun >= 3058 && iRun <= 3184)
    SnA = 124;

  STPC     = (Bool_t)atoi(gSystem->Getenv("STPC"));
  BigRIPS  = (Bool_t)atoi(gSystem->Getenv("BIGRIPS"));
  KyotoArry= (Bool_t)atoi(gSystem->Getenv("KYOTOARRY"));
  KATANA   = (Bool_t)atoi(gSystem->Getenv("KATANA"));
  NeuLAND  = (Bool_t)atoi(gSystem->Getenv("NEULAND"));
  
  std::cout << "Included data sets -> " ; 
  if(STPC)      std::cout << " SpiRIT-TPC"   ;
  if(BigRIPS)   std::cout << " & BigRIPS "   ;
  if(KyotoArry) std::cout << " & KyotoArry " ;
  if(KATANA)    std::cout << " & KATANA "    ;
  if(NeuLAND)   std::cout << " & NeuLAND "   ;
  std::cout << std::endl;

  if(STPC) {
    TString bbfitter = gSystem->Getenv("STBBFITTER");
    auto fitFile = new TFile(bbfitter);
    auto fit = (TF1 *) fitFile -> Get("fit_proton");

    fitterPara[0] = fit -> GetParameter(0);
    fitterPara[1] = fit -> GetParameter(1);

    if( fitFile != NULL)
      cout << "BetheBloch fitter is loaded. "
	   << " para 0 " << fitterPara[0]
	   << " 1 " << fitterPara[1]
	   << endl; 

    fitFile->Close();
    delete fitFile;
  }
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


  Int_t nEntry = 0;

  //----- TPC data ------------------------
  Long64_t nEvtTPC = 0;
  if( STPC ) {
    SetTPC();

    if(fChain) nEvtTPC = fChain -> GetEntries();
    std::cout << "Number of events in TPC: " << nEvtTPC << std::endl;
    if(nEvtTPC == 0) return;

    nEntry = nEvtTPC;
  }

  //----- BigRIPS data --------------------
  Long64_t nEvents = 0;
  if(BigRIPS) {
    SetBigRIPS();
    if( ribfChain) {
      nEvents   = ribfChain->GetEntries();
      Int_t nEventsBDC = bdcChain->GetEntries();

      if(nEvents != nEventsBDC) {
	std::cout << "Inconsistent event number in bigRIS RIDF (quit)" << std::endl;
	exit(0);
      }
    }
    std::cout << "Number of events in RIDF: " << nEvents << std::endl;

    if( nEntry == 0 ) nEntry = nEvents;

  }
  //----- KATANA Array --------------------  
  Long64_t nEvtKTN = 0;
  if(KATANA)   {
    SetKATANARoot_bt();
    nEvtKTN = kChain -> GetEntries();
    std::cout << "Number of events in KATANA: " << nEvtKTN << std::endl;

    if( nEntry == 0 ) nEntry = nEvtKTN;
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
    std::cout << "Number of events in KyotoArray: " << nEvtKyt << std::endl;

    if( nEntry == 0 ) nEntry = nEvtKyt;
  }

  //---- NeuLAND ---------------------------
  Long64_t nEvtNL = 0;
  if(NeuLAND) {
    if(SetNeuLANDRoot()) {
      nEvtNL = nlChain -> GetEntries();
      std::cout << "Number of events in NeuLAND: "<< nEvtNL << std::endl;

      LoadNeuLANDPID();
    }
    else{
      std::cout << "No NeuLAND data is found. " << std::endl;
      NeuLAND=kFALSE;
    }

    if( nEntry == 0 ) nEntry = nEvtNL;
  }

  //-------------------- event number check --------------------
  if( nEntry == 0 ) {
    std::cout << "No data wwas found so stopped...." << std::endl;
    return;
  }

  if( nEvtTPC != 0 && nEntry > nEvtTPC )
    nEntry = nEvents;

  if( nEvtKTN != 0 && nEntry > nEvtKTN)
    nEntry = nEvtKTN;
  
  if( nEvtKyt != 0 && nEntry > nEvtKyt)
    nEntry = nEvtKyt;

  if( nEvtNL  != 0 && nEntry > nEvtNL)
    nEntry = nEvtNL;


  std::cout << "Number of events will be processed ---->  " << nEntry << std::endl;

  //  nEntry = 100;
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


    // --------------- TPC ---------------
    if(STPC)
      fChain -> GetEntry(i);

    // --------------- RIDF --------------
    if(BigRIPS) {
      ribfChain->GetEntry(i);
      bdcChain ->GetEntry(i);
    }

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
    
    // --------------- NeuLAND ---------- 
    if(NeuLAND && i < nEvtNL)
      nlChain -> GetEntry(i);

    // --- BDC offset --------------------
    SetBeamOnTarget();
    
    // --- Beam on the target  --------------------
    bmpid = GetBeamPID();
    


    //-------------------- User Analysis --------------------
    //

    //------- TPC ------
    if(STPC) {
      Int_t numTracksFromArray = trackArray -> GetEntries();
      ntrack[0] = numTracksFromArray;

      UInt_t mtrack = 0;
      tpcParticle->Clear();

      TIter next(trackArray);
      STRecoTrack *trackFromArray = NULL;

    
      while( (trackFromArray = (STRecoTrack*)next()) ) {
      
	TClonesArray &ptpcParticle = *tpcParticle;

	auto parentvid = trackFromArray->GetVertexID();


	STVertex* vertex = NULL;
	if (parentvid > -1) {
	  vertex = (STVertex *) vertexArray -> At(parentvid);
	
	  STParticle *aParticle = new STParticle();
	  aParticle->SetRecoTrack(trackFromArray);
	  
	  //--- Set BetheBloch mass 
	  aParticle->SetBetheBlochMass(fitterPara);

	  //--- Rotate tracks along beam direction ---;
	  if(ProjA > -1000 && ProjB > -1000)
	    aParticle->RotateAlongBeamDirection(ProjA/1000., ProjB/1000.);


	  //--- check origin of the track ---;
	  aParticle->SetVertex(vertex); 

	  if( aParticle->GetMomentumAtTarget().Mag() == 0)
	    aParticle->SetMaxMomentumFlag(0);

	  else if( CheckVertex(aParticle) )   {


	    aParticle->SetBestTrackFlag(1);
	    ntrack[2]++;


	    //--- Set track quality flag ---;
	    if( aParticle->GetDistanceAtVertex() > 20 )
	      aParticle->SetDistanceAtVertexFlag(0);

	    if( aParticle->GetNDF() < 30)
	      aParticle->SetNDFFlag(0);
	    
	    if( aParticle->GetP() > 2500 )
	      aParticle->SetMaxMomentumFlag(0);
	    
	    if( aParticle->GetRotatedMomentum().Theta() >= 0.8 )
	      aParticle->SetMaxThetaFlag(0);

	    if( aParticle->GetdEdx() > 1000 )
	      aParticle->SetMaxdEdxFlag(0);

	  }

	  if( aParticle->GetBestTrackFlag() ) {
	    ntrack[3]++;
	  }


	  aParticle->SetTrackID(mtrack);      
	  new(ptpcParticle[mtrack]) STParticle(*aParticle);      
	  mtrack++;

	}
    
      }
      ntrack[1] = mtrack;
    }
    ////  ---- end of TPC  ----

    // ----- NeuLAND -----
    if(NeuLAND) {//NeuLAND        

      Int_t numnlCluster = nlcluster->GetEntries();
      nhitnl[0] = numnlCluster;
      //      cout << " number of cluster "<< numnlCluster << endl;
      TIter nlnext(nlcluster);
      STNeuLANDCluster *nlFromCluster = NULL;
    
      UInt_t veto_torq = 0;
      if (SnA == 108)
	veto_torq = 1;

      while( (nlFromCluster = (STNeuLANDCluster*)nlnext() ) ){
 	if( nlFromCluster->GetVetoHitAll(veto_torq) == 0) 
	  nhitnl[1]++;
	else
	  nhitnl[5]++;	  
	

	if( nlFromCluster->GetVetoHitOne(veto_torq) == 0) 
	  nhitnl[2]++;
	else
	  nhitnl[6]++;

	if( nlFromCluster->GetVetoHitMid(veto_torq) == 0) 
	  nhitnl[3]++;
	else
	  nhitnl[7]++;
	  
	
	auto nlPID = GetNeuLANDPID(nlFromCluster->GetEdep(), nlFromCluster->GetTOF(), nlFromCluster->GetVetoHitLoose(veto_torq));
	if(nlPID > 0){
	  nlFromCluster->SetMass(nlPID);

	  if( nlPID == 2112) 
	    nhitnl[4]++;
	  else 
	    nhitnl[8]++;
	}
	else{ // if NeuLAND PID is not found

	  if( nlFromCluster->GetVetoHitLoose(veto_torq) == 0)
	    nhitnl[4]++;
	  else
	    nhitnl[8]++;

	}
      }
    }
    //// ---- endof NeuLad -----


    //-------------------- end of track LOOP User Analysis --------------------
    if(ntrack[1] > 0 || nhitnl[0] > 0)
      flw->Fill();
  }

  std::cout << " Writing " << std::endl;

  flw->Write();
  
  std::cout << " Output root file is : " << fout->GetName() << std::endl;

  if(gROOT->IsBatch()) {
    fout->Close();
    
    std::cout << " is closed " << std::endl;
    exit(0);
  }
}

//##################################################//
void OutputTree(Int_t nmax)
{
  TString sdeb = ".s";
  if(nmax < 0)  sdeb = "";

  
  TString foutname = "data/run"+sRun+"_f0.v"+sVer+sdeb+".root";

  if( !STPC && NeuLAND ) 
    foutname = "data/run"+sRun+"_nl.v"+sVer+sdeb+".root";


  fout = new TFile(foutname,"recreate");
  flw  = new TTree("flw","Beam and TPC track");
  
  BeamonTarget = new TVector2();

  std::cout << "Output file is " << foutname << std::endl;

  //-- output
  flw->Branch("irun",&iRun,"irun/I");

  if(BigRIPS){
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
  }

  if(STPC) {
    tpcParticle = new TClonesArray("STParticle",120);

    flw->Branch("STVertex",&vertexArray);  
    flw->Branch("STParticle",&tpcParticle);
    flw->Branch("ntrack",ntrack,"ntrack[7]/I");
  }

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

  if(NeuLAND) {//NeuLAND
    flw->Branch("nhitnl",nhitnl,"nhitnl[9]/I");
    flw->Branch("STNeuLANDHit",   &nlhit);
    flw->Branch("STNeuLANDCluster",&nlcluster);
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

  for (Int_t m = 0; m < 9; m++) nhitnl[m] = 0;
}


Int_t GetBeamPID(){
  Int_t pid = 0;
  if(g132Sn){
    if(g132Sn->IsInside(aoq,z) && z < 50.536)
      pid = 132;
  }

  if(g108Sn){
    if(g108Sn->IsInside(aoq,z))
      pid = 108;
  }

  if(g124Sn){
    if(g124Sn->IsInside(aoq,z))
      pid = 124;
  }

  if(g112Sn){
    if(g112Sn->IsInside(aoq,z))
      pid = 112;
  }

  return pid;
}

void BeamPID()
{

  if(SnA == 132){
    auto gcutFile = new TFile("data/gcut132Sn.ROOT");
    g132Sn=(TCutG*)gcutFile->Get("sigma20");
    gcutFile->Close();

    if(g132Sn != NULL) std::cout << "gcut132Sn with sigma20 from data/gcut132Sn.ROOT" <<std::endl;
    else std::cout << "data/gcut132Sn.ROOT was not found" << std::endl;

  }

  else if(SnA == 108){
    auto gcutFile = new TFile("data/gcut108Sn.ROOT");
    g108Sn=(TCutG*)gcutFile->Get("sigma20");
    gcutFile->Close();

    if(g108Sn != NULL) std::cout << "gcut108Sn with sigma20 from data/gcut108Sn.ROOT" <<std::endl;
    else std::cout << "data/gcut108Sn.ROOT was not found" << std::endl;

  }
  else if(SnA == 124){
    auto gcutFile = new TFile("data/gcut124Sn.ROOT");
    g124Sn=(TCutG*)gcutFile->Get("sigma20");
    gcutFile->Close();

    if(g124Sn != NULL) std::cout << "gcut124Sn with sigma20 from data/gcut124Sn.ROOT" <<std::endl;
    else std::cout << "data/gcut124Sn.ROOT was not found" << std::endl;
  }
  else if(SnA == 112){
    auto gcutFile = new TFile("data/gcut112Sn.ROOT");
    g112Sn=(TCutG*)gcutFile->Get("sigma20");
    gcutFile->Close();
    
    if(g112Sn != NULL) std::cout << "gcut112Sn with sigma20 from data/gcut112Sn.ROOT" <<std::endl;
    else std::cout << "data/gcut112Sn.ROOT was not found" << std::endl;
  }    
}

Int_t GetNeuLANDPID(Double_t x, Double_t y, Int_t vhit)
{
  if( gcutNLNeutron == NULL ) return -1;

  Int_t pid = 0;
  if( gcutNLNeutron->IsInside(x,y) && vhit == 0)
    pid = 2112;
  else if( gcutNLProton->IsInside(x,y) )
    pid = 2212;
  else if( gcutNLDeuteron->IsInside(x,y) )
    pid = 1000010020;
  else if( gcutNLTrition->IsInside(x,y) )
    pid = 1000010030;
  else
    pid = 2000000000;

  return pid;
}

void LoadNeuLANDPID()
{
  TString gfname = "data/gcutNLNeutron.root";

  TFile gcutFilen(gfname);
  if(gcutFilen.IsOpen()) {
    std::cout << gfname << " is loaded." << std::endl;
    gcutNLNeutron  = (TCutG*)gcutFilen.Get("gcutNLNeutron");
    gcutFilen.Close();
  }
  else
    std::cout << gfname << " is not found." << std::endl;
    
  
  
  gfname = "data/gcutNLProton.root";
  TFile gcutFilep(gfname);
  if(gcutFilep.IsOpen()) {
    std::cout << gfname << " is loaded." << std::endl;
    gcutNLProton   = (TCutG*)gcutFilep.Get("gcutNLProton");
    gcutFilep.Close();
  }
  else
    std::cout << gfname << " is not found." << std::endl;


  gfname = "data/gcutNLDeuteron.root";
  TFile gcutFiled(gfname);
  if(gcutFiled.IsOpen()) {
    std::cout << gfname << " is loaded." << std::endl;
    gcutNLDeuteron = (TCutG*)gcutFiled.Get("gcutNLDeuteron");
    gcutFiled.Close();
  }
  else
    std::cout << gfname << " is not found." << std::endl;



  gfname = "data/gcutNLTriton.root";
  TFile gcutFilet(gfname);
  if(gcutFilet.IsOpen()) {
    std::cout << gfname << " is loaded." << std::endl;
    gcutNLTrition  = (TCutG*)gcutFilet.Get("gcutNLTriton");
    gcutFilet.Close();
  }
  else
    std::cout << gfname << " is not found." << std::endl;

}

void SetDataDirectory()
{

  if( gSystem -> Getenv("ST132DIR") != NULL ) {
    if(SnA == 132)
      rootDir = gSystem -> Getenv("ST132DIR");
    else if(SnA == 108)
      rootDir = gSystem -> Getenv("ST108DIR");
    else if(SnA == 124)
      rootDir = gSystem -> Getenv("ST124DIR");
    else if(SnA == 112)
      rootDir = gSystem -> Getenv("ST112DIR");
  }

  else
  //for v6.
    rootDir = gSystem -> Getenv("STTPCDIR");

}

void SetKATANADirectory()
{
  ktnrootDir = gSystem -> Getenv("STKATANADIR");
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

Bool_t CheckVertex(STParticle *aPart)
{   
  auto vec = aPart->GetVertex(); 

  if( abs( vec.Z() + 12.9 ) <= 3. &&
      abs( vec.X() ) <= 15.       &&
      abs( vec.Y() + 226.06 ) <= 20. )

    aPart->SetVertexAtTargetFlag(1);
  else
    aPart->SetVertexAtTargetFlag(0);
  

  return (Bool_t)aPart->GetVertexAtTargetFlag();
}
				   

Bool_t CheckBDCvsVertexCorrelation(TVector2 vxy)
{

  Double_t diffx = BeamonTarget->X() - vxy.X();
  Double_t diffy = BeamonTarget->Y() - vxy.Y() - beamVy_offset;


  if((TMath::Abs(diffx) < beamVx_nsig*beamVx_sigma) &&
     (TMath::Abs(diffy) < beamVy_nsig*beamVy_sigma)) {
    return kTRUE;
  }
  else
    return kFALSE;
}

Bool_t CheckBDCvsTrackCorrelation(TVector3 trackatTarget)
{
  Double_t diffx = BeamonTarget->X() - trackatTarget.X();
  Double_t diffy = BeamonTarget->Y() - trackatTarget.Y() + trackVy_offset;


  if(TMath::Abs(diffx) < trackVx_nsig*trackVx_sigma &&
     TMath::Abs(diffy) < trackVy_nsig*trackVy_sigma)
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

  TString kytDir = gSystem -> Getenv("STKYOTODIR");
  TString kytFile = "run"+sRun+".kyotopos.root";
  if( !gSystem->FindFile(kytDir, kytFile)) 
    return kFALSE;

  kaChain-> Add(kytFile);

  kaChain-> SetBranchAddress("multiplicity",&kynHit);
  kaChain-> SetBranchAddress("hitch",&kyhitch,&bkyhitch);
  kaChain-> SetBranchAddress("hitxpos",&kyhitx,&bkyhitx);
  kaChain-> SetBranchAddress("hitzpos",&kyhitz,&bkyhitz);
  std::cout << "Set Kyoto Array " << kytFile << std::endl;

  return kTRUE;
}

Bool_t SetKyotoMultiplicity()
{
  kaChain = new TChain("tree");

  TString kytDir = gSystem -> Getenv("STKYMLTDIR");
  TString kytFile = "run"+sRun+".mult.root";
  if( !gSystem->FindFile(kytDir, kytFile) ) 
    return kFALSE;

  kaChain-> Add(kytFile);

  kaChain-> SetBranchAddress("mult_ctcut",&kynHit);
  std::cout << "Set Kyoto Multiplicity " << kytFile << std::endl;
  
  return kTRUE;
}

Bool_t SetNeuLANDRoot()
{
  nlChain = new TChain("nl");
  
  TString nlDir  = gSystem -> Getenv("STNLDIR");
  TString nlFile = "neuland_run"+sRun+".root";

  if( !gSystem->FindFile(nlDir, nlFile) )
    return kFALSE;
      
  nlChain-> Add(nlFile);
    
  nlChain->SetBranchAddress("STNeuLANDHit",&nlhit);
  nlChain->SetBranchAddress("STNeuLANDCluster",&nlcluster);

  std::cout << "Set NeuLAND " << nlFile << std::endl;
  return kTRUE;
}

void SetBigRIPS()
{
  //----- BigRIPS data --------------------

  TString beamDir;
  if(SnA == 132) 
    beamDir = gSystem -> Getenv("STBEAM132");
  else if(SnA == 108) 
    beamDir = gSystem -> Getenv("STBEAM108");
  else if(SnA == 124)
    beamDir = gSystem -> Getenv("STBEAM124");
  else if(SnA == 112)
    beamDir = gSystem -> Getenv("STBEAM112");

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
  TString fileversion = gSystem->Getenv("STVERSION");
   
  Int_t i = 0;
  while(kTRUE && fileversion != ""){

    TString recoFile = Form("run"+sRun+"_s%d.reco."+fileversion+".root",i);
    std::cout << " recoFile " << rootDir+recoFile << std::endl;
    
    if(gSystem->FindFile(rootDir,recoFile)){
      fChain -> Add(recoFile);
    }
    else
      break;
    i++;

  }

  fChain -> SetBranchAddress("STRecoTrack", &trackArray);
  fChain -> SetBranchAddress("STVertex"   , &vertexArray);
}

Bool_t DefineVersion()
{
  Bool_t bfound = kFALSE;

  Ssiz_t end = sVer.First(".");

  std::cout << " end " << end << std::endl;
  if( end == -1) {

    iVer = atoi(sVer);

    bfound = kTRUE;
  }
  
  
  if(!bfound)
    std::cout << " missing version number v# " << std::endl;

  return bfound;

}

