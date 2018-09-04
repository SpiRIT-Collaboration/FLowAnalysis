#include "mc_process1.h"
//----------------------------------------------------------------------
// Macro: mc_process1
//  (c) Mizuki Kurata-Nishimura
//   2018 Aug. 23 
//----------------------------------------
//
//----------------------------------------------------------------------
void Setup()
{
  sVer = gSystem -> Getenv("VER");

  if( sRun=="" || sVer=="" || !DefineVersion() ) {
    std::cout << "Plase type " << std::endl;
    std::cout << "$ RUN=#### VER=# root flw_process1.C" << std::endl;
    exit(0);
  }

}

void mc_process1(Int_t nevt = -1)
{
  //////////////////////////////////////////////////////////
  // The calibrated beam_run#.ridf.root can be loaded.
  //////////////////////////////////////////////////////////

  Setup();

  //Reset ROOT and connect tree file
  gROOT->Reset();


  Int_t nEntry = 0;

  //----- TPC data ------------------------
  Long64_t nEvtTPC = 0;
  SetTPC();
  
  if(fChain) nEvtTPC = fChain -> GetEntries();
  std::cout << "Number of events in TPC: " << nEvtTPC << std::endl;
  if(nEvtTPC == 0) return;
  
  nEntry = nEvtTPC;


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

    //-------------------- User Analysis --------------------
    //

    //------- TPC ------
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
      STMCEventHeader* mcheader = NULL;

      if (parentvid > -1) {
	vertex   = (STVertex *) vertexArray -> At(parentvid);
	TVector2 beamVector = (STMCEventHeader*) mceventHeader -> GetBeamAngle();
	
	
	STParticle *aParticle = new STParticle();
	aParticle->SetRecoTrack(trackFromArray);

	//--- Rotate tracks along beam direction ---;
	  aParticle->RotateAlongBeamDirection(ProjA/1000., ProjB/1000.);


	//--- check origin of the track ---;
	aParticle->SetTrackAtTarget(vertex->GetPos()); 

	if( aParticle->GetMomentumAtTarget().Mag() == 0)
	  aParticle->SetMaxMomentumFlag(0);

	else if( CheckVertex(aParticle) )   {

	  aParticle->SetBestTrackFlag(1);
	  ntrack[2]++;


	  //--- Set track quality flag ---;
	  if( aParticle->GetDistanceAtVertex() >= 5 )
	    aParticle->SetDistanceAtVertexFlag(0);

	  if( aParticle->GetNDF() <= 30)
	    aParticle->SetNDFFlag(0);
	    
	  if( aParticle->GetP() >= 2500 )
	    aParticle->SetMaxMomentumFlag(0);
	    
	  if( aParticle->GetRotatedMomentum().Theta() >= 0.8 )
	    aParticle->SetMaxThetaFlag(0);

	  if( aParticle->GetdEdx() > 1000 )
	    aParticle->SetMaxdEdxFlag(0);

	}

	if( aParticle->GetBestTrackFlag() )
	  ntrack[3]++;


	aParticle->SetTrackID(mtrack);      
	new(ptpcParticle[mtrack]) STParticle(*aParticle);      
	mtrack++;

      }
    
    }
    ntrack[1] = mtrack;

    ////  ---- end of TPC  ----

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

  
  TString foutname = "data/mc_"+recoFile+"_f0.v"+sVer+sdeb+".root";


  fout = new TFile(foutname,"recreate");
  flw  = new TTree("flw","Beam and TPC track");

  std::cout << "Output file is " << foutname << std::endl;

  //-- output
  flw->Branch("irun",&iRun,"irun/I");

  flw->Branch("aoq",&aoq,"aoq/D");
  flw->Branch("z",&z,"z/D");
  flw->Branch("snbm",&bmpid,"snbm/I");
  flw->Branch("rpangle",&rpangle,"rpangle/D");

  tpcParticle = new TClonesArray("STParticle",120);

  flw->Branch("STVertex",&vertexArray);  
  flw->Branch("STParticle",&tpcParticle);
  flw->Branch("ntrack",ntrack,"ntrack[7]/I");

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
  auto vec = aPart->GetTrackAtTarget(); 

  if( vec.Z() < vrt_Zmin ||  vec.Z() > vrt_Zmax )
    aPart->SetVertexZAtTargetFlag(0);


  if( (vec.X() < trktgt_right || vec.X() > trktgt_left) ||
      (vec.Y() < trktgt_btm   || vec.Y() > trktgt_top) )	 
    aPart->SetVertexAtTargetFlag(0);
    
  if( aPart->GetVertexAtTargetFlag() * aPart->GetVertexZAtTargetFlag() )
    return kTRUE;
  else 
    return kFALSE;
}
				   

void SetTPC()
{
  fChain = new TChain("cbmsim");

  TString rootDir     = gSystem->Getenv("MCDIR");
  TString fileversion = gSystem->Getenv("STVERSION");
  TString mcsystem    = gSystem->Getenv("MCSYSTEM");
   
  Int_t i = 0;
  while(kTRUE && fileversion != ""){

    recoFile = Form( mcsystem + "_%d.reco.develop."+fileversion,i);
    std::cout << " recoFile " << rootDir+recoFile << std::endl;
    
    if(gSystem->FindFile(rootDir,recoFile+".root")){
      fChain -> Add(recoFile);
    }
    else
      break;
    i++;

  }

  fChain -> SetBranchAddress("STRecoTrack", &trackArray);
  fChain -> SetBranchAddress("STVertex"   , &vertexArray);
  fChain -> SetBranchAddress("STMCEventHeader", &mceventHeader);
  fChain -> SetBranchAddress("STMCTriggerResponse",&mctriggerResponse);
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

