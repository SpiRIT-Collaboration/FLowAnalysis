#include "STSpiRITTPCTask.hh"

ClassImp(STSpiRITTPCTask);

STSpiRITTPCTask::STSpiRITTPCTask() 
  : fIsPersistence(kTRUE),
    fIsFlowAnalysis(kTRUE),
    fIsFlowCorrection(kTRUE),
    fIsBootStrap(kFALSE),
    fIsSubeventAnalysis(kTRUE),
    selReactionPlanef(10),
    fEventID(-1),
    massFitter(new STMassFunction())
{

  fLogger = FairLogger::GetLogger();
  SetVerbose(1);

  fRunAna = FairRunAna::Instance();
  iRun  = fRunAna->getRunId();


  //------------------------//
  // initial setup
  fChain = nullptr;
}

STSpiRITTPCTask::~STSpiRITTPCTask()
{

  delete db;
  delete fChain;  //!
  delete eventHeader    ;  //!
  delete trackArray     ;  //!
  delete trackVAArray   ;  //!
  delete vertexArray    ;  //!
  delete vertexVAArray  ;  //!
  delete vertexBDCArray ;  //!

  delete tpcParticle;  //!

  //--- mass fitter
  delete massFitter;

  delete flowAnalysis; //!
  delete fflowinfo;    //!

  for(UInt_t i = 0; i < 2; i++) {
    if( aflowcorrArray[i] != nullptr )
      delete aflowcorrArray[i];
  }
}


void STSpiRITTPCTask::SetPersistence(Bool_t value) { fIsPersistence = value; }


InitStatus STSpiRITTPCTask::Init()
{
  LOG(INFO) << "STSpiRITTPCTask::InitTask() is called " << FairLogger::endl;

  beginTime.Copy(dtime);

  fRootManager = FairRootManager::Instance();
  if ( fRootManager == nullptr) {
    LOG(ERROR) << "Cannot find RootManager!" << FairLogger::endl;
    return kERROR;
  }

  if( fRunAna == nullptr ) {
    LOG(ERROR) << "STSpiRITTPCTask ==> Not defined. "  << FairLogger::endl;
    return kERROR;
  }

  auto fBDCArray = (TClonesArray *)fRootManager->GetObject("STBDC");
  if( fBDCArray == nullptr ) {
    LOG(ERROR) << "Register STBDC in advance! " << FairLogger::endl;
    return kERROR;
  }
  else
    LOG(INFO) << "STBDC is found. " << FairLogger::endl;
  tpcParticle = new TClonesArray("STParticle",100);
  fRootManager -> Register("STParticle","flow",tpcParticle, fIsPersistence);

  if( fIsFlowAnalysis ) {
    flowAnalysis = new TClonesArray("STFlowInfo",1);
    fRootManager -> Register("STFlow","flow",flowAnalysis, fIsPersistence);

    fflowinfo = new STFlowInfo();
    fflowtask = new STFlowTask(kTRUE, kTRUE, kFALSE); //(flattening, subevent, bootstrap)
    fIsFlowAnalysis = fflowtask->Init(iRun, sVer);
  }


  if( !SetupInputDataFile() ) {
    LOG(ERROR) << "STSpiRITTPCTask:: Cannot open input files" << FairLogger::endl;
    return kERROR;
  }

  if( !SetupParameters() ) {
    LOG(ERROR) << "STSpiRITTPCTask:: Cannot open parameter files" << FairLogger::endl;
    return kERROR;
  }

  return kSUCCESS;
}


Bool_t STSpiRITTPCTask::SetupParameters()
{
  Bool_t fstatus = kTRUE;

  //--- angle dependent mass fitter // 2019/02/19
  

  //--- angle dependeing mass fitter //
  fstatus = massFitter->SetFunction("/cache/scr/spirit/mizuki/Bethe-Bloch_Fitter/mk_NewFitter_20190111/",
				    "BBFitter.root",
				    "f1IterBBProtonRotate_" + STRunToBeamA::GetBeamSnA(iRun) );
  if( !fstatus ) LOG(ERROR) << " BB Mass fittter function was not loaded. " << FairLogger::endl;

  fstatus = massFitter->SetMassFitFunction("/cache/scr/spirit/mizuki/Bethe-Bloch_Fitter/mk_NewFitter_20190111/",
					   "MassFitter.root",
					   "f1IterMassFitRotate_" + STRunToBeamA::GetBeamSnA(iRun) );
  if( ! fstatus ) LOG(ERROR) << " BB Mass gaus fit function was not loaded. " << FairLogger::endl;

  //------------------------------

  if( kFALSE ) {
    //--- Single mass fitter
    TString bbfitter = gSystem->Getenv("STBBFITTER");
    auto fitFile = new TFile(bbfitter);
    auto fit = (TF1 *) fitFile -> Get("fit_proton");

    Double_t fitterPara[2];    
    if( fit ) {

      fitterPara[0] = fit -> GetParameter(0);
      fitterPara[1] = fit -> GetParameter(1);

      LOG(INFO) << "single BetheBloch fitter is loaded from a file. "
		<< " para 0 " << fitterPara[0]
		<< " 1 " << fitterPara[1]
		<< FairLogger::endl;
    }
    else {
      fitterPara[0] = -0.0040786335;
      fitterPara[1] = 10007.475;

      LOG(INFO) << "single BetheBloch fitter is setup. "
		<< " para 0 " << fitterPara[0]
		<< " 1 " << fitterPara[1]
		<< FairLogger::endl;
    }

    fitFile->Close();
    delete fitFile;
    //------------------------------

    //--- Load theoretical cluster number
    string dir = gSystem->Getenv("SPIRITROOT");
    //  string dir1 = "../../"+dir+"/ana/Momentum.config";
    string dir1 = "../../"+dir+"/ana/Momentum_tb_edge_ellipsoid_cut_clusternum_DB_4GeV_theta90_phi180.config";

    LOG(INFO) << "db file " << dir1 << endl;

    db = new ST_ClusterNum_DB();
    db->Initial_Config( dir1 );
    //  dir1 = "../../"+dir+"/ana/f1_DB_ClusterNum.root";
    dir1 ="../../"+ dir+"/ana/f1_tb_edge_ellipsoid_cut_clusternum_DB_theta90_phi180.root";
    db->ReadDB( dir1 );

    double Momentum_Range_Plus[2]  = {50,4000};
    double Momentum_Range_Minus[2] = {50,4000};
    db->Set_MomentumRange_Plus(Momentum_Range_Plus);
    db->Set_MomentumRange_Minus(Momentum_Range_Minus);
    //------------------------------
  }

  return fstatus;
}

Bool_t STSpiRITTPCTask::SetupInputDataFile() 
{
  LOG(INFO) << "STSpiRITTPCTask::SetupInputDataFile() is called " << FairLogger::endl;

  fChain = new TChain("cbmsim");

  UInt_t i = 0;
  while(kTRUE){

    TString recoFile = Form("run%04d_s%d.reco."+tVer+".root",iRun,i);


    if(gSystem->FindFile(rootDir,recoFile))
      fChain -> Add(recoFile);

    else 
      break;
    
    LOG(INFO) << i << " recoFile " << rootDir+recoFile << FairLogger::endl;
    i++;

  }

  if( i == 0 ) {
    LOG(ERROR) << "Cannot find TPC data! "  <<  FairLogger::endl;
    return kFALSE;
  }

  nEntry = fChain->GetEntries();
  LOG(INFO) << i << " files were chained. " << nEntry << " events" << FairLogger::endl;
  
  prcEntry = nEntry;


  fChain -> SetBranchAddress("STRecoTrack",   &trackArray);
  fChain -> SetBranchAddress("VATracks",      &trackVAArray);
  fChain -> SetBranchAddress("STVertex"   ,   &vertexArray);
  fChain -> SetBranchAddress("VAVertex"   ,   &vertexVAArray);
  fChain -> SetBranchAddress("BDCVertex"  ,   &vertexBDCArray);

  //  fChain->Print();

  return kTRUE;
}


void STSpiRITTPCTask::Exec(Option_t *opt)
{

  LOG(DEBUG) << "STSpiRITTPCTask::Exec() is called " << opt << FairLogger::endl;  
  Clear();

  fEventID++;
  
  UInt_t pdev = prcEntry/200;
  if(prcEntry < 200) pdev = 1;

  if(fEventID%pdev == 0) {
    dtime.Set();
    Int_t ptime = dtime.Get() - beginTime.Get();
    LOG(INFO) << "Process " 
	      << setw(4) << ((Double_t)fEventID/(Double_t)prcEntry)*100. << " % = "
	      << setw(8) << fEventID << "/"<< prcEntry 
	      << "--->"
	      << dtime.AsString() << " ---- "
	      << (Int_t)ptime/60 << " [min] "
	      << FairLogger::endl;
  }

  if( fEventID < prcEntry ) {

    SetupEventInfo();
    ProceedEvent();
    FinishEvent();
  }
  else {
    LOG(INFO) << "STSpiRITTPC:: Finishing analysis. " << FairLogger::endl;
    Finish();
  }
}

void STSpiRITTPCTask::FinishEvent()
{
  LOG(DEBUG) << "STSpiRITTPCTask::FinishEvent is called. " << FairLogger::endl;

  if( fIsFlowAnalysis ) {

    fflowtask->FinishEvent();
    fflowinfo = fflowtask->GetFlowInfo();
    
    TClonesArray &aflow = *flowAnalysis;
    new( aflow[0] ) STFlowInfo( *fflowinfo );
  }

  auto anaRun = FairRunAna::Instance();
  if( ntrack[2] == 0) 
    anaRun->MarkFill(kFALSE);
  else
    anaRun->MarkFill(kTRUE);

}


void STSpiRITTPCTask::SetupEventInfo()
{
  auto fBDCArray = (TClonesArray *)fRootManager->GetObject("STBDC");
  auto fBDC = (STBDC*)fBDCArray->At(0);
  
  ProjA   = fBDC->GetProjA();
  ProjB   = fBDC->GetProjB();  

  if( fIsFlowAnalysis && fflowinfo != nullptr) 
    fflowtask->SetupEventInfo(fEventID, fBDC);

}



void STSpiRITTPCTask::SetupTrackQualityFlag(STParticle *apart) 
{
  if( apart->GetP() == 0 || apart->GetP() > 3100 )
    apart->SetMomentumFlag(0);

  if( apart->GetdEdx() <= 0. )
    apart->SetdEdxFlag(0);


  if( apart->GetDistanceAtVertex() > 20 )
    apart->SetDistanceAtVertexFlag(0);

  if( apart->GetDistanceAtVertex() > 20 )
    apart->SetDistanceAtVertexFlag(0);

  if( abs( apart->GetVertex().Z() + 13.1 ) > 1.7*3. )
    apart->SetVertexAtTargetFlag(0);

  if( apart->GetNDF() < 20) 
    apart->SetNDFFlag(0);

}

void STSpiRITTPCTask::Clear()
{
  for (Int_t m = 0; m < 7; m++) ntrack[m] = 0;

  tpcParticle->Clear();

  if( fIsFlowAnalysis ) {
    fflowinfo->Clear();
    flowAnalysis->Clear();
  }
}

void STSpiRITTPCTask::ProceedEvent()
{

  fChain -> GetEntry(fEventID);  

  ntrack[0] = trackArray -> GetEntries();

  LOG(DEBUG) << "nttack[0] -------------------- > " << ntrack[0] << FairLogger::endl;
                 
  TIter next(trackVAArray);

  STRecoTrack *trackFromArray = NULL;
  TClonesArray &ptpcParticle = *tpcParticle;

  while( (trackFromArray = (STRecoTrack*)next()) ) {

    auto parentvid = trackFromArray->GetVertexID();

    STVertex* vertex = NULL;
    if (parentvid > -1) {

      ntrack[1]++;

      vertex = (STVertex *) vertexArray -> At(parentvid);

      STParticle *aParticle = new STParticle();
      aParticle->SetRecoTrack(trackFromArray);

      //--- Set vertex of the track ---;  
      aParticle->SetVertex(vertex);

      SetupTrackQualityFlag( aParticle );

      //      if( aParticle->GetGoodTrackFlag()%2 == 0) continue;

      //--- Rotate tracks along beam direction ---;                    
      if(ProjA > -1000 && ProjB > -1000)
	aParticle->RotateAlongBeamDirection(ProjA/1000., ProjB/1000.);

      Int_t    Charge   = aParticle->GetCharge();
      TVector3 VMom     = aParticle->GetRotatedMomentum();
      Double_t dEdx     = aParticle->GetdEdx();

      //--- Set MassFitter      
      Double_t massH = 0.;
      Double_t massHe = 0.;
      Int_t    pid_tight  = 0;
      Int_t    pid_loose  = 0;      
      massFitter->GetBBMass(VMom, dEdx, Charge, massH, massHe, pid_tight, pid_loose); 
      aParticle->SetBBMass(massH);      
      aParticle->SetBBMassHe(massHe);
      aParticle->SetPID(pid_tight);
      aParticle->SetPIDLoose(pid_loose);

      LOG(DEBUG) << " mass "  << massHe << " & pid " << pid_loose << FairLogger::endl;

      //--- number cluster ratio
      if( kFALSE ) {
	Double_t clustNum = VMom.Mag()>4000.?
	  db->GetClusterNum(Charge, VMom.Theta(), VMom.Phi(), 4000.):
	  db->GetClusterNum(Charge, VMom.Theta(), VMom.Phi(), VMom.Mag());
	//--- Set theoretical number of cluster
	aParticle->SetExpectedClusterNumber(clustNum);
      }
      
      if( aParticle->GetGoodTrackFlag() > 1000 )
	ntrack[3]++;

      aParticle->SetTrackID(ntrack[2]);
      new(ptpcParticle[ntrack[2]]) STParticle(*aParticle);
      ntrack[2]++;

    }
  
  }

  //--- Set up for flow ---;
  if( fIsFlowAnalysis && fflowinfo != nullptr) {
    fflowtask->SetNTrack(ntrack);
    fflowtask->SetFlowTask( ptpcParticle  );
  }
}
