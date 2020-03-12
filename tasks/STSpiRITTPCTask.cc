#include "STSpiRITTPCTask.hh"

ClassImp(STSpiRITTPCTask);

STSpiRITTPCTask::STSpiRITTPCTask() 
  : fIsPersistence(kTRUE),
    fIsFlowAnalysis(kTRUE),
    fIsFlowCorrection(kTRUE),
    fIsBootStrap(kFALSE),
    fIsSubeventAnalysis(kTRUE),
    selReactionPlanef(100),
    fEventID(0),
    massCal(new STMassCalculator())
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
  //  delete massFitter;
  delete massCal;

  delete flowAnalysis; //!
  //  delete fflowinfo;    //!

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
  //  massCal->SetParameter("db/BBFitter.root");

  UInt_t bmA = STRunToBeamA::GetBeamA(iRun);

  if( bmA == 124 )
    bmA = 132;  

  //  massCal->LoadCalibrationParameters("db/PIDCalib.root",bmA);
  massCal->LoadCalibrationParameters("db/FlattenPID.root",bmA);

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

    TString recoFile = Form("run%04d_s%02d.reco."+tVer+".root",iRun,i);

    if(gSystem->FindFile(rootDir,recoFile)) 
      fChain -> Add(recoFile);

    else if(i < 10) {
      recoFile = Form("run%04d_s%d.reco."+tVer+".root",iRun,i);
      if(gSystem->FindFile(rootDir,recoFile)) 
	fChain -> Add(recoFile);
      else
	break;
    }

    else 
      break;
    
    LOG(INFO) << i << " recoFile " << recoFile << FairLogger::endl;
    i++;

  }

  if( i == 0 ) {
    LOG(ERROR) << "Cannot find TPC data! "  <<  FairLogger::endl;
    return kFALSE;
  }

  nEntry = fChain->GetEntries();
  LOG(INFO) << i << " files were chained. " << nEntry << " events" << FairLogger::endl;
  
  prcEntry = nEntry;

  fChain -> SetBranchAddress("STEventHeader", &eventHeader);
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
    
  fChain -> GetEntry(fEventID);  

  ShowProcessTime();

  if(SetupEventInfo()) {
    if(ProceedEvent()) 
      FinishEvent();

    fEventID++;
  }

  if( fEventID >= prcEntry )  {
      LOG(INFO) << "STSpiRITTPC:: Finishing analysis. " << FairLogger::endl;
      Finish();
  }

}

void STSpiRITTPCTask::ShowProcessTime()
{
  UInt_t pdev = prcEntry/200;
  if(prcEntry < 200) pdev = 1;

  if(fEventID%pdev == 0) {
    dtime.Set();
    Int_t ptime = dtime.Get() - beginTime.Get();
    LOG(INFO) << "Process " 
	      << setw(4) << Int_t((Double_t)fEventID/(Double_t)prcEntry*100.) << " % = "
	      << setw(8) << fEventID << "/"<< prcEntry 
	      << "--->"
	      << dtime.AsString() << " ---- "
	      << (Int_t)ptime/60 << " [min] "
	      << FairLogger::endl;
  }
}

void STSpiRITTPCTask::FinishEvent()
{
  LOG(DEBUG) << "STSpiRITTPCTask::FinishEvent is called. " << FairLogger::endl;

  auto anaRun = FairRunAna::Instance();  
  // if( ntrack[2] == 0 || BeamPID == 0) {
  //   anaRun->MarkFill(kFALSE);
  //   return;
  // }

  if( fIsFlowAnalysis ) {

    fflowtask->SetupEventInfo(rEventID, BeamPID);
    fflowtask->FinishEvent();
    fflowinfo = fflowtask->GetFlowInfo();
    
    TClonesArray &aflow = *flowAnalysis;
    new( aflow[0] ) STFlowInfo( *fflowinfo );
  }

  anaRun->MarkFill(kTRUE);

}


Bool_t STSpiRITTPCTask::SetupEventInfo()
{
  rEventID = eventHeader -> GetEventID() - 1;
  
  auto fBDCArray = (TClonesArray *)fRootManager->GetObject("STBDC");
  auto fBDC = (STBDC*)fBDCArray->At(0);

  LOG(DEBUG) << " TPC vs BDC " << rEventID << " vs " << fBDC->GetEventID() << FairLogger::endl;

  if( rEventID != fBDC->GetEventID() ) {
    LOG(ERROR) << " not match TPC and BDC event ID " << rEventID << " vs " << fBDC->GetEventID() << FairLogger::endl;
    return kFALSE;
  }

  ProjA   = fBDC->ProjA;
  ProjB   = fBDC->ProjB;  
  ProjX   = fBDC->ProjX;
  ProjY   = fBDC->ProjY;
  BeamPID = fBDC->GetBeamPID();

  return kTRUE;
}



Bool_t STSpiRITTPCTask::GetVertexQuality(TVector3 vert) 
{
  if( ProjA < -99 || ProjA > 100) return kFALSE;

  if( ProjB < -99 ) return kFALSE;

  if( abs(ProjX) > 20 || abs(ProjY) > 20 ) return kFALSE;

  auto BeamIndex = STRunToBeamA::GetSystemID(iRun);
  if(BeamIndex >= 4) return kFALSE;

  if( abs( vert.Z() - VtxMean[BeamIndex].Z() ) > 2.*VtxSigm[BeamIndex] ||
      abs( vert.X() - VtxMean[BeamIndex].X() ) > 15. ||
      abs( vert.Y() - VtxMean[BeamIndex].Y() ) > 20. )
    return kFALSE;;

  return kTRUE;
}

void STSpiRITTPCTask::SetupTrackQualityFlag(STParticle *apart) 
{
  if( apart->GetDistanceAtVertex() > 10 )
    apart->SetDistanceAtVertexFlag(0);

  if( apart->GetNDF() < 15) 
    apart->SetNDFFlag(0);
}

void STSpiRITTPCTask::Clear()
{
  rEventID = 0;
  BeamPID = 0;
  ProjA = 0.;
  ProjB = 0;
  ProjX = -999.;
  ProjY = -999.;

  for (Int_t m = 0; m < 7; m++) ntrack[m] = 0;

  tpcParticle->Clear();

  if( fIsFlowAnalysis ) {
    fflowtask->Clear();
    flowAnalysis->Clear();
  }
}

Bool_t STSpiRITTPCTask::ProceedEvent()
{

  ntrack[0] = trackArray -> GetEntries();

  LOG(DEBUG) << "nttack[0] -------------------- > " << ntrack[0] << FairLogger::endl;
                 
  auto vertex = (STVertex *) vertexArray -> At(0);

  Bool_t vflag = kFALSE;
  if( vertex != NULL ) 
    vflag = GetVertexQuality(vertex->GetPos());
  else
    return kFALSE;


  TIter next(trackArray);
  STRecoTrack *trackFromArray = NULL;
  TClonesArray &ptpcParticle = *tpcParticle;

  while( (trackFromArray = (STRecoTrack*)next()) ) {

    ntrack[1]++;

    STParticle *aParticle = new STParticle();
    aParticle->SetRecoTrack(trackFromArray);

    //--- Set event and track quality ---;  
    aParticle->SetVertex(vertex);
    aParticle->SetVertexAtTargetFlag((Int_t)vflag);
    SetupTrackQualityFlag( aParticle );

    //--- Rotate tracks along beam direction ---;                    
    if( vflag ) 
      aParticle->RotateAlongBeamDirection(ProjA/1000., ProjB/1000.);
    else
      aParticle->SetBeamonTargetFlag(0);


    Int_t    Charge   = aParticle->GetCharge();
    TVector3 VMom     = aParticle->GetRotatedMomentum();
    Double_t dEdx     = aParticle->GetdEdx();


    //--- Set MassFitter      
    Double_t mass[2] = {0.,0.};
    if( dEdx > -1 ){
      mass[0]  = massCal->CalcMass(0, 1., VMom, dEdx);  // proton fitted
      
      if( mass[0] > 1500. )  // deuteron fitted
	mass[0]  = massCal->CalcMass(1, 1., VMom, dEdx);  
      
      mass[1]  = massCal->CalcMass(1, 2., VMom, dEdx);
    }

    aParticle->SetBBMass(mass[0]);      
    aParticle->SetBBMassHe(mass[1]);
      
    Int_t    pid_tight  = GetPIDTight(mass, dEdx);
    Int_t    pid_normal = GetPIDNorm(mass, dEdx);
    Int_t    pid_loose  = GetPIDLoose(mass, dEdx);
    
    //      massFitter->GetBBMass(VMom, dEdx, Charge, massH, massHe, pid_tight, pid_loose); 
    aParticle->SetPID(pid_tight);
    aParticle->SetPIDTight(pid_tight);
    aParticle->SetPIDNorm(pid_normal);
    aParticle->SetPIDLoose(pid_loose);
    

    LOG(DEBUG) << " mass H "  << mass[0]  << " & pid " << pid_loose << " : " << pid_tight << ": " << dEdx << FairLogger::endl;
    LOG(DEBUG) << " mass He"  << mass[1] << " & pid " << pid_loose << " : " << pid_tight  << ": " << dEdx << FairLogger::endl;

    //--- number cluster ratio
    if( kFALSE ) {
      Double_t clustNum = VMom.Mag()>4000.?
	db->GetClusterNum(Charge, VMom.Theta(), VMom.Phi(), 4000.):
	db->GetClusterNum(Charge, VMom.Theta(), VMom.Phi(), VMom.Mag());
      //--- Set theoretical number of cluster
      aParticle->SetExpectedClusterNumber(clustNum);
    }
      
    if( aParticle->GetGoodTrackFlag() >= 1 )
      ntrack[3]++;
    
    aParticle->SetTrackID(ntrack[2]);
    new(ptpcParticle[ntrack[2]]) STParticle(*aParticle);
    ntrack[2]++;
    
  }
  
  //--- Set up for flow ---;
  if( fIsFlowAnalysis ) {
    fflowtask->SetGoodEventFlag((UInt_t)vflag);
    fflowtask->SetNTrack(ntrack);
    fflowtask->SetFlowTask( ptpcParticle  );
  }

  return kTRUE;
}

Int_t STSpiRITTPCTask::GetPID(Double_t mass[2], Double_t dedx)
{
  // p, d, t                                                                                                                               
  if( mass[0] == 0 )
    return 0;

  else if( mass[1] < 2500 && mass[0] > 0 ) {

    for(UInt_t i = 0; i < 4; i++) {
      Double_t mass_low = MassRegion[i][0]-MassRegion[i][1]*MassRegion[i][2] ;
      Double_t mass_up  = MassRegion[i][0]+MassRegion[i][1]*MassRegion[i][3] ;

      if( mass[0] >= mass_low && mass[0] <= mass_up ) {
	STPID::PID pid = static_cast<STPID::PID>(i);
        return STPID::GetPDG(pid);
      }
    }
  }

  // He3, He4, He6                                                                                                                          
  else if( mass[0] >= 3100 && dedx <= 700) {
    for( UInt_t i = 4; i < 7; i++ ){
      Double_t mass_low = MassRegion[i][0]-MassRegion[i][1]*MassRegion[i][2] ;
      Double_t mass_up  = MassRegion[i][0]+MassRegion[i][1]*MassRegion[i][3] ;

      if( mass[1] >= mass_low && mass[1] <= mass_up ) {
	STPID::PID pid = static_cast<STPID::PID>(i);
        return STPID::GetPDG(pid);
      }
    }
  }
  return 0;
}

Int_t STSpiRITTPCTask::GetPIDLoose(Double_t mass[2], Double_t dedx)
{
  if( mass[0] == 0 )
    return 0;

  // p, d, t   
  else if( mass[1] < MassRegionLU_L[4][0] ) {
    for(UInt_t i = 0; i < 4; i++) {
      if( mass[0] >= MassRegionLU_L[i][0] && mass[0] < MassRegionLU_L[i][1] ) {
	STPID::PID pid = static_cast<STPID::PID>(i);
        return STPID::GetPDG(pid);
      }
    }
  }
  // He3, He4, He6
  else if( mass[0] >= 3100 && dedx <= 700) {
    for( UInt_t i = 4; i < 7; i++ ){
      if( mass[1] >= MassRegionLU_L[i][0] && mass[1] < MassRegionLU_L[i][1] ) {
	STPID::PID pid = static_cast<STPID::PID>(i);
        return STPID::GetPDG(pid);
      }
    }
  }
  return 0;
}
Int_t STSpiRITTPCTask::GetPIDTight(Double_t mass[2], Double_t dedx)
{
  if( mass[0] == 0 )
    return 0;

  // p, d, t   
  else if( mass[1] < MassRegionLU_T[4][0]  ) {
    for(UInt_t i = 0; i < 4; i++) {
      if( mass[0] >= MassRegionLU_T[i][0] && mass[0] < MassRegionLU_T[i][1] ) {
	STPID::PID pid = static_cast<STPID::PID>(i);
        return STPID::GetPDG(pid);
      }
    }
  }
  // He3, He4, He6
  else if( mass[0] >= 3100 && dedx <= 700) {
    for( UInt_t i = 4; i < 7; i++ ){
      if( mass[1] >= MassRegionLU_T[i][0] && mass[1] < MassRegionLU_T[i][1] ) {
	STPID::PID pid = static_cast<STPID::PID>(i);
        return STPID::GetPDG(pid);
      }
    }
  }
  return 0;
}
Int_t STSpiRITTPCTask::GetPIDNorm(Double_t mass[2], Double_t dedx)
{
  if( mass[0] == 0 )
    return 0;

  // p, d, t   
  else if( mass[1] < MassRegionLU_N[4][0] ) {
    for(UInt_t i = 0; i < 4; i++) {
      if( mass[0] >= MassRegionLU_N[i][0] && mass[0] < MassRegionLU_N[i][1] ) {
	STPID::PID pid = static_cast<STPID::PID>(i);
        return STPID::GetPDG(pid);
      }
    }
  }
  // He3, He4, He6
  else if( mass[0] >= 3100 && dedx <= 700) {
    for( UInt_t i = 4; i < 7; i++ ){
      if( mass[1] >= MassRegionLU_N[i][0] && mass[1] < MassRegionLU_N[i][1] ) {
	STPID::PID pid = static_cast<STPID::PID>(i);
        return STPID::GetPDG(pid);
      }
    }
  }
  return 0;
}
