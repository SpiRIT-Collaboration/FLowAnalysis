#include "STSpiRITTPCTask.hh"

ClassImp(STSpiRITTPCTask);

STSpiRITTPCTask::STSpiRITTPCTask() 
  : fIsPersistence(kTRUE),
    fIsFlowAnalysis(kFALSE),
    fIsFlowCorrection(kFALSE),
    fIsBootStrap(kFALSE),
    fIsSubeventAnalysis(kFALSE),
    selReactionPlanef(100),
    fEventID(0),
    massCal(),
    massCalH(new STMassCalSimpleBB("EmpiricalBB")),
    massCalHe(new STMassCalSimpleBB("EmpiricalBB")),
    massGateFile(NULL)
{

  fLogger = FairLogger::GetLogger();
  SetVerbose(1);

  fRunAna = FairRunAna::Instance();
  iRun  = fRunAna->getRunId();

  fBeam         = new STBDC();
  bmA   = STRunToBeamA::GetBeamA(iRun);
  fBeam -> SetRun(iRun);
  
  if( bmA == 100 )
    BeamPID = 100;
  else {
    if( !SetupBeamCut(bmA) ) {
      LOG(ERROR) << " Beam cut file is not found. " << FairLogger::endl; 
      exit(0);
    }
  }

  //------------------------//
  // initial setup
  fChain = nullptr;

  vertexVAArray = new TClonesArray("STVertex",1);
  vertexArray   = new TClonesArray("STVertex",1);
}

Bool_t STSpiRITTPCTask::SetupBeamCut(UInt_t SnA)
{
  TString gcutFileName;

  if(SnA == 132)
    gcutFileName = "gcut132Sn.ROOT";
  else if(SnA == 108)
    gcutFileName = "gcut108Sn.ROOT";
  else if(SnA == 124)
    gcutFileName = "gcut124Sn.ROOT";
  else if(SnA == 112)
    gcutFileName = "gcut112Sn.ROOT";

  if( gcutFileName == "" )
    return kFALSE;


  auto gcutFile = new TFile( "data/"+gcutFileName );
  gBeamCut = (TCutG*)gcutFile->Get("sigma20");
  gcutFile->Close();

  if(gBeamCut == NULL) {
    LOG(ERROR) << " Beam Cut " << gcutFileName << " is not opened. " <<FairLogger::endl;
    return kFALSE;
  }

  return kTRUE;
}



STSpiRITTPCTask::~STSpiRITTPCTask()
{
  delete fBeam;           //!

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
  delete massCalH;
  delete massCalHe;

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

  fRootManager -> RegisterAny("BeamInfo", fBeam, kTRUE);

  tpcParticle = new TClonesArray("STKParticle",100);
  fRootManager -> Register("STKParticle","flow",tpcParticle, fIsPersistence);


  if( fIsFlowAnalysis ) {
    //    flowAnalysis = new TClonesArray("STFlowInfo",1);
    fRootManager -> RegisterAny("STFlow", fflowinfo, fIsPersistence);

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

  if( bmA == 124 )
    bmA = 132;
  
  if( bmA != 100 ) {

    auto calFile = TFile::Open(Form("db/PIDCalib_%dSn.root",bmA));
    TH2D *h2ParamH[7], *h2ParamHe[7];
    if( calFile  ) {
      for(auto i: ROOT::TSeqI(7)){
	calFile->GetObject(Form("h2InterpolateNM_%dSn_Par%d"  ,bmA,i),h2ParamH[i]);
	calFile->GetObject(Form("h2InterpolateHeNM_%dSn_Par%d",bmA,i),h2ParamHe[i]);
      }
      massCalH  -> AddParameters(h2ParamH);
      massCalHe -> AddParameters(h2ParamHe);
    }
    else {
      LOG(ERROR) << " Mass calibration file " << Form("PIDCalib_%dSn.root",bmA) << " is not found." << FairLogger::endl;
      fstatus =  kFALSE;
    }
  }  

  fstatus *= SetupPIDFit();

  return fstatus;
}

Bool_t STSpiRITTPCTask::SetupInputDataFile() 
{
  LOG(INFO) << "STSpiRITTPCTask::SetupInputDataFile() is called " << FairLogger::endl;

  fChain = new TChain("cbmsim");
  rootDir += Form("/Sn%d",bmA);

  UInt_t i = 0;
  while(kTRUE){


    TString recoFile = Form("run%04d_s%02d.reco."+tVer+".root",iRun,i);
    if( bmA == 100 )
      recoFile = Form("mizuki_%06d_s%02d.reco.v1.04.root",iRun,i);

    LOG(INFO) << i << " recoFile " << rootDir << " / " << recoFile << FairLogger::endl;


    if(gSystem->FindFile(rootDir,recoFile)) 
      fChain -> Add(recoFile);

    else if(i < 10) {
      recoFile = Form("run%04d_s%d.reco."+tVer+".root",iRun,i);
      if( bmA == 100 )
	recoFile = Form("mizuki_%06d_s%d.reco.v1.04.root",iRun,i);

      if(gSystem->FindFile(rootDir,recoFile)) 
	fChain -> Add(recoFile);
      else
	break;
    }

    else 
      break;
    
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

  fChain -> SetBranchAddress("VAVertex"   ,   &vertexVAArray);
  fChain -> SetBranchAddress("BDCVertex"  ,   &vertexBDCArray);

  fChain -> SetBranchAddress("STBeamInfo" ,   &beamInfo);

  if( bmA == 100 ) {
    fChain -> SetBranchAddress("STVertex_1" ,   &vertexArray);
    fChain -> SetBranchAddress("STRecoTrack",   &trackVAArray);
  }
  else
    fChain -> SetBranchAddress("STVertex"   ,   &vertexArray);

  //  fChain->Print();

  return kTRUE;
}


void STSpiRITTPCTask::Exec(Option_t *opt)
{

  LOG(DEBUG) << "STSpiRITTPCTask::Exec() is called " << fEventID << FairLogger::endl;  
  Clear();

  auto anaRun = FairRunAna::Instance();  
    
  fChain -> GetEntry(fEventID);  

  ShowProcessTime();
  fEventID++;

  Bool_t bfill = SetupEventInfo();
  if( bfill ) {

    if( ProceedEvent() ) 
      FinishEvent();
    else
      bfill = kFALSE;
  }

  anaRun->MarkFill(bfill);

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
  LOG(DEBUG) << "STSpiRITTPCTask::FinishEvent is called. " << BeamPID << FairLogger::endl;


  if( fIsFlowAnalysis ) {

    fflowtask->SetupEventInfo(rEventID, BeamPID);
    fflowtask->FinishEvent();

    if( bmA != 100 )
      fflowinfo = fflowtask->GetFlowInfo();
    
    // TClonesArray &aflow = *flowAnalysis;
    // new( aflow[0] ) STFlowInfo( *fflowinfo );
  }
}


Bool_t STSpiRITTPCTask::SetupEventInfo()
{
  if( bmA == 100 ) {
    BeamPID = 100;

    LOG(DEBUG) << "STSpiRITTPCTask::SetupEventInfo() " << FairLogger::endl;
    return kTRUE;
  }

  rEventID = eventHeader -> GetEventID() - 1;

  LOG(DEBUG) << " eventHeader " << rEventID << FairLogger::endl;

  ProjA   = beamInfo->fRotationAngleATargetPlane;//fBDC->ProjA;
  ProjB   = beamInfo->fRotationAngleBTargetPlane;//fBDC->ProjB;  
  ProjX   = beamInfo->fXTargetPlane;
  ProjY   = beamInfo->fYTargetPlane;

  fBeam->evt   = rEventID;
  fBeam->ProjA = ProjA;
  fBeam->ProjB = ProjA;
  fBeam->ProjX = ProjX;
  fBeam->ProjY = ProjY;
  fBeam->aoq   = beamInfo->fBeamAoQ; 
  fBeam->z     = beamInfo->fBeamZ; 

  BeamPID      = GetBeamPID(bmA, fBeam->aoq, fBeam->z);
  fBeam->SnA   = BeamPID;

  if( BeamPID > 0 )
    return kTRUE;

  return kFALSE;

}

Int_t STSpiRITTPCTask::GetBeamPID(Int_t SnA, Double_t aoq, Double_t z)
{

  Int_t pid = 0;
  if( SnA == 132){
    if( gBeamCut -> IsInside(aoq,z) && z < 50.536)
      pid = 132;
  }
  else {
    if( gBeamCut -> IsInside(aoq,z) )
      pid = SnA;
    else 
      pid = 0;
  }

  return pid;
}



Bool_t STSpiRITTPCTask::GetVertexQuality(TVector3 vert) 
{
  if( bmA == 100 ) return kTRUE; 

  if( ProjA < -99 || ProjA > 0) return kFALSE;
  
  if( ProjB < -99 ) return kFALSE;
  
  if( abs(ProjX) > 20 || abs(ProjY) > 20 ) return kFALSE;


  auto BeamIndex = STRunToBeamA::GetSystemID(iRun);
  if(BeamIndex >= 4) return kFALSE;

  if( abs( vert.Z() - VtxMean[BeamIndex].Z() ) > 3.*VtxSigm[BeamIndex] ||
      abs( vert.X() - VtxMean[BeamIndex].X() ) > 15. ||
      abs( vert.Y() - VtxMean[BeamIndex].Y() ) > 20. )
    return kFALSE;;
  
  return kTRUE;
}

void STSpiRITTPCTask::SetupTrackQualityFlag(STKParticle *apart) 
{
  if( bmA == 100 ) {
    cout << " vDistance " << apart->GetDistanceAtVertex() << " NDF " << apart->GetNDF() << endl;
    return ;
  }
  
  //   if( apart->GetDistanceAtVertex() > 25 )
  // if( apart->GetDistanceAtVertex() > 20 ) //since v54
  //    apart->SetDistanceAtVertexFlag(0);
  
  if( apart->GetNCL() < 15 )  // since v55
    apart->SetNCLFlag(0);

   //   if( apart->GetNDF() < 10) 
   // if( apart->GetNDF() < 15) // since v54 
   //   apart->SetNDFFlag(0);


   UInt_t pids = apart->GetPID_seq();

   if( pids < 7 ) {
     if( apart->GetP() < momRange[pids][0] ) {
       apart->SetMomentumFlag(0);
     }
   }

}

void STSpiRITTPCTask::Clear()
{
  rEventID = 0;
  BeamPID = 0;
  ProjA = 0.;
  ProjB = 0;
  ProjX = -999.;
  ProjY = -999.;
  fBeam -> Clear();
  
  for (Int_t m = 0; m < 7; m++) ntrack[m] = 0;

  tpcParticle->Clear("C");

  if( fIsFlowAnalysis ) {
    fflowtask->Clear();
    //    flowAnalysis->Clear();
  }
}


Bool_t STSpiRITTPCTask::ProceedEvent()
{

  Bool_t bprint = 1;
  if( bprint ) {
    LOG(DEBUG) << " ProceedEvent ... " << fEventID-1 << FairLogger::endl;
    bprint = kTRUE;
  }

  ntrack[0] = trackArray -> GetEntries(); // after 20191214

  auto vertex = (STVertex *) vertexArray -> At(0);
  if( !vertex ) return kFALSE;
  TVector3 vert  = vertex->GetPos();

  if( !GetVertexQuality( vert ) ) return kFALSE;

  Bool_t vflag = kFALSE;

  TIter nextReco(trackArray);
  TIter nextVA(trackVAArray);
  STRecoTrack *atrackReco = NULL;
  STRecoTrack *atrackVA = NULL;


  UInt_t ii = 0;
  while( (atrackReco = (STRecoTrack*)nextReco() ) ) {
    ii++;

    TVector3 dist = atrackReco->GetPOCAVertex() - vert;
    if( dist.Mag() <= 20 ) 
      ntrack[1]++;
  }

   LOG(DEBUG) << "-----------------" << FairLogger::endl;

  STKParticle *aParticle    = new STKParticle();

  //  for( auto it : ROOT::TSeqL( ntrack[0] ) ) {
  nextReco.Reset();
  while( (atrackReco = (STRecoTrack*)nextReco() ) ) {
    
    auto helixID = atrackReco->GetHelixID();

    LOG(DEBUG) << " reco id " << helixID;
    
    if( !(atrackVA = (STRecoTrack*)nextVA() ) ) break;
    auto helixIDVA = atrackVA->GetHelixID();

    LOG(DEBUG) << " vaid " << helixIDVA;

    Int_t diffID = helixIDVA - helixID;

    LOG(DEBUG) << " diff " << diffID ;
    if( diffID > 0  ) {
      while( 1 ) {

	if( !(atrackReco = (STRecoTrack*)nextReco() ) ) break;
	helixID = atrackReco->GetHelixID(); 

	LOG(DEBUG) << " helixid again " << helixID;
	
	diffID = helixIDVA - helixID;

	LOG(DEBUG) << " diff again " << diffID <<" ->> " ;
	if( diffID <= 0 ) break;
      }
    }
  

    aParticle->SetRecoTrack(atrackReco, atrackVA);

    //--- Rotate tracks along beam direction ---;                    
    aParticle->RotateAlongBeamDirection(ProjA/1000., ProjB/1000.);

    //--- Set event and track quality ---;  
    aParticle->SetVertex(vert);
      
    Int_t    Charge   = aParticle->GetCharge();
    TVector3 VMom     = aParticle->GetMomentumAtTarget();
    Double_t dEdx     = aParticle->GetdEdx();



    //--- Set MassFitter      
    Double_t mass[2] = {-1.,-1.};
    UInt_t   nPID[2] = {0, 0};

    // updated for 20191214
    mass[0] = massCalH -> CalcMass(1., VMom, dEdx, kTRUE);
    mass[1] = massCalHe-> CalcMass(2., VMom, dEdx, kTRUE);

    if( mass[0] > 0 ) {
      aParticle->SetBBMass(mass[0]);      
      nPID[0] = GetPIDFit(1, mass[0], VMom, ntrack[1]);
    }

    if( mass[1] > 0 ) {
      aParticle->SetBBMassHe(mass[1]);      
      nPID[1] = GetPIDFit(2, mass[1], VMom, ntrack[1]);
    }

    Int_t  pid_loose  = GetPIDLoose(mass, VMom.Mag(), dEdx);

    // if( nPID[0] != 0 && nPID[1] != 0 )
    //   aParticle->SetDoubleFlag(0);

    if( nPID[1] != 0 ) {
      aParticle->SetPID(nPID[1]);
      aParticle->SetBBMass(mass[1]);
    }
    else if( nPID[0] != 0 ) 
      aParticle->SetPID(nPID[0]);

    else if( pid_loose == 211 ){
      aParticle->SetPID(pid_loose);
      aParticle->SetMassFlag(0);
    }
    else if( nPID[0] == 0 && nPID[1] == 0 && pid_loose > 0 ) {
      aParticle->SetPIDLoose(pid_loose);
      aParticle->SetMass(pid_loose);
      aParticle->SetMassFlag(0);
    }
    else
      aParticle->SetPID(0);


    if( 0 ) 
      LOG(INFO) 
	<< "ntrack[1] "<< setw(3) << ntrack[1]
	<< " dEdx = "  << setw(8)  << dEdx
	<< " p =   "   << setw(10) << VMom.Mag()
	<< " mass H "  << setw(8)  << mass[0] 
	<< " nPID[0] " << setw(8)  << nPID[0]
	<< " mass He " << setw(8)  << mass[1]
	<< " nPID[1] " << setw(8)  << nPID[1]
	<< " dist "    << setw(5)  << aParticle->GetDistanceAtVertex()
	<< FairLogger::endl;

    LOG(DEBUG) << " dEdx " << dEdx 
	 << " P " << VMom.Mag()
	 << " rDist " << aParticle->GetDistanceAtVertex() 
	 << " PID " << nPID[0] << " : " << nPID[1] << " : " << pid_loose << FairLogger::endl;;


    SetupTrackQualityFlag( aParticle );

    //    if( aParticle->GetGoodTrackFlag() >= 0 ) {
    if( aParticle->GetGoodTrackFlag() >= 100000 ) {
      aParticle->SetTrackID(ntrack[2]);
      new( (*tpcParticle)[ntrack[2]] ) STKParticle(*aParticle);
      ntrack[2]++;    
    
      if( aParticle->GetGoodTrackFlag() == 111111 ){
	ntrack[3]++;
	vflag = kTRUE;
      }

      LOG(DEBUG) << " good " << aParticle->GetTrackID()<< " :: " << ntrack[2] << " dist " << aParticle->GetDistanceAtVertex() << FairLogger::endl; 
    }
    else
      LOG(DEBUG) << aParticle->GetTrackID() << " : " << ntrack[2] << " / " << ntrack[0] << " " << aParticle->GetDistanceAtVertex()	<< FairLogger::endl;
  }


  LOG(DEBUG) << " event : " << fEventID-1 << " mtrack1 " << ntrack[1] << " mtrack2 " << ntrack[2] << FairLogger::endl;
  
  //--- Set up for flow ---;
  if( fIsFlowAnalysis ) {
    fflowtask->SetGoodEventFlag((UInt_t)vflag);
    fflowtask->SetFlowTask( *tpcParticle  );
    for( auto i : {0,1,2,3} )
      fflowtask->SetNTrack(i,ntrack[i]);
  }

  return kTRUE;
}

Int_t STSpiRITTPCTask::GetPID(Double_t mass[2], Double_t fMom,  Double_t dedx)
{
  // p, d, t 
  if( mass[0] == 0 )
    return 0;

  else if( mass[1] < 2500 && mass[0] > 0 ) {

    for(UInt_t i = 0; i < 4; i++) {
      Double_t mass_low = MassRegion[i][0]-MassRegion[i][1]*MassRegion[i][2] ;
      Double_t mass_up  = MassRegion[i][0]+MassRegion[i][1]*MassRegion[i][3] ;

      if( mass[0] >= mass_low && mass[0] <= mass_up ) {

	if( i == 1 && fMom > protonMaxMomentum ) continue;
 
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

Int_t STSpiRITTPCTask::GetPIDLoose(Double_t mass[2], Double_t fMom, Double_t dedx)
{
  if( mass[0] == 0 )
    return 0;

  // p, d, t   
  else if( mass[1] < MassRegionLU_L[4][0] ) {
    for(UInt_t i = 0; i < 4; i++) {
      if( mass[0] >= MassRegionLU_L[i][0] && mass[0] < MassRegionLU_L[i][1] ) {

	if( i == 1 && fMom > protonMaxMomentum ) continue;

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


Bool_t STSpiRITTPCTask::SetupPIDFit()
{
  if( bmA == 100 ) return kTRUE;

  TString pidName[]={"Proton","Deuteron","Triton","Helium3","Alpha"}; 

  std::vector< TString > label_rl = {"left","right"};

  UInt_t k = 0;
  for( auto il : label_rl ) {
    TString massFitName = Form("db/MassFit_%dSn.%s.root",bmA,il.Data());
    massGateFile = TFile::Open(massFitName);
    if( massGateFile == NULL ) {
      LOG(ERROR) << massFitName << " is not found. " << FairLogger::endl;
      return kFALSE;
    }    

    LOG(INFO) << massFitName << " is loaded. " << FairLogger::endl;

    for(auto pid: ROOT::MakeSeq(5)) 
      for(auto mbin: ROOT::MakeSeq(4)) 
	for(auto i: ROOT::TSeqI(2)) 
	  for(auto j: ROOT::TSeqI(4)) { 
	    TString f1name = Form("f1MassGate_%dSn_"+pidName[pid]+"_mbin%d_%d_%d",bmA,mbin,i,j);
	    massGateFile->GetObject(f1name,  f1MassGate[pid][mbin][i][j][k]);

	    if( f1MassGate[pid][mbin][i][j][k] == NULL ) {
	      LOG(ERROR) << f1name << " is not found. " << FairLogger::endl;
	      return kFALSE;
	    }
	  }
    massGateFile->Close();
    k++;
  }
  return kTRUE;
}

Int_t STSpiRITTPCTask::GetPIDFit(Int_t z, Double_t mass, TVector3 vMom, Int_t mult)
//Double_t mass[2], Double_t fMom, Int_t mbin, Int_t phibin)
{
  // coded by Kaneko 20200418
  // pid: p=0, d=1, t=2, 3he=3, 4he=4. 
  // mass[0]: calibrated mass of track-> p,d,t: mass[1]:helium mass
  // fMom: rigidity magnitude of track
  // mbin: centrality bin, mbin=0: M>=56, mbin=1: 50<=M<56, mbin=2: 40<=M<50, mbin=3: no selection

  Int_t fpid = -1;
  if( mass <= 0 ) return 0;

  Double_t fMom = vMom.Mag();
  if( fMom < 100) return 0;

  Double_t phi  = vMom.Phi()*TMath::RadToDeg();
  UInt_t phibin = 1;
  if( phi <= 20 && phi >= -30 ) phibin = 0;

  UInt_t mbin = 0;
  if( mult >= 56 ) mbin = 0;
  else if( mult >= 50 && mult < 56 ) mbin = 1;
  else if( mult >= 40 && mult < 50 ) mbin = 2;
  else mbin = 3;

  Bool_t bfind = kFALSE;

  if( z == 1 ) {
    for(auto i : {0,1,2} ) {

      Bool_t roughCut =  mass >= MassRange_Fit[i][0] && mass <= MassRange_Fit[i][1]; 

      if( roughCut ) {
    
	Int_t fitSigmaID =  3 ;
	Bool_t fitCut = mass >= f1MassGate[i][mbin][0][fitSigmaID][phibin]->Eval(fMom) && 
	  mass <= f1MassGate[i][mbin][1][fitSigmaID][phibin]->Eval(fMom);
    
	if( fitCut ) {
	  fpid = i;
	  //	  cout << " fpid " << fpid << " mass " << mass << endl;
	  bfind = kTRUE;
	  Break;
	}
	else if ( 0 ) {
	  cout << " mass " << mass 
	       << " p " << fMom
	       << " fitSigmaID " << fitSigmaID
	       << " mbin " << mbin
	       << " mult " << mult
	       << " " << f1MassGate[i][mbin][0][fitSigmaID][phibin]->Eval(fMom)
	       << " " << f1MassGate[i][mbin][1][fitSigmaID][phibin]->Eval(fMom)
	       << " phi " << phi
	       << " phibin " << phibin
	       << endl;
	}

      }
    }
  }
  else if( z == 2 ) {
    fMom *= 2.;
    for(auto i : {3,4} ) {
      Bool_t roughCut =  mass >= MassRange_Fit[i][0] && mass <= MassRange_Fit[i][1];
      if( roughCut ) {
    
	Int_t fitSigmaID =  2 ;
	Bool_t fitCut = mass >= f1MassGate[i][mbin][0][fitSigmaID][phibin]->Eval(fMom) && mass <= f1MassGate[i][mbin][1][fitSigmaID][phibin]->Eval(fMom);

	if( i == 3 )
	  fitCut *= mass >= 2850.-0.35*fMom/2.;

    
	if( fitCut ) {
	  fpid = i;
	  

	  bfind = kTRUE;
	  Break;
	}
      }
    }
  }



  STPID::PID pid = static_cast<STPID::PID>(fpid+1);
  Int_t rpid = STPID::GetPDG(pid) * bfind;

  LOG(DEBUG) << " GetPIDFit " << mass  << " fpid " << pid << " " << rpid <<FairLogger::endl;
        
  return rpid;
}

Int_t STSpiRITTPCTask::GetPIDTight(Double_t mass[2], Double_t fMom, Double_t dedx)
{
  if( mass[0] == 0 )
    return 0;

  // p, d, t   
  else if( mass[1] < MassRegionLU_T[4][0]  ) {
    for(Int_t i = 3; i >=0; i--) {
    //    for(UInt_t i = 0; i < 4; i++) {
      if( mass[0] >= MassRegionLU_T[i][0] && mass[0] < MassRegionLU_T[i][1] ) {

	if( i == 1 && fMom > protonMaxMomentum ) continue; 

	STPID::PID pid = static_cast<STPID::PID>(i);
	return STPID::GetPDG(pid);
      }
    }
  }

  // He3, He4, He6
  else if( mass[0] >= 3100 && dedx <= 700) {
    //    for( UInt_t i = 4; i < 7; i++ ){
    for( Int_t i = 6; i >= 4; i-- ){
      if( mass[1] >= MassRegionLU_T[i][0] && mass[1] < MassRegionLU_T[i][1] ) {
	STPID::PID pid = static_cast<STPID::PID>(i);
        return STPID::GetPDG(pid);
      }
    }
  }

  return 0;
}

