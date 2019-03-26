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
    fLink(new FairLink()),
    massFitter(new STMassFunction())
{
  fLogger = FairLogger::GetLogger();
  SetVerbose(1);

  fRunAna = FairRunAna::Instance();
  iRun  = fRunAna->getRunId();

  dtime.Set();
  rnd.SetSeed(dtime.GetSecond());  

  //------------------------//
  // initial setup
  fChain = nullptr;

  // output tree file

  fflowinfo = new STFlowInfo();


  if( fIsBootStrap ) {
    bs_unitP = new STBootStrap(10);
    LOG(INFO) << "BootStrap is called. " << FairLogger::endl;
  }

}

STSpiRITTPCTask::~STSpiRITTPCTask()
{

  delete db;
  delete fLink;

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

  if( bs_unitP != nullptr )
    delete bs_unitP;


}


void STSpiRITTPCTask::SetPersistence(Bool_t value) { fIsPersistence = value; }


InitStatus STSpiRITTPCTask::Init()
{
  beginTime.Copy(dtime);

  fRootManager = FairRootManager::Instance();
  if ( fRootManager == nullptr) {
    LOG(ERROR) << "Cannot find RootManager!" << FairLogger::endl;
    return kERROR;
  }

  LOG(INFO) << "STSpiRITTPCTask::InitTask() is called " << FairLogger::endl;

  auto fBDCArray = (TClonesArray *)fRootManager->GetObject("STBDC");
  if( fBDCArray == nullptr ) {
    LOG(ERROR) << "Register STBDC in advance! " << FairLogger::endl;
    return kERROR;
  }
  else
    LOG(INFO) << "STBDC is found. " << FairLogger::endl;

  if( fRunAna == nullptr ) {
    LOG(ERROR) << "STSpiRITTPCTask ==> Not defined. "  << FairLogger::endl;
    return kERROR;
  }

  if( !SetupInputDataFile() ) {
    LOG(ERROR) << "STSpiRITTPCTask:: Cannot open input files" << FairLogger::endl;
    return kERROR;
  }

  if( !SetupParameters() ) {
    LOG(ERROR) << "STSpiRITTPCTask:: Cannot open parameter files" << FairLogger::endl;
    return kERROR;
  }


  if( fIsFlowAnalysis && fIsFlowCorrection) {
    if( !SetupFlowDataBase() ) {
      LOG(ERROR) << "STSpiRITTPCTask:: Flow database cannot be setup" << FairLogger::endl;
      return kERROR;
    }
    else
      LOG(INFO) << "STSpiRITTPCTask:: Flow database are ready. " << FairLogger::endl;
  }

  tpcParticle = new TClonesArray("STParticle",100);
  fRootManager -> Register("STParticle","flow",tpcParticle, fIsPersistence);

  if( fIsFlowAnalysis ) {
    flowAnalysis = new TClonesArray("STFlowInfo",1);
    fRootManager -> Register("STFlow","flow",flowAnalysis, fIsPersistence);
  }

  return kSUCCESS;
}


Bool_t STSpiRITTPCTask::SetupParameters()
{
  Bool_t fstatus = kTRUE;

  //--- angle dependent mass fitter // 2019/02/19
  

  //--- angle dependeing mass fitter //
  //  fstatus = massFitter->SetFunction("/cache/scr/spirit/mizuki/Bethe-Bloch_Fitter/mk_NewFitter_20190111/","LVBBFitter.root");
  // fstatus = massFitter->SetFunction("/cache/scr/spirit/mizuki/Bethe-Bloch_Fitter/mk_NewFitter_20190111/",
  // 				    "IterBBFitter.root",
  // 				    "f1IterBBProton_" + STRunToBeamA::GetBeamSnA(iRun) );
  fstatus = massFitter->SetFunction("/cache/scr/spirit/mizuki/Bethe-Bloch_Fitter/mk_NewFitter_20190111/",
				    "BBFitter.root",
				    "f1IterBBProtonRotate_" + STRunToBeamA::GetBeamSnA(iRun) );
  if( !fstatus ) LOG(ERROR) << " BB Mass fittter function was not loaded. " << FairLogger::endl;


  // fstatus = massFitter->SetMassFitFunction("/cache/scr/spirit/mizuki/Bethe-Bloch_Fitter/mk_NewFitter_20190111/",
  // 					   "MassFitter.root",
  // 					   "f1MassFit_" + STRunToBeamA::GetBeamSnA(iRun) );
  fstatus = massFitter->SetMassFitFunction("/cache/scr/spirit/mizuki/Bethe-Bloch_Fitter/mk_NewFitter_20190111/",
					   "MassFitter.root",
					   "f1IterMassFitRotate_" + STRunToBeamA::GetBeamSnA(iRun) );
  if( ! fstatus ) LOG(ERROR) << " BB Mass gaus fit function was not loaded. " << FairLogger::endl;

  //------------------------------

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


  return fstatus;
}

Bool_t STSpiRITTPCTask::SetupInputDataFile() 
{
  LOG(INFO) << "STSpiRITTPCTask::SetupInputDataFile() is called " << FairLogger::endl;

  fChain = new TChain("cbmsim");

  UInt_t i = 0;
  while(kTRUE){

    TString recoFile = Form("run%04d_s%d.reco."+sVer+".root",iRun,i);


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

  //  fRootManager->SetInChain(fChain);

  fLink->SetLink(-1, -1, -1, 1, 1.);
  auto tempclone = fRootManager->GetCloneOfTClonesArray(*fLink);
  if( tempclone != nullptr)
    tempclone->Print();

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

  //  prcEntry = 100;

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

  if( fIsFlowAnalysis && fflowinfo != nullptr) {

    fflowinfo->SetNTrack(ntrack);
    fflowinfo->unitP   = unitP;

    if( fIsSubeventAnalysis ) {
      DoSubeventAnalysis();
      
      fflowinfo->unitP_1 = unitP_1;
      fflowinfo->unitP_2 = unitP_2;

      fflowinfo->mtrack_1 = mtrack_1;
      fflowinfo->mtrack_2 = mtrack_2;
    }



    DoIndividualReactionPlaneAnalysis();


    if( fIsFlowCorrection && aflowcorrArray[0] != NULL ) {

      fflowinfo->unitP_fc = Psi_FlatteningCorrection( 0, ntrack[4], TVector3(unitP.X(), unitP.Y(), 0.));
      fflowinfo->unitP_rc = Psi_ReCenteringCorrection(0, ntrack[4], TVector3(unitP.X(), unitP.Y(), 0.));

    }

    
    TClonesArray &aflow = *flowAnalysis;
    new( aflow[0] ) STFlowInfo( *fflowinfo );


  }

  auto anaRun = FairRunAna::Instance();
  if( ntrack[2] == 0) 
    anaRun->MarkFill(kFALSE);
  else
    anaRun->MarkFill(kTRUE);

}

void STSpiRITTPCTask::DoSubeventAnalysis()
{
  if( ntrack[4] == 0 ) return;
  
  STBootStrap* bs_unitP_1 = new STBootStrap(1000);
  STBootStrap* bs_unitP_2 = new STBootStrap(1000);

  UInt_t np = ntrack[4];
  if(np%2 == 1) np++;
  const UInt_t npart = np;
  UInt_t *rndArray = new UInt_t[npart];
  rndArray = RandomDivide2(npart);

  TIter next(tpcParticle);
  STParticle *apart = NULL;

  UInt_t itrack = 0;
  while( (apart = (STParticle*)next() ) ) {
    
    if( apart->GetReactionPlaneFlag() >= selReactionPlanef ) {
      Double_t wt = apart->GetRPWeight();
      TVector2 ptr= apart->GetRotatedPt();

      if( rndArray[itrack] == 0 ) {

        unitP_1+= wt * ptr.Unit();
	LOG(DEBUG) << " sub 1 " << unitP_1.X()  
		   << " + "     << wt * ptr.X()
		   << " : "     << apart->GetReactionPlaneFlag()
		   << FairLogger::endl;

        if( fIsBootStrap )
          bs_unitP_1->Add(wt * ptr.Unit());

        TVector2 ptpt = wt * ptr.Unit();
        apart->AddReactionPlaneFlag(100);
        mtrack_1++;
      }
      else  {
        unitP_2+= wt * ptr.Unit();
	LOG(DEBUG) << " sub 2    " << unitP_2.X()  
		   << " + "     << wt * ptr.X()
		   << " : "     << apart->GetReactionPlaneFlag()
		   << FairLogger::endl;
        if( fIsBootStrap )
          bs_unitP_2->Add(wt * ptr.Unit());

        apart->AddReactionPlaneFlag(200);
        mtrack_2++;
      }

      itrack++;
      if( itrack > npart ) break;

    }
  }

  if(fIsBootStrap && mtrack_1 > 0 && mtrack_2 > 0 ) {
    bs_unitP_1->BootStrapping();
    bs_unitP_2->BootStrapping();

    bsPhi_1[0] = bs_unitP_1->GetMean();
    bsPhi_1[1] = bs_unitP_1->GetStdDev();
    bsPhi_1[2] = bs_unitP_1->GetMod();

    bsPhi_2[0] = bs_unitP_2->GetMean();
    bsPhi_2[1] = bs_unitP_2->GetStdDev();
    bsPhi_2[2] = bs_unitP_2->GetMod();

  }

  delete bs_unitP_1;
  delete bs_unitP_2;
  
}
UInt_t *STSpiRITTPCTask::RandomDivide2(const UInt_t npart)
{
  UInt_t  *rndarray = new UInt_t[npart];

  UInt_t c1 = 0;
  UInt_t c2 = 0;
  UInt_t count = 0;
  while( count < npart ){

    Float_t rrd = rnd.Rndm() ;

    if( rrd < 0.5 ) {
      if( c1 < npart/2 ) {
        rndarray[count] = 0;
        c1++;
        count++;
      }
    }
    else if( rrd >= 0.5 ) {
      if( c2 < npart/2 ) {
        rndarray[count] = 1;
        c2++;
        count++;
      }
    }
  }


  return rndarray;
}


void STSpiRITTPCTask::SetupEventInfo()
{
  auto fBDCArray = (TClonesArray *)fRootManager->GetObject("STBDC");
  auto fBDC = (STBDC*)fBDCArray->At(0);
  
  ProjA   = fBDC->GetProjA();
  ProjB   = fBDC->GetProjB();  

  if( fIsFlowAnalysis && fflowinfo != nullptr) {

    fflowinfo->SetEventID( fEventID );
    fflowinfo->SetRun( iRun );
    fflowinfo->SnA     = fBDC->SnA;
    fflowinfo->beamPID = fBDC->beamPID;
    LOG(DEBUG) << "flowinfo " << fEventID << " * " << fflowinfo <<FairLogger::endl;
  }
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
}

void STSpiRITTPCTask::SetupTrackExtraQualityFlag(STParticle *apart)
{

  if( apart->GetNDF() < 20) 
    apart->SetNDFFlag(0);

  // if( apart->GetClusterRatio() < 0.7 || apart->GetClusterRatio() > 2 ) 
  //   apart->SetClusterRatioFlag(0);
}


void STSpiRITTPCTask::SetupFlow(STParticle &apart)
{
  // Setup for flow analysis

  auto pid    =  apart.GetPIDLoose();
  if( pid == 211 )
    apart.SetReactionPlaneFlag(1);

  else if( pid > 2000 &&
           apart.GetGoodTrackFlag()     > 0 &&
           apart.GetdEdxFlag()          > 0 &&
           apart.GetMomentumFlag()      > 0
           ) {
    apart.SetReactionPlaneFlag(1000);

    if(apart.GetNDFFlag())
      apart.SetReactionPlaneFlag(2000);
  }
  else
    apart.SetReactionPlaneFlag(0);


  // Pt weight
  TLorentzVector lrnzVec =  apart.GetLorentzVector();
  
  TVector3 boostVec = STLorentzBoostVector::GetBoostVector(4); //4: p+p

  lrnzVec.Boost(-boostVec);

  auto rapiditycm = lrnzVec.Rapidity();

  apart.SetRapiditycm(rapiditycm);

  if( rapiditycm  <  0 )
    apart.SetRPWeight(-1);
  else
    apart.SetRPWeight(1);

}

void STSpiRITTPCTask::DoFlowAnalysis(STParticle &apart)
{
  
  SetupFlow( apart );

  if( apart.GetReactionPlaneFlag() >= selReactionPlanef ){
    ntrack[4]++;

    unitP += apart.GetRPWeight() * apart.GetRotatedPt().Unit();

    if( fIsBootStrap )
      bs_unitP->Add(apart.GetRPWeight() * apart.GetRotatedPt().Unit());
  }

  if( apart.GetReactionPlaneFlag() == 20 )
    ntrack[5]++;  
}



void STSpiRITTPCTask::DoIndividualReactionPlaneAnalysis()
{
  TIter orgnext(tpcParticle);
  STParticle *apart = NULL;

  while( (apart = (STParticle*)orgnext()) ) {

    Double_t wt = apart->GetRPWeight();
    TVector2 pt = apart->GetRotatedPt();

   // individual reaction plane
    UInt_t itraex = 0;
    TVector2 mExcRP(0.,0.);
    TIter next(tpcParticle);
    STParticle *restpart = NULL;

    while( (restpart = (STParticle*)next() ) ){
      LOG(DEBUG) << " piched id " << apart->GetTrackID()  << " and " << restpart->GetTrackID() << " @" << restpart << FairLogger::endl;

      if( restpart->GetTrackID() != apart->GetTrackID() && restpart->GetReactionPlaneFlag() > 1000 ) {
	Double_t wt_rp = restpart->GetRPWeight();
	TVector2 pt_rp = restpart->GetRotatedPt();

	mExcRP += wt_rp * pt_rp.Unit();
	itraex++;
      }
      else 
	LOG(DEBUG) << "rejected  track id @" << restpart << " " << restpart->GetTrackID() << FairLogger::endl; 
    }

    TVector3 rp_recv = TVector3(-999.,-999.,-999.);
    if(itraex > 0 && aflowcorrArray[0] != NULL) 
      rp_recv = Psi_FlatteningCorrection(0, ntrack[4] , TVector3(mExcRP.X(), mExcRP.Y(), 0.));
    else
      rp_recv =  TVector3(mExcRP.X(), mExcRP.Y(), 0.);
      
    apart->SetIndividualRPAngle( (Double_t)TVector2::Phi_mpi_pi( rp_recv.Phi() ));
    apart->SetAzmAngle_wrt_RP( (Double_t)TVector2::Phi_mpi_pi( apart->GetRotatedPt().Phi() - rp_recv.Phi()));
  }
}


void STSpiRITTPCTask::Clear()
{
  for (Int_t m = 0; m < 7; m++) ntrack[m] = 0;

  tpcParticle->Clear();

  fflowinfo->Clear();

  flowAnalysis->Clear();

  unitP_fc = TVector3(0.,0.,0.);
  unitP_rc = TVector3(0.,0.,0.);

  unitP    = TVector2(0.,0.);
  unitP_1  = TVector2(0.,0.);
  unitP_2  = TVector2(0.,0.);

  mtrack_1 = 0;
  mtrack_2 = 0;


  // bootstrap parameter
  for(UInt_t i = 0; i < 3; i++){
    bsPhi[i]   = -999.;
    bsPhi_1[i] = -999.;
    bsPhi_2[i] = -999.;
  }
}

TVector3 STSpiRITTPCTask::Psi_FlatteningCorrection(UInt_t isel, Int_t ival, TVector3 Pvec)
{
  Int_t    iBIN = GetCorrectionIndex(isel, ival, Pvec.Theta());

  TVector3 Psi_cf;
  if(iBIN >= 0 && aflowcorrArray[isel] != NULL){
    auto flowcorr = (STFlowCorrection*)(aflowcorrArray[isel]->At(iBIN));
    Psi_cf = flowcorr->ReCenteringFourierCorrection(Pvec);
    //    Psi_cf.SetPhi(flowcorr->GetCorrection(Pvec.Phi())); 
  }

  return Psi_cf;
}

TVector3 STSpiRITTPCTask::Psi_ReCenteringCorrection(UInt_t isel, Int_t ival, TVector3 Pvec)
{
  Int_t    iBIN = GetCorrectionIndex(isel, ival, Pvec.Theta());

  TVector3 Psi_cf;
  if(iBIN >= 0 && aflowcorrArray[isel] != NULL) {
    auto flowcorr = (STFlowCorrection*)aflowcorrArray[isel]->At(iBIN);

    Psi_cf = flowcorr->ReCentering(Pvec);
  }

  return Psi_cf;
}


Int_t STSpiRITTPCTask::GetCorrectionIndex(UInt_t isel, UInt_t ival, Double_t fval)
{
  Int_t index = GetMultiplicityCorretionIndex(isel, ival);
  index =  GetThetaCorrectionIndex(isel, index, fval);
  return index;
}

Int_t STSpiRITTPCTask::GetMultiplicityCorretionIndex(UInt_t isel, UInt_t ival)
{
  std::vector< std::pair<Double_t, Double_t> >::iterator itr;

  UInt_t ink = mtkbin[isel].size()-1;

  for(itr = pbinmin[isel].end()-1; itr != pbinmin[isel].begin()-1; itr-= mtkbin[isel].at(ink), ink--){

    if(isel == 0)
      LOG(DEBUG) << "mult " << itr->first<< " : " << itr->second << " -- " << mtkbin[isel].at(ink) << " at " 
		 << itr - pbinmin[isel].begin()   << FairLogger::endl;

    if(ival >= itr->first) {

      if(isel == 0)
	LOG(DEBUG) << " ntrack  " << ival << " : " << itr - pbinmin[isel].begin() << FairLogger::endl;

      return itr - pbinmin[isel].begin();
    }
  }

  return -1;
}

Int_t STSpiRITTPCTask::GetThetaCorrectionIndex(UInt_t isel, Int_t ival, Double_t fval)
{
  std::vector< std::pair<Double_t, Double_t> >::iterator itr;

  for(itr = pbinmin[isel].begin()+ival; itr != pbinmin[isel].begin()-1; itr--){

    if( fval >= itr->second) {
      return itr - pbinmin[isel].begin();
    }
  }
  return -1;
}

Bool_t STSpiRITTPCTask::SetupFlowDataBase()
{
  UInt_t version = (UInt_t)atoi(gSystem -> Getenv("VER"));
  Int_t ncount = 1;
  TString  fname[2];

  //  TString fname_psi    = Form("db/132Sn.v%d.psi.m%dn%d.txt",version, imtk, ihmm);
  //  TString fname_subpsi = Form("db/132Sn.v%d.subpsi.m%dn%d.txt",version, imtk, ihmm);

  //  fname[0] = Form("db/%3dSn.v%d.psi.",fflowinfo->SnA,version);
  //fname[0] = "run3062_rf.v11.0.Psi2rtcv0.";;
  fname[0] = "132Sn.v18.psi.";

  if( fIsSubeventAnalysis ) {
    //    fname[1] = Form("db/%03Sn.v%d.subpsi.",fflowinfo->SnA,version);
    //fname[1] = "run3062_rf.v11.0.Psis1rcv0.";
    fname[1] = "132Sn.v18.subpsi1.";
    ncount++;
  }


  for(UInt_t i = 0; i < TMath::Min(ncount, 2); i++){
    LOG(INFO) << " Database name is " << fname[i]  << FairLogger::endl;

    UInt_t ihmsum = 0;
    UInt_t imtk = 0;
    while(1){

      UInt_t ihm = 0;
      UInt_t ihmm = 0;
      while(1){

        TString ffname = fname[i] +  Form("m%dn%d.txt",imtk,ihmm);
        if( gSystem->FindFile("db",ffname) ){
          vfname[i].push_back(ffname);
          ihm++;  ihmm++;
	  LOG(INFO) << " Databse " << ffname << " is loaded. " << FairLogger::endl;
        }
        else
          break;
      }

      if(ihm > 0 ) {
        mtkbin[i].push_back(ihm);
        ihmsum += ihm;
      }
      else if(ihm == 0)
        break;

      imtk++;
    }


    auto nBin = (UInt_t)vfname[i].size();
    if(nBin == 0) return nBin;

    aflowcorrArray[i] = new TClonesArray("STFlowCorrection",30);
    TClonesArray &arr = *aflowcorrArray[i];

    LOG(INFO) << " nBin " << nBin << FairLogger::endl;

    std::vector<TString>::iterator itb;
    pbinmin[i].clear();
    UInt_t ihm = 0;
    
    STFlowCorrection *flowcorr;
    for(itb = vfname[i].begin(); itb != vfname[i].end(); itb++, ihm++){

      new(arr[ihm]) STFlowCorrection();

      flowcorr = (STFlowCorrection*)arr.At(ihm);
      flowcorr->SetRealOrMix(1);
      //      flowcorr->SetRealOrMix((UInt_t)bMix);
      flowcorr->SetFileName(vfname[i].at(ihm));
      flowcorr->LoadCorrectionFactor();

      pbinmin[i].push_back(make_pair(flowcorr->GetBin_min(0),flowcorr->GetBin_min(1)));
      LOG(INFO) << " $$$$$$$$$$$ ---->  nbin " << ihm  << " " << pbinmin[i].at(ihm).first << " : "<< pbinmin[i].at(ihm).second << FairLogger::endl;
    }

    //    binpara   = flowcorr->GetBinParameter(1);
  }

  auto nBin = (UInt_t)vfname[0].size();
  
  if( nBin == 0 ) return kFALSE;

  return kTRUE;
}



void STSpiRITTPCTask::ProceedEvent()
{

  fChain -> GetEntry(fEventID);  

  ntrack[0] = trackArray -> GetEntries();

  LOG(DEBUG) << "nttack[0] -------------------- > " << ntrack[0] << FairLogger::endl;
                 
  TIter next(trackVAArray);

  STRecoTrack *trackFromArray = NULL;

  while( (trackFromArray = (STRecoTrack*)next()) ) {

    TClonesArray &ptpcParticle = *tpcParticle;

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
      Double_t clustNum = VMom.Mag()>4000.?
	db->GetClusterNum(Charge, VMom.Theta(), VMom.Phi(), 4000.):
	db->GetClusterNum(Charge, VMom.Theta(), VMom.Phi(), VMom.Mag());
      //--- Set theoretical number of cluster
      aParticle->SetExpectedClusterNumber(clustNum);

      
      //--- Set extra quality flag ---; 
      SetupTrackExtraQualityFlag( aParticle );

      //--- Set up for flow ---;
      if( fIsFlowAnalysis && fflowinfo != nullptr) 
	DoFlowAnalysis( *aParticle );

      if( aParticle->GetGoodTrackFlag() > 1000 )
	ntrack[3]++;

      aParticle->SetTrackID(ntrack[2]);
      new(ptpcParticle[ntrack[2]]) STParticle(*aParticle);
      ntrack[2]++;

    }
  }
}
