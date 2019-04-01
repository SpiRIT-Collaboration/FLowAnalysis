#include "STFlowTask.hh"

STFlowTask::STFlowTask(Bool_t bfltn, Bool_t bsub, Bool_t bbst) :
  fIsFlowCorrection(bfltn),
  fIsSubeventAnalysis(bsub),
  fIsBootStrap(bbst)
{

  TDatime dtime;
  rnd.SetSeed(dtime.GetSecond());  

  //------------------------//
  fflowinfo = new STFlowInfo();

  if( fIsBootStrap ) {
    bs_unitP = new STBootStrap(10);
    LOG(INFO) << "BootStrap is called. " << FairLogger::endl;
  }
}

STFlowTask::~STFlowTask()
{

  delete tpcParticle;  //!
  delete fflowinfo;    //!

  for(UInt_t i = 0; i < 2; i++) {
    if( aflowcorrArray[i] != nullptr )
      delete aflowcorrArray[i];
  }

  if( bs_unitP != nullptr )
    delete bs_unitP;

}


void STFlowTask::SetFlowTask( TClonesArray &atpcParticle )
{
  Clear();

  tpcParticle = &atpcParticle;
  TIter next(tpcParticle);

  STParticle* aParticle = NULL;
  
  while( (aParticle = (STParticle*)next() ) ) {
    
    SetupFlow( *aParticle );
    DoFlowAnalysis( *aParticle );

  }
}


Bool_t STFlowTask::Init(UInt_t irun, TString sver)
{
  sVer = sver;

  fflowinfo->SetRun( irun );  
  
  LOG(INFO) << "STFlowTask::Init() >> " << irun << " Ver. " << sver << FairLogger::endl;

  if( fIsFlowCorrection) {
    if( !SetupFlowDataBase() ) {
      LOG(ERROR) << "STFlowTask:: Flow database cannot be setup" << FairLogger::endl;
      fIsFlowCorrection = kFALSE;
    }
    else
      LOG(INFO) << "STFlowTask:: Flow database are ready. " << FairLogger::endl;
  }

  return kTRUE;
}

void STFlowTask::SetFlowInfo(STFlowInfo *aflowinfo) 
{
  fflowinfo = aflowinfo;

  ntrack[0] = fflowinfo->mtrack0;
  ntrack[1] = fflowinfo->mtrack1;
  ntrack[2] = fflowinfo->mtrack2;
  ntrack[3] = fflowinfo->mtrack3;
  ntrack[4] = fflowinfo->mtrack4;
  ntrack[5] = fflowinfo->mtrack5;
  ntrack[6] = fflowinfo->mtrack6;
}

void   STFlowTask::SetNTrack(UInt_t *nval)
{
  for(UInt_t i = 0; i < 7; i++) 
    ntrack[i] = *(nval + i);
}

Bool_t STFlowTask::SetupParameters()
{
  Bool_t fstatus = kTRUE;

  return fstatus;
}

void STFlowTask::FinishEvent()
{
  LOG(DEBUG) << "STFlowTask::FinishEvent is called. " << ntrack[4] << FairLogger::endl;

  fflowinfo->SetNTrack(ntrack);

  DoIndividualReactionPlaneAnalysis();

  if( fIsSubeventAnalysis ) 
    DoSubeventAnalysis();

  if( fIsFlowCorrection ) {
    if( !DoFlattening() )
      LOG(ERROR) << " Fail flatteing " << FairLogger::endl;
  }
}

Bool_t STFlowTask::DoFlattening()
{
  if( aflowcorrArray[0] != NULL ) {
    
    fflowinfo->unitP_fc = Psi_FlatteningCorrection( 0, ntrack[4], fflowinfo->unitP);
    fflowinfo->unitP_rc = Psi_ReCenteringCorrection(0, ntrack[4], fflowinfo->unitP);
    
    return kTRUE;
  }
 
  return kFALSE;
}

Bool_t STFlowTask::DoFlatteningSub()
{
  if( aflowcorrArray[1] != NULL ) {

    fflowinfo->unitP_1fc = Psi_FlatteningCorrection( 1, fflowinfo->mtrack_1, fflowinfo->unitP_1);
    fflowinfo->unitP_2fc = Psi_FlatteningCorrection( 1, fflowinfo->mtrack_2, fflowinfo->unitP_2);
    
    return kTRUE;
  }
  
  return kFALSE;
}

Bool_t STFlowTask::DoFlattening(STParticle &part)
{
  if( aflowcorrArray[0] != NULL ) {
    
    TVector3 rcvec =  Psi_FlatteningCorrection( 0, ntrack[4], part.GetIndividualRPVector() );

    part.SetIndividualRPAngle( rcvec.Phi() );
    part.SetAzmAngle_wrt_RP( (Double_t)TVector2::Phi_mpi_pi( part.GetRotatedPt().Phi() - rcvec.Phi()));
  
    return kTRUE;
  }

  return kFALSE;
}


void STFlowTask::DoSubeventAnalysis()
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
      TVector2 ptr= apart->GetRotatedPt().Unit();

      if( rndArray[itrack] == 0 ) {

        fflowinfo->unitP_1 += wt * TVector3(ptr.X(), ptr.Y(), 0.);
	LOG(DEBUG) << " sub 1 " << fflowinfo->unitP_1.X()  
		   << " + "     << wt * ptr.X()
		   << " : "     << apart->GetReactionPlaneFlag()
		   << FairLogger::endl;

        if( fIsBootStrap )
          bs_unitP_1->Add(wt * ptr);

        apart->AddReactionPlaneFlag(100);
        fflowinfo->mtrack_1++;
      }
      else  {
        fflowinfo->unitP_2+= wt * TVector3(ptr.X(), ptr.Y(), 0.);
	LOG(DEBUG) << " sub 2    " << fflowinfo->unitP_2.X()  
		   << " + "     << wt * ptr.X()
		   << " : "     << apart->GetReactionPlaneFlag()
		   << FairLogger::endl;
        if( fIsBootStrap )
          bs_unitP_2->Add(wt * ptr);

        apart->AddReactionPlaneFlag(200);
        fflowinfo->mtrack_2++;
      }

      itrack++;
      if( itrack > npart ) break;

    }
  }

  if(fIsBootStrap && fflowinfo->mtrack_1 > 0 && fflowinfo->mtrack_2 > 0 ) {
    bs_unitP_1->BootStrapping();
    bs_unitP_2->BootStrapping();
  }

  DoFlatteningSub();


  delete bs_unitP_1;
  delete bs_unitP_2;
  
}
UInt_t *STFlowTask::RandomDivide2(const UInt_t npart)
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


void STFlowTask::SetupEventInfo(Long64_t eventid, STBDC *aBDC)
{
  if( fflowinfo != nullptr) {

    fflowinfo->SetEventID( eventid );
    fflowinfo->SnA     = aBDC->SnA;
    fflowinfo->beamPID = aBDC->beamPID;

    LOG(DEBUG) << "flowinfo " << eventid << " * " << fflowinfo <<FairLogger::endl;
  }
}

void STFlowTask::SetupTrackExtraQualityFlag(STParticle *apart)
{

  // if( apart->GetClusterRatio() < 0.7 || apart->GetClusterRatio() > 2 ) 
  //   apart->SetClusterRatioFlag(0);
}


void STFlowTask::SetupFlow(STParticle &apart)
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

void STFlowTask::DoFlowAnalysis(STParticle &apart)
{
  
  SetupFlow( apart );

  if( apart.GetReactionPlaneFlag() >= selReactionPlanef ){
    ntrack[4]++;

    TVector2 upt = apart.GetRotatedPt().Unit();
    fflowinfo->unitP += apart.GetRPWeight() * TVector3( upt.X(), upt.Y(), 0.);

    if( fIsBootStrap )
      bs_unitP->Add(apart.GetRPWeight() * apart.GetRotatedPt().Unit());
  }

  if( apart.GetReactionPlaneFlag() == 20 )
    ntrack[5]++;  
}



void STFlowTask::DoIndividualReactionPlaneAnalysis()
{
  TIter orgnext(tpcParticle);
  STParticle *apart = NULL;

  while( (apart = (STParticle*)orgnext()) ) {

    Double_t wt = apart->GetRPWeight();
    TVector2 pt = apart->GetRotatedPt();

   // individual reaction plane
    UInt_t itraex = 0;
    TVector3 mExcRP(0.,0.,0.);
    TIter next(tpcParticle);
    STParticle *restpart = NULL;

    while( (restpart = (STParticle*)next() ) ){
      LOG(DEBUG) << " piched id " << apart->GetTrackID()  << " and " << restpart->GetTrackID() << " @" << restpart << FairLogger::endl;

      if( restpart->GetTrackID() != apart->GetTrackID() && restpart->GetReactionPlaneFlag() > 1000 ) {
	Double_t wt_rp = restpart->GetRPWeight();
	TVector2 pt_rp = restpart->GetRotatedPt().Unit();

	mExcRP += wt_rp * TVector3( pt_rp.X(), pt_rp.Y(), 0.);
	itraex++;
      }
      else 
	LOG(DEBUG) << "rejected  track id @" << restpart << " " << restpart->GetTrackID() << FairLogger::endl; 
    }

    apart->SetIndividualRPVector( mExcRP );

    if(itraex > 0 )
      DoFlattening( *apart );
      
  }
}

void STFlowTask::Clear()
{
  ntrack[4] = 0;
  ntrack[5] = 0;

  fflowinfo->Clear();
}

TVector3 STFlowTask::Psi_FlatteningCorrection(UInt_t isel, Int_t ival, TVector3 Pvec)
{
  Int_t    iBIN = GetMultiplicityCorretionIndex(isel, ival);

  TVector3 Psi_cf;
  if(iBIN >= 0 && aflowcorrArray[isel] != NULL){
    auto flowcorr = (STFlowCorrection*)(aflowcorrArray[isel]->At(iBIN));
    Psi_cf = flowcorr->ReCenteringFourierCorrection(Pvec);
  }

  LOG(DEBUG) << Pvec.X() << " iBIN : " << iBIN << " pvec.theta " << Pvec.Theta() << " Psi_cf " << Psi_cf.X() << FairLogger::endl;

  return Psi_cf;
}

TVector3 STFlowTask::Psi_ReCenteringCorrection(UInt_t isel, Int_t ival, TVector3 Pvec)
{
  Int_t    iBIN = GetMultiplicityCorretionIndex(isel, ival);

  TVector3 Psi_cf;
  if(iBIN >= 0 && aflowcorrArray[isel] != NULL) {
    auto flowcorr = (STFlowCorrection*)aflowcorrArray[isel]->At(iBIN);

    Psi_cf = flowcorr->ReCentering(Pvec);
  }

  return Psi_cf;
}

Int_t STFlowTask::GetMultiplicityCorretionIndex(UInt_t isel, UInt_t ival)
{
  // isel : 0 full Psi
  // isel : 1 Subevetn psi
  std::vector< std::pair<Double_t, Double_t> >::iterator itr;

  UInt_t ink = mtkbin[isel].size()-1;

  for(itr = pbinmin[isel].end()-1; itr != pbinmin[isel].begin()-1; itr--, ink--){

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

Int_t STFlowTask::GetThetaCorrectionIndex(UInt_t isel, Int_t ival, Double_t fval)
{
  std::vector< std::pair<Double_t, Double_t> >::iterator itr;

  for(itr = pbinmin[isel].begin()+ival; itr != pbinmin[isel].begin()-1; itr--){

    if( fval >= itr->second) {
      return itr - pbinmin[isel].begin();
    }
  }
  return -1;
}


Bool_t STFlowTask::SetupFlowDataBase()
{
  TString  fname[2];
  Int_t    ncount = 1;
  fname[0] = "132Sn.v"+sVer+".psi.";

  if( fIsSubeventAnalysis ) {
    fname[1] = "132Sn.v"+sVer+".subpsi1.";;
    ncount++;
  }

  

  for(UInt_t i = 0; i < TMath::Min(ncount, 2); i++){
    LOG(INFO) << " Database name is " << fname[i]  << " / " << ncount << " ( " << fIsSubeventAnalysis << FairLogger::endl;

    UInt_t ihmsum = 0;
    UInt_t imtk = 0;
    while(1){

      TString ffname = fname[i] +  Form("m%d.txt",imtk);
      if( gSystem->FindFile("db",ffname) ){
	vfname[i].push_back(ffname);
	LOG(INFO) << " Databse " << ffname << " is loaded. " << FairLogger::endl;
      }
      else
	break;

      mtkbin[i].push_back(imtk);

      imtk++;
    }


    auto nBin = (UInt_t)vfname[i].size();
    if(nBin == 0) return kFALSE;

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
      flowcorr->SetFileName(vfname[i].at(ihm));
      flowcorr->LoadCorrectionFactor();

      pbinmin[i].push_back(make_pair(flowcorr->GetBin_min(0),flowcorr->GetBin_min(1)));
      LOG(INFO) << " $$$$$$$$$$$ ---->  nbin " << ihm  << " " 
		<< pbinmin[i].at(ihm).first << " : "
		<< pbinmin[i].at(ihm).second << FairLogger::endl;
    }
  }

  auto nBin = (UInt_t)vfname[0].size();
  
  if( nBin == 0 ) return kFALSE;

  return kTRUE;
}



