#include "STFlowTask.hh"

STFlowTask::STFlowTask(Bool_t bfltn, Bool_t bsub, Bool_t bbst) :
  fIsFlowCorrection(bfltn),
  fIsSubeventAnalysis(bsub),
  fIsBootStrap(bbst),
  selReactionPlanef(1000)
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

  //  delete tpcParticle;  //!
  // if( fflowinfo != nullptr )
  //   delete fflowinfo;    //!

  // for(UInt_t i = 0; i < 2; i++) {
  //   if( aflowcorrArray[i] != nullptr )
  //     delete aflowcorrArray[i];
  // }

  // if( bs_unitP != nullptr )
  //   delete bs_unitP;

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
  
  Double_t N = (Double_t)(ntrack[4]-1);
  Double_t rpsigma = 1./(2.* N) * (sum_omg2/N / pow(sum_omg/N,2) );
  if( std::isinf( rpsigma ) ) rpsigma = 0;

  fflowinfo->SetRPSigma( sqrt(rpsigma) );
}


Bool_t STFlowTask::Init(UInt_t irun, TString sver)
{
  sVer = sver;
  iRun = irun;

  iSystem = 4;
  if(iRun >= 2841 && iRun <= 3039)
    iSystem = 0; // 132            
  else if(iRun >= 2261 && iRun <= 2509)
    iSystem = 1; // 108
  else if(iRun >= 3059 && iRun <= 3184)
    iSystem = 2; // 124
  else if(iRun >= 2520 && iRun <= 2653)
    iSystem = 3; // 112 


  fflowinfo->SetRun( irun );  
  
  if( fIsFlowCorrection) {
    if( !SetupFlowDataBase() ) {
      LOG(ERROR) << "STFlowTask:: Flow database cannot be setup" << FairLogger::endl;
      fIsFlowCorrection = kFALSE;
      // return kFALSE;
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
  //DoSubeventAnalysisFixedMultiplicity(20);

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

TVector3 STFlowTask::DoFlattening(TVector3 mvec, UInt_t ntrk)
{
  if( aflowcorrArray[0] == NULL ) 
    return TVector3(-999.,0.,0.);
      
  TVector3 rcvec =  Psi_FlatteningCorrection( 0, ntrk, mvec );
  
  return rcvec;
}


void STFlowTask::DoSubeventAnalysis()
{
  if( ntrack[4] == 0 ) return;
  
  STBootStrap* bs_unitP_1 = new STBootStrap(1);
  STBootStrap* bs_unitP_2 = new STBootStrap(1);

  UInt_t np = ntrack[4];
  if(np%2 == 1) np++;
  const UInt_t npart = np;
  UInt_t *rndArray = new UInt_t[npart];
  rndArray = RandomDivide2(npart);

  TIter next(tpcParticle);
  STParticle *apart = NULL;

  UInt_t itrack = 0;

  LOG(DEBUG) << " ntrack[4] " << ntrack[4] << FairLogger::endl;

  while( (apart = (STParticle*)next() ) ) {
    
    if( apart->GetReactionPlaneFlag() %2 == 1 ) {
      Double_t wt = apart->GetRPWeight();
      TVector2 ptr= apart->GetRotatedPt().Unit();

      if( rndArray[itrack] == 0 ) {
        apart->AddReactionPlaneFlag(100);

        fflowinfo->unitP_1 += wt * TVector3(ptr.X(), ptr.Y(), 0.);
	LOG(DEBUG) << " sub 1 " << fflowinfo->unitP_1.X()  
		   << " + "     << wt * ptr.X()
		   << " : "     << apart->GetReactionPlaneFlag()
		   << FairLogger::endl;

        if( fIsBootStrap )
          bs_unitP_1->Add(wt * ptr);

        fflowinfo->mtrack_1++;
      }
      else  {
        apart->AddReactionPlaneFlag(500);

        fflowinfo->unitP_2+= wt * TVector3(ptr.X(), ptr.Y(), 0.);
	LOG(DEBUG) << " sub 2    " << fflowinfo->unitP_2.X()  
		   << " + "     << wt * ptr.X()
		   << " : "     << apart->GetReactionPlaneFlag()
		   << FairLogger::endl;
        if( fIsBootStrap )
          bs_unitP_2->Add(wt * ptr);

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

void STFlowTask::DoSubeventAnalysisFixedMultiplicity(UInt_t val)
{
  if( ntrack[4] < val+5 ) return;
  
  STBootStrap* bs_unitP_1 = new STBootStrap(1);
  STBootStrap* bs_unitP_2 = new STBootStrap(1);


  UInt_t totaltrack = tpcParticle->GetEntries();
  const UInt_t npart = val;

  UInt_t *index = RandumPickUp(val, totaltrack);

  UInt_t *rndArray = new UInt_t[npart];
  rndArray = RandomDivide2(npart);

  STParticle *apart = NULL;

  UInt_t itrack = 0;
  
  for( UInt_t i = 0; i < (UInt_t)totaltrack; i++ ) {
  
    apart = (STParticle*)tpcParticle->At( *(index+i) );
    
    if( apart->GetReactionPlaneFlag() %2 == 1 ) {

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

UInt_t *STFlowTask::RandumPickUp(const UInt_t val, const UInt_t npart)
{
  Float_t vsort[npart];
  
  rnd.RndmArray(npart, vsort);
  
  std::vector< std::pair< Float_t, UInt_t > > psort;

  for( UInt_t i = 0; i < npart; i++ )    
    psort.push_back( std::make_pair(vsort[i], i) );

  std::sort( psort.begin(), psort.end() );

  std::vector< UInt_t > idxsort;
  for( UInt_t i = 0; i < npart; i++ ) 
    idxsort.push_back( psort[i].second ); 

  std::sort( idxsort.begin(), idxsort.end() );

  UInt_t *index = new UInt_t[npart];
  for( UInt_t i = 0; i < npart; i++ ) 
    *(index + i) = idxsort[i];
  
  return index;
}

UInt_t *STFlowTask::RandomDivide2(const UInt_t npart)
{
  UInt_t  *rndarray = new UInt_t[npart];
  UInt_t nmax = npart;

  if( npart%2 == 1 ) {
    nmax = npart -1 ;
    rndarray[npart-1] = 0;
  }

  UInt_t c1 = 0;
  UInt_t c2 = 0;
  UInt_t count = 0;
  while( count < nmax ){

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


void STFlowTask::SetupEventInfo(Long64_t eventid, UInt_t val)
{
  if( fflowinfo != nullptr) {

    fflowinfo->SetEventID( eventid );
    fflowinfo->beamPID = val;
    fflowinfo->SnA     = val;

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
    apart.SetReactionPlaneFlag(10);

  else if( pid > 2000 &&  apart.GetFromTargetFlag() ) { //fTargetf = fVatTargetf*fdistanceatvertexf;
    apart.SetReactionPlaneFlag(1000);

    if(apart.GetNDFFlag())
      apart.AddReactionPlaneFlag(10000);
  }
  else
    apart.SetReactionPlaneFlag(0);


  // Pt weight
  TLorentzVector lrnzVec =  apart.GetLorentzVector();
  
  TVector3 boostVec = STLorentzBoostVector::GetBoostVector(iSystem); 

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

    apart.AddReactionPlaneFlag(1);

    TVector2 upt = apart.GetRotatedPt().Unit();
    fflowinfo->unitP += apart.GetRPWeight() * TVector3( upt.X(), upt.Y(), 0.);


    sum_omg2 += pow(apart.GetRPWeight(), 2);
    sum_omg  += apart.GetRPWeight();


    if( fIsBootStrap )
      bs_unitP->Add(apart.GetRPWeight() * apart.GetRotatedPt().Unit());
  }

  if( apart.GetReactionPlaneFlag() >= selReactionPlanef && apart.GetReactionPlaneFlag()%2==1 )
    ntrack[5]++;  
}



void STFlowTask::DoIndividualReactionPlaneAnalysis( )
{
  TIter orgnext(tpcParticle);
  STParticle *apart = NULL;

  while( (apart = (STParticle*)orgnext()) ) {

    SetIndividualReactionPlane( *apart );

  }
}

void STFlowTask::SetIndividualReactionPlane( STParticle &apart )
{
  UInt_t itraex = 0;
  TVector3 mExcRP(0.,0.,0.);
  TIter next(tpcParticle);
  STParticle *restpart = NULL;

  auto befv = apart.GetIndividualRPAngle();

  while( (restpart = (STParticle*)next() ) ){

    if( restpart->GetTrackID() != apart.GetTrackID() && restpart->GetReactionPlaneFlag()%2 == 1 ) {
      Double_t wt_rp = restpart->GetRPWeight();
      TVector2 pt_rp = restpart->GetRotatedPt().Unit();
      
      mExcRP += wt_rp * TVector3( pt_rp.X(), pt_rp.Y(), 0.);

      LOG(DEBUG) << itraex << " x " << wt_rp*pt_rp.X() << " y " << wt_rp*pt_rp.Y() << " flag " << restpart->GetReactionPlaneFlag() << FairLogger::endl;

      itraex++;
    }
    else 
      LOG(DEBUG) << "rejected  track id @" << restpart << " " << restpart->GetTrackID() << FairLogger::endl; 
  }
  
  LOG(DEBUG) << " RP x " << mExcRP.X() << " <<- " << befv  << " num " << tpcParticle->GetEntries() << FairLogger::endl;

  if(itraex > 0 ) { 
    auto rcvec = DoFlattening( mExcRP, itraex );
    apart.SetIndividualRPVector( rcvec );
    apart.SetIndividualRPAngle( rcvec.Phi() );
    apart.SetAzmAngle_wrt_RP( (Double_t)TVector2::Phi_mpi_pi( apart.GetRotatedMomentum().Phi() - rcvec.Phi()));
  }
  else {
    apart.SetIndividualRPVector( TVector3(-999.,0.,0.) );
    apart.SetIndividualRPAngle( -9. );
    apart.SetAzmAngle_wrt_RP( -9. );
  }
}


void STFlowTask::Clear()
{
  ntrack[4] = 0;
  ntrack[5] = 0;

  fflowinfo->Clear();

  sum_omg2 = 0;
  sum_omg  = 0;
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

Bool_t STFlowTask::SetupFlowDataBase()
{
  TString  fname[2];
  Int_t    ncount = 1;
  TString  SNA = STRunToBeamA::GetBeamSnA(iRun);
  //  fname[0] = "132Sn.v"+sVer+".psi.";
  fname[0] = SNA + ".v"+sVer+".psi.";

  if( fIsSubeventAnalysis ) {
    //    fname[1] = "132Sn.v"+sVer+".subpsi1.";;
    fname[1] = SNA + ".v"+sVer+".subpsi1.";;
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



