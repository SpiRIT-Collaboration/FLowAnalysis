#include "STFlowTask.hh"

STFlowTask::STFlowTask(Bool_t bfltn, Bool_t bsub) 
{
  STFlowTask(bfltn, bsub, kFALSE);
}


STFlowTask::STFlowTask(Bool_t bfltn, Bool_t bsub, Bool_t bbst) :
  fIsFlowCorrection(bfltn),
  fIsSubeventAnalysis(bsub),
  fIsBootStrap(bbst),
  selReactionPlanef(1000),
  fRPMidCut(0.),
  PID_sel(0)
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
  tpcParticle = &atpcParticle;

  fflowinfo->AllClear();
  SetFlowTask();
}

void STFlowTask::SetFlowTask()
{
  sum_omg2 = 0;
  sum_omg  = 0;  
  ntrack[4] = 0;
  ntrack[5] = 0;
  ntrack[6] = 0;

  fflowinfo->Clear();
  fflowinfo->SetPIDSelection(PID_sel);

  TIter next(tpcParticle);

  STKParticle* aParticle = NULL;
  
  while( (aParticle = (STKParticle*)next() ) ) {

    // if( aParticle->GetDistanceAtVertex() <= 20 && aParticle->GetNumCluster() >=15 )
    //   ntrack[6]++;
    
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


  iSystem =  STRunToBeamA::GetSystemID(irun);
  LOG(INFO) << " Irun " << irun << " system " << iSystem << FairLogger::endl;

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
  
  fflowinfo->SetRPMidCut(fRPMidCut);
  //  fflowinfo->SetPIDSelection(PID_sel);

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

  LOG(DEBUG) << "STFlowTask::SetFlowInfo is called. ntrack[3]=" << ntrack[3] << FairLogger::endl;

}

void   STFlowTask::SetNTrack(UInt_t i, UInt_t val)
{
  fflowinfo->SetNTrack(i, val);
}
Bool_t STFlowTask::SetupParameters()
{
  Bool_t fstatus = kTRUE;

  return fstatus;
}

void STFlowTask::FinishEvent()
{
  for( auto i : {4,5,6}  )
    fflowinfo->SetNTrack(i,ntrack[i]);

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
 
    //fflowinfo->unitP_fc = Psi_FlatteningCorrection( 0, ntrack[2], fflowinfo->unitP);
    //fflowinfo->unitP_rc = Psi_ReCenteringCorrection(0, ntrack[2], fflowinfo->unitP);
    fflowinfo->unitP_fc = Psi_FlatteningCorrection( 0, ntrack[4], fflowinfo->unitP);
    fflowinfo->unitP_rc = Psi_ReCenteringCorrection(0, ntrack[4], fflowinfo->unitP);

    auto psic = Psi_FlatteningCorrection( 2, ntrack[4], fflowinfo->unit2P);
	//fflowinfo->unit2P_fc = Psi_FlatteningCorrection( 2, ntrack[2], fflowinfo->unit2P);

    return kTRUE;
  }
 
  return kFALSE;
}

Bool_t STFlowTask::DoFlatteningSub()
{
  Bool_t bcorrcted = kFALSE;

  if( aflowcorrArray[1] != NULL ) {

    // original
    fflowinfo->unitP_1fc = Psi_FlatteningCorrection( 1, fflowinfo->mtrack_1, fflowinfo->unitP_1);
    fflowinfo->unitP_2fc = Psi_FlatteningCorrection( 1, fflowinfo->mtrack_2, fflowinfo->unitP_2);
    //fflowinfo->unitP_1fc = Psi_FlatteningCorrection( 1, ntrack[2]/2, fflowinfo->unitP_1);
    //fflowinfo->unitP_2fc = Psi_FlatteningCorrection( 1, ntrack[2]/2, fflowinfo->unitP_2);

    fflowinfo->cosdPsi  = cos(fflowinfo->unitP_1fc.Phi() - fflowinfo->unitP_2fc.Phi());

    bcorrcted = kTRUE;
  }

  if( aflowcorrArray[3] != NULL ) {

    fflowinfo->unit2P_1fc = Psi_FlatteningCorrection( 3, fflowinfo->mtrack_1, fflowinfo->unit2P_1);
    fflowinfo->unit2P_2fc = Psi_FlatteningCorrection( 3, fflowinfo->mtrack_2, fflowinfo->unit2P_2);
    
    fflowinfo->cos2dPsi = cos(2.*( fflowinfo->unit2P_1fc.Phi()/2. - fflowinfo->unit2P_2fc.Phi()/2. ));
    
    bcorrcted = kTRUE;
  }
  
  return bcorrcted;
}

// TVector3 STFlowTask::DoFlattening(TVector3 mvec, UInt_t ntrk)
// {
//   if( aflowcorrArray[0] == NULL ) 
//     return TVector3(-999.,0.,0.);
      
//   TVector3 rcvec =  Psi_FlatteningCorrection( 0, ntrk, mvec );
  
//   return rcvec;
// }

TVector3 STFlowTask::DoFlattening(UInt_t isel,  TVector3 mvec, UInt_t ntrk)
{
  if( aflowcorrArray[isel] == NULL ) 
    return TVector3(-999.,0.,0.);
      
  TVector3 rcvec =  Psi_FlatteningCorrection( isel, ntrk, mvec );
  
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
  STKParticle *apart = NULL;

  UInt_t itrack = 0;

  LOG(DEBUG) << " ntrack[4] " << ntrack[4] << FairLogger::endl;

  while( (apart = (STKParticle*)next() ) ) {
    
    if( apart->GetReactionPlaneFlag() %2 == 1 ) {
      Double_t wgt = apart->GetRPWeight();
      TVector3 rMom= apart->GetRotatedMomentum();

      if( rndArray[itrack] == 0 ) {
        apart->AddReactionPlaneFlag(100);

        fflowinfo->unitP_1  += wgt * TVector3(cos(   rMom.Phi()), sin(   rMom.Phi()), 0.);
	//        fflowinfo->unit2P_1 += wgt * TVector3(cos(2.*rMom.Phi()), sin(2.*rMom.Phi()), 0.);
        fflowinfo->unit2P_1 +=  TVector3(cos(2.*rMom.Phi()), sin(2.*rMom.Phi()), 0.);
	LOG(DEBUG) << " sub 1 " << fflowinfo->unitP_1.X()  
		   << " + "     << wgt * rMom.X()
		   << " : "     << apart->GetReactionPlaneFlag()
		   << FairLogger::endl;

        if( fIsBootStrap )
          bs_unitP_1->Add(wgt * rMom.Pt());

        fflowinfo->mtrack_1++;
      }
      else  {
        apart->AddReactionPlaneFlag(500);

        fflowinfo->unitP_2  += wgt * TVector3(cos(   rMom.Phi()), sin(   rMom.Phi()), 0.);
	//        fflowinfo->unit2P_2 += wgt * TVector3(cos(2.*rMom.Phi()), sin(2.*rMom.Phi()), 0.);
        fflowinfo->unit2P_2 +=  TVector3(cos(2.*rMom.Phi()), sin(2.*rMom.Phi()), 0.);

	LOG(DEBUG) << " sub 2    " << fflowinfo->unitP_2.X()  
		   << " + "     << wgt * rMom.X()
		   << " : "     << apart->GetReactionPlaneFlag()
		   << FairLogger::endl;
        if( fIsBootStrap )
          bs_unitP_2->Add(wgt * rMom.Pt());

        fflowinfo->mtrack_2++;
      }

      itrack++;
      if( itrack > npart ) break;
    }
  }

  auto pphi = fflowinfo->unit2P_1.Phi();
  //  fflowinfo->unit2P_1.SetPhi( pphi/2. );
  pphi = fflowinfo->unit2P_2.Phi(); 
  //  fflowinfo->unit2P_2.SetPhi( pphi/2. ); 


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

  STKParticle *apart = NULL;

  UInt_t itrack = 0;
  
  for( UInt_t i = 0; i < (UInt_t)totaltrack; i++ ) {
  
    apart = (STKParticle*)tpcParticle->At( *(index+i) );
    
    if( apart->GetReactionPlaneFlag() %2 == 1 ) {

      Double_t wgt = apart->GetRPWeight();
      TVector3 rMom = apart->GetRotatedMomentum();
      
      if( rndArray[itrack] == 0 ) {

	fflowinfo->unitP_1 += wgt * TVector3( cos(rMom.Phi()), sin(rMom.Phi()), 0.);
	if( fIsBootStrap )
	  bs_unitP_1->Add(wgt * rMom.Pt());
	
	apart->AddReactionPlaneFlag(100);
	fflowinfo->mtrack_1++;
      }
      else  {
	fflowinfo->unitP_2+= wgt * TVector3( cos(rMom.Phi()), sin(rMom.Phi()), 0.);
	if( fIsBootStrap )
	  bs_unitP_2->Add(wgt * rMom.Pt());
	
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

void STFlowTask::SetupTrackExtraQualityFlag(STKParticle *apart)
{

  // if( apart->GetClusterRatio() < 0.7 || apart->GetClusterRatio() > 2 ) 
  //   apart->SetClusterRatioFlag(0);
}


void STFlowTask::SetupFlow(STKParticle &apart)
{
  // Setup for flow analysis
  double yNN[]   = {0.3696, 0.3697, 0.3705, 0.3706};
  auto pid    =  apart.GetPID();

  
  if( pid == 211 )
    apart.SetReactionPlaneFlag(10);

  ///@@@@@@check
  else if( iSystem == 5 || // for simulation
	   // ( apart.GetBBMass() > 500 &&
	   //   apart.GetMass() > 0 ) ) {

	   ( apart.GetBBMass() > 500 &&
	     apart.GetMass() > 0 &&
	     apart.GetGoodTrackFlag() >= 100000 ) ){

    apart.SetReactionPlaneFlag(1001);
    ntrack[6]++;    

    // Pt weight
    auto rapidity =  apart.GetRapidity();

    auto rapiditycm = rapidity/yNN[iSystem]-1.;
    apart.SetRapiditycm( rapiditycm );
    
    if( rapiditycm  <  0 )
      apart.SetRPWeight(-1);
    else
      apart.SetRPWeight(1);
    
    if( abs( rapiditycm ) < fRPMidCut ) 
      apart.AddReactionPlaneFlag(3);

  }
  else
    apart.SetReactionPlaneFlag(0);

}

void STFlowTask::DoFlowAnalysis(STKParticle &apart)
{
  
  //  SetupFlow( apart );

  if( apart.GetReactionPlaneFlag()%2 == 1 ){
    ntrack[4]++;

    TVector3 rMom = apart.GetRotatedMomentum();
    fflowinfo->unitP  += apart.GetRPWeight() * TVector3( cos(   rMom.Phi()), sin(   rMom.Phi()), 0.);
    //    fflowinfo->unit2P += apart.GetRPWeight() * TVector3( cos(2.*rMom.Phi()), sin(2.*rMom.Phi()), 0.);
    fflowinfo->unit2P +=  TVector3( cos(2.*rMom.Phi()), sin(2.*rMom.Phi()), 0.);

    sum_omg2 += pow(apart.GetRPWeight(), 2);
    sum_omg  += apart.GetRPWeight();

    TVector2 ptv(apart.GetRotatedMomentum().X(), apart.GetRotatedMomentum().Y());
    if( fIsBootStrap )
      bs_unitP->Add(apart.GetRPWeight() * ptv.Unit());
  }
  else
    ntrack[5]++;
}



void STFlowTask::DoIndividualReactionPlaneAnalysis( )
{
  TIter orgnext(tpcParticle);
  STKParticle *apart = NULL;

  while( (apart = (STKParticle*)orgnext()) ) {

    //    SetIndividualReactionPlane_recal( *apart );
    SetIndividualReactionPlane( *apart );

  }
}
 
void STFlowTask::SetIndividualReactionPlane( STKParticle &apart )
{
  Double_t wgt = apart.GetRPWeight();
  TVector3 rMom = apart.GetRotatedMomentum();
      
  TVector3 uvec  = wgt * TVector3( cos(   rMom.Phi()), sin(   rMom.Phi()), 0.);
  //  TVector3 uvec2 = wgt * TVector3( cos(2.*rMom.Phi()), sin(2.*rMom.Phi()), 0.);
  TVector3 uvec2 = TVector3( cos(2.*rMom.Phi()), sin(2.*rMom.Phi()), 0.);

  TVector3 mExcRP  = fflowinfo->unitP;
  TVector3 mExc2RP = fflowinfo->unit2P;

  if( apart.GetReactionPlaneFlag()%2 == 1 ) {  // remove auto correlation
    mExcRP =  fflowinfo->unitP  - uvec;
    mExc2RP = fflowinfo->unit2P - uvec2;  // Psi_2
  }
  
  LOG(DEBUG) << " bef : Psi = " << mExcRP.Phi() << " " ;

  TVector3 rcvec = fIsFlowCorrection == kTRUE ?  DoFlattening(0, mExcRP, fflowinfo->mtrack4 ) : mExcRP ;        

  apart.SetIndividualRPVector( rcvec );
  apart.SetIndividualRPAngle( rcvec.Phi() );
  //  apart.SetIndividualRPAngle( (fflowinfo->unitP).Phi() );
  apart.SetAzmAngle_wrt_RP( (Double_t)TVector2::Phi_mpi_pi( apart.GetRotatedMomentum().Phi() - rcvec.Phi()));

  LOG(DEBUG) << " Psi = " << rcvec.Phi() 
	    << " relative " << (Double_t)TVector2::Phi_mpi_pi( apart.GetRotatedMomentum().Phi() - rcvec.Phi())
	    << FairLogger::endl;

  rcvec = fIsFlowCorrection == kTRUE ?  DoFlattening(2, mExc2RP, fflowinfo->mtrack4 ) : mExc2RP ;
  apart.SetIndividualRPVector2( rcvec );
  apart.SetIndividualRPAngle2 ( rcvec.Phi() );
  apart.SetAzmAngle2_wrt_RP( (Double_t)TVector2::Phi_mpi_pi( apart.GetRotatedMomentum().Phi() - rcvec.Phi()/2.) );
  
}

void STFlowTask::SetIndividualReactionPlane_recal( STKParticle &apart )
{
  UInt_t itraex = 0;
  TVector3 mExcRP(0.,0.,0.);
  TVector3 mExc2RP(0.,0.,0.);
  TIter next(tpcParticle);
  STKParticle *restpart = NULL;

  auto befv = apart.GetIndividualRPAngle();

  while( (restpart = (STKParticle*)next() ) ){

    //    if( restpart->GetTrackID() != apart.GetTrackID() && restpart->GetReactionPlaneFlag()%2 == 1 ) {
    if( restpart->GetReactionPlaneFlag()%2 == 1 ) {
      Double_t wt_rp = restpart->GetRPWeight();
      TVector2 pt_rp(restpart->GetRotatedMomentum().X(), restpart->GetRotatedMomentum().Y());
      pt_rp = pt_rp.Unit();
      

      if( restpart->GetTrackID() != apart.GetTrackID() ) {
	mExcRP  += wt_rp * TVector3( pt_rp.X(), pt_rp.Y(), 0.);
	mExc2RP += wt_rp * TVector3( cos(2.*pt_rp.Phi()), sin(2.*pt_rp.Y()), 0.);

	itraex++;
      }	
      else
	LOG(INFO) << ntrack[4] << " x " << wt_rp*pt_rp.X() << " y " << wt_rp*pt_rp.Y() << " flag " << restpart->GetReactionPlaneFlag() << FairLogger::endl;


    }
    else 
      LOG(DEBUG) << "rejected  track id @" << restpart << " " << restpart->GetTrackID() << FairLogger::endl; 
  }
  
  LOG(DEBUG) << " RP x " << mExcRP.X() << " <<- " << befv  << " num " << tpcParticle->GetEntries() << FairLogger::endl;

  if(itraex > 0 ) {
    TVector3 rcvec = fIsFlowCorrection == kTRUE ?  DoFlattening(0, mExcRP, itraex + 1 ) : mExcRP ;

    apart.SetIndividualRPVector( rcvec );
    apart.SetIndividualRPAngle( rcvec.Phi() );
    apart.SetAzmAngle_wrt_RP( (Double_t)TVector2::Phi_mpi_pi( apart.GetRotatedMomentum().Phi() - rcvec.Phi()));

    //v2
    TVector3 rvec = fIsFlowCorrection == kTRUE ?  DoFlattening(0, mExc2RP, itraex + 1 ) : mExc2RP ;
    //    rcvec.SetY(abs(rvec.Y()));
    
    apart.SetIndividualRPVector2( rcvec );
    apart.SetIndividualRPAngle2( rcvec.Phi() );
    apart.SetAzmAngle2_wrt_RP( (Double_t)TVector2::Phi_mpi_pi( 2.*(apart.GetRotatedMomentum().Phi() - rcvec.Phi())));

    LOG(DEBUG) << " phi " << rcvec.Phi() << " <<- " << mExcRP.Phi() << " " << fIsFlowCorrection << FairLogger::endl;
  }
  else {
    apart.SetIndividualRPVector( TVector3(-999.,0.,0.) );
    apart.SetIndividualRPAngle( -9. );
    apart.SetAzmAngle_wrt_RP( -9. );

    apart.SetIndividualRPVector2( TVector3(-999.,0.,0.) );
    apart.SetIndividualRPAngle2( -9. );
    apart.SetAzmAngle2_wrt_RP( -9. );
  }
}


void STFlowTask::Clear()
{
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

  // if( isel == 2 || isel == 3 ){
  //   Psi_cf.SetY( abs(Psi_cf.Y()) );
  // }


  if( isel == 3 )
    LOG(DEBUG) << Pvec.Phi() << " iBIN : " << iBIN << " pvec.theta " << Pvec.Theta() << " Psi_cf " << Psi_cf.X() << FairLogger::endl;

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
  // isel : 0 full Psi
  // isel : 1 Subevent subpsi
  // isel : 2 Subevent 2psi
  // isel : 3 Subevent sub2psi
  TString  fname[4];
  Int_t    ncount = 1;
  TString  SNA = STRunToBeamA::GetBeamSnA(iRun);
  fname[0] = SNA + ".v"+sVer+".psi.";
  ncount++;
  fname[2] = SNA + ".v"+sVer+".2psi.";;
  ncount++;

  if( fIsSubeventAnalysis ) {
    //    fname[1] = "132Sn.v"+sVer+".subpsi1.";;
    fname[1] = SNA + ".v"+sVer+".subpsi1.";;
    ncount++;
    fname[3] = SNA + ".v"+sVer+".sub2psi.";;
    ncount++;
  }
  
  for(UInt_t i = 0; i < TMath::Min(ncount, 4); i++){
    LOG(INFO) << " Database name is " << fname[i]  << " / " << ncount << " (" << fIsSubeventAnalysis << ")" <<FairLogger::endl;

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

    aflowcorrArray[i] = new TClonesArray("STFlowCorrection",22);
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



