#include "STBigRIPSTask.hh"

ClassImp(STBigRIPSTask);

STBigRIPSTask::STBigRIPSTask()
{
  auto anaRun = FairRunAna::Instance();
  iRun  = anaRun->getRunId();

  SetupBeamA(iRun);

  fEventID = -1;
  fIsPersistence = kTRUE;

  fLogger = FairLogger::GetLogger();

}

void STBigRIPSTask::Clear()
{
  fBDCArray->Clear();
}

InitStatus STBigRIPSTask::Init()
{

  fRootManager = FairRootManager::Instance();

  if ( fRootManager == nullptr) {
    LOG(ERROR) << "Cannot find RootManager!" << FairLogger::endl;
    return kERROR;
  }

  LOG(INFO) << "STBigRIPSTask::InitTask() is called " << FairLogger::endl;
  if( !SetupInputDataFile() ) {
    LOG(ERROR) << "STBigRIPSTask:: Cannot open input files" << FairLogger::endl;
    return kERROR;
  }


  if( !SetupBeamCut() ) {
    LOG(ERROR) << " Beam PID cut file is not opened. " << FairLogger::endl;
    return kERROR;
  }

  fBDCArray = new TClonesArray("STBDC");
  fRootManager -> Register("STBDC","flow",fBDCArray, 1);

  return kSUCCESS;
}

void STBigRIPSTask::SetupBeamA(UInt_t vRun)
{
  if(vRun >= 2174 && vRun <= 2509)
    SnA = 108;
  else if( vRun >= 2520 && vRun <= 2653)
    SnA = 112;
  else if( vRun >= 2836 && vRun <= 3039)
    SnA = 132;
  else if( vRun >= 3058 && vRun <= 3184)
    SnA = 124;
}


Bool_t STBigRIPSTask::SetupInputDataFile() 
{

  ribfChain = new TChain("TBeam");
  bdcChain  = new TChain("TBDC");

  TString beamDir;
  if(SnA == 132)
    beamDir = gSystem -> Getenv("STBEAM132");
  else if(SnA == 108)
    beamDir = gSystem -> Getenv("STBEAM108");
  else if(SnA == 124)
    beamDir = gSystem -> Getenv("STBEAM124");
  else if(SnA == 112)
    beamDir = gSystem -> Getenv("STBEAM112");

  TString beamFile = Form("beam_run%04d.ridf.root",iRun);


  if( gSystem->FindFile(beamDir, beamFile) ) 
    LOG(INFO) << " BigRIPS data : " << beamFile << " is opened." << FairLogger::endl;

  else {
    LOG(ERROR) << "RIDF data was not found. " << FairLogger::endl;
    return kFALSE;
  }

  ribfChain-> Add(beamFile);
  bdcChain -> Add(beamFile);

  //----- Set branch addresses.       

  ribfChain->SetBranchAddress("neve",&neve);
  ribfChain->SetBranchAddress("z",&z);
  ribfChain->SetBranchAddress("aoq",&aoq);
  ribfChain->SetBranchAddress("beta",&beta);
  ribfChain->SetBranchAddress("brho",&brho);
  ribfChain->SetBranchAddress("isGood",&isGood);
  ribfChain->SetBranchAddress("intZ",&intZ);
  ribfChain->SetBranchAddress("intA",&intA);

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


  return kTRUE;
}

void STBigRIPSTask::ProceedEvent()
{
  ribfChain->GetEntry(fEventID);
  bdcChain->GetEntry(fEventID);

  
  fBDC = new STBDC();
  
  fBDC->evt    = (Long64_t)neve;
  fBDC->SnA    = SnA;
  fBDC->aoq    = aoq ;
  fBDC->z      = z ;
  fBDC->tof    = tof ;
  fBDC->beta   = beta ;
  fBDC->brho   = brho ;
  fBDC->isGood = isGood ;
  fBDC->intZ   = intZ ;
  fBDC->intA   = intZ ;

  fBDC->bdcax  = bdcax ;
  fBDC->bdcby  = bdcby ;
  fBDC->ProjX  = ProjX ;
  fBDC->ProjY  = ProjY ;
  fBDC->ProjZ  = ProjZ ;
  fBDC->ProjP  = ProjP ;
  fBDC->ProjPX = ProjPX ;
  fBDC->ProjPY = ProjPY ;
  fBDC->ProjPZ = ProjPZ ;
  fBDC->ProjA  = ProjA ;
  fBDC->ProjB  = ProjB ;

  fBDC->SetRun( iRun ); 
  fBDC->SetEventID( neve );

  beamPID = GetBeamPID();
  fBDC->SetBeamPID( beamPID );

  TClonesArray &sBDC = *fBDCArray;
  new( sBDC[0] ) STBDC(*fBDC);

}


Bool_t STBigRIPSTask::SetupBeamCut()
{
  if(SnA == 132){
    auto gcutFile = new TFile("data/gcut132Sn.ROOT");
    g132Sn=(TCutG*)gcutFile->Get("sigma20");
    gcutFile->Close();

    if(g132Sn != NULL) LOG(INFO) << "gcut132Sn with sigma20 from data/gcut132Sn.ROOT" <<FairLogger::endl;
    else {
      LOG(ERROR) << "data/gcut132Sn.ROOT was not found" << FairLogger::endl;
      return kFALSE;
    }
  }

  else if(SnA == 108){
    auto gcutFile = new TFile("data/gcut108Sn.ROOT");
    g108Sn=(TCutG*)gcutFile->Get("sigma20");
    gcutFile->Close();

    if(g108Sn != NULL) LOG(INFO) << "gcut108Sn with sigma20 from data/gcut108Sn.ROOT" <<FairLogger::endl;
    else {
      LOG(ERROR) << "data/gcut108Sn.ROOT was not found" << FairLogger::endl;
      return kFALSE;
    }
  }
  else if(SnA == 124){
    auto gcutFile = new TFile("data/gcut124Sn.ROOT");
    g124Sn=(TCutG*)gcutFile->Get("sigma20");
    gcutFile->Close();

    if(g124Sn != NULL) LOG(INFO) << "gcut124Sn with sigma20 from data/gcut124Sn.ROOT" <<FairLogger::endl;
    else {
      LOG(ERROR) << "data/gcut124Sn.ROOT was not found" << FairLogger::endl;
      return kFALSE;
    }
  }
  else if(SnA == 112){
    auto gcutFile = new TFile("data/gcut112Sn.ROOT");
    g112Sn=(TCutG*)gcutFile->Get("sigma20");
    gcutFile->Close();

    if(g112Sn != NULL) LOG(INFO) << "gcut112Sn with sigma20 from data/gcut112Sn.ROOT" <<FairLogger::endl;
    else {
      LOG(ERROR) << "data/gcut112Sn.ROOT was not found" << FairLogger::endl;
      return kFALSE;
    }
  }

  return kTRUE;
}

Int_t STBigRIPSTask::GetBeamPID()
{

  Int_t pid = 0;
  if(g132Sn != NULL && SnA == 132){
    if(g132Sn->IsInside(aoq,z) && z < 50.536)
      pid = 132;
  }

  if(g108Sn != NULL && SnA == 108){
    if(g108Sn->IsInside(aoq,z))
      pid = 108;
  }

  if(g124Sn != NULL && SnA == 124){
    if(g124Sn->IsInside(aoq,z))
      pid = 124;
  }

  if(g112Sn != NULL && SnA == 112){
    if(g112Sn->IsInside(aoq,z))
      pid = 112;
  }

  return pid;

}


void STBigRIPSTask::Exec(Option_t *opt)
{
  LOG(DEBUG) << "STBigRIPSTask::Exec is called " << FairLogger::endl;
  fEventID++;

  ProceedEvent();


  //  Clear();
}
