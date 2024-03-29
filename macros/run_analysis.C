//----------------------------------------------------------------------
// Macro: run_analysis
//  (c) Mizuki Kurata-Nishimura
//   31 Jan. 2019
//----------------------------------------
// Assemble TPC reconstructed data, Kyoto,  Katana, BigRIPS, and NeuLAND
//----------------------------------------------------------------------
void run_analysis(Int_t nevt = -1)
{
  //Reset ROOT and connect tree file
  gROOT->Reset();
  gROOT->ProcessLine("gErrorIgnoreLevel = kBreak");
  gROOT->ProcessLine("gErrorIgnoreLevel = kSysError");
  gROOT->ProcessLine("gErrorIgnoreLevel = kFatal");

  TString sRun = gSystem -> Getenv("RUN");
  TString sVer = gSystem -> Getenv("VER");
  TString dVer = gSystem -> Getenv("DBVER");
  TString tDir = gSystem -> Getenv("TPCDIR");
  TString tVer = gSystem -> Getenv("RCVER");
  TString bDir = gSystem -> Getenv("BDCDIR");
  TString ssfx = gSystem -> Getenv("SUFX");
  TString mxevt= gSystem -> Getenv("MXEVT");



  FairLogger* fLogger = FairLogger::GetLogger();
  if( sRun=="" || sVer=="" || tDir=="" ) {
    LOG(ERROR) << "Plase type " << FairLogger::endl;
    LOG(ERROR) << "$ RUN=#### VER=# TPCDIR= RCVER= root run_analysis.C" << FairLogger::endl;
    exit(0);
  }


  FairRunAna* anaRun = new FairRunAna();
  UInt_t vRun = atoi(sRun);
  anaRun->SetRunId(vRun);

  TString foutname = "data/run" + sRun + "_" + ssfx + ".v" + sVer + ".root";
  anaRun->SetOutputFile(foutname);
    
  FairRootManager* fman = FairRootManager::Instance();
  
  FairLogger *logger = FairLogger::GetLogger();
  logger ->SetLogToScreen(true);

  cout << sRun << endl;

  // UInt_t SnA = STRunToBeamA::GetBeamA(iRun);
  // tDir += Form("/Sn%d",SnA);

  // auto BDCTask     = new STBigRIPSTask();  // This should be called earlier than TPC
  // anaRun->AddTask(BDCTask);

  auto TPCTask     = new STSpiRITTPCTask();
  anaRun->AddTask(TPCTask);

  TPCTask->SetRunInfo(tDir, tVer, dVer);
  TPCTask->SetFlowAnalysis(1);  // Flow analysis is activated.
  //  TPCTask->SetFlowAnalysis(0);

  anaRun->Init();
  
  Long64_t maxevt = TPCTask->GetEntries();    
  if( mxevt != "" ) {
    nevt = atoi( mxevt );
    maxevt = nevt;
  }
  LOG(INFO) << " maxevent " << maxevt << FairLogger::endl;

  TPCTask->SetProcessingNumberOfEvent(maxevt);

  anaRun->Run(0, maxevt);

  cout << " num " << anaRun->GetNTasks() << " " << foutname << endl;
  gApplication -> Terminate();


  delete anaRun;
}

