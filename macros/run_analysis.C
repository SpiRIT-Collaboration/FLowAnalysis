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
  TString tDir = gSystem -> Getenv("TPCDIR");
  TString tVer = gSystem -> Getenv("RCVER");
  TString bDir = gSystem -> Getenv("BDCDIR");

  if( sRun=="" || sVer=="" || tDir=="" ) {
    std::cout << "Plase type " << std::endl;
    std::cout << "$ RUN=#### VER=# TPCDIR= RCVER= root run_analysis.C" << std::endl;
    exit(0);
  }


  FairRunAna* anaRun = new FairRunAna();
  anaRun->SetRunId(atoi(sRun));

  FairLogger *logger = FairLogger::GetLogger();
  logger ->SetLogToScreen(true);

  TString sOut;
  auto BDCTask     = new STBigRIPSTask();  // This should be called earlier than TPC
  if( BDCTask != NULL ) sOut += "B";

  auto TPCTask     = new STSpiRITTPCTask();
  TPCTask->SetRunInfo(tDir, tVer);
  if( TPCTask!=NULL ) sOut += "T";
  
  TPCTask->SetFlowAnalysis(kTRUE);  // Flow analysis is activated.


  TString foutname = "data/run"+sRun+"_"+sOut+".v"+sVer+".root";
  anaRun->SetOutputFile(foutname);


  anaRun->AddTask(BDCTask);
  anaRun->AddTask(TPCTask);
  
  anaRun->Init();

  anaRun->Run(0, 1e+10);

  cout << " num " << anaRun->GetNTasks() << " " << foutname << endl;
  gApplication -> Terminate();

}

void temp()
{

  //  auto BigRIPSTask = new STBigRIPSTask(sRun, sDir, SnA);
  //  anaManager->AddTask(BigRIPSTask);




  //  anaManager->Run(0);

  //  gApplication->Teminate();

  //  ST_ClusterNum_DB* db = new ST_ClusterNum_DB();
}

