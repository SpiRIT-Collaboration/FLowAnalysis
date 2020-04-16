#include "openRunAna.h"
#include  "../analysisformat/STRunToBeamA.hh"

FairRootManager* fman = FairRootManager::Instance();


void openRunAna(Bool_t bopen=kTRUE)
{
  gROOT->Reset();
  gStyle->SetOptStat(0);

  FairLogger* fLogger = FairLogger::GetLogger();
    
  sRun   = gSystem -> Getenv("RUN");         // RUN number
  sSuf   = gSystem -> Getenv("SUFX");
  sVer   = gSystem -> Getenv("IVER");         // input version
  oVer   = gSystem -> Getenv("OVER");        // output version
  dVer   = gSystem -> Getenv("DBVER");       // database version


  if( sRun != "" || sSuf != "" || sVer != "" || dVer != ""){ 
    LOG(INFO) << " RUN -> " << sRun 
	      << " : Suffix -> " << sSuf 
	      << " : Ver --> " << sVer
	      << " : DataBaser Version --> " << dVer
	      << FairLogger::endl;

    // Set RUN number
    iRun    =  atoi(sRun);
    isys    =  STRunToBeamA::GetSystemID(iRun);

    if( bopen )
      OpenChain();

  }

}


void OpenChain()
{
  
  cout << "RUN.length " << sRun.Length() << endl;

  Int_t nrun = (sRun.Length())/5;

  printHeader = "run"+sRun(1,4) + "_" + sSuf + ".v" + sVer; 

  LOG(INFO)   << " total " << nrun 
	      << " Print Header Name " << printHeader << FairLogger::endl;;


  vector<Int_t> lrun;
  Int_t ist = 1;
  while(ist < sRun.Length() - 4) {
    TString prun = sRun(ist,4);
    lrun.push_back( atoi(prun) );
    

    ist+=5;
  }


  isys    =  STRunToBeamA::GetSystemID(lrun.at(0));
  sysName =  STRunToBeamA::GetSystemName(lrun.at(0)); 
  cout << " system ID >> " << isys << " systemName >> " << sysName << endl;
 
  // set tree
  TString treename = "cbmsim";  

  // loading file
  TString fform = "_" +sSuf+ ".v" +sVer+ ".root";
 
  UInt_t ifnd = 0;
  rChain = new TChain(treename);

  for(Int_t i = 0; i < (Int_t)lrun.size(); i++){
    TString rootdir = "data/"; 

    TString fname = Form("run%04d",lrun.at(i))+fform;
    cout << fname << endl;

    if(gSystem->FindFile(rootdir,fname)) {
      rChain->Add(fname);
      ifnd++;
    }
    else 
      std::cout << " File is not found " << fname << std::endl;
  }
  

  if(rChain->GetListOfFiles()->GetEntries() == 0){
    std::cout << " No files are loaded. " << endl;    
    exit(0);
  }


  if( rChain!=NULL ){
    std::cout << "Runs are linked to rChain" 
	      << " " << rChain->GetListOfFiles()->GetEntries() << " / " << nrun
	      << std::endl;
    

    rChain->SetName("rChain");
    std::cout << " Entries :" << rChain->GetEntries() << std::endl;

  }
}


void SaveCanvas(TString fopt = "", Int_t isel=-1)
{
  if(isel > -1)
    gROOT->GetListOfCanvases()->At(isel)->SaveAs(printHeader+fopt+Form("_%d",isel)+".png");

  else {
    
    TString mdir = sSuf+fopt+gSystem->Now().AsString();
    gSystem->mkdir(mdir);
    gSystem->cd(mdir);
    
    Int_t iCanvas = gROOT->GetListOfCanvases()->GetEntries();  
    for(Int_t i = 0; i < iCanvas; i++)
      gROOT->GetListOfCanvases()->At(i)->SaveAs(sSuf+fopt+Form("_%d",i)+".png");

    gSystem->cd("..");
    std::cout << " Figures were saved in " << mdir << std::endl;

  }
}


Long64_t SetBranch()
{
  if(rChain == NULL) {
    std::cout << " no file is loaded " << std::endl;
    return 0;
  }

  if(aArray != NULL)
    aArray->Clear();

  rChain->SetBranchAddress("STParticle",&aArray);
  rChain->SetBranchAddress("STFlow"    ,&aFlowArray);

  if( isys == 5 )
    rChain->SetBranchAddress("RPPsi",&RPPsi);


  Long64_t totalevent = rChain->GetEntries() ;

  Long64_t maxevt = totalevent;
  TString smaxevt = gSystem -> Getenv("MXEVT");
  if( smaxevt != "" )
    maxevt = atoi(smaxevt);

  nEntry = maxevt < totalevent ? maxevt : totalevent;

  LOG(INFO) << " NOTICE !!!! Process to " << nEntry << " (" << totalevent << ")" << FairLogger::endl;

  return nEntry;
}


void ShowProcess(Long64_t ievt)
{
  TDatime dtime;
  static TDatime btime(dtime);

  Long64_t eprint = nEntry/10;

  if(ievt%eprint == 0) {
    dtime.Set();
    Int_t ptime = dtime.Get() - btime.Get();
    std::cout << "Process " << std::setw(8) << ievt << "/"<< nEntry << " = " 
	      << ((Double_t)(ievt)/(Double_t)nEntry)*100. << " % --->"
	      << dtime.AsString() << " ---- "
	      << std::setw(5) << (Int_t)ptime/60 << " [min] "
	      << std::endl;
  }
}
