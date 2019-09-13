#include "openRunAna.h"

void OpenChain();

void openRunAna()
{
  gROOT->Reset();
  gStyle->SetOptStat(0);
    
  aRun = gSystem -> Getenv("RUN");
  sSuf = gSystem -> Getenv("SUFX");
  sVer = gSystem -> Getenv("VER");  

  if( aRun != "" || sSuf != "" || sVer != ""){ 
    std::cout << " system " << sysName
	      << " RUN -> " << aRun 
	      << " : Suffix -> " << sSuf 
	      << " : Ver --> " << sVer
	      << std::endl;
    OpenChain();
  }
}

void GetSystem(UInt_t ival)
{
  // system selection
  // 0: 132Sn + 124Sn : 2841 ~ 3039 
  // 1: 108Sn + 112Sn : 2261 ~ 2509
  // 2: 124Sn + 112Sn : 2520 - 2653
  // 3: 112Sn + 124Sn

  if(ival >= 2841 && ival <= 3039){
    isys = 0; // 132
    sysName = "132Sn";
  }    
  else if(ival >= 2261 && ival <= 2509){
    isys = 1; // 108
    sysName = "108Sn";
  }
  else if(ival >= 3059 && ival <= 3184){
    isys = 2; // 124
    sysName = "124Sn";
  }
  else if(ival >= 2520 && ival <= 2653){
    isys = 3; // 112
    sysName = "112Sn";
  }
}

void OpenChain()
{
  
  cout << "aRUN.length " << aRun.Length() << endl;

  Int_t nrun = (aRun.Length())/5;

  printHeader = "run"+aRun(1,4) + "_" + sSuf + ".v" + sVer; 

  cout << " RUN -> " << aRun 
       << " total " << nrun 
       << " Print Header Name " << printHeader << endl;;


  vector<Int_t> lrun;
  Int_t ist = 1;
  while(ist < aRun.Length() - 4) {
    TString prun = aRun(ist,4);
    lrun.push_back( atoi(prun) );
    

    ist+=5;
  }

  GetSystem(lrun.at(0));
  cout << " system " << isys << " : " << sysName << endl;

  
  // set tree
  TString treename = "cbmsim";  


  // loading file
  TString fform = "_" +sSuf+ ".v" +sVer+ ".root";
 
  UInt_t ifnd = 0;
  rChain = new TChain(treename);

  for(Int_t i = 0; i < (Int_t)lrun.size(); i++){
    TString rootdir = "data/"; 

    TString fname = Form("run%d",lrun.at(i))+fform;
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


UInt_t SetBranch()
{
  if(rChain == NULL) {
    std::cout << " no file is loaded " << std::endl;
    return 0;
  }

  std::cout << " Nentry ->  " << rChain->GetEntries() << std::endl;

  if(aArray != NULL)
    aArray->Clear();

  rChain->SetBranchAddress("STParticle",&aArray);
  rChain->SetBranchAddress("STFlow"    ,&aFlowArray);

  return rChain->GetEntries();
}
