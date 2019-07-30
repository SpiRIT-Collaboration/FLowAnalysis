#include "openFlw.h"

static UInt_t id = 0;
 

UInt_t OpenChain(UInt_t m = 0);

void openFlw()
{
  gROOT->Reset();
  gStyle->SetOptStat(0);
    
  for(UInt_t i = 0 ; i < nconfig; i++){
    TString form = Form("RUN%d",i);
    aRun[i] = gSystem -> Getenv(form);

    form = Form("DB%d",i);
    sDB[i] = gSystem -> Getenv(form);
    if( i > 0 && sDB[i] == "")
      sDB[i] = sDB[0];
  }


  UInt_t ichain = 0;

  for(UInt_t i = 0; i < nconfig; i++){
    if( aRun[i] != "" && sDB[i] != ""){ 
      //      std::cout << " RUN" << i << "-> " << aRun[i] << " : DB " << i << "-> " << sDB[i] << std::endl;

      isys[ichain] = OpenChain(ichain);
      if(isys[ichain] < 10) ichain++;
    }
  }

}

UInt_t GetSystem(UInt_t ival)
{
  if(ival > 4)
    return 9;

  return isys[ival];
}

UInt_t OpenChain(UInt_t m)
{
  if( m > 4) exit(0);
  
  cout << "aRUN.length " << aRun[m].Length() << endl;

  Int_t nrun = (aRun[m].Length())/5;
  cout << aRun[m] << "-> nrun " << nrun << endl;;

  printHeader = "FlwRUN"+aRun[0](1,4)+ Form("m%d",nrun) + sDB[0]; 
  TString printName ;
  vector<Int_t> lrun;
  Int_t ist = 1;
  while(ist < aRun[m].Length() - 4) {
    TString prun = aRun[m](ist,4);
    lrun.push_back( atoi(prun) );

    ist+=5;
  }


  cout << "sDB[0] " << sDB[0] << " " << sDB[0](10,1) << endl;

  // set tree
  TString treename = "cflw";  
  if(sDB[0]( sDB[0].First("c"), 2 ) == "cv") // output of flw_process4.C  run2844_rf.v0.0.0.root  
      treename = "cflw";
  
  else if(sDB[0](1,2) == "f0")       // output of flw_process1.C is run2844_f0.v0.root 
    treename = "flw";

  else if(sDB[0](1,2) == "rf")  // output of flw_prosess2.C run2841_rf_v4.1.0.root
    treename = "rflw";

  else if(sDB[0](1,2) == "mf")  // output of flw_process2.C run2841_mf_v4.1.0.root
    treename = "mflw";


  else if(sDB[0](1,2) == "nl")
    treename = "flw";



  // loading file
  TString fform = sDB[m] + ".root";
 
  UInt_t ifnd = 0;
  rChain[m] = new TChain(treename);

  for(Int_t i = 0; i < (Int_t)lrun.size(); i++){
    TString rootdir = "data/"; 

    TString fname = Form("run%d",lrun.at(i))+fform;
    cout << fname << endl;

    if(gSystem->FindFile(rootdir,fname)) {
      rChain[m]->Add(fname);
      ifnd++;
    }
    else 
      std::cout << " File is not found " << fname << std::endl;
  }
  

  if(rChain[m]->GetListOfFiles()->GetEntries() == 0){
    std::cout << " No files are loaded. " << endl;    
    exit(0);
  }


  if( rChain[m]!=NULL ){
    std::cout << "Runs are linked to rChain"<< m  
	      << " " << rChain[m]->GetListOfFiles()->GetEntries() << " / " << nrun
	      << std::endl;
    

    rChain[m]->SetName(Form("rChain%d",m));
    std::cout << " Entries :" << rChain[m]->GetEntries() << std::endl;

    m_end = m+1;
  }


  // system selection
  // 0: 132Sn + 124Sn : 2841 ~ 3039 
  // 1: 108Sn + 112Sn : 2261 ~ 2509
  // 2: 124Sn + 112Sn : 2520 - 2653
  // 3: 112Sn + 124Sn
  UInt_t system = 10;
  if(lrun.at(0) >= 2841 && lrun.at(0) <= 3039){
    system = 0; // 132
    sysName[0] = "132Sn";
  }    
  else if(lrun.at(0) >= 2261 && lrun.at(0) <= 2509){
    system = 1; // 108
    sysName[1] = "108Sn";
  }
  else if(lrun.at(0) >= 3059 && lrun.at(0) <= 3184){
    system = 2; // 124
    sysName[2] = "124Sn";
  }
  else if(lrun.at(0) >= 2520 && lrun.at(0) <= 2653){
    system = 3; // 112
    sysName[3] = "112Sn";
  }

  return system;
}


void SaveCanvas(TString fopt = "", Int_t isel=-1)
{
  if(isel > -1)
    gROOT->GetListOfCanvases()->At(isel)->SaveAs(printHeader+fopt+Form("_%d",isel)+".png");

  else {
    
    TString mdir = sDB[0](4,sDB[0].Sizeof())+fopt+gSystem->Now().AsString();
    gSystem->mkdir(mdir);
    gSystem->cd(mdir);
    
    Int_t iCanvas = gROOT->GetListOfCanvases()->GetEntries();  
    for(Int_t i = 0; i < iCanvas; i++)
      gROOT->GetListOfCanvases()->At(i)->SaveAs(sDB[0](4,sDB[0].Sizeof())+fopt+Form("_%d",i)+".png");

    gSystem->cd("..");
    std::cout << " Figures were saved in " << mdir << std::endl;

  }
}


