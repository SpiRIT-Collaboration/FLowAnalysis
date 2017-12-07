#include "openFlw.h"

static UInt_t id = 0;
 

UInt_t OpenChain(UInt_t m = 0);

void openFlw()
{
  gROOT->Reset();
    
  for(UInt_t i = 0 ; i < nconfig; i++){
    TString form = Form("RUN%d",i);
    aRun[i] = gSystem -> Getenv(form);
    form = Form("DB%d",i);
    sDB[i] = gSystem -> Getenv(form);
  }


  UInt_t ichain = 0;

  for(UInt_t i = 0; i < nconfig; i++){
    if( aRun[i] != "" && sDB[i] != ""){ 
      std::cout << " RUN" << i << "-> " << aRun[i] << " : DB " << i << "-> " << sDB[i] << std::endl;

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
  if( m > 2) exit(0);
  
  Int_t nrun = (aRun[m].Length()-1)/5;
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

  for(Int_t i = 0; i < nrun; i++)
    cout << " lrun = " << lrun.at(i) << endl;

  
  TString treename = "cflw";
  
  if(sDB[0](1,2) == "f0")       // output of flw_process1.C is run2844_f0.v0.root 
    treename = "flw";

  else if(sDB[0](10,3) == "_db" ) // output of flw_process4.C  run2844_rf.v0.0.0_db2844.v0.0.cv0.root  
      treename = "cflw";

  else if(sDB[0](1,2) == "rf")  // output of flw_prosess2.C run2841_rf_v4.1.0.root
    treename = "rflw";

  else if(sDB[0](1,2) == "mf")  // output of flw_process2.C run2841_mf_v4.1.0.root
    treename = "mflw";



  TString fform = sDB[m] + ".root";
 
  LChain[m] = new TChain(treename);

  for(Int_t i = 0; i < (Int_t)lrun.size(); i++){
    TString rootdir = "data/"; 

    TString fname = Form("run%d",lrun.at(i))+fform;
    cout << fname << endl;

    if(gSystem->FindFile(rootdir,fname))
      LChain[m]->Add(fname);
    else
      std::cout << " File is not found " << fname << std::endl;
  }
  

  if(LChain[m]->GetListOfFiles()->GetEntries() == 0){
    std::cout << " No files are loaded. " << endl;    
    exit(0);
  }


  if( LChain[m]!=NULL ){
    std::cout << "Runs are linked to rChain"<< m  << std::endl;
    LChain[m]->SetName(Form("rChain%d",m));
    std::cout << " Entries :" << LChain[m]->GetEntries() << std::endl;
  }


  // system selection
  // 0: 132Sn + 124Sn : 2841 ~ 3039 
  // 1: 108Sn + 112Sn : 2261 ~ 2509
  // 2: 124Sn + 112Sn : 2520 - 2653
  // 3: 112Sn + 124Sn
  UInt_t system = 10;
  if(lrun.at(0) >= 2841 && lrun.at(0) <= 3039)
    system = 0; // 132
  else if(lrun.at(0) >= 2261 && lrun.at(0) <= 2509)
    system = 1; // 108
  else if(lrun.at(0) >= 3059 && lrun.at(0) <= 3184)
    system = 2; // 124
  else if(lrun.at(0) >= 2520 && lrun.at(0) <= 2653)
    system = 3; // 112

  return system;
}


void SaveCanvas(TString fopt = "")
{
  Int_t iCanvas = gROOT->GetListOfCanvases()->GetEntries();  
  for(Int_t i = 0; i < iCanvas; i++)
    gROOT->GetListOfCanvases()->At(i)->SaveAs(printHeader+fopt+Form("_%d",i)+".png");
}


