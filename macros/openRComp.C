TChain *LChain[4];
static UInt_t id = 0;
 
TString printHeader="";
TString aRun[2];
TString sDB[2];
UInt_t  isys[4]={9, 9, 9, 9};


UInt_t OpenChain(UInt_t m = 0);

void openRComp()
{
  gROOT->Reset();
    
  aRun[0] = gSystem -> Getenv("RUN0");
  
  sDB[0] = gSystem -> Getenv("DB0");

  aRun[1] = gSystem -> Getenv("RUN1");
  
  sDB[1] = gSystem -> Getenv("DB1");

  UInt_t ichain = 0;

  if( aRun[0] != "" && sDB[0] != ""){ 
    isys[ichain] = OpenChain(ichain);
    if(isys[ichain] < 10) ichain++;
  }
  else
    std::cout << " RUN0 -> " << aRun[0] << " : DB0 -> " << sDB[0] << std::endl;

  if( aRun[1] != "" && sDB[1] != "") 
    isys[ichain] = OpenChain(ichain);
  
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
  
  if(sDB[0](1,4) == "flw_") // output of Assemble_flwv4.C run2841_flw_v4.1.root 
    treename = "flw";

  else if(sDB[0](14,4) == "crdb" ) // output of getFalatten2DCorrected 
      treename = "cflw";

  else if(sDB[0](1,6) == "rdflw_")  // output of AsmFlw_getEvent_v2.C run2841_rdflw_v4.1.0.root
    treename = "rflw";

  else if(sDB[0](1,6) == "mxflw_") // output of AsmFlw_getEvent_v2.C run2841_mxflw_v4.1.0.root
    treename = "mflw";




  TString fform = sDB[m] + ".root";
 
  LChain[m] = new TChain(treename);

  for(Int_t i = 0; i < (Int_t)lrun.size(); i++){
    TString rootdir = "../data"; 

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
    system = 0;
  else if(lrun.at(0) >= 2261 && lrun.at(0) <= 2509)
    system = 1;
  else if(lrun.at(0) >= 2520 && lrun.at(0) <= 2653)
    system = 2;


  return system;
}


void SaveCanvas(TString fopt = "")
{
  Int_t iCanvas = gROOT->GetListOfCanvases()->GetEntries();  
  for(Int_t i = 0; i < iCanvas; i++)
    gROOT->GetListOfCanvases()->At(i)->SaveAs(printHeader+fopt+Form("_%d",i)+".png");
}


