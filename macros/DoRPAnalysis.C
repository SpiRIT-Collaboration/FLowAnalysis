#include "openRunAna.C"
#include "DoRPAnalysis.h"
#include  "../analysisformat/STRunToBeamA.hh"
//----------------------------------------------------------------------
// DoFlowAnalysis.C 
// input: 
// output:
// (c) Mizuki Kurata-Nishimura 
//----------------------------------------
//----------------------------------------------------------------------
void DoRPAnalysis()
{
  openRunAna(0);

  SetEnvironment();

  Open();

  // I/O
  Bool_t bCorr = kTRUE;
  if( bRedo ) bCorr = kFALSE;

  if( bCorr )
    std::cout << " Flattening will be done. " << std::endl;

  auto fflowtask = new STFlowTask(bCorr, kTRUE, kFALSE);
  auto fIsFlowAnalyais = fflowtask->Init(iRun, dVer);

  if( !fIsFlowAnalyais ) {
    std::cout << " A flow correction dbase is not fond. " << std::endl;
    exit(0);
  }
  

  OutputTree();

  Long64_t nEntry = SetBranch();

  for(Long64_t ievt = 0; ievt < nEntry; ievt++){
    
    ShowProcess(ievt);
    rChain->GetEntry(ievt);
    
    //------ Event selection
    aFlowInfo = (STFlowInfo*)aFlowArray->At(0);
    if( aFlowInfo == NULL ) continue;
    beamPID = aFlowInfo->beamPID;

    if( iVer[0] > 38 ) {
      if( aFlowInfo->goodEventf == 0 || beamPID == 0 )
	continue;
    }
    else
      if( beamPID == 0 ) continue;
    
  
    //------ end of event selection
    fflowtask->SetFlowInfo( aFlowInfo );
    fflowtask->SetRPMidRapidityCut(dMct);
    
    fflowtask->SetParticleArray( *aArray );

    if( !bRedo ) { 

      fflowtask->DoFlattening();
      fflowtask->DoFlatteningSub();
      fflowtask->DoIndividualReactionPlaneAnalysis();

    }
    else {  // Redo Reaction plane calculation

      fflowtask->SetPIDSelection(0);
      fflowtask->SetFlowTask();
      fflowtask->FinishEvent();
    }
  
    aFlowInfo = fflowtask->GetFlowInfo();

    TClonesArray &aflow = *anewFlow;
    new( aflow[0] ) STFlowInfo( *aFlowInfo );

    mflw->Fill();
  }


  
  fout->cd();
  fout->Write("",TObject::kWriteDelete);
  if( fv1y != NULL ) fv1y->Write();
  if( fv2y != NULL ) fv2y->Write();
  std::cout << fout->GetName() << std::endl;

  if(gROOT->IsBatch()) {
    fout->Close();

    delete fflowtask;
    exit(0);
  }
}

void SetEnvironment()
{

  TString sMct   = gSystem -> Getenv("MDCUT"); // mid-rapidity cut abs(y_cm)< MDCUT
  if( sMct != "")
    dMct = atof(sMct);

  TString sRedo  = gSystem -> Getenv("REDO");
  
  if(sRun =="" || sSuf == "" || sVer == "" || !DefineVersion()) {
    cout << " Please type " 
	 << " sRun = " << sRun
	 << " sSuf = " << sSuf
	 << " sVer = " << sVer
	 << endl;
    cout << "Check settings of environments " << endl;
    exit(0);
  }

  finname  = "run"+sRun+"_"+sSuf+".v"+sVer + ".root";
  foutname = "run"+sRun+"_"+sSuf+".v"+oVer + ".root"; 

  // set re calculation flag
  bRedo = sRedo=="1" ? kTRUE : kFALSE;

  // Print Run configuration 
  cout << " ---------- CONFIGURATION ----------  " << endl;
  cout << " Input file  -> " << finname 
       << " Output file -> " << foutname
       << " flow correction database version -> v." << dVer
       << " ReDo ? " << bRedo
    
       << endl;
}

Bool_t DefineVersion()
{
  Bool_t bfound = kFALSE;

  TString ver = sVer + ".";
  
  for ( Int_t i = 0; i < 2; i++) {
    if( ver.First(".") > 0 ) {

      Ssiz_t end = ver.First(".")  ;
      TString ver1 = ver(0, end);
      
      ver = ver(end+1, ver.Length());
      
      iVer[i] = atoi(ver1);


      if(i==1) bfound = kTRUE;
      break;
    }
  }
   
  if ( !bfound ) {
    iVer[0] = atoi( sVer );
    iVer[1] = 0;
    std::cout << " Input ersion number : v" << iVer[0] << "  " << iVer[1]  << std::endl;

    if( iVer[0] >= 0 ) bfound = kTRUE;
  }

  return bfound;
}


void OutputTree()
{
  //@@@
  // if( bRedo ) 
  //  fname = Form("run"+sRun+"_"+sSuf+"R.v%d",iVer[0]);
  //fname += Form(".%d.root",iVer[1]);
  TString fo = foutname;

  if( !gROOT->IsBatch() && gSystem->FindFile("data/",fo) ) {
    std::cout << fo << " is existing. Do you like to overwrite? (y/n)" << std::endl;
    TString sAns;
    std::cin >> sAns;
    if(sAns == "y" || sAns == "Y")
      fout  = new TFile(fo,"recreate");
    else {
      cout << " Retry" << endl;
      exit(0);
    }      
  }
  else {
    std::cout << "Output file is " << foutname << endl;
    foutname = "data/" + foutname;
    
    fout  = new TFile(foutname,"recreate");
  }

  mflw = new TTree("cbmsim","Flow corrected");


  //-- output                                                                                                              
  if( aBDC != NULL && beamPID != 100)  mflw->Branch("STBDC",&aBDC);
  if( aArray != NULL)   mflw->Branch("STParticle",&aArray);
  if( anewFlow != NULL) mflw->Branch("STFlow",&anewFlow);
  if( isys  == 5 )      mflw->Branch("RPPsi" ,&RPPsi);
    
 
}

void Open()
{ 
  //  fname = Form("run"+sRun+"_"+sSuf+".v%d",iVer[0]);

  TString fn = finname;
  std::cout << fn << std::endl;

  if( !gSystem->FindFile("data/", fn) ) {
    std::cout << finname << " : Input file desn't exist. "  << std::endl;
    exit(0);
  }

  if( finname == foutname ) {
    std::cout << "[NOTICE] An input and an ouput files is the sane name." << std::endl;
    
    gSystem->Rename(fn, "data/tmp_"+finname);
    fn = "data/tmp_"+finname;
  }

  rChain = new TChain("cbmsim");
  rChain->Add(fn);
  
  if(!rChain) {
    cout << " No data was found " << fn << endl;
    exit(0);
  }
  cout << "Input file is " << fn << endl;

  aBDC     = new TClonesArray("STBDC", 1);
  aArray   = new TClonesArray("STParticle",100);
  anewFlow = new TClonesArray("STFlowInfo",1);

  
  rChain->SetBranchAddress("STParticle",&aArray);
  rChain->SetBranchAddress("STFlow",    &aFlowArray);

  if( aArray == NULL ) 
    LOG(ERROR) << " Particle is not active " << FairLogger::endl;

  
  if( isys < 4)
    rChain->SetBranchAddress("STBDC"     ,&aBDC);
  else if( isys == 5 ) { 
    rChain->SetBranchAddress("RPPsi"     ,&RPPsi);
    auto fin = TFile::Open(fn);
    fv1y = (TF1*)fin->Get("fv1y");
    fv2y = (TF1*)fin->Get("fv2y");
    fin->Close();
  }
}

