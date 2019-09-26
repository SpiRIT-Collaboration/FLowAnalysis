#include "DoRPAnalysis.h"
//----------------------------------------------------------------------
// DoFlowAnalysis.C 
// input: 
// output:
// (c) Mizuki Kurata-Nishimura 
//----------------------------------------

//
//----------------------------------------------------------------------
void DoRPAnalysis(Long64_t nmax = -1)
{
  SetEnvironment();

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

  Open();

  OutputTree();

  nEntry = rChain->GetEntries();
  for(Long64_t ievt = 0; ievt < nEntry; ievt++){

    PrintProcess(ievt);

    rChain->GetEntry(ievt);

    //------ Event selection
    auto aFlowInfo = (STFlowInfo*)aFlowArray->At(0);
    if( aFlowInfo == NULL ) continue;

    if( iVer[0] > 38 ) {
      if( aFlowInfo->goodEventf == 0 || aFlowInfo->beamPID == 0 )
	continue;
    }
    else
      if( aFlowInfo->beamPID == 0 ) continue;
    
    if( aFlowInfo->beamPID == 124 ){
      if( (aFlowInfo->mtrack1*0.8-20) > aFlowInfo->mtrack4 ) continue;
    }
    //------ end of event selection

    fflowtask->SetFlowInfo( aFlowInfo );
    fflowtask->SetParticleArray( *aParticleArray );

    if( !bRedo ) { 

      fflowtask->DoFlattening();
      fflowtask->DoFlatteningSub();

      TIter next(aParticleArray);
      STParticle *apart = NULL; 
      UInt_t idx = 0;
      while( (apart = (STParticle*)next()) ) {

	fflowtask->SetIndividualReactionPlane( *apart );
	TVector3 rpvec = apart->GetIndividualRPVector();
	apart->SetAzmAngle_wrt_RP( TVector2::Phi_mpi_pi( apart->GetRotatedMomentum().Phi() - rpvec.Phi() ) );

	TClonesArray &arr = *aParticleArray;
	apart = (STParticle*)arr.At(idx);
	
	//	cout << "aft1 ->" << apart->GetIndividualRPVector().Phi() << endl;

	idx++ ;
      }
    }

    else {  // Redo Reaction plane calculation

      fflowtask->SetFlowTask();
      fflowtask->FinishEvent();
    }
  
    aFlowInfo = fflowtask->GetFlowInfo();

    TClonesArray &aflow = *anewFlow;
    new( aflow[0] ) STFlowInfo( *aFlowInfo );


    mflw->Fill();

    if( ievt > nmax && nmax != -1 ) break;

  }

  fout->cd();
  fout->Write();
  std::cout << fout->GetName() << std::endl;

  if(gROOT->IsBatch()) {
    fout->Close();

    delete fflowtask;
    exit(0);
  }
}

void SetEnvironment()
{

  sRun   = gSystem -> Getenv("RUN");  // RUN number
  sSuf   = gSystem -> Getenv("SUFX");
  sVer   = gSystem -> Getenv("VER");  // Version ID
  dVer   = gSystem -> Getenv("DBVER");
  TString sRedo  = gSystem -> Getenv("REDO");
  
  if(sRun =="" || sSuf == "" || sVer == "" || !DefineVersion()) {
    cout << " Please type " << endl;
    cout << "$ RUN=#### VER=#.# SUFX=BTt root DoFlow_Analysis.C " << endl;
    exit(0);
  }

  finname  = "run"+sRun+"_"+sSuf+".v"+sVer + ".root";
  foutname = "run"+sRun+"_"+sSuf+".v"+dVer + ".root"; 

  // Set RUN number
  iRun = atoi(sRun);

  // set re calculation flag
  bRedo = sRedo=="1" ? kTRUE : kFALSE;

  // Print Run configuration 
  cout << " ---------- CONFIGURATION ----------  " << endl;
  cout << "RUN = "      << sRun   << " with ver "  << sVer  
       << "Flow correction ver "  << dVer   << " databases " 
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



void PrintProcess(Int_t ievt)
{
  TDatime dtime;
  static TDatime btime(dtime);

  UInt_t eprint = nEntry/10;

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

  aBDC = new TClonesArray("STBDC", 1);
  aParticleArray = new TClonesArray("STParticle",100);

  rChain->SetBranchAddress("STBDC"     ,&aBDC);
  rChain->SetBranchAddress("STParticle",&aParticleArray);
  rChain->SetBranchAddress("STFlow",    &aFlowArray);

}
void OutputTree()
{
  //@@@
  // if( bRedo ) 
  //  fname = Form("run"+sRun+"_"+sSuf+"R.v%d",iVer[0]);
  //fname += Form(".%d.root",iVer[1]);
  TString fo = foutname;

  if( !gROOT->IsBatch() && gSystem->FindFile("data/",fo) ) {
    std::cout << fo << " is existing. Do you recreate? (y/n)" << std::endl;
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

  anewFlow = new TClonesArray("STFlowInfo",1);

  //-- output                                                                                                              
  if( aBDC != NULL )          mflw->Branch("STBDC",&aBDC);
  if( aParticleArray != NULL) mflw->Branch("STParticle",&aParticleArray);
  if( anewFlow != NULL)       mflw->Branch("STFlow",&anewFlow);

}


