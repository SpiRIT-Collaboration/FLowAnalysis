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
  auto fflowtask = new STFlowTask(kTRUE, kTRUE, kFALSE);
  auto fIsFlowAnalyais = fflowtask->Init(iRun, dVer);

  Open();

  OutputTree();

  nEntry = rChain->GetEntries();
  for(Long64_t ievt = 0; ievt < nEntry; ievt++){

    PrintProcess(ievt);

    rChain->GetEntry(ievt);

    auto abeamInfo = (STBDC*)aBDC->At(0);
    if( !abeamInfo->GetBeamPID() ) continue;

    auto aFlowInfo = (STFlowInfo*)aFlowArray->At(0);
    

    fflowtask->SetFlowInfo( aFlowInfo );
    fflowtask->DoFlattening();
    fflowtask->DoFlatteningSub();

    aFlowInfo = fflowtask->GetFlowInfo();
      
    TClonesArray &aflow = *anewFlow;
    new( aflow[0] ) STFlowInfo( *aFlowInfo );


    fflowtask->SetParticleArray( *aParticleArray );
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
  
  if(sRun =="" || sSuf == "" || sVer == ""|| !DefineVersion()) {
    cout << " Please type " << endl;
    cout << "$ RUN=#### VER=#.# SUFX=BTt root DoFlow_Analysis.C " << endl;
    exit(0);
  }

  // Set RUN number
  iRun = atoi(sRun);

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
    if( ver.First(".") < 0 ) break;

    Ssiz_t end = ver.First(".")  ;
    TString ver1 = ver(0, end);

    ver = ver(end+1, ver.Length());

    iVer[i] = atoi(ver1);

    if(i==1) bfound = kTRUE;

  }
  
  if(!bfound)
    std::cout << " missing version number : v" << iVer[0] << "." << iVer[1]  << std::endl;

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
  fname = Form("run"+sRun+"_"+sSuf+".v%d",iVer[0]);
  TString fn = fname + ".root";
  std::cout << fn << std::endl;

  if( !gSystem->FindFile("data/", fn) ) {
    std::cout << fname << " : Input file desn't exist. "  << std::endl;
    exit(0);
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
  fname += Form(".%d.root",iVer[1]);
  TString fo = fname;

  if( !gROOT->IsBatch() && gSystem->FindFile("data/",fo) ) {
    cout << fname << " is existing. Do you recreate? (y/n)" << endl;
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
    std::cout << "Output file is " << fname << endl;
    fname = "data/" + fname;
    fout  = new TFile(fname,"recreate");
  }

  mflw = new TTree("cbmsim","Flow corrected");

  anewFlow = new TClonesArray("STFlowInfo",1);

  //-- output                                                                                                              
  if( aBDC != NULL )          mflw->Branch("STBDC",&aBDC);
  if( aParticleArray != NULL) mflw->Branch("STParticle",&aParticleArray);
  if( anewFlow != NULL)       mflw->Branch("STFlow",&anewFlow);
}


