#include "flw_process2.h"
//----------------------------------------------------------------------
//  Flow analysis proces 2
// (c) Mizuki Kurata-Nishimura
//    23 Nov 2017 
//----------------------------------------
//   Particle selection and mixed events will be done.
//
//----------------------------------------------------------------------
void flw_process2(Long64_t nmax = -1)
{
  //  nmax = 30;
  SetEnvironment();

  TDatime dtime;
  rnd.SetSeed(dtime.GetSecond());

  Int_t itime = dtime.Get();
  gRandom = new TRandom3(0);
  gRandom->SetSeed(itime);


  // I/O
  Open();

  OutputTree(nmax);

  SetNumberOfProcess(nmax);

  //  auto bs_unitP = new STBootStrap(100);

  for(Long64_t ievt = 0; ievt < maxProc; ievt++){

    PrintProcess(ievt);

    Initialize();

    if(ievt == 0)
      DefineLorentzBoostVector();

    //    bs_unitP->Clear();
    
    Long64_t nTrack = 0;

    //if 
    if(bMix) {
      nTrack = (Long64_t)GetRandomNumTrack(); // get number of track 

      // std::cout << " Number of tracks in mixed event is " << nTrack << std::endl;
    }
    else {
      fTree->GetEntry(ievt);

      if( (snbm != 132 && snbm != 108 && snbm != 124 && snbm != 112) ||
	  (ProjA < -1000 || ProjB < -1000) ) 
	continue;

      nTrack = aParticleArray->GetEntries();
    }

    TClonesArray &npar = *nParticleArray;

    if(ntrack[2] > 0 || bMix) {
     
      Long64_t nLoop = 0;

      while(kTRUE) {

	if(  bMix && numGoodTrack > nTrack - 1 ) break;  // mGoodTrack-1

	if( !bMix && nLoop >=  nTrack ) break;

	STParticle *aPart1;
	Long64_t    mixEvt = (Long64_t)ievt;
	
	Int_t ntrk = nTrack;

	if( bMix ) {
	  aPart1 = GetMixedTrack(&mixEvt, &ntrk); // Get a track randomly
	}
	else 
	  aPart1 = GetRealTrack(nLoop);

	// Good track and fragments are selected
	if(  CheckParticle(aPart1) ){
	  
	  SetFlowFlag(aPart1);

	  SetPtWeight(aPart1);
	    
	  if( bMix ) {
	    aPart1->SetMixedEventID(mixEvt);
	    event.push_back(mixEvt);
	  }
	  else 
	    event.push_back(ievt);
	  
	  aPart1->SetMixedNtrack(ntrk);

	  trackID.push_back(numGoodTrack);
	    
	  new(npar[numGoodTrack]) STParticle( *aPart1 );	  
	    
	  numGoodTrack++;

	  //	  cout << "flag " << aPart1->GetReactionPlaneFlag() << " >= " << selReactionPlanef ;;
	  if( aPart1->GetReactionPlaneFlag() >= selReactionPlanef ){
	    mtrack++;

	    unitP_ave  += aPart1->GetRotatedMomentum().Unit(); 
	    unitP_rot  += aPart1->GetRPWeight() * aPart1->GetRotatedMomentum().Unit();

	    unitP2_ave += aPart1->GetRotatedPt().Unit(); 
	    unitP2_rot += aPart1->GetRPWeight() * aPart1->GetRotatedPt().Unit();

	    //	    bs_unitP->Add(aPart1->GetRPWeight() * aPart1->GetRotatedPt().Unit());
	  }

	  if( aPart1->GetReactionPlaneFlag() == 20 )
	    ntrack[5]++;
	}
	nLoop++;
      }


      ntrack[4] = mtrack;
    }

    
    if(ntrack[4] > 0) {
      // bs_unitP->BootStrapping();
      // bsPhi[0] = bs_unitP->GetMean();
      // bsPhi[1] = bs_unitP->GetStdDev();
      // bsPhi[2] = bs_unitP->GetMod();

      SetSubEvent(npar, ntrack[4]);
    }
    else {
      unitP_ave.SetX(-999.);  unitP_ave.SetY(-999.);
      unitP_rot.SetX(-999.);  unitP_rot.SetY(-999.);
      unitP2_ave.SetX(-999.); unitP2_ave.SetY(-999.);
      unitP2_rot.SetX(-999.); unitP2_rot.SetY(-999.);
    }
    ntrack[6] = trackID.size();

    if(ntrack[3] > 0) {
      mflw->Fill();

      if(!bMix)
	hgtc->Fill(ntrack[3]);
    }
    else if(bMix)
      ntrack[3] = trackID.size();
  }

  if(!bMix && mhfile != NULL) {
    mhfile->cd();
    hgtc->Write();
  }


  fout->cd();
  fout->Write();
  std::cout << fout->GetName() << std::endl;

  if(gROOT->IsBatch()) {
    fout->Close();
    exit(0);
  }

  //  delete bs_unitP;
  delete gRandom;
}

void SetSubEvent(TClonesArray &pararray, const UInt_t npart)
{
  auto bs_unitP_1 = new STBootStrap(1000);
  auto bs_unitP_2 = new STBootStrap(1000);

  UInt_t *rndArray = new UInt_t[npart];
  
  rndArray = RanndomDivide2(npart);

  TIter next(&pararray);
  STParticle *aPart1 = NULL;

  mtrack_1 = 0;
  mtrack_2 = 0;

  UInt_t itrack = 0;
  while( (aPart1 = (STParticle*)next()) ){

    if(aPart1->GetReactionPlaneFlag() >= selReactionPlanef){
      Double_t wt = aPart1->GetRPWeight();
      TVector2 ptr= aPart1->GetRotatedPt();

      if( rndArray[itrack] == 0 ) {

	unitP_1r+= wt * ptr.Unit();
	bs_unitP_1->Add(wt * ptr.Unit());

	TVector2 ptpt = wt * ptr.Unit();
	//	std::cout << ptpt.Phi() << ", "; 

	aPart1->AddReactionPlaneFlag(100);
        mtrack_1++;
      }
      else  {
        unitP_2r+= wt * ptr.Unit();
	bs_unitP_2->Add(wt * ptr.Unit());

	aPart1->AddReactionPlaneFlag(200);
        mtrack_2++;
      }

      itrack++;
      if( itrack > npart ) break;

    }
  }

  //  std::cout << endl;

  // if( mtrack_1 > 0 && mtrack_2 > 0 ) {
  //   bs_unitP_1->BootStrapping();
  //   bs_unitP_2->BootStrapping();
    
  //   bsPhi_1[0] = bs_unitP_1->GetMean();
  //   bsPhi_1[1] = bs_unitP_1->GetStdDev();
  //   bsPhi_1[2] = bs_unitP_1->GetMod();
    
  //   bsPhi_2[0] = bs_unitP_2->GetMean();
  //   bsPhi_2[1] = bs_unitP_2->GetStdDev();
  //   bsPhi_2[2] = bs_unitP_2->GetMod();
  // }

  delete bs_unitP_1;
  delete bs_unitP_2;
}

UInt_t *RanndomDivide2(UInt_t npart)
{
  if(npart%2 == 1) npart++;

  const UInt_t npart2 = npart; 
  UInt_t  *rndarray = new UInt_t[npart2];

  UInt_t c1 = 0;
  UInt_t c2 = 0;
  UInt_t count = 0;
  while( count < npart2 ){
    
    Float_t rrd = rnd.Rndm() ;
    //    cout << rrd << "  ";

    if( rrd < 0.5 ) {
      if( c1 < npart2/2 ) {
	rndarray[count] = 0;
	c1++;
	count++;
      }
    }
    else if( rrd >= 0.5 ) {
      if( c2 < npart2/2 ) {
	rndarray[count] = 1;
	c2++;
	count++;
      }
    }
  }


  return rndarray;
}

Bool_t CheckParticle(STParticle *apart)
{
  Bool_t bsel = kFALSE;


  if( apart == NULL ) return bsel;


  if( !apart->GetBestTrackFlag() ) return bsel;

  //  ResetPID(apart); // graphical cut PID


  auto pid    =  apart -> GetPID();
  auto charge =  apart -> GetCharge();
  
  if( pid == 211 ) 
    bsel = kTRUE;
  else if( pid > 2000 && charge > 0 ) 
    bsel = kTRUE;
  
  return bsel;
} 

void SetFlowFlag(STParticle *apart)
{
  auto pid    =  apart -> GetPID();

  if( pid == 211 )
    apart->SetReactionPlaneFlag(1);
  
  else if( pid > 2000 && 
	   apart->GetBestTrackFlag()     > 0 &&
	   apart->GetMaxdEdxFlag()       > 0 &&
	   apart->GetMaxMomentumFlag()   > 0
	   ) {
    apart->SetReactionPlaneFlag(1000);

    if(apart->GetNDFFlag())
      apart->SetReactionPlaneFlag(2000);
  }
  else
    apart->SetReactionPlaneFlag(0);
}


void SetEnvironment()
{

  sRun = gSystem -> Getenv("RUN");  // RUN number
  sVer = gSystem -> Getenv("VER");  // Version ID
  sMix = gSystem -> Getenv("MIX");


  if(sRun =="" || sVer == "" ||sMix == "" ||!DefineVersion()) {
    cout << " Please type " << endl;
    cout << "$ RUN=#### VER=#.# MIX=0(real) or 1(mix) root flw_process2.C(Number of event) " 
	 << endl;
    exit(0);
  }


  // Print Run configuration 
  cout << " ---------- CONFIGURATION ----------  " << endl;
  cout << "RUN = "      << sRun << " Assemble v" << sAsm 
       << " with ver "  << sVer  << " mix = " << sMix 
       << " Flatten ; " << nBin
       << endl;

  cout << "sBinp = " << sBinp << " Flatten RUN = " << sbRun << endl;

  // Real or mixed event 
  if (sMix == "1") bMix = kTRUE;
  else bMix = kFALSE;


  // Set RUN number
  iRun = atoi(sRun);

  if( sbRun == "0" )
    nBin = 0;

}

void PrintProcess(Int_t ievt)
{
  TDatime dtime;
  static TDatime btime(dtime);

  UInt_t eprint = 10000;
  if(bMix) eprint = 100;

  if(ievt%eprint == 0) {
    dtime.Set();
    Int_t ptime = dtime.Get() - btime.Get();
    std::cout << "Process " << setw(8) << ievt << "/"<< maxProc << " = " 
	      << ((Double_t)(ievt)/(Double_t)maxProc)*100. << " % --->"
	      << dtime.AsString() << " ---- "
	      << setw(5) << (Int_t)ptime/60 << " [min] "
	      << std::endl;
  }
}

void SetNumberOfProcess(Int_t nmax)
{
  Long64_t mEvt = GetMultiplicityDistribution();
  nEntry = fTree->GetEntries();

  if( bMix ) {
    if( nmax == -1 )  maxProc = mEvt;
    else if( nmax == -2 ) maxProc = mEvt*10;
    else
      maxProc = nmax;
  }
  else{
    if( nmax == -1 || nmax > nEntry )  maxProc = nEntry;
    else   maxProc = nmax;
  }
      


  std::cout << " Number of Events -->" << nEntry   << std::endl;
  std::cout << " Maximum Events   -->" << maxProc  << std::endl;
  std::cout << " Set     Events   -->" << nmax  << std::endl;

}

void Open()
{
  //  LoadPIDFile();

  TString fnameBase = "run%d_f0.v";

  Ssiz_t end = sVer.First(".");

  TString fn = Form(fnameBase + sVer(0,end)+".root",iRun);
  if( !gSystem->FindFile("data/", fn) ) {
    

    Int_t iver = atoi((TString)sVer(end+1,1));

    while( iver > -1 ){
      TString ss = Form(".%d",iver);
      fn = Form(fnameBase + sVer(0,end)+ss+".root",iRun);

      if( !gSystem->FindFile("data/", fn) ) 
	iver--;
      else 
	break;
    }

  }
  
  auto fFile = new TFile(fn);
  fTree = (TTree*)fFile->Get("flw");

  if(!fTree) {
    cout << " No data was found " << fn << endl;
    exit(0);
  }
  cout << "Input file is " << fn << endl;


  aParticleArray   = new TClonesArray("STParticle",150);

  fTree->SetBranchAddress("STParticle",&aParticleArray);
  fTree->SetBranchAddress("ntrack",ntrack);
  //  fTree->SetBranchAddress("kymult",&kymult);

  fTree->SetBranchAddress("aoq",&aoq);
  fTree->SetBranchAddress("z",&z);
  fTree->SetBranchAddress("ProjA",&ProjA);
  fTree->SetBranchAddress("ProjB",&ProjB);
  fTree->SetBranchAddress("snbm",&snbm);

  fTree->SetBranchAddress("STNeuLANDCluster",&aNLCluster);
  //  p_rot  = new TClonesArray("TVector3",150);
  //  p_org  = new TClonesArray("TVector3",150);
  
}

void Initialize()
{

  trackID.clear();

  mtrack = 0;
  mxntrk = 0;  
  numGoodTrack = 0;

  //  for(Int_t i = 0; i< 7; i++) ntrack[i] = 0;
  
  event.clear();

  aParticleArray->Clear();
  nParticleArray->Clear();

  unitP_ave  = TVector3(0.,0.,0.);
  unitP_rot  = TVector3(0.,0.,0.);
  unitP2_ave = TVector2(0.,0.);
  unitP2_rot = TVector2(0.,0.);
  unitP_1r   = TVector2(0.,0.);
  unitP_2r   = TVector2(0.,0.);

  for(UInt_t i = 0; i < 3; i++){
    bsPhi[i]   = -999.;
    bsPhi_1[i] = -999.;
    bsPhi_2[i] = -999.;
  }

}


Long64_t GetMultiplicityDistribution()
{
  // Multplicity histgram
  TString mHistFile = Form("MultRoot/run%d.ngt_v"+sVer(0,1)+".root",iRun);
  TString hf = mHistFile;

  Long64_t nEvt = 0;

  if(bMix) { 
    mhfile = new TFile(mHistFile);
    if( mhfile == NULL ) {
      cout << mHistFile << " is not found." << endl;
      exit(0);
    }
    cout << "Multiplicity : "  << mHistFile << " is loaded." << endl;


    if( (TH1I*)mhfile->FindObject("hgtc") ) 
      exit(0);
    
    histGT_r =  (TH1I*)mhfile->Get("hgtc");
    histGT_r -> SetName("hgtc_r");
    histGT_r -> SetDirectory(gROOT);
    nEvt = (Long64_t)histGT_r -> GetEntries();
  }
  else {
    //    if( !gSystem->FindFile(".",hf) ) {
    mhfile = new TFile(mHistFile,"recreate");
    cout << mHistFile << " is created." << endl;
    //    }
    // else
    //   cout << mHistFile << " is existing." << endl;
  }

  return nEvt;

}


void ResetPID(STParticle *apart)
{

  auto p    = apart->GetP();
  auto dedx = apart->GetdEdx();

  if(gProton && gDeuteron && gTriton && gPip && gPim){
    Int_t gpid  = 0;

    if(gProton->IsInside(p, dedx))
      gpid = 2212;
    else if(gDeuteron->IsInside(p, dedx))
      gpid = 1000010020;
    else if(gTriton->IsInside(p, dedx))
      gpid = 1000010030;
    else if(gPip->IsInside(p, dedx))
      gpid = 211;
    else if(gPim->IsInside(p, dedx))
      gpid = 211;

    else
      gpid = 12212; // temporal

    apart->SetPID(gpid);
  }
}


void LoadPIDFile()
{
  TString fname = "gcutPID132Sn.ROOT";
  TFile *gcutFile;

  if( !gSystem->FindFile("data",fname) )
    std::cout << " Graphical PID selection is not determined " << std::endl;
  else {
    std::cout << " Graphical PID selection, " << fname << ", is loaded " << std::endl;

    gcutFile = new TFile(fname);

    gProton   = (TCutG*)gcutFile->Get("gcutProton132Sn2");
    gDeuteron = (TCutG*)gcutFile->Get("gcutDeutron132Sn");
    gTriton   = (TCutG*)gcutFile->Get("gcutTriton132Sn");
    gPip      = (TCutG*)gcutFile->Get("gcutPip132Sn");
    gPim      = (TCutG*)gcutFile->Get("gcutPim132Sn");
    

    gcutFile->Close();
  }
}



void OutputTree(Long64_t nmax)
{

  // ROOT out
  TString sdeb = ".s";
  if(nmax < 0)  sdeb = "";
  TString foutname = "data/run"+sRun+"_";

  if( bMix ) {
    foutname += "mf";
    foutname += ".v"+sVer+".root";
  }
  else {
    foutname += "rf";
    foutname += ".v"+sVer+sdeb+".root";
  }
  
  TString fo = foutname;

  if( !gROOT->IsBatch() && gSystem->FindFile(".",fo) ) {
    cout << foutname << " is existing. Do you recreate? (y/n)" << endl;
    TString sAns;
    cin >> sAns;
    if(sAns == "y" || sAns == "Y")
      fout  = new TFile(foutname,"recreate");
    else {
      cout << " Retry" << endl;
      exit(0);
    }      
  }
  else
    fout  = new TFile(foutname,"recreate");



  if( bMix ) 
    mflw = new TTree("mflw","FLOW analysis track mixing");    
  else
    mflw = new TTree("rflw","FLOW analysis");

  cout << "Output file is " << foutname << endl;

  //-- output                                                                                                              
  mflw->Branch("irun",&iRun,"irun/I");
  mflw->Branch("aoq",&aoq,"aoq/D");
  mflw->Branch("z",&z,"z/D");
  mflw->Branch("snbm", &snbm, "snbm/I");
  mflw->Branch("ProjA",&ProjA,"ProjA/D");
  mflw->Branch("ProjB",&ProjB,"ProjB/D");

  nParticleArray = new TClonesArray("STParticle",100);
  mflw->Branch("STParticle",&nParticleArray);

  mflw->Branch("ntrack"    ,ntrack,"ntrack[7]/I");
  mflw->Branch("event"     ,&event);
  mflw->Branch("mxntrk"    ,&mxntrk);
  mflw->Branch("unitP_ave" ,&unitP_ave);
  mflw->Branch("unitP_rot" ,&unitP_rot);
  mflw->Branch("unitP2_ave",&unitP2_ave);
  mflw->Branch("unitP2_rot",&unitP2_rot);

  mflw->Branch("unitP_1r",&unitP_1r);
  mflw->Branch("unitP_2r",&unitP_2r);
  mflw->Branch("bsPhi"     ,bsPhi   ,"bsPhi[3]/D");
  mflw->Branch("bsPhi_1"   ,bsPhi_1 ,"bsPhi_1[3]/D");
  mflw->Branch("bsPhi_2"   ,bsPhi_2 ,"bsPhi_2[3]/D");

  mflw->Branch("mtrack_1",&mtrack_1);
  mflw->Branch("mtrack_2",&mtrack_2);
  


  if(aNLCluster != NULL)
    mflw->Branch("STNeuLANDCluster",&aNLCluster);
      

  hgtc         = new TH1I("hgtc","number of good fragment",100,0,100);
}

Long64_t GetRandomNumTrack()
{
  // histGT_r is filled with ntrack[3]

  if(!histGT_r){
    cout << " no histgram was found for randowm number " << endl;
    return 0;
  }

  Long64_t bcont = (Long64_t)histGT_r->GetRandom();


  return bcont;
}

STParticle *GetRealTrack(Long64_t ival)
{
  STParticle *realPart = NULL;
  realPart = (STParticle*)aParticleArray -> At( ival );
  
  //  cout << " real part " << realPart << endl;
  return realPart;
}


STParticle *GetMixedTrack(Long64_t *ival, Int_t *kval)
{
  STParticle  *mixPart = NULL;
  static Long64_t mevt = *ival;

  // ---------- shfiting -- search similar multiplcity events
  while(kTRUE){	
    mevt += 2;
    if(mevt > nEntry) 	mevt = 0;
           
    UInt_t kLoop = 0;
    while( kLoop < 5 ){
      fTree->GetEntry(mevt);

      mevt++;
      if(mevt > nEntry) {
	mevt = 0;
	kLoop++;
	if(kLoop > 1) cout << " Too much loops " << kLoop << " kval " << *kval << " mevt " << mevt << " / " << nEntry <<  endl;
      }

      if( snbm != 132 && snbm != 108 && snbm != 124 && snbm != 112) continue;

      if( std::abs(ntrack[3] - *kval) < ntr_diff ) {
	//	cout << "ntrack[3] " << ntrack[3] << " kval " << *kval  << " - " << mevt << endl;
	break;
      }
    }


    Int_t nmixTrack = aParticleArray->GetEntries();

    Int_t imixTrack = (Int_t)pran.Uniform(0,nmixTrack);

    if( imixTrack < nmixTrack ) {
      mixPart = (STParticle*)aParticleArray->At(imixTrack);	        

      if( mixPart->GetBestTrackFlag()) {
	mixPart->SetMixedNtrack(nmixTrack);
	mixPart->SetMixTrackID(imixTrack);

	*kval = ntrack[3];

	*ival = mevt;
	return mixPart;
      }

    }

  }

  *ival = -1;
  return mixPart;
}


void  DefineLorentzBoostVector()
{
  Double_t amu  = 931.4940954; //MeV/c2                                                                                                     
  Double_t c    = 299792458.; //m/s                                                                                                         

  TString system[5];
  Double_t AB[5];
  Double_t mB[5];
  Double_t eB_lb[5];
  Double_t mT[5];

  system[0] =  "(132Sn + 124Sn)";
  AB[0]     =  132.;
  mB[0]     =  131.8906 ; //amu                                                                                                             
  eB_lb[0]  =  268.9 ;    //MeV/amu; incident energy                                                                                        
  mT[0]     =  123.8773895;   //amu Target mass                                                                                             

  system[1] = "(108Sn + 112Sn)";
  AB[1]     =  108.;
  mB[1]     =  107.8844964; //amu                                                                                                           
  eB_lb[1]  =  268.9;
  mT[1]     =  111.8773895;;

  system[2] = "(124Sn + 112Sn)";
  AB[2]     =  124.;
  mB[2]     =  123.8778449; //amu                                                                                                           
  eB_lb[2]  =  270.2;
  mT[2]     =  111.8773895;;

  system[3] = "(112Sn + 124Sn)";
  AB[3]     =  112.;
  mB[3]     =  111.8773895; //amu                                                                                                           
  eB_lb[3]  =  270.2;
  mT[3]     =  123.8773895;

  system[4] = "(p + p)";
  AB[4]     = 1.;
  mB[4]     = 1.00727646688;
  eB_lb[4]  = 268.9;
  mT[4]     = 1.00727646688;

  UInt_t sysid = 4;
  switch(snbm){
  case 132:
    sysid = 0;
    break;
  case 108:
    sysid = 1;
    break;
  case 124:
    sysid = 2;
    break;
  case 112:
    sysid = 3;
    break;
  }
  

  sysid = 4;

  std::cout << " Lorentz Boost with " << system[sysid] << std::endl;

  Double_t EkB_lb    = eB_lb[sysid]  * mB[sysid];
  mB[sysid] *= amu;
  mT[sysid] *= amu;

  // Beam                                                                                                                                   
  Double_t EB_lb  = EkB_lb + mB[sysid];
  Double_t PB_lb  = sqrt(EB_lb*EB_lb - mB[sysid]*mB[sysid]);

  auto bmVec = new TLorentzVector( TVector3(0., 0., PB_lb), EB_lb );
  auto tgVec = new TLorentzVector( TVector3(0., 0., 0.), mT[sysid] );
  auto totalVec = *bmVec + *tgVec;
  
  boostVec = totalVec.BoostVector();

}


void SetPtWeight(STParticle *apart)
{
  apart->SetRapidity();

  Double_t Etot = sqrt(apart->GetMomentum().Mag2() + pow(apart->GetMass(),2) );
  TLorentzVector lrnzVec( apart->GetMomentum(), Etot);
  
  lrnzVec.Boost(-boostVec);

  auto PZcm = lrnzVec.Z();
  auto Ecm  = lrnzVec.E();
  auto rapidity = 0.5*log( (Ecm + PZcm)/(Ecm - PZcm) );

  rapidity = 0.5*log( (Ecm + PZcm)/(Ecm - PZcm) );

  
  apart->SetRapiditycm(rapidity);
  //  cout << " rap cm " << apart->GetRapiditycm() << endl;
  

  if( rapidity  <  0 ) 
    apart->SetRPWeight(-1);
  else
    apart->SetRPWeight(1);

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
    cout << " missing version number : " << iVer[0] << "." << iVer[1]  << endl;

  return bfound;

}

