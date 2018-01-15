#include "flw_process4.h"
//----------------------------------------------------------------------
// flw_process4.C 
// input: 
// output:
// (c) Mizuki Kurata-Nishimura 
//----------------------------------------

//
//----------------------------------------------------------------------
void flw_process4(Long64_t nmax = -1)
{
  SetEnvironment();

  // Set number of bins divinding rapidity.
  if(sbRun != "0")
    nBin = SetDatabaseFiles();

  if(nBin == 0) {
    std::cout << " Correction database was not found. " << std::endl;
    exit(0);
  }

  // Print Run configuration 
  cout << " ---------- CONFIGURATION ----------  " << endl;
  cout << "RUN = "      << sRun << " Assemble v" << sAsm 
       << " with ver "  << sVer  << " mix = " << sMix 
       << endl;


  // I/O
  Open();

  OutputTree();

  //--------------------start
  //
  SetNumberOfProcess(nmax);

  for(Long64_t ievt = 0; ievt < maxProc; ievt++){

    PrintProcess(ievt);

    Initialize();
   

    fTree->GetEntry(ievt);
    Long64_t  nGoodTrack = aParticleArray->GetEntries();


    Int_t mtkBIN = -1;
    seltrack = ntrack[seltrackID];

    if(ntrack[2] > 0 ) {
     
      mtkBIN = GetMultiplicityCorretionIndex(seltrack);

      Long64_t nLoop = 0;

      TIter next(aParticleArray);
      STParticle *aPart1 = NULL;

      while( (aPart1 = (STParticle*)next()) ) {

	FlatteningCorrection(aPart1,mtkBIN);
	  
	// Particle selection
	if ( aPart1->GetReactionPlaneFlag() == selReactionPlanef ){  // p/d/t/He

	  unitP += aPart1->GetFlattenMomentum().Unit();
	
	  unitP_lang += aPart1->GetRPWeight() * (aPart1->GetFlattenPt()).Unit();
	
	  unitP_rot  += aPart1->GetRPWeight() * (aPart1->GetRotatedPt()).Unit();
	}
      }

      if(seltrack > 2) {
	SubEventAnalysis();
	AzmAngleRPTReactionPlane();
      }

      mflw->Fill();
    }
  }

  fout->cd();
  fout->Write();
  std::cout << fout->GetName() << std::endl;

  if(gROOT->IsBatch()) {
    fout->Close();
    exit(0);
  }

}

void SetEnvironment()
{

  sRun = gSystem -> Getenv("RUN");  // RUN number
  sVer = gSystem -> Getenv("VER");  // Version ID
  sbRun= gSystem -> Getenv("BRUN"); // If BRUN=0, flattening is not done.
  sbVer= gSystem -> Getenv("BVER"); // BRUN version
  //scVer= gSystem -> Getenv("CVER"); // corrected version

  sMix = gSystem -> Getenv("MIX");
  
  if(sRun =="" || sVer == "" ||sMix == "" || sbRun == "" || sbVer == ""|| !DefineVersion()) {
    cout << " Please type " << endl;
    cout << "$ RUN=#### VER=#.#.# MIX=0(real) or 1(mix) BRUN=# BVER=#.#.cv# root flw_process4.C(Number of event) " << endl;
    exit(0);
  }



  // Real or mixed event 
  if (sMix == "1") bMix = kTRUE;
  else bMix = kFALSE;


  // Set RUN number
  iRun = atoi(sRun);


}

void PrintProcess(Int_t ievt)
{
  TDatime dtime;
  static TDatime btime(dtime);

  UInt_t eprint = 1000;
  if(bMix) eprint = 1000;

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
  nEntry = fTree->GetEntries();

  maxProc = nEntry;
  std::cout << " Number of Events -->" << nEntry   << std::endl;


}

void Open()
{ 
  foutname = Form("run%d_rf.v%d.%d",iRun,iVer[0],iVer[1]);

  if( bMix )
    foutname = Form("run%d_mf.v%d.%d",iRun,iVer[0],iVer[1]);

  cout << foutname << endl;

  TString fn = foutname + ".root";;

  if( !gSystem->FindFile("data/", fn) ) {
    std::cout << fn << " : Input file desn't exist. "  << std::endl;
    exit(0);
  }

  
  auto fFile = new TFile(fn);
  if( bMix )
    fTree = (TTree*)fFile->Get("mflw");
  else
    fTree = (TTree*)fFile->Get("rflw");

  if(!fTree) {
    cout << " No data was found " << fn << endl;
    exit(0);
  }
  cout << "Input file is " << fn << endl;


  aParticleArray   = new TClonesArray("STParticle",150);

  fTree->SetBranchAddress("STParticle",&aParticleArray);
  fTree->SetBranchAddress("ntrack",ntrack);
  fTree->SetBranchAddress("mtrack",&mtrack);

  fTree->SetBranchAddress("aoq",&aoq);
  fTree->SetBranchAddress("z",&z);
  fTree->SetBranchAddress("snbm",&snbm);
}

void Initialize()
{

  trackID.clear();

  unitP   = TVector3(0,0,0);
  unitP_lang  = TVector2(0,0);
  unitP_rot   = TVector2(0,0);
  unitP_1 = TVector2(0.,0.);
  unitP_2 = TVector2(0.,0.);
  unitP_1r= TVector2(0.,0.);
  unitP_2r= TVector2(0.,0.);

  mtrack_1 = 0;
  mtrack_2 = 0;
	  
  numGoodTrack = 0;

  for(UInt_t i = 0; i < 2; i++){
    bsPhi[i]   = -999.;
    bsPhi_1[i] = -999.;
    bsPhi_2[i] = -999.;
  }
  for(UInt_t i = 0; i < 3; i++){
    bsPhi_ex[i]= -999.;
  }
  

  aParticleArray->Clear();

  TDatime dtime;
  rnd.SetSeed(dtime.Get());
}


Int_t GetPID(Double_t valx, Double_t valy)
{
  if(gProton){
    if(gProton->IsInside(valx,valy))
      return 12212;
  }
  return 0;
}

void LoadPIDFile()
{
  auto gcutFile = new TFile("gcutProton.root");
  gProton=(TCutG*)gcutFile->Get("gProton");
  gcutFile->Close();
}



void OutputTree()
{
  //@@@
  foutname += Form(".%d",iVer[2]);
  foutname += "_db"+sbRun+".v"+sbVer+".root";
  // Set like  -> run2843_rf.v0.0.0_db2843.v0.0.cv0.root

  TString fo = foutname;

  if( !gROOT->IsBatch() && gSystem->FindFile("data/",fo) ) {
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
    fout  = new TFile("data/"+foutname,"recreate");




  mflw = new TTree("cflw","FLOW analysis");

  cout << "Output file is " << foutname << endl;

  //-- output                                                                                                              
  mflw->Branch("irun",&iRun,"irun/I");
  mflw->Branch("aoq",&aoq,"aoq/D");
  mflw->Branch("z",&z,"z/D");
  mflw->Branch("snbm",&snbm,"snbm/I");

  mflw->Branch("STParticle",&aParticleArray);

  mflw->Branch("ntrack",ntrack,"ntrack[7]/I");
  mflw->Branch("numGoodTrack",&numGoodTrack);
  mflw->Branch("mtrack",&mtrack,"mtrack/I");
  mflw->Branch("mtrack_1",&mtrack_1,"mtrack_1/I");
  mflw->Branch("mtrack_2",&mtrack_2,"mtrack_1/I");

  mflw->Branch("unitP"     ,&unitP);
  mflw->Branch("unitP_1"   ,&unitP_1);
  mflw->Branch("unitP_2"   ,&unitP_2);
  mflw->Branch("unitP_lang",&unitP_lang);
  mflw->Branch("unitP_rot" ,&unitP_rot);
  mflw->Branch("unitP_1r"  ,&unitP_1r);
  mflw->Branch("unitP_2r"  ,&unitP_2r);

  mflw->Branch("bsPhi"     ,bsPhi   ,"bsPhi[2]/D");
  mflw->Branch("bsPhi_1"   ,bsPhi_1 ,"bsPhi_1[2]/D");
  mflw->Branch("bsPhi_2"   ,bsPhi_2 ,"bsPhi_2[2]/D");
  mflw->Branch("bsPhi_ex"  ,bsPhi_ex,"bsPhi_ex[3]/D");

}

Bool_t DefineVersion()
{
  Bool_t bfound = kFALSE;

  TString ver = sVer + ".";
  
  for ( Int_t i = 0; i < 3; i++) {
    if( ver.First(".") < 0 ) break;

    Ssiz_t end = ver.First(".")  ;
    TString ver1 = ver(0, end);

    ver = ver(end+1, ver.Length());

    iVer[i] = atoi(ver1);

    if(i==2) bfound = kTRUE;

  }
  
  if(!bfound)
    std::cout << " missing version number : v" << iVer[0] << "." << iVer[1] << "." << iVer[2] << std::endl;

  return bfound;

}

UInt_t SetDatabaseFiles()
{

  //  std::vector<TString> vfname;
  if(sbRun != "0") {
    TString fname = "run"+sbRun;
    if(bMix)
      fname += "_mf.v";
    else
      fname += "_rf.v";

    fname += sbVer + ".";

        cout << " fname " << fname  << endl;

    UInt_t ihmsum = 0;
    UInt_t imtk = 0;
    while(1){

      UInt_t ihm = 0;
      UInt_t ihmm = 0;
      while(1){

	TString ffname = fname +  Form("m%dn%d.txt",imtk,ihmm);
		cout << ffname << endl;
	
	if( gSystem->FindFile("db",ffname) ){
	  vfname.push_back(ffname);
	  ihm++;  ihmm++;
	}
	else 
	  break;
      }

      if(ihm > 0 ) {
	mtkbin.push_back(ihm);
	ihmsum += ihm;
	///	cout << " ihm " << ihm  << " <- " << imtk << " sum " << ihmsum << endl;
      }
      else if(ihm == 0)
	break;

      imtk++;
    }
  }

  
  nBin = (UInt_t)vfname.size();
  if(nBin == 0) return nBin;
  
  aflowcorrArray = new TClonesArray("STFlowCorrection",300);  
  TClonesArray &arr = *aflowcorrArray;

  cout << " nBin " << nBin << endl;

  std::vector<TString>::iterator itb;
  pbinmin.clear();
  UInt_t ihm = 0;
 
  for(itb = vfname.begin(); itb != vfname.end(); itb++, ihm++){

    new(arr[ihm]) STFlowCorrection();

    flowcorr = (STFlowCorrection*)arr.At(ihm);
    flowcorr->SetRealOrMix(1);
    flowcorr->SetRealOrMix((UInt_t)bMix);
    flowcorr->SetFileName(vfname.at(ihm));
    flowcorr->GetCorrectionFactor();
    
    pbinmin.push_back(make_pair(flowcorr->GetBin_min(0),flowcorr->GetBin_min(1)));
    cout << " $$$$$$$$$$$ ---->  nbin " << ihm  << " " << pbinmin.at(ihm).first << " : "<< pbinmin.at(ihm).second << endl;    
    //    flowcorr->ShowParameters();
  }

  binpara   = flowcorr->GetBinParameter(1);

  return nBin;    
}


void SetPtWeight(STParticle *apart)
{

  Double_t rpd   = apart->GetRapidity() - ycm;


  if( rpd < 0 )
    apart->SetRPWeight(-1);
  else
    apart->SetRPWeight(1);

}


Int_t GetMultiplicityCorretionIndex(UInt_t ival)
{
  std::vector< std::pair<Double_t, Double_t> >::iterator itr;

  UInt_t ink = mtkbin.size()-1;

  for(itr = pbinmin.end()-1; itr != pbinmin.begin(); itr-= mtkbin.at(ink), ink--){
    
    //    cout << "mult " << itr->first<< " : " << itr->second << " -- " << mtkbin.at(ink) << " at " << itr - pbinmin.begin()<< endl;

    if(ival >= itr->first) {
      return itr - pbinmin.begin();
    }
  }

  return -1;
}

Int_t GetThetaCorrectionIndex(Double_t fval, Int_t ival)
{

  std::vector< std::pair<Double_t, Double_t> >::iterator itr;
  // itr =  pbinmin.begin();
  // itr += ival;
  
  for(itr = pbinmin.begin()+ival; itr != pbinmin.begin()-1; itr--){
    
    if( fval >= itr->second) {      
      return itr - pbinmin.begin();
    }
  }

  return -1;
}

void CheckPlot(UInt_t ival)
{
  if(aflowcorrArray == NULL) return;

  auto cplot = new TCanvas();
  cplot->Divide(2,2);

  std::cout << "Checking histgram " << ival << std::endl;
  flowcorr = (STFlowCorrection*)aflowcorrArray->At(ival);

  flowcorr->PrintRange();

  hvphi  = new TH1D("hvphi"  ,"phi"   ,100,-3.2,3.2);
  hvthet = new TH1D("hvtheta","theta" ,100,0.,1.4);
  hvmtk  = new TH1I("hvmtk"  ,"mtrack", 60,0,60);

  std::vector<Double_t>::iterator itr;
  std::vector<Int_t>::iterator   iitr;

  std::vector<Double_t> vec1 =  flowcorr->GetOriginalPhi();
  for(itr=vec1.begin(); itr!=vec1.end(); itr++)
    hvphi->Fill(*itr);
  vec1.clear();

  cplot->cd(1);
  hvphi->Draw();
  

  vec1 =  flowcorr->GetTheta();
  for(itr=vec1.begin(); itr!=vec1.end(); itr++)
    hvthet->Fill(*itr);

  cplot->cd(2);
  hvthet->Draw();

  std::vector<Int_t> vec2 =  flowcorr->GetMTrack();
  for(iitr = vec2.begin(); iitr != vec2.end(); iitr++)
    hvmtk->Fill(*iitr);

  cplot->cd(3);
  hvmtk->Draw();

}

void SubEventAnalysis()
{
  //subevent analysis
  TIter next(aParticleArray);
  STParticle *aPart1 = NULL;

  //  UInt_t np = aParticleArray->GetEntries();
  UInt_t np = seltrack;
  Float_t arr[np];
  rnd.RndmArray(np, arr);
  UInt_t isel = 0;


  std::vector< TVector2 > elem1;
  std::vector< TVector2 > elem2;

  while( (aPart1 = (STParticle*)next()) ) {
    
    if(aPart1->GetReactionPlaneFlag() == selReactionPlanef){
      Double_t wt = aPart1->GetRPWeight();
      TVector2 pt = aPart1->GetCorrectedPt();
      TVector2 ptr= aPart1->GetRotatedPt();

      
      if( (UInt_t)(arr[isel]*np)%2 ==0 && mtrack_1 < seltrack/2 ) {
	unitP_1 += wt * pt.Unit();
	unitP_1r+= wt * ptr.Unit();

	elem1.push_back( wt * pt.Unit() );

	aPart1->AddReactionPlaneFlag(100);
	mtrack_1++;
      }
      else if( mtrack_2 < seltrack/2 ) {
	unitP_2 += wt * pt.Unit();
	unitP_2r+= wt * ptr.Unit();

	elem2.push_back( wt * pt.Unit() );

	aPart1->AddReactionPlaneFlag(200);
	mtrack_2++;
      }
      else{
	unitP_1 += wt * pt.Unit();
	unitP_1r+= wt * ptr.Unit();

	elem1.push_back( wt * pt.Unit() );

	aPart1->AddReactionPlaneFlag(100);
	mtrack_1++;
      }

      isel++;
      //cout << " _2 " << mtrack_2 << " / " << itra << endl;     
      if( mtrack_1 + mtrack_2 > seltrack ) break; 

    }
  }


  UInt_t nboot = 500;
  // auto btsp1 = new STBootStrap(nboot, &elem1); 
  // bsPhi_1[0] = btsp1->GetMean(2);
  // bsPhi_1[1] = btsp1->GetStdDev(2);

  // auto btsp2 = new STBootStrap(nboot, &elem2); 
  // bsPhi_2[0] = btsp2->GetMean(2);
  // bsPhi_2[1] = btsp2->GetStdDev(2);
  
}

void AzmAngleRPTReactionPlane()
{

  Int_t itra = 0;

  TIter next(aParticleArray);
  STParticle *aPart1 = NULL;


  while( (aPart1 = (STParticle*)next()) ) {

    if(aPart1->GetReactionPlaneFlag() > 100){
      Double_t wt = aPart1->GetRPWeight();
      TVector2 pt = aPart1->GetCorrectedPt();

      std::vector < TVector2 > exVec;
      TVector2 mExcRP(0.,0.);

      TIter rest(aParticleArray);
      STParticle *restPart = NULL;

      itra++;
      
      UInt_t itraex = 0;
      while( (restPart = (STParticle*)rest()) ) {

	if( aPart1 != restPart && restPart->GetReactionPlaneFlag() > 100 ) {

	  Double_t wt_rp = restPart->GetRPWeight();
	  TVector2 pt_rp = restPart->GetCorrectedPt();

	  mExcRP += wt_rp * pt_rp.Unit();
	  exVec.push_back( wt_rp * pt_rp.Unit() );


	  itraex++;
	}
      }

      // bootstrap method
      // auto btsp  = new STBootStrap(500, &exVec);
      
      // bsPhi_ex[0] = btsp->GetMean(2);
      // bsPhi_ex[1] = btsp->GetStdDev(2);
      // bsPhi_ex[2] = btsp->GetPhiOriginal();
      
      // aPart1->SetAzmAngle_wrt_RP  ( (Double_t)TVector2::Phi_mpi_pi( pt.Phi() - bsPhi_ex[0]));
      // aPart1->SetIndividualRPAngle( (Double_t)TVector2::Phi_mpi_pi( bsPhi_ex[0] ));

      aPart1->SetAzmAngle_wrt_RP( (Double_t)TVector2::Phi_mpi_pi( pt.Phi()-mExcRP.Phi()));
      aPart1->SetIndividualRPAngle( (Double_t)TVector2::Phi_mpi_pi( mExcRP.Phi() ));
      //      cout << " rest part " << itraex << " : " << itra << " @ "<< mExcRP.Phi() << " = " << aPart1->GetIndividualRPAngle()
      //   << endl;          
    }
  }
}


void FlatteningCorrection(STParticle *apart, Int_t ival)
{
  
  Int_t    mtkBIN = ival;
  UInt_t   iBin   = 999;
  Double_t binParameter;

  if( binpara == "pz") 
    binParameter = (apart->GetRotatedMomentum()).Z();
  else if( binpara == "theta" || binpara == "mtktheta")
    binParameter = (apart->GetRotatedMomentum()).Theta();
  else {
    std::cout << " This DB is not allowed " << binpara << std::endl;
    exit(0);
  }

  
  Int_t iBIN = GetThetaCorrectionIndex(binParameter, mtkBIN);

  //  cout << "binpara " << binParameter << " mtkBIN " << mtkBIN << " ibin " << iBIN <<endl;

  if(iBIN >= 0){
    flowcorr = (STFlowCorrection*)aflowcorrArray->At(iBIN);

    //    cout << flowcorr->GetBin_min(0) << " < " << flowcorr->GetBin_max(0) << " : " 
    //	 << flowcorr->GetBin_min(1) << " < " << flowcorr->GetBin_max(1) << endl; 

    //cout << " before " << apart->GetRotatedMomentum().Phi() << endl;
    Double_t phi =  flowcorr->GetCorrection((apart->GetRotatedMomentum()).Phi());
    apart->Flattening( phi );
    //    cout << " after " << apart->GetFlattenMomentum().Phi() << endl;
    

  }
  else
    std::cout << " A correction file is not found. " << binParameter << " with " << ntrack[3] << std::endl;

}
