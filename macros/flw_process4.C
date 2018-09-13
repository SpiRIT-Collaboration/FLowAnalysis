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

  nBin = SetDatabaseFiles();

  if(nBin == 0) {
    std::cout << " Correction database was not found. " << std::endl;
    exit(0);
  }

  // Print Run configuration 
  cout << " ---------- CONFIGURATION ----------  " << endl;
  cout << "RUN = "      << sRun   << " with ver "  << sVer  
       << "Flow correction ver "  << sbVer  << nBin << " databases " 
       << " mix = " << sMix 
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

      unitP_fc = Psi_FlatteningCorrection( 0, seltrack, TVector3(unitP2_rot->X(), unitP2_rot->Y(), 0.));
      unitP_rc = Psi_ReCenteringCorrection(0, seltrack, TVector3(unitP2_rot->X(), unitP2_rot->Y(), 0.));

      SubEventAnalysis();
      AzmAngleWRTReactionPlane();

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
  sbVer= gSystem -> Getenv("FLC"); // BRUN version
  sMix = gSystem -> Getenv("MIX");
  ssVer= gSystem -> Getenv("FLCS"); // BRUN version
  
  if(sRun =="" || sVer == "" ||sMix == "" || sbVer == ""|| !DefineVersion()) {
    cout << " Please type " << endl;
    cout << "$ RUN=#### VER=#.#.# MIX=0(real) or 1(mix) FLC=run2900_rf.v4.0.Psicv0 root flw_process4.C(Number of event) " << endl;
    exit(0);
  }

  // Real or mixed event 
  if (sMix == "1") bMix = kTRUE;
  else bMix = kFALSE;

  // Set RUN number
  iRun = atoi(sRun);

  // Correction version
  Ssiz_t cv = sbVer.First("c");
  ssbVer = sbVer( cv, 3);
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
  fTree->SetBranchAddress("ntrack"    ,ntrack);
  fTree->SetBranchAddress("mtrack"    ,&mtrack);
  fTree->SetBranchAddress("unitP_ave" ,&unitP_ave,  &bunitP_ave);
  fTree->SetBranchAddress("unitP_rot" ,&unitP_rot,  &bunitP_rot);
  fTree->SetBranchAddress("unitP2_ave",&unitP2_ave, &bunitP2_ave);
  fTree->SetBranchAddress("unitP2_rot",&unitP2_rot, &bunitP2_rot);
  fTree->SetBranchAddress("mtrack_1"  ,&mtrack_1);
  fTree->SetBranchAddress("mtrack_2"  ,&mtrack_2);
  fTree->SetBranchAddress("unitP_1r"  ,&unitP_1r  , &bunitP_1r);
  fTree->SetBranchAddress("unitP_2r"  ,&unitP_2r  , &bunitP_2r);

  fTree->SetBranchAddress("aoq" ,&aoq);
  fTree->SetBranchAddress("z"   ,&z);
  fTree->SetBranchAddress("snbm",&snbm);
  fTree->SetBranchAddress("ProjA",&ProjA);
  fTree->SetBranchAddress("ProjB",&ProjB);

  fTree->SetBranchAddress("STNeuLANDCluster",&aNLCluster);
}

void Initialize()
{
  trackID.clear();

  unitP       = TVector3(0.,0.,0.);

  unitP_lang  = TVector2(0.,0.);
  unitP_1     = TVector2(0.,0.);
  unitP_2     = TVector2(0.,0.);
  unitP_fc    = TVector3(0.,0.,0.);
  unitP_rc    = TVector3(0.,0.,0.);

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
  foutname += "." + ssbVer + ".root";
  //  foutname += "_db"+sbVer+".root";
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
  mflw->Branch("ProjA",&ProjA,"ProjA/D");
  mflw->Branch("ProjB",&ProjA,"ProjB/D");

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
  mflw->Branch("unitP2_rot",&unitP2_rot);
  mflw->Branch("unitP_1r"  ,&unitP_1r);
  mflw->Branch("unitP_2r"  ,&unitP_2r);
  mflw->Branch("unitP_fc"  ,&unitP_fc);
  mflw->Branch("unitP_rc"  ,&unitP_rc);

  mflw->Branch("bsPhi"     ,bsPhi   ,"bsPhi[2]/D");
  mflw->Branch("bsPhi_1"   ,bsPhi_1 ,"bsPhi_1[2]/D");
  mflw->Branch("bsPhi_2"   ,bsPhi_2 ,"bsPhi_2[2]/D");
  mflw->Branch("bsPhi_ex"  ,bsPhi_ex,"bsPhi_ex[3]/D");

  if(aNLCluster != NULL)
    mflw->Branch("STNeuLANDCluster",&aNLCluster);

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

  TString  fname[2];
  fname[0] = sbVer + ".";
  
  UInt_t ncount = 1;
  if(ssVer != "") {
    fname[1] = ssVer + ".";
    ncount = 2;
  }


  for(UInt_t i = 0; i < ncount; i++){
    cout << " Database name is " << fname[i]  << endl;

    UInt_t ihmsum = 0;
    UInt_t imtk = 0;
    while(1){
      
      UInt_t ihm = 0;
      UInt_t ihmm = 0;
      while(1){
	
	TString ffname = fname[i] +  Form("m%dn%d.txt",imtk,ihmm);
	if( gSystem->FindFile("db",ffname) ){
	  vfname[i].push_back(ffname);
	  ihm++;  ihmm++;
	  std::cout << " Databse " << ffname << " is loaded. " << std::endl;
	}
	else 
	  break;
      }
      
      if(ihm > 0 ) {
	mtkbin[i].push_back(ihm);
	ihmsum += ihm;
      }
      else if(ihm == 0)
	break;
      
      imtk++;
    }
    
  
    nBin = (UInt_t)vfname[i].size();
    if(nBin == 0) return nBin;

    aflowcorrArray[i] = new TClonesArray("STFlowCorrection",30);    
    TClonesArray &arr = *aflowcorrArray[i];

    cout << " nBin " << nBin << endl;

    std::vector<TString>::iterator itb;
    pbinmin[i].clear();
    UInt_t ihm = 0;
    
    for(itb = vfname[i].begin(); itb != vfname[i].end(); itb++, ihm++){
      
      new(arr[ihm]) STFlowCorrection();

      flowcorr = (STFlowCorrection*)arr.At(ihm);
      flowcorr->SetRealOrMix(1);
      flowcorr->SetRealOrMix((UInt_t)bMix);
      flowcorr->SetFileName(vfname[i].at(ihm));
      flowcorr->LoadCorrectionFactor();
      
      pbinmin[i].push_back(make_pair(flowcorr->GetBin_min(0),flowcorr->GetBin_min(1)));
      cout << " $$$$$$$$$$$ ---->  nbin " << ihm  << " " << pbinmin[i].at(ihm).first << " : "<< pbinmin[i].at(ihm).second << endl;    
      //    flowcorr->ShowBinInformation();
    }

    binpara   = flowcorr->GetBinParameter(1);
  }

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

Int_t GetCorrectionIndex(UInt_t isel, UInt_t ival, Double_t fval)
{
  Int_t index = GetMultiplicityCorretionIndex(isel, ival);

  index =  GetThetaCorrectionIndex(isel, index, fval);

  return index;
}

Int_t GetThetaCorrectionIndex(UInt_t isel, Int_t ival, Double_t fval)
{

  std::vector< std::pair<Double_t, Double_t> >::iterator itr;
  
  for(itr = pbinmin[isel].begin()+ival; itr != pbinmin[isel].begin()-1; itr--){
    
    if( fval >= itr->second) {      
      return itr - pbinmin[isel].begin();
    }
  }

  return -1;
}
Int_t GetMultiplicityCorretionIndex(UInt_t isel, UInt_t ival)
{
  std::vector< std::pair<Double_t, Double_t> >::iterator itr;


  UInt_t ink = mtkbin[isel].size()-1;

  for(itr = pbinmin[isel].end()-1; itr != pbinmin[isel].begin()-1; itr-= mtkbin[isel].at(ink), ink--){
    
    if(isel == 9)
      cout << "mult " << itr->first<< " : " << itr->second << " -- " << mtkbin[isel].at(ink) << " at " << itr - pbinmin[isel].begin()<< endl;

    if(ival >= itr->first) {

      if(isel == 9)
      cout << " ntrack  " << ival << " : " << itr - pbinmin[isel].begin() << endl;

      return itr - pbinmin[isel].begin();
    }
  }

  return -1;
}

void CheckPlot(UInt_t ival)
{
  if(aflowcorrArray[0] == NULL) return;

  auto cplot = new TCanvas();
  cplot->Divide(2,2);

  std::cout << "Checking histgram " << ival << std::endl;
  flowcorr = (STFlowCorrection*)aflowcorrArray[0]->At(ival);

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
  if(mtrack_1 > 0 && mtrack_2 > 0){
    
    TVector3 uvec;
    uvec = Psi_FlatteningCorrection(1, mtrack_1, TVector3(unitP_1r->X(), unitP_1r->Y(), 0.));
    unitP_1.SetX(uvec.X()); unitP_1.SetY(uvec.Y());

    uvec = Psi_FlatteningCorrection(1, mtrack_2, TVector3(unitP_2r->X(), unitP_2r->Y(), 0.));
    unitP_2.SetX(uvec.X()); unitP_2.SetY(uvec.Y());

  }

  // UInt_t nboot = 50;
  // auto btsp1 = new STBootStrap(nboot, (UInt_t)unitP_1.Size(), unitP_1); 
  // bsPhi_1[0] = btsp1->GetMean(2);
  // bsPhi_1[1] = btsp1->GetStdDev(2);

  // auto btsp2 = new STBootStrap(nboot, (UInt_t)unitP_2.Size(), unitP_2); 
  // bsPhi_2[0] = btsp2->GetMean(2);
  // bsPhi_2[1] = btsp2->GetStdDev(2);
  
}

void AzmAngleWRTReactionPlane()
{

  Int_t itra = 0;

  TIter next(aParticleArray);
  STParticle *aPart1 = NULL;

  while( (aPart1 = (STParticle*)next()) ) {

    //    if( aPart1->GetReactionPlaneFlag() > 1000 || aPart1->GetPID() == 211){

      Double_t wt = aPart1->GetRPWeight();
      TVector2 pt = aPart1->GetRotatedPt();

      std::vector < TVector2 > exVec;
      TVector2 mExcRP(0.,0.);

      TIter rest(aParticleArray);
      STParticle *restPart = NULL;

      itra++;
      
      UInt_t itraex = 0;
      while( (restPart = (STParticle*)rest()) ) {

	if( aPart1 != restPart && restPart->GetReactionPlaneFlag() > 1000 ) {

	  Double_t wt_rp = restPart->GetRPWeight();
	  TVector2 pt_rp = restPart->GetRotatedPt();

	  mExcRP += wt_rp * pt_rp.Unit();
	  exVec.push_back( wt_rp * pt_rp.Unit() );

	  
	  itraex++;
	}
      }

      // ReCentering and Shifting correction
      TVector3 rp_recv = TVector3(-999.,-999.,-999.);
      if(itraex > 0)
	rp_recv = Psi_FlatteningCorrection(0, seltrack , TVector3(mExcRP.X(), mExcRP.Y(), 0.));


      // bootstrap method
      // auto btsp  = new STBootStrap(500, &exVec);
      
      // bsPhi_ex[0] = btsp->GetMean(2);
      // bsPhi_ex[1] = btsp->GetStdDev(2);
      // bsPhi_ex[2] = btsp->GetPhiOriginal();
      
      // aPart1->SetAzmAngle_wrt_RP  ( (Double_t)TVector2::Phi_mpi_pi( pt.Phi() - bsPhi_ex[0]));
      // aPart1->SetIndividualRPAngle( (Double_t)TVector2::Phi_mpi_pi( bsPhi_ex[0] ));

      // aPart1->SetAzmAngle_wrt_RP( (Double_t)TVector2::Phi_mpi_pi( pt.Phi()-mExcRP.Phi()));
      aPart1->SetIndividualRPAngle( (Double_t)TVector2::Phi_mpi_pi( rp_recv.Phi() ));

      aPart1->SetAzmAngle_wrt_RP( (Double_t)TVector2::Phi_mpi_pi( pt.Phi() - rp_recv.Phi()));
      //aPart1->SetAzmAngle_wrt_RP( (Double_t)TVector2::Phi_mpi_pi( unitP2_rot->Phi() ));
      //      aPart1->SetIndividualRPAngle( (Double_t)TVector2::Phi_mpi_pi( rp_rec.Phi() ));

      //}
      // else
      // aPart1->SetAzmAngle_wrt_RP(-10.);
  }
}


void FlatteningCorrection(UInt_t isel, STParticle *apart, Int_t ival)
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

  
  Int_t iBIN = GetThetaCorrectionIndex(isel, mtkBIN, binParameter);

  if(iBIN >= 0){
    flowcorr = (STFlowCorrection*)aflowcorrArray[0]->At(iBIN);

    Double_t phi =  flowcorr->GetCorrection((apart->GetRotatedMomentum()).Phi());
    apart->Flattening( phi );

  }
  else
    std::cout << " A correction file is not found. " << binParameter << " with " << ival << std::endl;

}

TVector3 Psi_FlatteningCorrection(UInt_t isel, Int_t ival, TVector3 Pvec)
{
  Int_t    iBIN = GetCorrectionIndex(isel, ival, Pvec.Theta());


  TVector3 Psi_cf;
  if(iBIN >= 0){
    flowcorr = (STFlowCorrection*)aflowcorrArray[isel]->At(iBIN);
    Psi_cf = flowcorr->ReCenteringFourierCorrection(Pvec);
    //    Psi_cf.SetPhi(flowcorr->GetCorrection(Pvec.Phi()));

    if( isel == 9 && ival == 1 )  {
       flowcorr->ShowBinInformation();
      cout << "---##---> ntrack " << ival
	   << " iBIN " << iBIN << " Phi " << Pvec.Phi() << " -> " << Psi_cf.Phi() 
     	   << " X : Y : Z " << Pvec.X() << " " << Pvec.Y() << " " << Pvec.Z() 
     	   << endl;
    } 
  }

  return Psi_cf;
}

TVector3 Psi_ReCenteringCorrection(UInt_t isel, Int_t ival, TVector3 Pvec)
{
  Int_t    iBIN = GetCorrectionIndex(isel, ival, Pvec.Theta());

  TVector3 Psi_cf;
  if(iBIN >= 0 && aflowcorrArray[isel] != NULL) {
    flowcorr = (STFlowCorrection*)aflowcorrArray[isel]->At(iBIN);
    
    Psi_cf = flowcorr->ReCentering(Pvec);
  }

  return Psi_cf;
}
