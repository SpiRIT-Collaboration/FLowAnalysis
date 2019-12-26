#include "./thermalFunction.h"
#include "../tasks/STLorentzBoostVector.hh"
#include "ProcessInfo.C"

Double_t y_cm[]  = { 0.382453, 0.364873, 0.390302, 0.354066};
Double_t y_bm[]  = { 0.360199, 0.377779, 0.354065, 0.390301};
auto mp   =    938.2720813;
Bool_t bRandMult = kTRUE;

void RPSim(UInt_t irun=5, UInt_t maxevt=100)
{
  UInt_t mevent = maxevt;

  //@@@@@
  //UInt_t mlt = 100;
  UInt_t mlt = 40;
  //  UInt_t mlt = 80;
  //UInt_t mlt = 10;

  auto fmult  = TFile::Open("data/mlt108Sn.v41.1.root");
  auto hmult = (TH1I*)fmult->Get("hmult");

  if( hmult != NULL ) {
    std::cout << " Multiplicity distribution is taken from data " << std::endl;
    hmult->Print();
  }
  else
    bRandMult = kFALSE; 


  Bool_t bSingle = kFALSE;
  Bool_t bRPRand = kTRUE;   // 1: Reaction plane is romdom direction

  TDatime dtime;
  TRandom3 grnd(dtime.GetSecond());
  gRandom->SetSeed(dtime.GetSecond());

  //@@-->
  // fgcut = TFile::Open("data/RPSim_AnglegCut.root");
  // auto gcang = (TCutG*)fgcut->Get("gCutAngle");
  // if( gcang == NULL )
  //   std::cout << " Angle cut data is not opned. " << std::endl;
  // fgcut->Close();

  //@@
  auto fgcut = TFile::Open("data/gThetaPhiCut.root");
  auto gcang = (TCutG*)fgcut->Get("gThetaPhiCut");
  if( gcang == NULL )
    std::cout << " Angle cut data is not opned. " << std::endl;
  fgcut->Close();
  //@@<--

  fgcut = TFile::Open("data/hyawpitch_prt_v41.2.10.root");
  auto hyp = (TH2D*)fgcut->Get("hyawpitch");
  TH1D *hyp_x;
  TH1D *hyp_y;
  if( hyp != NULL ) {
    hyp_x = hyp->ProjectionX();
    hyp_y = hyp->ProjectionY();
  }    
  fgcut->Close();

  TFile *fout = new TFile(Form("data/run%04d_rpsim.v0.root",irun),"recreate");
  TTree *ftree= new TTree("cbmsim","Flow simulation");

  UInt_t ic = -1;
  UInt_t id =  1;
  TCanvas *cc;
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1200,1400);
  cc->Divide(3,4);
  
  //  auto vboost = LorentzBoost(0);
  auto vboost = STLorentzBoostVector::GetBoostVector(5);

  //---- 
  auto fdndy  = new TF1("fdndy" ,"gaus",-1.,1.);
  fdndy->SetParameter(0,1.);
  fdndy->SetParameter(1,0.);
  fdndy->SetParameter(2,0.78);

  auto fdndut = new TF1("fdndut","expo", 0.,1.1);
  fdndut->SetParameter(0,0.);
  fdndut->SetParameter(1,-0.1);

  thermalFunction *aTherm = new thermalFunction();
  aTherm->m = mp/1000.;
  aTherm->y = 0.1;
  auto fdndpt = new TF1("fdndpt",aTherm,&thermalFunction::Blastwave,0.,1.,3,"thermalFunction","Blastwave");
  fdndpt->SetParameter(0,2.e+12);
  fdndpt->SetParameter(1,0.04);
  fdndpt->SetParameter(2,0.25);

  if( bSingle ) { 
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic)); }
  else { 
    cc->cd(id); id++;}
  fdndpt->Draw();

  auto fdndp = new TF1("dndp","gaus",0.,1500.);
  fdndp->SetParameter(0,1.);
  fdndp->SetParameter(1,555.);
  fdndp->SetParameter(2,250.);

  //-----------------------------------------


  auto hThetaPhi = new TH2D("hThetaPhi",";#Theta;#Phi",100,0.,TMath::Pi()/2.,100, -TMath::Pi(), TMath::Pi());
  auto hRapPt    = new TH2D("hRapPt"   ,";y_{cm}; P_t",100,-1.,1.,100,0.,500.);
  auto hRapPx    = new TH2D("hRapPx"   ,";y_{cm}; P_x",100,-1.,1.,100,-300.,300.);
  
  auto hThetaPhi_cm = new TH2D("hThetaPhi_cm",";#Theta_{lab};#Phi",100, 0.,TMath::Pi()/2.,100, -TMath::Pi(), TMath::Pi());
  auto hsub_phi      = new TH1D("hsub_phi"," RP; #Psi_{A}" ,100,-TMath::Pi(), TMath::Pi());
  auto hsub1_phi     = new TH1D("hsub1_phi","Sub event A; #Psi_{A}",100,-TMath::Pi(), TMath::Pi());
  auto hsub_dphi     = new TH1D("hsub_dphi","Sub event diff; #Psi_{A}-#Psi_{B}",100,-TMath::Pi(), TMath::Pi());
  auto hsub_12phi    = new TH2D("hsub_12phi","; #Psi_{A};#Psi_{B}",100,-TMath::Pi(), TMath::Pi(),100,-TMath::Pi(), TMath::Pi());
  auto hphi_yb       = new TH1D("hphi_yb"   ,"y<(0.5,1.0); #phi",100,-TMath::Pi(), TMath::Pi());
  auto hphi_yt       = new TH1D("hphi_yt"   ,"y<(-1,-0.5); #phi",100,-TMath::Pi(), TMath::Pi());
  auto hphi_ym       = new TH1D("hphi_ym"   ,"y<(-0.2,0.2); #phi",100,-TMath::Pi(), TMath::Pi());
  auto hcos_yb       = new TH1D("hcos_yb"   ,"y<( 0.5,1.0); #phi",100,-1., 1.);
  auto hcos_yt       = new TH1D("hcos_yt"   ,"y<(-1.,-0.5); #phi",100,-1., 1.);
  auto hcos_ym       = new TH1D("hcos_ym"   ,"y<(-0.2,0.2); #phi",100,-1., 1.);
  auto hyawpitch     = new TH2D("hyawpitch" ,"; pitch; yaw",100,-TMath::Pi()/2., TMath::Pi()/2.,100,-TMath::Pi()/2., TMath::Pi()/2.);

  //---- Flow parameter -----
  //-----------------------------------------
  auto fazim = new TF1("fazim","1+2.*[0]*cos(x-[2]) + 2.*[1]*cos(2.* (x-[2]))",-TMath::Pi(), TMath::Pi());
  fazim->SetParameter(2, 0.);

  auto fv1y = new TF1("fv1y","[0]+[1]*x+[2]*x^3"  ,-1.,1.);
  fv1y->SetParameter(0,0);
  fv1y->SetParameter(1,5.18056e-01);
  fv1y->SetParameter(1,5.18056e-01);
  fv1y->SetParameter(2,-1.84025e-01);

  //@@@@@
  auto fv2y = new TF1("fv2y","[0]+[1]*x^2+[2]*x^4",-1.,1.);
  fv2y->SetParameter(0, -0.08);
  fv2y->SetParameter(1,  0.1);
  fv2y->SetParameter(2, -0.02);

  // fv2y->SetParameter(0, -0.08);
  // fv2y->SetParameter(1,  0.);
  // fv2y->SetParameter(2,  0.);

  if( bSingle ) {
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic)); }
  else { 
    cc->cd(id); id++;}
  
  fv1y->Draw("");

  if( bSingle ) {
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic)); }
  else { 
    cc->cd(id); id++;}
  fv2y->Draw("");

  TClonesArray    aPartArray("STParticle",80);
  TClonesArray    aFlow("STFlowInfo",1);
  STParticle*     aParticle = new STParticle();
  STFlowInfo*     aFlowInfo = new STFlowInfo();
  STFlowTask*     aFlowTask = new STFlowTask(kFALSE, kTRUE, kFALSE);
  aFlowTask -> Init(irun, "0");
  
  Double_t        RPPsi = 0.;

  ftree->Branch("STParticle",&aPartArray);
  ftree->Branch("STFlow"    ,&aFlow);
  ftree->Branch("RPPsi"     ,&RPPsi);

  for( UInt_t j = 0; j < mevent; j++) {

    ProcessInfo((Long64_t)mevent, (Long64_t)j);

    aPartArray.Clear();
    aFlow.Clear();
    aFlowInfo->Clear();

    // Random reaction plane angle
    RPPsi = bRPRand == kTRUE ? 2.*TMath::Pi()*(grnd.Rndm() - 0.5) : 0.;

    if( bRandMult ) { 
      mlt = (Int_t)hmult->GetRandom();
      if( mlt < 10 ) continue;
    }


    UInt_t i = 0;
    while( i < mlt ){

      aParticle->Clean();

      // cm frame      
      Double_t aRapcm= fdndy->GetRandom();    

      Double_t aV1   = fv1y->Eval(aRapcm);
      Double_t aV2   = fv2y->Eval(aRapcm);
      fazim->SetParameter(0, aV1);
      fazim->SetParameter(1, aV2);
      fazim->SetParameter(2, RPPsi);
      Double_t aPhi = fazim->GetRandom();

      //@@@@
      // if( i >= mlt*0.9 )
      //  	aPhi = 2.*TMath::Pi()*(grnd.Rndm() - 0.5);

      //Lab frame
      Double_t aPt    = fdndpt->GetRandom()*1000.;
      Double_t aRap   = (aRapcm + 1.0)*y_cm[0];  
      Double_t aEovPz = (exp(2.*aRap) + 1.)/(exp(2.*aRap) - 1.);
      Double_t aTan   = sqrt( aPt*aPt * ( aEovPz*aEovPz - 1.)/(mp*mp + aPt*aPt) );
      Double_t aTheta = TMath::ATan(aTan);

      Double_t aPx  = aPt*cos(aPhi);
      Double_t aPy  = aPt*sin(aPhi);
      Double_t aPz  = aPt/aTan;

      aParticle->SetTrackID(i);
      aParticle->SetPID(2212);
      aParticle->SetPIDLoose(2212);
      aParticle->SetMass(2212);      
      aParticle->SetRotatedMomentum( TVector3(aPx, aPy, aPz ) );

      auto yaw   = aParticle->GetYawAngle();
      auto pitch = aParticle->GetPitchAngle();
      auto vecP   = aParticle->GetRotatedMomentum();

      //@@@@@
      //      if( !gcang->IsInside( yaw, pitch ) )continue; //v5
      //   if( !gcang->IsInside( pitch, yaw ) )continue; //v7, v8 OK
      if( !gcang->IsInside(vecP.Theta() , vecP.Phi() )) continue; //v7, v8 OK
      
      hyawpitch->Fill(yaw, pitch);
      hThetaPhi->Fill(aTheta, aPhi);
      hRapPx->Fill(aRapcm, aPx); 
      hRapPt->Fill(aRapcm, aPt);

      aParticle->SetFromTargetFlag(1);
      aParticle->SetDistanceAtVertex(0.);
      aParticle->SetGoodTrackFlag(1);
      aParticle->SetReactionPlaneFlag(1111);

      TLorentzVector aLorentz(aParticle->GetLorentzVector());
      aLorentz.Boost(-vboost);

      hThetaPhi_cm->Fill( aLorentz.Theta(), aLorentz.Phi() );

      new( aPartArray[i] ) STParticle(*aParticle);

      if( aRapcm >= 0.5 && aRapcm <= 1 ) {
	hphi_yb->Fill(aParticle->GetRotatedMomentum().Phi());
	hcos_yb->Fill(cos(aParticle->GetRotatedMomentum().Phi()));
      }

      if( aRapcm >= -1 && aRapcm <= -0.5 ) {
	hphi_yt->Fill(aParticle->GetRotatedMomentum().Phi());
	hcos_yt->Fill(cos(aParticle->GetRotatedMomentum().Phi()));
      }

      if( abs(aRapcm) <= 0.2 ) {
	hphi_ym->Fill(aParticle->GetRotatedMomentum().Phi());
	hcos_ym->Fill(cos(2.*aParticle->GetRotatedMomentum().Phi()));
      }      

      i++;
    }

    aFlowTask->SetFlowTask( aPartArray );
    aFlowTask->FinishEvent();
    aFlowInfo = aFlowTask->GetFlowInfo();
    aFlowInfo->goodEventf = 1;
    aFlowInfo->SetBeamPID(100);
    aFlowInfo->evt = j;

    hsub_phi ->Fill(aFlowInfo->unitP.Phi());
    hsub1_phi->Fill(aFlowInfo->unitP_1.Phi());
    hsub_dphi->Fill(TVector2::Phi_mpi_pi(aFlowInfo->unitP_1.Phi()-aFlowInfo->unitP_2.Phi()));

    hsub_12phi->Fill(aFlowInfo->unitP_1.Phi(),aFlowInfo->unitP_2.Phi());

    new( aFlow[0] ) STFlowInfo(*aFlowInfo); 

    ftree->Fill();
  }

  std::cout << " RP angle " << RPPsi << std::endl;
  std::cout << " beam   rapidity <cos>  = " << hcos_yb->GetMean() << " +- " << hcos_yb->GetStdDev()/sqrt((Double_t)hcos_yb->GetEntries()) << std::endl;
  std::cout << " Target rapidity <cos>  = " << hcos_yt->GetMean() << " +- " << hcos_yt->GetStdDev()/sqrt((Double_t)hcos_yt->GetEntries()) << std::endl;
  std::cout << " Mid    rapidity <cos2> = " << hcos_ym->GetMean() << " +- " << hcos_ym->GetStdDev()/sqrt((Double_t)hcos_ym->GetEntries()) << std::endl;



  if( bSingle ) { 
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic)); }
  else { 
    cc->cd(id); id++;}
  hphi_yb->SetNormFactor(100);
  hphi_yb->SetLineColor(3);
  hphi_yb->Draw("");
  hphi_ym->SetNormFactor(100);
  hphi_ym->SetLineColor(2);
  hphi_ym->Draw("same");
  hphi_yt->SetNormFactor(100);
  hphi_yt->SetLineColor(4);
  hphi_yt->Draw("same");

  if( bSingle ) { 
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic)); }
  else { 
    cc->cd(id); id++;}
  hsub_phi->Draw();

  if( bSingle ) { 
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic)); }
  else { 
    cc->cd(id); id++;}
  hsub_dphi->Draw();

  if( bSingle ) { 
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic)); }
  else { 
    cc->cd(id); id++;}
  hsub_12phi->Draw("colz");

  if( bSingle ) { 
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic)); }
  else { 
    cc->cd(id); id++;}
  hThetaPhi->Draw("colz");

  if( bSingle ) {
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic)); }
  else { 
    cc->cd(id); id++;}
  hRapPt->Draw("colz");

  if( bSingle ) {
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic)); }
  else { 
    cc->cd(id); id++;}
  hRapPx->Draw("colz");  

  if( bSingle ) {
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic)); }
  else { 
    cc->cd(id); id++;}

  hyawpitch->Draw("colz");
  
  ftree->Write("",TObject::kWriteDelete);
  fv1y->Write();
  fv2y->Write();
  hThetaPhi->Write();
  hyawpitch->Write();

  std::cout << fout->GetName() << " is created. " << std::endl;

}
