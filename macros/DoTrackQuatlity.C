#include "DoRPRes.C" //#include "openRunAna.C" #include "DoFlow.h" #include "SimFunction.C" 

void DoTrackQuatlity()
{
  // Command
  //->   RUN={$RNF132} VER=43.0 root DoTrackQuality.C 
  //**************************************************

  gROOT->Reset();
  
  openRunAna();
  if(rChain != NULL)     
    LOG(INFO) << " DoFLow_adv: System " << isys << "  -> " << sysName << FairLogger::endl; 

  else
    exit(0);

  // Configuration ================================
  TString su = gSystem -> Getenv("UC");
  if( su != "" ) {
    Ucent = (UInt_t)atoi(su);
    if( Ucent < 0 )
      Ucent = 80;
  }

  su = gSystem -> Getenv("LC");
  if( su != "" ) {
    Lcent = (UInt_t)atoi(su);
    if( Lcent > Ucent || Lcent < 0)
      Lcent = 0;
  }

  LOG(INFO) << " Multiplicity :: " << Lcent << " to " << Ucent << FairLogger::endl;



  TH2D *hac20 = new TH2D("hac20","Proton        ; y_{cm}/y_{beam}; Pt[MeV/c]"  ,200,-0.5,1.1,200,0.,800);
  TH2D *hac21 = new TH2D("hac21","Proton m0to45 ; y_{cm}/y_{beam}; Pt[MeV/c]"  ,200,-0.5,1.1,200,0.,800);
  TH2D *hac22 = new TH2D("hac22","Proton m20to40; y_{cm}/y_{beam}; Pt[MeV/c]"  ,200,-0.5,1.1,200,0.,800);
  TH2D *hac30 = new TH2D("hac30","Deuteron      ; y_{cm}/y_{beam}; Pt[MeV/c]"  ,200,-0.5,1.1,200,0.,800);
  TH2D *hac31 = new TH2D("hac31","Deuteron m0to45 ; y_{cm}/y_{beam}; Pt[MeV/c]",200,-0.5,1.1,200,0.,800);
  TH2D *hac32 = new TH2D("hac32","Deuteron m20to40; y_{cm}/y_{beam}; Pt[MeV/c]",200,-0.5,1.1,200,0.,800);
  TH2D *hac40 = new TH2D("hac40","Triton        ; y_{cm}/y_{beam}; Pt[MeV/c]"  ,200,-0.5,1.1,200,0.,800);
  TH2D *hac41 = new TH2D("hac41","Triton m0to45 ; y_{cm}/y_{beam}; Pt[MeV/c]"  ,200,-0.5,1.1,200,0.,800);
  TH2D *hac42 = new TH2D("hac42","Triton m20to40; y_{cm}/y_{beam}; Pt[MeV/c]"  ,200,-0.5,1.1,200,0.,800);
  TH2D *hac50 = new TH2D("hac50","3He        ; y_{cm}/y_{beam}; Pt[MeV/c]"     ,200,-0.5,1.1,200,0.,800);
  TH2D *hac51 = new TH2D("hac51","3He m0to45 ; y_{cm}/y_{beam}; Pt[MeV/c]"     ,200,-0.5,1.1,200,0.,800);
  TH2D *hac52 = new TH2D("hac52","3He m20to40; y_{cm}/y_{beam}; Pt[MeV/c]"     ,200,-0.5,1.1,200,0.,800);


  TH1D *hdy20 = new TH1D("hdy20","Proton     m50    ; y_{cm}"    ,200,-0.5,1.1);
  TH1D *hdy21 = new TH1D("hdy21","Proton     m0to45 ; y_{cm}"    ,200,-0.5,1.1);
  TH1D *hdy22 = new TH1D("hdy22","Proton     m20to40; y_{cm}"    ,200,-0.5,1.1);
  TH1D *hdy30 = new TH1D("hdy30","Deuteron   m50    ; y_{cm}"    ,200,-0.5,1.1);
  TH1D *hdy31 = new TH1D("hdy31","Deuteron   m0to45 ; y_{cm}"    ,200,-0.5,1.1);
  TH1D *hdy32 = new TH1D("hdy32","Deuteron   m20to40; y_{cm}"    ,200,-0.5,1.1);
  TH1D *hdy40 = new TH1D("hdy40","Triton     m50    ; y_{cm}"    ,200,-0.5,1.1);
  TH1D *hdy41 = new TH1D("hdy41","Triton     m0to45 ; y_{cm}"    ,200,-0.5,1.1);
  TH1D *hdy42 = new TH1D("hdy42","Triton     m20to40; y_{cm}"    ,200,-0.5,1.1);
  TH1D *hdy50 = new TH1D("hdy50","3He        m50    ; y_{cm}"    ,200,-0.5,1.1);
  TH1D *hdy51 = new TH1D("hdy51","3He        m0to45 ; y_{cm}"    ,200,-0.5,1.1);
  TH1D *hdy52 = new TH1D("hdy52","3He        m20to40; y_{cm}"    ,200,-0.5,1.1);


  //-------------------------------------------------- 
  //--- Event Loop
  //--------------------------------------------------  
  TDatime beginTime;
  TDatime dtime;

  Int_t nevt = SetBranch();
  Int_t nevt_begin = 0;
  Int_t maxevt = nevt_begin+nevt;
  nevt = maxevt < nevt ? maxevt : nevt;
  LOG(INFO) << " NOTICE !!!! " << nevt_begin << " to " << nevt << FairLogger::endl;

  for(Int_t i = nevt_begin; i < nevt; i++){

    rChain->GetEntry(i);    

    beginTime.Copy(dtime);
    if(i%(UInt_t)(nevt/50) == 0) {
      dtime.Set();
      Int_t ptime = dtime.Get() - beginTime.Get();

      LOG(INFO) << "Processing .... " 
		<< setw(4) << Int_t(((Double_t)i/(Double_t)nevt)*100.) << " % = "
		<< setw(8) << i << "/"<< nevt
		<< "--->"
		<< dtime.AsString() << " ---- "
		<< FairLogger::endl;
    }


    STFlowInfo *aflow = (STFlowInfo*)aFlowArray->At(0);
    //    if(aflow->mtrack2 > Ucent || aflow->mtrack2 <= Lcent || aflow->mtrack4 < 6) continue;
    Bool_t bmt[3] = {kFALSE, kFALSE, kFALSE};
    if( aflow->mtrack2 >= 50 ) bmt[0] = kTRUE;
    if( aflow->mtrack2 <= 45 && aflow->mtrack2 >  0 ) bmt[1] = kTRUE; 
    if( aflow->mtrack2 <= 40 && aflow->mtrack2 > 20 ) bmt[2] = kTRUE; 


    TIter next(aArray);
    STParticle *aPart = NULL;
    //----------------------------------------
    //----- Main loop 
    //--------------------------------------------------
    UInt_t mtk = 0;
    while( (aPart = (STParticle*)next()) ) {


      auto pt    = aPart->GetRotatedMomentum().Pt();
      auto rapid  = aPart->GetRapiditycm();;	

      auto pid       = aPart->GetPIDTight();  
      auto pid_los   = aPart->GetPIDLoose();  
      auto pid_nrm   = aPart->GetPIDNorm();  

      auto yaw   = aPart->GetYawAngle();
      auto pitch = aPart->GetPitchAngle();

      //*********************************************
      // ---- track quality selection ---
      if( aPart->GetGoodTrackFlag() != 1111 ) continue;
      //------------------------------
      //      if( (pitch/yaw) > 1. || yaw < 0 ) continue;
      //if( (pitch/yaw) > 1. || yaw > 0 ) continue;

      //if ( abs( pitch / yaw ) > 1. ) continue;
      //if( abs( pitch / yaw ) <= 1. ) continue;
      //*********************************************

      
      if( bmt[0] ) {
	if( pid == 2212)       hac20->Fill(rapid,pt);
	if( pid == 1000010020) hac30->Fill(rapid,pt);
	if( pid == 1000010030) hac40->Fill(rapid,pt);
	if( pid == 1000020030) hac50->Fill(rapid,pt);

	if( pid == 2212)       hdy20->Fill(rapid);
	if( pid == 1000010020) hdy30->Fill(rapid);
	if( pid == 1000010030) hdy40->Fill(rapid);
	if( pid == 1000020030) hdy50->Fill(rapid);
      }


      if( bmt[1] ) {
	if( pid == 2212) hac21->Fill(rapid,pt);
	if( pid == 1000010020) hac31->Fill(rapid,pt);
	if( pid == 1000010030) hac41->Fill(rapid,pt);
	if( pid == 1000020030) hac51->Fill(rapid,pt);

	if( pid == 2212)       hdy21->Fill(rapid);
	if( pid == 1000010020) hdy31->Fill(rapid);
	if( pid == 1000010030) hdy41->Fill(rapid);
	if( pid == 1000020030) hdy51->Fill(rapid);
      }

      if( bmt[2] ) {
	if( pid == 2212) hac22->Fill(rapid,pt);
	if( pid == 1000010020) hac32->Fill(rapid,pt);
	if( pid == 1000010030) hac42->Fill(rapid,pt);
	if( pid == 1000020030) hac52->Fill(rapid,pt);

	if( pid == 2212)       hdy22->Fill(rapid);
	if( pid == 1000010020) hdy32->Fill(rapid);
	if( pid == 1000010030) hdy42->Fill(rapid);
	if( pid == 1000020030) hdy52->Fill(rapid);
      }

    }
  }

  TCanvas *cvx1 = new TCanvas("cvx1","cvx1",1200,1500);
  cvx1->Divide(3,4);
  UInt_t idd = 1;

  cvx1->cd(idd); idd++;   hac20->Draw("colz");
  cvx1->cd(idd); idd++;   hac21->Draw("colz");
  cvx1->cd(idd); idd++;   hac22->Draw("colz");

  cvx1->cd(idd); idd++;   hac30->Draw("colz");
  cvx1->cd(idd); idd++;   hac31->Draw("colz");
  cvx1->cd(idd); idd++;   hac32->Draw("colz");

  cvx1->cd(idd); idd++;   hac40->Draw("colz");
  cvx1->cd(idd); idd++;   hac41->Draw("colz");
  cvx1->cd(idd); idd++;   hac42->Draw("colz");

  cvx1->cd(idd); idd++;   hac50->Draw("colz");
  cvx1->cd(idd); idd++;   hac51->Draw("colz");
  cvx1->cd(idd); idd++;   hac52->Draw("colz");

  TCanvas *cvx2 = new TCanvas("cvx2","cvx2",1200,1500);
  cvx2->Divide(3,4);
  idd = 1;

  cvx2->cd(idd); idd++;   hdy20->Draw();
  cvx2->cd(idd); idd++;   hdy21->Draw();
  cvx2->cd(idd); idd++;   hdy22->Draw();

  cvx2->cd(idd); idd++;   hdy30->Draw();
  cvx2->cd(idd); idd++;   hdy31->Draw();
  cvx2->cd(idd); idd++;   hdy32->Draw();

  cvx2->cd(idd); idd++;   hdy40->Draw();
  cvx2->cd(idd); idd++;   hdy41->Draw();
  cvx2->cd(idd); idd++;   hdy42->Draw();

  cvx2->cd(idd); idd++;   hdy50->Draw();
  cvx2->cd(idd); idd++;   hdy51->Draw();
  cvx2->cd(idd); idd++;   hdy52->Draw();


}
