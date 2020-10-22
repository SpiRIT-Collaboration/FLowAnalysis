#include "DoRPRes.C" //#include "openRunAna.C" #include "DoFlow.h" #include "SimFunction.C"

//drawing

UInt_t selReactionPlanef = 10000;

TFile* fefffile;

Bool_t bccPsi;


// functions
void     PlotPtDependence(UInt_t selid);
TString  SetupOutputFile(TString fopt);
Bool_t   SetupEffCorrection();

Bool_t   LoadAcceptanceCorrection(UInt_t ipid);
Double_t *GetPsiRPResolution (Double_t *rpcor, UInt_t ival);
Double_t *Get2PsiRPResolution(Double_t *rpcor, UInt_t ival);
Double_t *GetMultRPResolution(Double_t *rpcor, UInt_t vn, UInt_t mult);

UInt_t   GetV1RapidityIndex(Double_t y);
UInt_t   GetV2RapidityIndex(Double_t y);
UInt_t   GetV1cmRapidityIndex(Double_t y);
UInt_t   GetV2cmRapidityIndex(Double_t y);
UInt_t   GetV1PtIndex(Double_t val);
UInt_t   GetV2PtIndex(Double_t val);

void     DrawUt(TH2D *h2, UInt_t sel, TString opt="");
TString  GetPsiRPLoadFileName();
void     GetdNdydPt(Double_t ycm, Double_t pt);
void     GetdNdydUt( Double_t ycm, Double_t ut, TH2D* h2c);
void     GetdNdydUt( Int_t ih, Int_t iy, Int_t ipt, TH2D *h2c);
void     GetYUtCorrection(TH2D &h2, TH2D &h2c);

TFile*        GraphSave;
TH2D*         hAcpCorr;
TH2D*         hAcpYUtCorr;

Double_t      fmass;
Double_t      *acpcorr = new Double_t[2];
//-------------------//
void DoFlow_fin(Int_t isel = 0) 
{
  gROOT->Reset();

  openRunAna();

  if(rChain != NULL) 
    LOG(INFO) << " DoFLow_adv: System " << isys << "  -> " << sysName << FairLogger::endl; 

  else
    exit(0);


  gROOT->ProcessLine(".! grep -i void DoFlow_adv.C ");


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

  LOG(INFO) << "Multiplicity :: " << Lcent << " to " << Ucent << FairLogger::endl;

  //--------------------------------------------------
  TString rpBase = gSystem -> Getenv("RPBS");
  RPBase = rpBase != "" ? atoi(rpBase): 0;

  LOG(INFO) << "Reaction plane base is " << RPBase << FairLogger::endl;

  TString ccPsi = gSystem -> Getenv("CCPSI");
  bccPsi = ccPsi == "" ? kTRUE: kFALSE;
  LOG(INFO) << "Angle dependent correction is ";
  if( bccPsi )
    LOG(INFO) << " !!! ON !!!"  << FairLogger::endl;
  else
    LOG(INFO) << " OFF " << FairLogger::endl;


  //==================================================
  acpcorr[0] = 1.;
  acpcorr[1] = 0.;

  dpt1 = ut_max[isel-2]/(Double_t)ptbin1;
  dpt2 = ut_max[isel-2]/(Double_t)ptbin2;


  if( isel > 0 ) 
    PlotPtDependence((UInt_t)isel);  
  else if( isel == -1 )
    PlotPtDependence(0);  
}


//@@/--- pdt ------------------------------------------------------- 
//******************** Main ****************************************
void PlotPtDependence(UInt_t selid = 2)       //%% Executable :
{

  TDatime beginTime;
  TDatime dtime;

  gStyle->SetOptStat(0);


  ///------>>>> Acceptance Correction Option: --------
  ///$$$$$////
  //Bool_t bAcc_corr = kFALSE;
  Bool_t bAcc_corr = kTRUE;

  if( bAcc_corr && LoadAcceptanceCorrection(selid) ) {
    LOG(INFO) << " Acceptance correction is found. " << FairLogger::endl;
  }
  else {
    LOG(INFO) << " Acceptance correction is OFF." << FairLogger::endl;
    bAcc_corr = kFALSE;
  }


  LOG(INFO) << "PlotPtDependence(" << selid << ")" << FairLogger::endl;

  //-- Define Output file name 
  TString fHeader = "finYPt_"+ sysName + "_" + partname[selid]+".v"+sVer+".";
  if( oVer != "" )
    fHeader = "finYPt_"+ sysName + "_" + partname[selid]+".v"+oVer+".";

  auto fName = SetupOutputFile( fHeader );
  
  Bool_t bEffCorr = kFALSE;
  if( selid == 2 || selid == 4 )
    auto bEffCorr = SetupEffCorrection();

  if( bEffCorr )
    LOG(INFO) << fefffile->GetName() << FairLogger::endl;

  GraphSave  = new TFile(fName,"recreate");


  LOG(INFO) << " Rapidity binning " << ybin1 << FairLogger::endl;

  //------------------------------
  LOG(INFO) << " v1 BIN y " << ybin1 << " pt " << ptbin1 << " mult " << mbin << FairLogger::endl;
  LOG(INFO) << " v2 BIN y " << ybin2 << " pt " << ptbin2 << " mult " << mbin << FairLogger::endl;

  auto hpid       = new TH2D("hpid",   "All particles    ; P/Z[MeV/c]; dEdx[ADC/mm]",400,-800.,3000.,300,0.,1500);
  auto hpidsel    = new TH2D("hpidsel","selected particle; P/Z[MeV/c]; dEdx[ADC/mm]",400,-800.,3000.,300,0.,1500);
  auto hyawpitch  = new TH2D("hyawpitch","; Yaw angle; Pitch angle",200,-1.5,1.5,200,-1.5,1.5);
  auto hmass      = new TH2D("hmass",   ";P/Q; Mass [MeV/c^2]"     ,200,  0.,2500., 200, 0.,7000);
  TString hlabel  = (TString)Form("mtrack1 %2d to %2d ; Y_{cm}; Pt [MeV/c]",Lcent,Ucent);

  auto hypteff    = new TH2D("hypteff", hlabel , 40, -2., 2., 50, 0., 2.5);
  auto hyuteff    = new TH2D("hyuteff", hlabel , 40, -2., 2., 50, 0., 2.5);
  auto hdndydut   = new TH2D("hdndydut",";y;ut", 40, -2., 2., 50, 0., 2.5);
  // auto hyutacp    = new TH2D("hyutacp"   , hlabel ,200, -1., 1.4, 100, 0., 1.5);  // w/o corr. y-ut
  // auto hyutacpcrr = new TH2D("hyutacpcrr", hlabel ,200, -1., 1.4, 100, 0., 1.5);  // corrected y-ut
  auto hyutacp    = new TH2D("hyutacp"   , hlabel , 40, -2., 2., 50, 0., 2.5);
  auto hyutacpcrr = new TH2D("hyutacpcrr", hlabel , 40, -2., 2., 50, 0., 2.5);

  auto hyptacp    = new TH2D("hyptacp"   , hlabel ,200, -1., 1.4, 200, 0., 1100);
  auto huy        = new TH2D("huy","; y_{cm}/y_{beam}; u_{t0} ", 200,-1., 1.8, 200,0.,2.5);
  auto hpsi       = new TH2D("hpsi",";#Psi;",15, 0., 15., 200,0.,1.);
  auto hpsi1      = new TH1D("hpsi1",";#Psi",500,0.3,0.8);
  auto hpsindx    = new TH2D("hpsindx",";#Psi ;index ",100,-TMath::Pi(),TMath::Pi(),20,0.,20.);
  auto hiphi      = new TH1I("hiphi","hiphi",13,0,13);
  auto hPsi       = new TH1D("hPsi"   ,";RP #Psi"  ,100,-TMath::Pi(),TMath::Pi());
  auto hPsinc     = new TH1D("hPsinc" ,"w/o Cut;RP #Psi"  ,100,-TMath::Pi(),TMath::Pi());
  auto hRPPsi     = new TH1D("hRPPsi",";Indiv. #Psi"  ,100,-TMath::Pi(),TMath::Pi());
  auto hRPPsipsi  = new TH2D("hRPPsipsi",";RP #Psi; Indiv #Psi"  ,100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi());
  auto hphi       = new TH1D("hphi",     ";#phi"  ,100,-TMath::Pi(),TMath::Pi());

  TH1I *hmult = new TH1I("hmult",";Multiplicity",100,0,100);
  TH1D *hdphi0_180   = new TH1D("hdphi0_180"  , " 0to180",100,0.,3.2);
  TH1D *hdphi90_180  = new TH1D("hdphi90_180" , "90to180",100,0.,3.2);

  TH2D *hirap1 = new TH2D("hirap1",";Rapidity; irap1",100,-1.0,1.5, 12,0,12);  
  TH2D *hirap2 = new TH2D("hirap2",";Rapidity; irap1",100,-1.0,1.5, 12,0,12);  


  TH1D *hutphi10[ybin1];
  TH1D *hdy1[ybin1];               // y_nrm
  TH1D *hdyv1[ybin1];
  TH1D *hdycos1[ybin1];              // cos vs y_nrm
  TH1D *hdydut1[ybin1][ptbin1];
  TH1D *hdydutcos1[ybin1][ptbin1];
  TH1D *hdydut1cut[ybin1][ptbin1];
  TH1D *hdydutcos1cut[ybin1][ptbin1];

  TH1D *hutphi20[ybin2];
  TH1D *hdy2[ybin2];
  TH1D *hdyv2[ybin2];
  TH1D *hdycos2[ybin2];
  TH1D *hdydut2[ybin2][ptbin2];
  TH1D *hdydutcos2[ybin2][ptbin2] ;
  TH1D *hdydut2cut[ybin2][ptbin2];
  TH1D *hdydutcos2cut[ybin2][ptbin2];

  //-----  booking
  TString rangeLabel1[ybin1];
  TString rangeLabel2[ybin2];

  // replace y_cm -> y_bm

  // booking for v1
  for( UInt_t iy = 0; iy < ybin1; iy++ ) {
    rangeLabel1[iy] = Form("%5.2f <= y < %5.2f"      ,yrange1nrm[iy]/y_cm[isys+6],yrange1nrm[iy+1]/y_cm[isys+6]);

    hutphi10[iy]   = new TH1D(Form("hutphi10+%d",iy)  ,rangeLabel1[iy]+";sub#Delta #Psi",100, -TMath::Pi(), TMath::Pi());
    hdy1[iy]       = new TH1D(Form("hdy1_y%d",iy)     ,rangeLabel1[iy]+";Rapidity", 500,-2.5,2.5);
    hdyv1[iy]      = new TH1D(Form("hdyv1_%d",iy)     ,rangeLabel1[iy]+"ut>0.4;<cos#phi>;",100,-1.,1.);
    hdycos1[iy]    = new TH1D(Form("hdycos1%dp",iy)   ,rangeLabel1[iy]+"; cos(#phi - #Psi)"  , 100, -1., 1.);

    for( UInt_t ipt = 0; ipt < ptbin1; ipt++ ) {
      hdydut1[iy][ipt]    = new TH1D(Form("hdydut1_y%dpt%d",iy,ipt),rangeLabel1[iy]+"; Pt", 500,0., 2.5);
      hdydutcos1[iy][ipt] = new TH1D(Form("hdydutcos1_y%dpt%d",iy,ipt),rangeLabel1[iy]+";u_{t0}", 100, -1., 1.);

      hdydut1cut[iy][ipt]    = new TH1D(Form("hdydut1cut_y%dpt%d",iy,ipt),rangeLabel1[iy]+"; Pt", 500,0., 2.5);
      hdydutcos1cut[iy][ipt] = new TH1D(Form("hdydutcos1cut_y%dpt%d",iy,ipt),rangeLabel1[iy]+";u_{t0}", 100, -1., 1.);
    }
  }
  

  for( UInt_t iy = 0; iy < ybin2; iy++ ) {
    rangeLabel2[iy] = Form("%5.2f <= y < %5.2f",yrange2nrm[iy]/y_cm[isys+6],yrange2nrm[iy+1]/y_cm[isys+6]);
    hutphi20[iy]    = new TH1D(Form("hutphi20+%d",iy)  ,rangeLabel1[iy]+";sub#Delta #Psi",100, -TMath::Pi(), TMath::Pi());
    hdy2[iy]        = new TH1D(Form("hdy2_y%d",iy),rangeLabel2[iy]+";Rapidity", 500,-2.5,2.5);
    hdyv2[iy]       = new TH1D(Form("hdyv2_y%d",iy),rangeLabel2[iy]+"ut>0.4;<cos#phi>;Rapidity", 500,-2.5,2.5);
    hdycos2[iy]     = new TH1D(Form("hdycos2%dp",iy),rangeLabel2[iy]+"; cos2(#phi - #Psi)" , 100, -1., 1.);

    for( UInt_t ipt = 0; ipt < ptbin2; ipt++ ) {
      hdydut2[iy][ipt]       = new TH1D(Form("hdydut2_y%dpt%d",iy,ipt),rangeLabel2[iy]+"; Pt", 500,0., 2.5);
      hdydutcos2[iy][ipt]    = new TH1D(Form("hdydutcos2_y%d_ut%d",iy,ipt),rangeLabel2[iy]+"; cos(#phi - #Psi)" , 100, -1., 1.);

      hdydut2cut[iy][ipt]    = new TH1D(Form("hdydut2cut_y%dpt%d",iy,ipt),rangeLabel2[iy]+"; Pt", 500,0., 2.5);
      hdydutcos2cut[iy][ipt] = new TH1D(Form("hdydutcos2cut_y%dpt%d",iy,ipt),rangeLabel2[iy]+";u_{t0}", 100, -1., 1.);
    }
  }

  //------------------------------
  Double_t u_p = 0.355151 * 1.06974; //p+p(268.9MeV/u)

  // Double_t *rppsires  = new Double_t[3];  // Psi dependent correction
  // Double_t *rpres     = new Double_t[4];  // 

  Long64_t nEntry = SetBranch();
  
  //--------------------------------------------------
  //--- Event Loop
  //--------------------------------------------------

  for(Long64_t i = 0; i < nEntry; i++){

    ShowProcess(i);

    rChain->GetEntry(i);
    STFlowInfo *aflow = (STFlowInfo*)aFlowArray->At(0);

    /// for reaction plane resolution
    Bool_t bFill = kFALSE;
    Bool_t bRes  = kFALSE;


    /// centrality selection
    if(aflow->mtrack2 > Ucent || aflow->mtrack2 <= Lcent || aflow->mtrack4 < 6) continue;

    hmult->Fill( aflow->mtrack2 );
    ///    auto RPangle = GetRPBaseAngle(aflow);
    auto RPangle = aflow->unitP_fc.Phi();
    hPsinc->Fill(RPangle);


    bRes = kTRUE; //@1

    Double_t subevt_phi = abs(TVector2::Phi_mpi_pi((aflow->unitP_1fc).Phi()-
						   (aflow->unitP_2fc).Phi()));


    TIter next(aArray);
    STParticle *aPart = NULL;

    //--------------------------------------------------
    //----- Main loop 
    //--------------------------------------------------
    UInt_t mtk = 0;

    while( (aPart = (STParticle*)next()) ) {

      mtk++;
      //&&&&&
      if( isys != 5 && 
	  //	  ((selid >= 2 && aPart->GetReactionPlaneFlag()%2 != 1) ||
	  //	   (selid <= 1 && aPart->GetGoodTrackFlag() != 1111)) )

	  (aPart->GetGoodTrackFlag() != 1111 || aPart->GetReactionPlaneFlag() == 0 ))
	  //	  (aPart->GetGoodTrackFlag() != 1111))
	continue;
      //------------------------------

      auto yaw   = aPart->GetYawAngle();
      auto pitch = aPart->GetPitchAngle();
      auto pt    = aPart->GetRotatedMomentum().Pt();
      auto mom   = aPart->GetRotatedMomentum().Mag();
      auto bmass = aPart->GetBBMass();
      auto phi   = aPart->GetRotatedMomentum().Phi();
      auto theta = aPart->GetRotatedMomentum().Theta();
      auto charge= aPart->GetCharge();

      auto pid   = aPart->GetPID();  
      //auto pid   = aPart->GetPIDTight();
      //  auto pid   = aPart->GetPIDLoose();  
      //  auto pid   = aPart->GetPIDNorm();  


      //----- Particle Selection -----
      Bool_t bpart = kFALSE;
      auto MassNumber = partA[selid];
      if( selid == 8 ) {   //H
	if( pid == partid[2] ) {
	  MassNumber = partA[2];
	  bpart = kTRUE;
	}
	else if( pid == partid[3] ){
	  MassNumber = partA[3];
	  bpart = kTRUE;
	}
	else if( pid == partid[4] ) {
	  MassNumber = partA[4];
	  bpart = kTRUE;
	}
      }
      else if( selid == 9 ){ 
	if( pid == partid[4] && mom > 1000  ) {
	  pid = partid[4];
	  MassNumber = partA[4];
	  bpart = kTRUE;
	}
      }
      else if( selid == 1 && charge == 1 && aPart->GetPIDTight() == partid[0] &&
	       aPart->GetBBMass() >= 100 && aPart->GetBBMass() <= 160.  ) //pi+
	bpart = kTRUE;

      else if( selid == 0 && charge == -1 && aPart->GetPIDTight() == partid[0] &&
	       aPart->GetBBMass() >= 100 && aPart->GetBBMass() <= 160. ) //pi-
	bpart = kTRUE;

      else if( selid < 8  && pid == partid[selid] && charge == 1 ){ //p,d,t,3He,4He
	MassNumber = partA[selid];
	bpart = kTRUE;
      }
      else if( isys == 5 ) {
	MassNumber = partA[2];
	bpart = kTRUE;
      }

      //&&&----------------------------------------
      //if( theta > 45.*TMath::DegToRad() ) continue;
      //if( theta > 40.*TMath::DegToRad() && theta < 50.*TMath::DegToRad()) continue;
      
      //      if( phi < 0 ) continue;
      //      if( phi > 0 ) continue;

      //if( yaw < 0 ) continue; 
      //if( yaw > 0 ) continue; 

      //if( abs( phi ) >  30.*TMath::DegToRad() ) continue;
      //if( abs( phi ) < 150.*TMath::DegToRad() ) continue;
      //*********************************************

      hpid->Fill(aPart->GetRotatedMomentum().Mag()*aPart->GetCharge(), aPart->GetdEdx());

      //-----------------
      if( !bpart && isys != 5 ) continue; //default
      //-----------------
      hpidsel->Fill(aPart->GetRotatedMomentum().Mag()*aPart->GetCharge(), aPart->GetdEdx());

      bFill = kTRUE;

      y_norm = y_cm[isys+6];

      auto rapidl = aPart->GetRapidity();
      auto rapid  = aPart->GetRapiditycm();;	
      auto rapidn = rapid / y_norm;
      auto rpphi  = aPart->GetIndividualRPAngle();
      auto dphi   = aPart->GetAzmAngle_wrt_RP();
      auto dphi2  = aPart->GetAzmAngle2_wrt_RP();
      fmass  = aPart->GetMass();
      Double_t u_t0  = aPart->GetRotatedMomentum().Pt()/fmass/u_p;
      Double_t ou_t0 = aPart->GetMomentumAtTarget().Pt()/fmass/u_p;

      //      if( isys == 5 ) 
      //	dphi = TVector2::Phi_mpi_pi(aPart->GetRotatedMomentum().Phi() - RPPsi);

      hPsi     ->Fill(RPangle);
      hRPPsi   ->Fill(rpphi);
      hphi     ->Fill(phi);
      hRPPsipsi->Fill(RPPsi, rpphi);
      hyptacp  ->Fill(rapid/y_norm, pt);
      hypteff  ->Fill(rapid/y_norm, pt/1000.);
      hyuteff  ->Fill(rapid/y_norm, u_t0);
      hyutacp  ->Fill(rapid/y_norm, u_t0);

      hyawpitch->Fill(yaw,   pitch);
      hmass    ->Fill(aPart->GetRotatedMomentum().Mag(), bmass);



      UInt_t irapid1 = GetV1cmRapidityIndex(rapidn);
      UInt_t ipt1    = GetV1PtIndex(u_t0);


      UInt_t irapid2 = GetV2cmRapidityIndex(rapidn);
      UInt_t ipt2    = GetV2PtIndex(u_t0);

      hirap1 -> Fill( rapid/y_norm , (Double_t)irapid1 );
      hirap2 -> Fill( rapid/y_norm , (Double_t)irapid2 );


      //@@@@@
      //v1 -----------
      hutphi10[irapid1] -> Fill( subevt_phi );

      hdy1[irapid1]    -> Fill( rapid/y_norm );
      hdycos1[irapid1] -> Fill( cos(dphi) );

      hdydut1[irapid1][ipt1]    -> Fill( u_t0 );
      hdydutcos1[irapid1][ipt1] -> Fill( cos(dphi) );

      //      if( (selid >= 2 && u_t0 > 0) || selid < 2 ) {

      //v2 ------------
      hutphi20[irapid2] -> Fill( subevt_phi );

      hdy2[irapid2]     -> Fill( rapid/y_norm );
      hdycos2[irapid2]  -> Fill( cos(2.*dphi) );

      hdydut2[irapid2][ipt2] -> Fill( u_t0 );
      hdydutcos2[irapid2][ipt2] -> Fill( cos(2.*dphi) );

      ///$$$$$////
      if( (selid >= 2 && u_t0 > 0.) || selid < 2 ) {
      //      if( (selid >= 2 && u_t0 > 0.4) || selid < 2 ) {

	huy->Fill( rapid / y_norm, u_t0 );

	hdyv1[irapid1]               -> Fill( cos(dphi) );
	hdydut1cut[irapid1][ipt1]    -> Fill( u_t0 );
	hdydutcos1cut[irapid1][ipt1] -> Fill( cos(dphi) );

       	hdyv2[irapid2]               -> Fill( cos(2.*dphi) );
	hdydut2cut[irapid2][ipt2]    -> Fill( u_t0 );
	hdydutcos2cut[irapid2][ipt2] -> Fill( cos(2.*dphi) );
	
      }
    }
    
    if( bFill ){
      hdphi0_180->Fill( subevt_phi );
      if( subevt_phi > TMath::Pi()/2. )
	hdphi90_180->Fill( subevt_phi );
    }
  }

  // acceptance correction

  if( bAcc_corr ) {
    GetYUtCorrection(*hyutacp, *hyutacpcrr);
    DrawUt(hyutacpcrr,1);
    DrawUt(huy,1,"same");

    DrawUt(hyutacpcrr,2);
    DrawUt(huy,2,"same");
  }

  //--------------------------------------------------
  //--- Enf of event Loop
  //--------------------------------------------------
  LOG(INFO) << " End of Event Loop " << FairLogger::endl;


  //--------------------------------------------------

  TGraphErrors *gy_v1 = new TGraphErrors();
  gy_v1->SetName("gy_v1");
  TGraphErrors *gy_v2 = new TGraphErrors();
  gy_v2->SetName("gy_v2");

  TGraphErrors *gu_v1 = new TGraphErrors();
  gu_v1->SetName("gu_v1"); 
  gu_v1->SetTitle("u_{t0}; y_{cm}/y_{cm}; v1");
  TGraphErrors *gu_v2 = new TGraphErrors();
  gu_v2->SetName("gu_v2");
  gu_v2->SetTitle("u_{t0}; y_{cm}/y_{cm}; v2");

  TGraphErrors *gUt_v1[ybin1];
  TGraphErrors *gUt_v2[ybin2];

  for(UInt_t kn = 0; kn < ybin1 ; kn++){
    gUt_v1[kn] = new TGraphErrors();
    gUt_v1[kn]->SetName((TString)Form("gUt_v1%d",kn));
    TString sname = partname[selid]+" "+rangeLabel1[kn]+"; u_{t0}; v1";
    gUt_v1[kn]->SetTitle(sname); 
  }

  for(UInt_t kn = 0; kn < ybin2 ; kn++){
    gUt_v2[kn] = new TGraphErrors();
    gUt_v2[kn]->SetName((TString)Form("gUt_v2%d",kn));
    TString sname = partname[selid]+" "+rangeLabel2[kn]+"; u_{t0}; v2";
    gUt_v2[kn]->SetTitle(sname);
  }


  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++@@
  // This function returen resolution from subevent correlation with overall events.
  Double_t *rpresall = new Double_t[2];
  GetRPResolutionwChi(rpresall, hdphi0_180, hdphi90_180, 1.);
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //--- v1 ------------------------------
  UInt_t iyv = 0;
  UInt_t iyy = 0;
  for( UInt_t iy = 0; iy < ybin1; iy++ ) {

    Double_t yn    = (Double_t)hdy1[iy]->GetEntries();    
    Double_t ymean; Double_t ystdv;
    if( yn <= 0 ) continue;

    ymean = hdy1[iy]->GetMean() ;
    ystdv = hdy1[iy]->GetStdDev()/sqrt(yn);

    auto mcos  = hdyv1[iy]->GetMean();   // ut>0.4
    auto mcose = hdyv1[iy]->GetStdDev()/sqrt(yn);
    auto mcosee= GetError(mcos, rpresall[0], mcose, rpresall[1]);
      
    auto cosc = mcos/rpresall[0];
    gy_v1->SetPoint(iyy, ymean, cosc);
    gy_v1->SetPointError(iyy, ystdv, mcosee );
    iyy++;

    
    UInt_t iptv    = 0;
    //--- Acceptance correction
    Double_t v1u_sum = 0.;
    Double_t v1u_n   = 0.;
    Double_t v1u_ste = 0.;

    for( UInt_t ipt = 0; ipt < ptbin1; ipt++ ) {  // ptbin
      
      Double_t utn    = (Double_t)hdydut1[iy][ipt]->GetEntries();
      if( utn == 0 ) continue;

      Double_t utmean = hdydut1[iy][ipt]->GetMean();
      Double_t ute    = hdydut1[iy][ipt]->GetStdDev()/sqrt(utn);

      Double_t v1u   = hdydutcos1[iy][ipt]->GetMean();
      Double_t v1ue  = hdydutcos1[iy][ipt]->GetStdDev()/sqrt(utn);

      //      hdydutcos1[iy][ipt]->Write();
	
      Double_t v1ut  = v1u / rpresall[0];
      Double_t v1ute = GetError(v1u, rpresall[0], v1ue, rpresall[1]);
	
      if( !std::isnan(utmean) && !std::isnan(v1ut) ) {
	gUt_v1[iy]->SetPoint(iptv, utmean, v1ut );
	gUt_v1[iy]->SetPointError(iptv, ute, v1ute);
	iptv++;
      }


      utn   = hdydut1cut[iy][ipt]->GetEntries();

      if( utn > 0 ) {

	///--------------------
	// acceptance correction
	if( bAcc_corr ) {
	  Double_t ptmean = utmean * fmass * u_p /1000.;

	  if( !hyutacpcrr -> IsZombie() ) { 
	    //	    GetdNdydPt(ymean, ptmean);
	    //	    GetdNdydUt(ymean, utmean,  huy);
	    //	    GetdNdydUt(ymean, utmean,  hyutacpcrr);
	    //	    GetdNdydUt(1, iy, ipt, huy);
	    GetdNdydUt(1, iy, ipt, hyutacpcrr);

	    ///$$$$$////
	    hdndydut->Fill(yrange1nrm[iy], dpt1*ipt, acpcorr[0]);
	  }
	}
	//---------------------

	v1u   = hdydutcos1cut[iy][ipt]->GetMean();
	v1ue  = hdydutcos1cut[iy][ipt]->GetStdDev()/sqrt(utn);
	v1ut  = v1u / rpresall[0];
	v1ute = GetError(v1u, rpresall[0], v1ue, rpresall[1]);
	v1ut *= acpcorr[0];
	v1ute = pow(v1ute,2) /acpcorr[0];

	v1u_sum += v1ut;	
	v1u_ste += v1ute ;
	v1u_n   += acpcorr[0];
      }
    }

    if( yn > 0 && v1u_n > 0) {
      Double_t v1u_ave = v1u_sum / v1u_n;
      Double_t v1u_err = sqrt(v1u_ste);
      
      gu_v1->SetPoint( iyv, ymean, v1u_ave);
      gu_v1->SetPointError( iyv, ystdv, v1u_err);
      
      iyv++;
    }

    gUt_v1[iy]->Write();
  }
  

  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //--- v2

  GetRPResolutionwChi(rpresall, hdphi0_180, hdphi90_180, 2.);

  iyv = 0;
  iyy = 0;
  for( UInt_t iy = 0; iy < ybin2; iy++ ) {

    Double_t yn    = (Double_t)hdy2[iy]->GetEntries();
    Double_t ymean; Double_t ystdv;
    if( yn <= 0 ) continue;

    ymean = hdy2[iy]->GetMean();
    ystdv = hdy2[iy]->GetStdDev()/sqrt(yn);
      
    auto mcos  = hdyv2[iy]->GetMean();
    auto mcose = hdyv2[iy]->GetStdDev()/sqrt(yn);
    auto mcosee= GetError(mcos, rpresall[0], mcose, rpresall[1]);
    
    auto cosc = mcos/rpresall[0];

    gy_v2->SetPoint(iyy, ymean, cosc);
    gy_v2->SetPointError(iyy, ystdv, mcosee );
    iyy++;
  

    UInt_t   iptv    = 0;
    //--- Acceptance correction
    Double_t v2u_sum = 0.;
    Double_t v2u_n   = 0.;
    Double_t v2u_ste = 0.;
    
    for( UInt_t ipt = 0; ipt < ptbin2; ipt++ ){
    
      Double_t utn    = (Double_t)hdydut2[iy][ipt]->GetEntries();
      if( utn == 0 ) continue;

      Double_t utmean = hdydut2[iy][ipt]->GetMean();
      Double_t ute    = hdydut2[iy][ipt]->GetStdDev()/sqrt(utn);

      Double_t v2u    = hdydutcos2[iy][ipt]->GetMean();
      Double_t v2ue   = hdydutcos2[iy][ipt]->GetStdDev()/sqrt(utn);


      Double_t v2ut  = v2u / rpresall[0];
      Double_t v2ute = GetError(v2u, rpresall[0], v2ue, rpresall[1] );

      if( !std::isnan(utmean) && !std::isnan(ute) ) {
      	gUt_v2[iy]->SetPoint(iptv, utmean, v2ut);
      	gUt_v2[iy]->SetPointError(iptv, ute, v2ute);
      	iptv++;
      }


      utn  = hdydut2cut[iy][ipt]->GetEntries();

      ///$$$$$////
      if( utn > 0 ) {
	///--------------------
	// acceptance correction   
	if( bAcc_corr ) {
	  Double_t ptmean = utmean * fmass * u_p /1000.;

	  if( !hyutacpcrr -> IsZombie() ) {
	    //	    GetdNdydPt(ymean, ptmean);
	    //      GetdNdydUt(ymean, utmean, hyutacpcrr);
	    //	    GetdNdydUt(ymean, utmean, huy);
	    //	    GetdNdydUt(iy, ipt, huy);
	    //	    GetdNdydUt(2, iy, ipt, huy);
	    GetdNdydUt(2, iy, ipt, hyutacpcrr);
	  }
	}

	v2u   = hdydutcos2cut[iy][ipt] -> GetMean();
	v2ue  = hdydutcos2cut[iy][ipt] -> GetStdDev()/ sqrt(utn);          
	v2ut  = v2u / rpresall[0];
	v2ute = GetError(v2u, rpresall[0], v2ue, rpresall[1]) ;	
	v2ut *= acpcorr[0];
	//	v2ute = GetNError(v2ut, acpcorr[0],  v2ute, acpcorr[1]);
	v2ute = pow(v2ute,2)/acpcorr[0] ;

	v2u_sum += v2ut;
	v2u_ste += pow(v2ute,2);
	v2u_n   += acpcorr[0];
      }

      hdydutcos2[iy][ipt]->Write();
      hdydutcos2cut[iy][ipt]->Write();
    }

    if( v2u_n > 0 ) {
               
      Double_t v2u_ave = v2u_sum / v2u_n;
      Double_t v2u_err = sqrt(v2u_ste);
      
      gu_v2->SetPoint( iyv, ymean, v2u_ave);
      gu_v2->SetPointError( iyv, ystdv, v2u_err );
      
      iyv++;
    }
    gUt_v2[iy]->Write();
  }

  ///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  hdphi0_180->Write();
  hdphi90_180->Write();
  gu_v1->Write();
  gu_v2->Write();
  gy_v1->Write();
  gy_v2->Write();
  hmult->Write();
  hyawpitch->Write();
  hmass    ->Write();
  hyptacp  ->Write();
  hypteff  ->Write();
  hyuteff  ->Write();
  hyutacp  ->Write();
  hyutacpcrr  ->Write();
  hdndydut ->Write();
  hirap1   ->Write();
  hirap2   ->Write();
  huy      ->Write();
  hpsi     ->Write();
  hpsi1    ->Write();
  hpsindx  ->Write();
  hiphi    ->Write();
  hPsi     ->Write();
  hPsinc   ->Write();
  hRPPsi    ->Write();
  hRPPsipsi ->Write();
  hphi      ->Write();
  hpid      ->Write();
  hpidsel   ->Write();

  if( hAcpCorr ) {
    hAcpYUtCorr->Write();
    hAcpCorr   ->Write();
  }
  //  hdyucos2 ->Write();
  //--------------------------------------------------
  //--- Plotting
  //--------------------------------------------------
  id=1;
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),2000,500);
  cc->Divide(10);
  for( UInt_t iy = 0; iy < ybin1-1; iy++ ) {
    cc->cd(id); id++;
    gUt_v1[iy] -> SetMarkerStyle(20);
    gUt_v1[iy] -> SetMarkerColor(4);
    gUt_v1[iy] -> Draw("ALP");
  }


  id=1;
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),2000,500);
  cc->Divide(7);
  for( UInt_t iy = 0; iy < ybin2-1; iy++ ) {
    cc->cd(id); id++;
    gUt_v2[iy] -> SetMarkerStyle(20);
    gUt_v2[iy] -> SetMarkerColor(4);
    gUt_v2[iy] -> Draw("ALP");
  }

   
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  gu_v1->SetMarkerStyle(20);
  gu_v1->SetMarkerColor(2);
  gu_v1->SetLineColor(2);
  gu_v1->Draw("ALP");
  if( isys == 5 ) {
    FlowFunction();
    fv1y->Draw("same");
  }
  gy_v1->SetLineColor(4);
  gy_v1->Draw("same");

  auto LineV = new TLine(0.,gu_v1->GetYaxis()->GetXmin(), 0., gu_v1->GetYaxis()->GetXmax());
  auto LineH = new TLine(gu_v1->GetXaxis()->GetXmin(),    0., gu_v1->GetXaxis()->GetXmax(), 0.);
  LineV->SetLineStyle(3);
  LineH->SetLineStyle(3);
  LineV->Draw();
  LineH->Draw();

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  cc->Divide(2,2);
  cc->cd(1);
  huy->Draw("colz");
  cc->cd(2);
  hdndydut->Draw("colz");

  if( bAcc_corr ) {
    cc->cd(3);
    hAcpCorr->Draw("colz");
    cc->cd(4);
    hyutacpcrr->Draw("colz");
  }
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  cc->Divide(2,2);
  cc->cd(1);
  hPsinc->Draw();
  cc->cd(2);
  hPsi->Draw();
  cc->cd(3);
  hRPPsipsi->Draw("colz");
  cc->cd(4);
  hphi->Draw();

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  gu_v2->SetMarkerStyle(20);
  gu_v2->SetMarkerColor(2);
  gu_v2->SetLineColor(2);
  gu_v2->Draw("ALP");

  if(isys == 5) {
    fv2y->Draw("same");
    gu_v2->Print();
  }

  gy_v2->SetLineColor(4);
  gy_v2->Draw("same");
  

  gSystem->cd("..");
}



void PtDistribution(UInt_t selid)
{
  TH1D *hut[ybin1][2];
  TString rangeLabel1[ybin1];
  Double_t u_p = 0.355151 * 1.06974; //p+p(268.9MeV/u)
  
  for( UInt_t iy = 0; iy < ybin1; iy++ ) {
    rangeLabel1[iy] = Form("%5.2f <= y < %5.2f"      ,yrange1nrm[iy],yrange1nrm[iy+1]);
    hut[iy][0]    = new TH1D(Form("hut0_%d",iy)      ,rangeLabel1[iy]+"_yaw>0;U_{t};",100,0., 2.5);
    hut[iy][1]    = new TH1D(Form("hut1_%d",iy)      ,rangeLabel1[iy]+"_yaw<0;U_{t};",100,0., 2.5);
  }

  Long64_t nevt = SetBranch();

  //--------------------------------------------------
  //--- Event Loop
  //--------------------------------------------------
  for(Long64_t i = 0; i < nevt; i++){

    rChain->GetEntry(i);

    TIter next(aArray);
    STParticle *aPart = NULL;

    //--------------------------------------------------
    //----- Main loop 
    //--------------------------------------------------
    while( (aPart = (STParticle*)next()) ) {

      auto pid   = aPart->GetPIDTight();  // = GetPID()
      if( pid != partid[selid] ) continue;

      auto phi   = aPart->GetRotatedMomentum().Phi();
      auto rapid = aPart->GetRapiditycm();;	

      Double_t u_t0  = aPart->GetRotatedMomentum().Pt()/aPart->GetMass()/u_p;
      UInt_t irapid1 = GetV1RapidityIndex(rapid);

      if( abs(phi) <= 30.*TMath::DegToRad() )
	hut[irapid1][0] -> Fill( u_t0 );
      else if( abs(phi) >= 150.*TMath::DegToRad() )
	hut[irapid1][1] -> Fill( u_t0 );

    }
  }

  id = 1;
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1200,600);
  cc->Divide(ybin1/2, 2);
  for( UInt_t i = 0; i < ybin1; i++ ){
    cc->cd(id); id++;
    hut[i][0]->SetNormFactor(100/2.5);
    hut[i][1]->SetNormFactor(100/2.5);
    hut[i][0]->SetLineColor(2);
    hut[i][1]->SetLineColor(4);
    hut[i][0]->Draw();
    hut[i][1]->Draw("same");
  }

}


///######################################################################
/// sub functions


//**************************************************
Bool_t LoadAcceptanceCorrection(UInt_t selid)
{
  //  TString fname = "LCPSpectra_b0.15.root";
  TString fname = "UnfoldedLCPSpectra.slim.root";

  if( !gSystem->FindFile("data",fname) ){
    LOG(ERROR) << fname << " is not found " << FairLogger::endl; 
    return kFALSE;
  }

  TFile *fopen = TFile::Open(fname);
  
  //  TString hname = "h2PtY_"+ sysName + "_" +Partname[selid];
  TString hname = "h2UtYCorr_108Sn_" +Partname[selid]+"_mbin0_iter1";
  hAcpCorr = (TH2D*)fopen->Get(hname);
  
  if( !hAcpCorr ) {
    LOG(ERROR) << hname << " is not found. " << FairLogger::endl;
    fopen->Close();
    return kFALSE;
  }

  hAcpCorr -> SetDirectory( 0 );

  hname = "h2UtYEff_108Sn_" +Partname[selid]+"_mbin0_iter1";
  hAcpYUtCorr = (TH2D*)fopen->Get(hname);

  if( !hAcpYUtCorr ) {
    LOG(ERROR) << hname << " is not found " << FairLogger::endl;
    fopen->Close();
    return kFALSE;
  }
  hAcpYUtCorr -> SetDirectory(0);

  fopen->Close();
  return kTRUE;
}
void DrawUt(TH2D *h2, UInt_t sel, TString opt)
{
  TString hname = h2->GetName();

  if( opt != "same" ) {
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
    h2->Draw("colz");
  }
    
  Double_t *yrange = &yrange1nrm[0];
  UInt_t ybin = ybin1;

  if( sel == 2 ) {
    yrange = &yrange2nrm[0];
    ybin = ybin2;
  }


  TH1D *hacp;
  if( opt != "same" ) {
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1200,600);
    cc->Divide(4, ybin/4);
  }

  for( UInt_t i = 0; i < ybin; i++ ) {
    auto firstbin = h2->GetXaxis()->FindBin(*(yrange+i));
    auto lastbin  = h2->GetXaxis()->FindBin(*(yrange+i+1));

    //    cout << " ybin " << i << " " << *(yrange+i) << " y " << firstbin << " ~ " << *(yrange+
    hacp = (TH1D*) h2->ProjectionY(hname+Form("_%d_%d",sel,i), firstbin, lastbin);
    hacp -> SetTitle(Form("%3.2f < y < %3.2f" , *(yrange+i), *(yrange+i+1))); 
    if( hacp->GetEntries() > 0)
      hacp -> SetNormFactor(50);
      //      hacp -> Scale(1./hacp->GetEntries() );
    cc->cd(i+1);

    if( opt == "same" )
      hacp->SetLineColor(2);

    hacp -> Draw(opt);
  }
}

void GetdNdydPt( Double_t ycm, Double_t pt)
{
  auto ybin = hAcpCorr->GetXaxis()->FindBin(ycm);
  auto xbin = hAcpCorr->GetYaxis()->FindBin(pt);

  acpcorr[0] = hAcpCorr->GetBinContent(ybin, xbin);
  acpcorr[1] = hAcpCorr->GetBinError(ybin, xbin);  
}

void GetYUtCorrection(TH2D &h2, TH2D &h2c)
{
  auto xnbins = h2.GetXaxis()->GetNbins();
  auto ynbins = h2.GetYaxis()->GetNbins();


  for( UInt_t i = 1; i <= xnbins; i++ ) {
      auto xcenter = h2.GetXaxis()->GetBinCenter(i); 
      auto xbin = hAcpYUtCorr->GetXaxis()->FindBin(xcenter);

      for( UInt_t j = 1; j <= ynbins; j++ ){
	auto ycenter = h2.GetYaxis()->GetBinCenter(j);
	auto ybin = hAcpYUtCorr->GetYaxis()->FindBin(ycenter);

	auto crf = hAcpYUtCorr->GetBinContent(xbin, ybin);
	auto cre = hAcpYUtCorr->GetBinError(xbin, ybin);

	auto binCont = h2.GetBinContent(i,j);
	
	if( crf != 0 )
	  h2c.Fill(xcenter, ycenter, binCont/crf);
      }
  }
}

void GetdNdydUt( Int_t ih, Int_t ixbin, Int_t iybin, TH2D *h2c)
{
  acpcorr[0] = 0.;
  acpcorr[1] = 0.;
  
  Double_t x[2] = {0.,0.};
  Double_t y[2] = {0.,0.};

  if( ih == 1 ) {
    if( ixbin < 0 || ixbin > ybin1 ) return;
    if( iybin < 0 || iybin > ptbin1 ) return;

    x[0] = yrange1nrm[ixbin];
    x[1] = yrange1nrm[ixbin+1];
    y[0] = dpt1*iybin;
    y[1] = y[0] + dpt1;

  }
  else if( ih == 2 ) {
    if( ixbin < 0 || ixbin > ybin2 ) return;
    if( iybin < 0 || iybin > ptbin2 ) return;

    x[0] = yrange2nrm[ixbin];
    x[1] = yrange2nrm[ixbin+1];
    y[0] = dpt2*iybin;
    y[1] = y[0] + dpt2;
  }
  else
    return;
    
  ///$$$$$////
  auto xbin_frst = h2c->GetXaxis()->FindBin(x[0]);
  auto xbin_last = h2c->GetXaxis()->FindBin(x[1])-1;
  auto ybin_frst = h2c->GetYaxis()->FindBin(y[0]);
  auto ybin_last = h2c->GetYaxis()->FindBin(y[1])-1;


  Double_t content = 0.;
  Double_t error = 0.;

  content = h2c->IntegralAndError(xbin_frst, xbin_last, ybin_frst, ybin_last, error, "");


  // cout << ih << " :: irapid = " << ixbin << " ipt = " << iybin << endl;	
  // cout << " rapid " << xbin_frst << " : " << x[0] <<" - " << xbin_last << " "<< x[1] 
  //      << " pt " << ybin_frst << " : " << y[0] << " - " << ybin_last << " " << y[1] << endl;
  // cout << " --> " << content << " +- " << error << endl;


  //test
  // TF1 *utslp = new TF1("utslp","exp(-x/1.)+1000.",0.,2.);
  // content = utslp->Eval(ut);
  // error   = 0.;

  acpcorr[0] = content;
  acpcorr[1] = error/sqrt(content);
  //  acpcorr[1] = sqrt(content);
}

void GetdNdydUt( Double_t ycm, Double_t ut, TH2D *h2c)
{
  Int_t wbin = 1;

  auto ybin     = h2c->GetXaxis()->FindBin(ycm);
  auto ybin_low = h2c->GetXaxis()->FindBin(ycm)-wbin;
  auto ybin_hig = h2c->GetXaxis()->FindBin(ycm)+wbin;

  wbin = 2;
  auto xbin     = h2c->GetYaxis()->FindBin(ut);
  auto xbin_low = h2c->GetYaxis()->FindBin(ut)-wbin;
  auto xbin_hig = h2c->GetYaxis()->FindBin(ut)+wbin;
  
  //  cout << " ycm " << ycm << " ut " << ut << endl;
  Double_t content = 0.;
  Double_t error = 0.;
  for( Int_t iy = ybin_low; iy <= ybin_hig; iy++ ) {
    for( Int_t ix = xbin_low; ix <= xbin_hig; ix++ ) {
      auto num = h2c->GetBinContent(iy, ix);
      content += num;

      // cout << " y = " << iy << " x = " << ix <<  " " << h2c->GetBinContent(iy, ix) 
      //  	   << " +- " << h2c->GetBinError(iy, ix) << " : " << sqrt(error)
      //  	   << endl;
    }
  }

  //test
  // TF1 *utslp = new TF1("utslp","exp(-x/1.)+1000.",0.,2.);
  // content = utslp->Eval(ut);
  // error   = 0.;

  acpcorr[0] = content;
  acpcorr[1] = sqrt(content);
}


UInt_t GetV1RapidityIndex(Double_t y)
{
  UInt_t irapid1 = ybin1-1;
  for( UInt_t i = 1; i < ybin1; i++){
    if( y < yrange1nrm[i]){
      irapid1 = i-1;
      break;
    }
  }
  return irapid1;
}
UInt_t GetV1cmRapidityIndex(Double_t y)
{
  UInt_t irapid1 = ybin1-1;
  for( UInt_t i = 1; i < ybin1; i++){
    if( y < yrange1[i]/y_norm ){
      irapid1 = i-1;
      break;
    }
  }
  return irapid1;
}

UInt_t GetV2RapidityIndex(Double_t y)
{
  UInt_t irapid2 = ybin2-1;
  for( UInt_t i = 1; i < ybin2; i++){
    if( y < yrange2nrm[i]){
      irapid2 = i-1;
      break;
    }
  }
  return irapid2;
}
UInt_t GetV2cmRapidityIndex(Double_t y)
{
  UInt_t irapid2 = ybin2-1;
  for( UInt_t i = 1; i < ybin2; i++){
    if( y < yrange2[i]/y_norm){
      irapid2 = i-1;
      break;
    }
  }
  return irapid2;
}

UInt_t GetV1PtIndex(Double_t val)
{
  return UInt_t(val/dpt1);



  UInt_t ipt = ptbin1 - 1;
  for(UInt_t i = 0; i < ptbin1; i++){
    if( val < dpt1*(i+1)) {
      ipt = i;
      break;
    }
  }
  return ipt;
}

UInt_t GetV2PtIndex(Double_t val)
{
  return UInt_t(val/dpt2);


  UInt_t ipt = ptbin2 - 1;
  for(UInt_t i = 0; i < ptbin2; i++){
    if( val < dpt2*(i+1)) {
      ipt = i;
      break;
    }
  }
  return ipt;
}

//--------------------------------------------------


TString SetupOutputFile(TString fopt)
{
  //--------------------------------------------------
  //----- SETUP OUTPUT FILE --------------------------
  //--------------------------------------------------
  gSystem->cd("data");
  TString oVerSub= gSystem -> Getenv("OSUBV");
  TString fName = fopt + oVerSub;

  if( oVerSub == "" ) {
    UInt_t kVer = 0; 
    TString checkName;

    while( kTRUE ) {
      checkName = Form(fName + "%d.root", kVer);
      if( !gSystem->FindFile(".",checkName) ) {
	break;
      }
      else
	kVer++;
    }
    fName = Form(fName + "%d.root", kVer);
  }
  else {
    TString checkName = fName + ".root";;
    if( !gROOT->IsBatch() && gSystem->FindFile(".",checkName) ) {
      LOG(INFO) << checkName << " is existing. Do you recreate? (y/n)" << FairLogger::endl;
      TString sAns;
      std::cin >> sAns;
      if(sAns == "N" || sAns == "n") {
	LOG(INFO) << " Retry" << FairLogger::endl;
	exit(0);
      }
    }
    fName += ".root";
  }

  LOG(INFO) << "File " << fName << " will be created. " << FairLogger::endl;  

  
  return fName;
}
//--------------------------------------------------

Bool_t SetupEffCorrection()
{
  fefffile = TFile::Open("UnfoldedLCPSpectra.slim.root");

  if( fefffile != NULL ) 
    return kTRUE;
  else
    return kFALSE;
}


//---------------------------------------------------

void DetectorBias()
{
  UInt_t m = 0;

  Long64_t nEntry =  SetBranch();

  Double_t  cosphi = 0.;
  Double_t  sinphi = 0.;
  Double_t  cos2phi= 0.;
  Double_t  sin2phi= 0.;

  Double_t count = 0.;
  for(Int_t i = 0; i < nEntry; i++){

    rChain->GetEntry(i);
    STFlowInfo *aflow = (STFlowInfo*)aFlowArray->At(0);

    
    if(aflow->mtrack1 > Ucent || aflow->mtrack1 <= Lcent ) continue;
    count++;

    cosphi += cos( (aflow->unitP_fc).Phi());      
    sinphi += sin( (aflow->unitP_fc).Phi());

    cos2phi += cos(2.*(aflow->unitP_fc).Phi());      
    sin2phi += sin(2.*(aflow->unitP_fc).Phi());

  }
  
  LOG(INFO) << " <cos Phi> " << cosphi/count 
	    << " <sin phi> " << sinphi/count << FairLogger::endl;

  LOG(INFO) << " <cos 2Phi> " << cos2phi/count 
	    << " <sin 2phi> " << sin2phi/count << FairLogger::endl;

}
//---------------------------------------------------


//---------------------------------------------------
void PlotSubEvent(Double_t ml=30, Double_t mu=80)   
{

  Double_t mlt[] = {ml, mu};
  //    Double_t mlt[2] = {0., 8.};
  //Double_t mlt[2] = {8., 16.};
  // Double_t mlt[2] = {16., 24.};
  // Double_t mlt[2] = {24., 32.};
  // Double_t mlt[2] = {32., 40.};
  // Double_t mlt[2] = {40., 100.};

  TCut mcrot = Form("aFlow->mtrack4>%f&&aFlow->mtrack4<%f",mlt[0]*2.,mlt[1]*2.);
  TCut mc1r  = Form("mtrack_1>%f&&mtrack_1<%f"  ,mlt[0],mlt[1]);
  TCut mc2r  = Form("mtrack_2>%f&&mtrack_2<%f"  ,mlt[0],mlt[1]);

  TString sname;
  
  sname = mcrot.GetTitle();
  auto *hrotx = new TH1D("hrotx","All   "+sname+ ";Qx" ,100,-12.,12.);
  auto *h1rx  = new TH1D("h1rx", "sub_1 "+sname+ ";Qx" ,100,-12.,12.);
  auto *h2rx  = new TH1D("h2rx", "sub_2 "+sname+ ";Qx" ,100,-12.,12.);
  auto *hroty = new TH1D("hroty","All   "+sname+ ";Qy" ,100,-12.,12.);
  auto *h1ry  = new TH1D("h1ry", "sub_1 "+sname+ ";Qy" ,100,-12.,12.);
  auto *h2ry  = new TH1D("h2ry", "sub_2 "+sname+ ";Qy" ,100,-12.,12.);
 

  rChain->Project("hrotx","unitP2_fc.X()",mcrot);
  rChain->Project("h1rx" ,"unitP_1fc.X()"  ,mc1r);
  rChain->Project("h2rx" ,"unitP_2fc.X()"  ,mc2r);
					      
  rChain->Project("hroty","unitP2_fc.Y()",mcrot);
  rChain->Project("h1ry" ,"unitP_1fc.Y()"  ,mc1r);
  rChain->Project("h2ry" ,"unitP_2fc.Y()"  ,mc2r);


  //----- Drawing       
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);

  hrotx->SetLineColor(2);
  hrotx->SetNormFactor(1);
  h1rx ->SetLineColor(4);
  h1rx ->SetNormFactor(1);
  h2rx ->SetLineColor(6);
  h2rx ->SetNormFactor(1);

  h1rx->Draw();
  h2rx->Draw("same");
  hrotx->Draw("same");

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);

  hroty->SetLineColor(2);
  hroty->SetNormFactor(1);
  h1ry ->SetLineColor(4);
  h1ry ->SetNormFactor(1);
  h2ry ->SetLineColor(6);
  h2ry ->SetNormFactor(1);

  h1ry->Draw();
  h2ry->Draw("same");
  hroty->Draw("same");
}

///--------------------------------------------------
///######################################################################


