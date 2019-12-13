#include "DoRPRes.C" //#include "openRunAna.C" #include "DoFlow.h" #include "SimFunction.C"

//drawing

UInt_t selReactionPlanef = 10000;

TFile* fefffile;

// functions
void     PlotPtDependence(UInt_t selid);
TString  SetupOutputFile(TString fopt);
Bool_t   SetupEffCorrection();

Bool_t   SetPsiRPResolution();
Bool_t   LoadRPResolution();
Double_t *GetPsiRPResolution(UInt_t ival);
Double_t *GetRPResolution(UInt_t vn, UInt_t mult);
UInt_t   GetV1RapidityIndex(Double_t y);
UInt_t   GetV2RapidityIndex(Double_t y);
UInt_t   GetV1cmRapidityIndex(Double_t y);
UInt_t   GetV2cmRapidityIndex(Double_t y);
UInt_t   GetV1PtIndex(Double_t val);
UInt_t   GetV2PtIndex(Double_t val);

TString   GetPsiRPLoadFileName();


//-------------------//
void DoFlow_adv(Int_t isel = 0) 
{
  gROOT->Reset();

  openRunAna();

  if(rChain != NULL)     
    std::cout << " System " << isys << "  -> " << sysName << std::endl; 

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

  std::cout << " Multiplicity :: " << Lcent << " to " << Ucent << std::endl;

  //--------------------------------------------------
  TString rpBase = gSystem -> Getenv("RPBS");
  RPBase = rpBase != "" ? atoi(rpBase): 0;

  std::cout << " Reaction plane base is " << RPBase << std::endl;

  //==================================================
  if( isel > -1 )
    std::cout << " Output Version v" << sVer << "." << gSystem->Getenv("OUTVER") << std::endl;


  if( isel > 0 ) 
    PlotPtDependence((UInt_t)isel);  
  else if( isel == -1 )
    PlotPtDependence(0);  
}


//pdt
void PlotPtDependence(UInt_t selid = 2)       //%% Executable :
{

  TDatime beginTime;
  TDatime dtime;

  gStyle->SetOptStat(0);

  if(!LoadRPResolution()) {
    std::cout << " RP resolution is not found " << std::endl;
    //    return;
  }

  if(!SetPsiRPResolution()) {
    std::cout << " RP resolution is not found " << std::endl;
    exit(0);
  }
   
  std::cout << " RP resolution is ready " << std::endl;

  std::cout << "PlotPtDependence(" << selid << ")" << std::endl;
  dpt1 = pt_max/(Double_t)ptbin1;
  dpt2 = pt_max/(Double_t)ptbin2;

  // auto cutfile = new TFile("db/RegionCut.root");
  // TCutG *goodThetaPhi = (TCutG*)cutfile->Get("goodThetaPhi");
  // cutfile->Close();    


  TString fHeader = "advYPt_"+ sysName + "_" + partname[selid]+".v"+sVer+".";
  auto fName = SetupOutputFile( fHeader );
  
  Bool_t bEffCorr = kFALSE;
  if( selid == 2 || selid == 4 )
    auto bEffCorr = SetupEffCorrection();

  if( bEffCorr )
    std::cout << "correctedPt.phi60.nomultcut.root is opened." << std::endl;

  auto GraphSave  = new TFile(fName,"recreate");

  std::cout << " Rapidity binning " << ybin1 << std::endl;

  //------------------------------
  std::cout << " v1 BIN y " << ybin1 << " pt " << ptbin1 << " mult " << mbin << std::endl;
  std::cout << " v2 BIN y " << ybin2 << " pt " << ptbin2 << " mult " << mbin << std::endl;


  auto hyawpitch  = new TH2D("hyawpitch","; Yaw angle; Pitch angle",200,-1.5,1.5,200,-1.5,1.5);
  auto hmass      = new TH2D("hmass",   ";P/Q; Mass [MeV/c^2]"     ,200,  0.,2500., 200, 0.,7000);
  TString hlabel  = (TString)Form("mtrack1 %2d to %2d ; Y_{cm}; Pt [MeV/c]",Lcent,Ucent);
  auto hyptacp    = new TH2D("hyptacp", hlabel ,200, -1., 1.4, 200, 0., 1100);
  auto huy        = new TH2D("huy","; y_{cm}/y_{beam}; u_{t0} ", 200,-1., 1.8, 200,0.,1.8);
  auto hpsi       = new TH2D("hpsi",";#Psi;",100, -3.14,3.14, 200,0.,1.);
  auto hpsi1      = new TH1D("hpsi1",";#Psi",200,0.,1.);
  auto hpsindx    = new TH2D("hpsindx",";#Psi ;index ",100,-TMath::Pi(),TMath::Pi(),20,0.,20.);
  auto hiphi      = new TH1I("hiphi","hiphi",13,0,13);
  auto hPsi       = new TH1D("hPsi" ,";RP #Psi"  ,100,-TMath::Pi(),TMath::Pi());
  auto hRPPsi     = new TH1D("hRPPsi",";Indiv. #Psi"  ,100,-TMath::Pi(),TMath::Pi());
  auto hRPPsipsi  = new TH2D("hRPPsipsi",";RP #Psi; Indiv #Psi"  ,100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi());
  auto hphi       = new TH1D("hphi",     ";#phi"  ,100,-TMath::Pi(),TMath::Pi());

  TH1I *hmult = new TH1I("hmult",";Multiplicity",100,0,100);
  TH1F *hdydptcos1[ybin1][ptbin1][mbin];
  TH1F *hdydptcos2[ybin1][ptbin1][mbin];
  TH1D *hdydpt1[ybin1][ptbin1];
  TH1D *hdy1[ybin1];
  TH1D *hdympt1[ybin1];
  TH1D *hdydpt2[ybin1][ptbin1];
  TH1D *hdy2[ybin1];

  TH1D *hdydut1[ybin1][ptbin1];
  TH1D *hdydut2[ybin1][ptbin1];
  TH1D *hdyucos1[ybin1][psibin];
  TH1D *hdyucos2[ybin2][psibin];
  TH1D *hdydutcos1[ybin1][ptbin1][psibin];
  TH1D *hdydutcos2[ybin2][ptbin2][psibin];
  TH1D *hdydutdphi1[ybin1][ptbin1][psibin];
  TH1D *hutphi10[ybin1];
  TH1D *hutphi20[ybin2];
  TH1D *hutphi190[ybin1];
  TH1D *hutphi290[ybin2];
  TH1D *hut[ybin1][2];

  TH1D *hphi0_180  = new TH1D("hphi0_180" , "0to180",100,0.,3.2);
  TH1D *hphi90_180 = new TH1D("hphi90_180","90to180",100,0.,3.2);
  TH1D *h2phi  = new TH1D("h2phi","2x(#Delta #phi) at mid-rapidity",100,-1.*TMath::Pi(), TMath::Pi());

  TH2D *hirap1 = new TH2D("hirap1",";Rapidity; irap1",100,-1.0,1.5, 12,0,12);  
  TH2D *hirap2 = new TH2D("hirap2",";Rapidity; irap1",100,-1.0,1.5, 12,0,12);  

  //-----  booking
  TString rangeLabel1[ybin1];
  TString rangeLabel2[ybin2];

  for( UInt_t iy = 0; iy < ybin1; iy++ ) {
    rangeLabel1[iy] = Form("%5.2f <= y < %5.2f"      ,yrange1[iy]/y_cm[isys],yrange1[iy+1]/y_cm[isys]);
    hdy1[iy]      = new TH1D(Form("hdy1_y%d",iy)     ,rangeLabel1[iy]+";Rapidity", 500,-2.5,2.5);
    hdympt1[iy]   = new TH1D(Form("hdympt1%d",iy)    ,rangeLabel1[iy]+";Rapidity; <Pt>/A", 500,-1000.,1000.);
    hutphi10[iy]  = new TH1D(Form("hutphi10_%d",iy)  ,rangeLabel1[iy]+";#Delta #Psi",100,0.,3.2);
    hutphi190[iy] = new TH1D(Form("hutphi190_%d",iy) ,rangeLabel1[iy]+";#Delta #Psi",100,0.,3.2);
    hut[iy][0]    = new TH1D(Form("hut0_%d",iy)      ,rangeLabel1[iy]+"_yaw>0;U_{t};",100,0., 2.5);
    hut[iy][1]    = new TH1D(Form("hut1_%d",iy)      ,rangeLabel1[iy]+"_yaw<0;U_{t};",100,0., 2.5);

    for( UInt_t ips = 0; ips < psibin; ips++) 
      hdyucos1[iy][ips]  = new TH1D(Form("hdyucos1%dp%d",iy,ips)   ,rangeLabel1[iy]+"; cos(#phi - #Psi)"  , 100, -1., 1.);

    for( UInt_t ipt = 0; ipt < ptbin1; ipt++ ) {
      hdydpt1[iy][ipt] = new TH1D(Form("hdydpt1_y%dpt%d",iy,ipt),rangeLabel1[iy]+";Pt"    , 500,0., 1300);
      hdydut1[iy][ipt] = new TH1D(Form("hdydut1_y%dpt%d",iy,ipt),rangeLabel1[iy]+";u_{t0}", 500,0., 2.5);

      for( UInt_t ips = 0; ips < psibin; ips++) {
	hdydutcos1[iy][ipt][ips] = new TH1D(Form("hdydutcos1_y%dut%dp%d",iy,ipt,ips),rangeLabel1[iy]+"; cos(#phi - #Psi)" , 100, -1., 1.);
	hdydutdphi1[iy][ipt][ips] = new TH1D(Form("hdydutdphi1_y%dut%dp%d",iy,ipt,ips),rangeLabel1[iy]+"; #delta #phi",100,-3.15,3.15);
      }

      for( UInt_t im = 0; im < mbin; im++ ) {
	hdydptcos1[iy][ipt][im] = new TH1F(Form("hdydptcos1_y%dpt%dm%d",iy,ipt,im),"; cos(#phi - #Psi)" , 100., -1., 1.);
      }
    }
  }

  for( UInt_t iy = 0; iy < ybin2; iy++ ) {
    rangeLabel2[iy] = Form("%5.2f <= y < %5.2f",yrange2[iy]/y_cm[isys],yrange2[iy+1]/y_cm[isys]);
    hdy2[iy] = new TH1D(Form("hdy2_y%d",iy),rangeLabel2[iy]+";Rapidity", 500,-2.5,2.5);
    hutphi20[iy]  = new TH1D(Form("hutphi20_%d",iy),rangeLabel2[iy]+";#Delta #Psi",100,0.,3.2);
    hutphi290[iy]= new TH1D(Form("hutphi290_%d",iy),rangeLabel2[iy]+";#Delta #Psi",100,0.,3.2);

    for( UInt_t ips = 0; ips < psibin; ips++) 
      hdyucos2[iy][ips]  = new TH1D(Form("hdyucos2%dp%d",iy,ips),rangeLabel2[iy]+"; cos2(#phi - #Psi)" , 100, -1., 1.);

    for( UInt_t ipt = 0; ipt < ptbin2; ipt++ ) {
      hdydpt2[iy][ipt] = new TH1D(Form("hdydpt2_y%dpt%d",iy,ipt),rangeLabel2[iy]+"; Pt", 500,0., 1300.);
      hdydut2[iy][ipt] = new TH1D(Form("hdydut2_y%dpt%d",iy,ipt),rangeLabel2[iy]+"; Pt", 500,0., 2.5);

      for( UInt_t ips = 0; ips < psibin; ips++) 
	hdydutcos2[iy][ipt][ips]  = new TH1D(Form("hdydutcos2_y%dut%dp%d",iy,ipt,ips),rangeLabel2[iy]+"; cos(#phi - #Psi)" , 100, -1., 1.);
      

      for( UInt_t im = 0; im < mbin; im++ ) {
	hdydptcos2[iy][ipt][im] = new TH1F(Form("hdydptcos2_y%dpt%dm%d",iy,ipt,im),"; cos2(#phi - #Psi)", 100.,-1.,1.);
      }
    }
  }

  //------------------------------
  Double_t u_p = 0.355151 * 1.06974; //p+p(268.9MeV/u)

  Double_t *rppsires  = new Double_t[6];  // Psi dependent correction
  Double_t *rpres     = new Double_t[4];  // 
  Double_t *rpresy    = new Double_t[4];  // rapidity dependent correction

  Int_t nevt = SetBranch();
  
  //--------------------------------------------------
  //--- Event Loop
  //--------------------------------------------------
  Int_t nevt_begin = 0;
  //  Int_t maxevt = nevt_begin+1e+7; 
  Int_t maxevt = nevt_begin+2000000; 

  nevt = maxevt < nevt ? maxevt : nevt;
  cout << " NOTICE !!!! " << nevt_begin << " to " << nevt << endl;

  for(Int_t i = nevt_begin; i < nevt; i++){

    rChain->GetEntry(i);
    STFlowInfo *aflow = (STFlowInfo*)aFlowArray->At(0);

    /// for reaction plane resolution
    Bool_t bFill = kFALSE;
    Bool_t bRes  = kFALSE;

    beginTime.Copy(dtime);
    if(i%(UInt_t)(nevt/50) == 0) {
      dtime.Set();
      Int_t ptime = dtime.Get() - beginTime.Get();

      std::cout << "Processing .... " 
		<< setw(4) << Int_t(((Double_t)i/(Double_t)nevt)*100.) << " % = "
		<< setw(8) << i << "/"<< nevt
		<< "--->"
		<< dtime.AsString() << " ---- "
		<< std::endl;
    }

    /// centrality selection
    if(aflow->mtrack2 > Ucent || aflow->mtrack2 <= Lcent || aflow->mtrack4 < 6) continue;

    hmult->Fill( aflow->mtrack2 );

    bRes = kTRUE; //@1

    Double_t subevt_phi = abs(TVector2::Phi_mpi_pi((aflow->unitP_1fc).Phi()-
						   (aflow->unitP_2fc).Phi()));
      
    ///    auto RPangle = GetRPBaseAngle(aflow);

    auto RPangle = aflow->unitP_fc.Phi();

    UInt_t ipsi = GetPsiRPIndex( TVector2::Phi_mpi_pi(RPangle) );

    if( ipsi > npsi ) continue;
    
    TIter next(aArray);
    STParticle *aPart = NULL;

    //--------------------------------------------------
    //----- Main loop 
    //--------------------------------------------------
    UInt_t mtk = 0;
    while( (aPart = (STParticle*)next()) ) {

      mtk++;
      if( isys == 4 && mtk > (UInt_t)aflow->mtrack4*0.6 ) break; 

      // ---- track quality selection ---
      if( aPart->GetGoodTrackFlag() != 11 ) continue;
      //------------------------------

      auto yaw   = aPart->GetYawAngle();
      auto pitch = aPart->GetPitchAngle();
      auto pt    = aPart->GetRotatedMomentum().Pt();
      auto mom   = aPart->GetRotatedMomentum().Mag();
      auto bmass = aPart->GetBBMass();
      auto phi   = aPart->GetRotatedMomentum().Phi();
      auto theta = aPart->GetRotatedMomentum().Theta();
      auto charge= aPart->GetCharge();

      // auto pid   = aPart->GetPIDTight();  // = GetPID()
      auto pid   = aPart->GetPIDLoose();  
      //auto pid   = aPart->GetPIDNorm();  

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
      else if( selid < 8 ) { //pi- pi+ p d t 3He 4He
	if(pid == partid[selid] ){
	  MassNumber = partA[selid];
	    bpart = kTRUE;

	  if(selid == 0 && charge == 1 )
	    bpart = kFALSE;
	  else if( selid == 1 && charge == -1 )
	    bpart = kFALSE;

	}
      }

      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      if( !bpart ) continue; //default
      //-----------------
      // if( abs( phi ) > 30.*TMath::DegToRad() ) continue;
      //if( abs( phi ) < 150.*TMath::DegToRad() ) continue;
      //      if( abs( phi ) > 30.*TMath::DegToRad() && abs( phi ) < 150.*TMath::DegToRad() ) continue;
      //if( theta < 0.05 ) continue;
      //if( yaw < 0 ) continue;
      //if( yaw > 0 ) continue;
      //if( pitch > 0 ) continue;
      //if( yaw == 0.0 ) continue;
      //if( abs( pitch / yaw ) > 1. ) continue;
      // if( abs( pitch / yaw ) <= 1. ) continue;
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      bFill = kTRUE;
      auto dphi   = aPart->GetAzmAngle_wrt_RP();
      auto rapid  = aPart->GetRapiditycm();;	
      auto rapidn = rapid / y_cm[isys];
      auto rpphi  = aPart->GetIndividualRPAngle();

      hPsi->Fill(RPPsi);
      hRPPsi->Fill(rpphi);
      hphi->Fill(phi);
      hRPPsipsi->Fill(RPPsi, rpphi);
      rppsires = GetPsiRPResolution(ipsi);
      hpsi1->Fill(rppsires[4]);
      hpsi->Fill(RPangle, rppsires[4]);
      hpsindx->Fill(RPangle, (Double_t)ipsi);

      Double_t u_t0  = aPart->GetRotatedMomentum().Pt()/aPart->GetMass()/u_p;
      Double_t ou_t0 = aPart->GetMomentumAtTarget().Pt()/aPart->GetMass()/u_p;
      
      huy->Fill( rapid / y_cm[isys], u_t0);

      hyawpitch->Fill(yaw,   pitch);
      hmass    ->Fill(aPart->GetRotatedMomentum().Mag(), bmass);
      hyptacp  ->Fill(rapid/y_cm[isys], pt);

      UInt_t irapid1 = GetV1cmRapidityIndex(rapidn);
      UInt_t ipt1    = GetV1PtIndex(pt);
      
      UInt_t irapid2 = GetV2cmRapidityIndex(rapidn);
      UInt_t ipt2    = GetV2PtIndex(pt);
      UInt_t im      = GetRPCorrIndex((Double_t)aflow->mtrack4);

      hdy1[irapid1]->Fill( rapid/y_cm[isys] );
      hdy2[irapid2]->Fill( rapid/y_cm[isys] );

      hirap1 -> Fill( rapid/y_cm[isys] , (Double_t)irapid1 );
      hirap2 -> Fill( rapid/y_cm[isys] , (Double_t)irapid2 );

      hutphi10[irapid1] -> Fill( subevt_phi );
      hutphi20[irapid2] -> Fill( subevt_phi );

      if( subevt_phi > TMath::Pi()/2. ) {
	hutphi190[irapid1] -> Fill( subevt_phi );
	hutphi290[irapid2] -> Fill( subevt_phi );
      }

      if(irapid1 == 0 )
	h2phi->Fill(TVector2::Phi_mpi_pi(2.*dphi));

      hdympt1[irapid1]->Fill( pt*cos(dphi)/MassNumber );
      hdydpt1[irapid1][ipt1] -> Fill( pt );
      hdydpt2[irapid2][ipt2] -> Fill( pt );
      
      hdydut1[irapid1][ipt1] -> Fill( u_t0 );
      hdydut2[irapid2][ipt2] -> Fill( u_t0 );

      hdydptcos1[irapid1][ipt1][im]->Fill( cos(dphi) );
      hdydptcos2[irapid2][ipt2][im]->Fill( cos(2.*dphi) );

      hdydutcos1[irapid1][ipt1][ipsi]->Fill( cos(dphi) );
      hdydutcos2[irapid2][ipt2][ipsi]->Fill( cos(2.*dphi) );
      hdydutdphi1[irapid1][ipt1][ipsi]->Fill( dphi );
      

      if( abs(phi) <= 30.*TMath::DegToRad() )
	hut[irapid1][0] -> Fill( u_t0 );
      else if( abs(phi) >= 150.*TMath::DegToRad() )
	hut[irapid1][1] -> Fill( u_t0 );

      //      if( u_t0 > 0. ) {
      if( u_t0 > 0.4 ) {
	hdyucos1[irapid1][ipsi]->Fill( cos(dphi) );
	hdyucos2[irapid2][ipsi]->Fill( cos(2.*dphi) );
      }
    }
    
    if( bFill ){
      hphi0_180->Fill( subevt_phi );
      if( subevt_phi > TMath::Pi()/2. )
	hphi90_180->Fill( subevt_phi );
    }
  }
  //--------------------------------------------------
  //--- Enf of event Loop
  //--------------------------------------------------
  LOG(INFO) << " End of Event Loop " << FairLogger::endl;


  TGraphErrors *gv_v1 = new TGraphErrors();
  gv_v1->SetName("gv_v1");
  TGraphErrors *gu_v1 = new TGraphErrors();
  gu_v1->SetName("gu_v1"); 
  gu_v1->SetTitle("u_{t0}; y_{cm}/y_{cm}; v1");

  TGraphErrors *gv_v2 = new TGraphErrors();
  gv_v2->SetName("gv_v2");
  TGraphErrors *gu_v2 = new TGraphErrors();
  gu_v2->SetName("gu_v2");
  gu_v2->SetTitle("u_{t0}; y_{cm}/y_{cm}; v2");

  TGraphErrors *gPt_v1[ybin1];
  TGraphErrors *gPt_v2[ybin2];
  TGraphErrors *gUt_v1[ybin1];
  TGraphErrors *gUt_v2[ybin2];

  TH1D *corrPt_v1[ybin1];
  TH1D *corrPt_v2[ybin2];

  for(UInt_t kn = 0; kn < ybin1 ; kn++){
    gPt_v1[kn] = new TGraphErrors();
    gPt_v1[kn]->SetName((TString)Form("gPt_v1%d",kn));
    TString sname = partname[selid]+" "+rangeLabel1[kn]+"; Pt [MeV/c]; v1";
    gPt_v1[kn]->SetTitle(sname);

    gUt_v1[kn] = new TGraphErrors();
    gUt_v1[kn]->SetName((TString)Form("gUt_v1%d",kn));
    sname = partname[selid]+" "+rangeLabel1[kn]+"; u_{t0}; v1";
    gUt_v1[kn]->SetTitle(sname);
  }

  for(UInt_t kn = 0; kn < ybin2 ; kn++){
    gPt_v2[kn] = new TGraphErrors();
    gPt_v2[kn]->SetName((TString)Form("gPt_v2%d",kn));
    TString sname = partname[selid]+" "+rangeLabel2[kn]+"; Pt [MeV/c]; v2";
    gPt_v2[kn]->SetTitle(sname);

    gUt_v2[kn] = new TGraphErrors();
    gUt_v2[kn]->SetName((TString)Form("gUt_v2%d",kn));
    sname = partname[selid]+" "+rangeLabel2[kn]+"; u_{t0}; v2";
    gUt_v2[kn]->SetTitle(sname);
  }

  TGraphErrors *gmpx = new TGraphErrors();
  gmpx->SetName("gmpx");
  gmpx->SetTitle("; Rapidity/y_{cm}; <px>");


  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //--- v1
  UInt_t iyv = 0;
  UInt_t impx = 0;
  for( UInt_t iy = 0; iy < ybin1; iy++ ) {
    Double_t v_y   = 0.;
    Double_t v_ye2 = 0.;
    UInt_t iptv    = 0;

    for( UInt_t ipt = 0; ipt < ptbin1; ipt++ ) {  // ptbin
      
      Double_t ptn    = (Double_t)hdydpt1[iy][ipt]->GetEntries();
      Double_t ptmean = hdydpt1[iy][ipt]->GetMean();
      Double_t pte    = hdydpt1[iy][ipt]->GetStdDev()/sqrt(ptn);

      Double_t utn    = (Double_t)hdydut1[iy][ipt]->GetEntries();
      Double_t utmean = hdydut1[iy][ipt]->GetMean();
      Double_t ute    = hdydut1[iy][ipt]->GetStdDev()/sqrt(utn);

      // multiplicity dependent correction -->
      Double_t v_pt    = 0;
      Double_t vn_pt   = 0;
      Double_t v_pte2  = 0.;

      for( UInt_t im = 0; im < mbin; im++ ){  
	
	rpres = GetRPResolution(1, im);

	Double_t vn  = (Double_t)hdydptcos1[iy][ipt][im]->GetEntries();
	Double_t vm  = hdydptcos1[iy][ipt][im]->GetMean();
	Double_t ve  = hdydptcos1[iy][ipt][im]->GetStdDev()/sqrt(vn);
	Double_t vee = GetError(vm, rpres[1], ve, rpres[2]);  

	if( rpres[1] > 0 && vn > 0 ) {
	  v_pt   += vm/rpres[1] * vn;
	  vn_pt  += vn;
	  v_pte2 += pow(vee*vn,2);
	}
      }

      Double_t v1  = v_pt / vn_pt;
      Double_t v1e = sqrt(v_pte2)/vn_pt;

      if( !std::isnan(ptmean) && !std::isnan(pte) ) {
	gPt_v1[iy]->SetPoint(iptv, ptmean, v1);
	gPt_v1[iy]->SetPointError(iptv, pte, v1e);

	v_y   += v1*ptn;
	v_ye2 += pow(v1e*ptn,2); 
      }
      // <--


      //  phi dependent correction  ->
      Double_t v1u_sum = 0.;
      Double_t v1u_n   = 0.;
      Double_t v1u_ste = 0.;

      for( UInt_t ips = 0; ips < psibin; ips++ ){  
	Double_t nv1u  = (Double_t)hdydutcos1[iy][ipt][ips]->GetEntries();

	if( nv1u > 0 ) {
	  rppsires = GetPsiRPResolution(ips);
	  Double_t v1u   = hdydutcos1[iy][ipt][ips]->GetMean();
	  Double_t v1ue  = hdydutcos1[iy][ipt][ips]->GetStdDev()/sqrt(nv1u);
	  Double_t v1uee = GetError(v1u, rppsires[1], v1ue, rppsires[2]);

	  v1u_sum += v1u / rppsires[1] * nv1u;
	  v1u_n   += nv1u;
	  v1u_ste += pow(v1uee*nv1u,2);
	}
      }


      Double_t v1u_ave = v1u_sum / v1u_n;
      Double_t v1u_err = sqrt(v1u_ste)/v1u_n;

      if( !std::isnan(ptmean) && !std::isnan(pte) ) {
	gUt_v1[iy]->SetPoint(iptv, utmean, v1u_ave );
	gUt_v1[iy]->SetPointError(iptv, ute, v1u_err);
	iptv++;
      }
      //<- phi dependent correction 

    } ///// Pt bin
    
    Double_t yn    = (Double_t)hdy1[iy]->GetEntries();    
    Double_t ymean = hdy1[iy]->GetMean() ;
    Double_t ystdv = hdy1[iy]->GetStdDev()/sqrt(yn);

    if( yn > 0) {
      Double_t v1u_sum = 0.;
      Double_t v1u_n   = 0.;
      Double_t v1u_ste = 0.;

      for( UInt_t ips = 0; ips < psibin; ips++ ){
	Double_t nv1u  = (Double_t)hdyucos1[iy][ips]->GetEntries();

	if( nv1u > 0 ) {
	  rppsires = GetPsiRPResolution(ips);
	  Double_t v1u   = hdyucos1[iy][ips]->GetMean();
	  Double_t v1ue  = hdyucos1[iy][ips]->GetStdDev()/sqrt( nv1u );
	  Double_t v1uee = GetError(v1u, rppsires[1], v1ue, rppsires[2]);
	  
	  v1u_sum += v1u / rppsires[1] * nv1u;
	  v1u_n   += nv1u;
	  v1u_ste += pow(v1uee*nv1u,2);
	}

      }
      Double_t v1_y  = v_y / yn;
      Double_t v1_ye = sqrt(v_ye2)/yn;

      gv_v1->SetPoint( iyv, ymean, v1_y);
      gv_v1->SetPointError( iyv, ystdv, v1_ye);
      
      Double_t v1u_ave = v1u_sum / v1u_n;
      Double_t v1u_err = sqrt(v1u_ste)/v1u_n;
      
      gu_v1->SetPoint( iyv, ymean, v1u_ave);
      gu_v1->SetPointError( iyv, ystdv, v1u_err);

      iyv++;
    }
  }


  //--- <px>
  Double_t *rpresall = new Double_t[4];
  rpresall = GetRPResolutionwChi(hphi0_180, hphi90_180);

  if( !std::isnan(rpresall[0]) ){
    for( UInt_t iy = 0; iy < ybin1; iy++ ) {

      Double_t yn    = (Double_t)hdy1[iy]->GetEntries();    
      Double_t ymean = hdy1[iy]->GetMean() ;
      Double_t ystdv = hdy1[iy]->GetStdDev()/sqrt(yn);

      Double_t mpxn = (Double_t)hdympt1[iy]->GetEntries();
      Double_t mpxc = hdympt1[iy]->GetMean() / rpresall[0] ;
      Double_t mpxe = hdympt1[iy]->GetStdDev()/sqrt(mpxn);

      Double_t mpxee = GetError(mpxc, rpresall[0], mpxe, rpresall[1]);
      if( mpxn > 0 ) {
	gmpx->SetPoint( impx, ymean, mpxc ); 
	gmpx->SetPointError( impx, ystdv, mpxe ); 
	impx++;
      }      
      gPt_v1[iy]->Write();
      gUt_v1[iy]->Write();
    }
  }
  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //--- v2
  iyv = 0;
  for( UInt_t iy = 0; iy < ybin2; iy++ ) {
    Double_t v_y   = 0.;
    Double_t v_ye2 = 0.;
    UInt_t iptv    = 0;

    for( UInt_t ipt = 0; ipt < ptbin2; ipt++ ){
    
      Double_t ptn    = (Double_t)hdydpt2[iy][ipt]->GetEntries();
      Double_t ptmean = hdydpt2[iy][ipt]->GetMean();
      Double_t pte    = hdydpt2[iy][ipt]->GetStdDev()/sqrt(ptn);

      Double_t utn    = (Double_t)hdydut2[iy][ipt]->GetEntries();
      Double_t utmean = hdydut2[iy][ipt]->GetMean();
      Double_t ute    = hdydut2[iy][ipt]->GetStdDev()/sqrt(utn);

      // multiplicity dependent correction -->
      Double_t v_pt    = 0;
      Double_t vn_pt   = 0;
      Double_t v_pte2  = 0.;

      for( UInt_t im = 0; im < mbin; im++ ){
	
	rpres = GetRPResolution(2, im);
	
	Double_t vn  = (Double_t)hdydptcos2[iy][ipt][im]->GetEntries();
	Double_t vm  = hdydptcos2[iy][ipt][im]->GetMean();
	Double_t ve  = hdydptcos2[iy][ipt][im]->GetStdDev()/sqrt(vn);
	Double_t vee = GetError( vm, rpres[1], ve, rpres[2] );

	if( rpres[1] > 0 && vn > 0 ) {
	  v_pt   += vm/rpres[1] * vn;
	  vn_pt  += vn;
	  v_pte2 += pow(vee*vn,2);
	}
      }

      Double_t v2  = v_pt / vn_pt;
      Double_t v2e = sqrt(v_pte2)/vn_pt;

      if( !std::isnan(ptmean) && !std::isnan(pte) ) {
	gPt_v2[iy]->SetPoint(iptv, ptmean, v2);
	gPt_v2[iy]->SetPointError(iptv, pte, v2e);

	v_y   += v2*ptn;
        v_ye2 += pow(v2e*ptn,2);
      }
      // <--


      //  phi dependent correction  ->
      Double_t v2u_sum = 0.;
      Double_t v2u_n   = 0.;
      Double_t v2u_ste = 0.;

      for( UInt_t ips = 0; ips < psibin; ips++ ){
	Double_t nv2u  = (Double_t)hdydutcos2[iy][ipt][ips]->GetEntries();

	if( nv2u > 0 ) {
	  rppsires = GetPsiRPResolution(ips);
	  Double_t v2u   = hdydutcos2[iy][ipt][ips]->GetMean();
	  Double_t v2ue  = hdydutcos2[iy][ipt][ips]->GetStdDev()/sqrt( nv2u );
	  Double_t v2uee = GetError(v2u, rppsires[4], v2ue, rppsires[5]);

	  v2u_sum += v2u / rppsires[4] * nv2u;
	  v2u_n   += nv2u;
	  v2u_ste += pow(v2uee*nv2u,2);

	  LOG(DEBUG) << " v2 corr " << ips << " res " << rppsires[4] << " +- " << rppsires[5] 
		     << " n " << nv2u << " v2u = " << v2u << " -> " << v2u/rppsires[4]
		     << FairLogger::endl;
	}
      }

      Double_t v2u_ave = v2u_sum / v2u_n;
      Double_t v2u_err = sqrt(v2u_ste)/v2u_n;

      if( !std::isnan(ptmean) && !std::isnan(pte) ) {
	gUt_v2[iy]->SetPoint(iptv, utmean, v2u_ave);
	gUt_v2[iy]->SetPointError(iptv, ute, v2u_err);
	iptv++;
      }
      //<- phi dependent correction 

    }  ///// Pt bin


    Double_t yn    = (Double_t)hdy2[iy]->GetEntries();
    if( yn > 0 ) {

      Double_t ymean = hdy2[iy]->GetMean();
      Double_t ystdv = hdy2[iy]->GetStdDev();

      Double_t v2u_sum = 0.;
      Double_t v2u_n   = 0.;
      Double_t v2u_ste = 0.;

      for( UInt_t ips = 0; ips < psibin; ips++ ){
	rppsires = GetPsiRPResolution(ips);
	Double_t nv2u  = (Double_t)hdyucos2[iy][ips]->GetEntries();

	
	if( nv2u > 0 ) {
	  Double_t v2u   = hdyucos2[iy][ips]->GetMean();
	  Double_t v2ue  = hdyucos2[iy][ips]->GetStdDev()/sqrt( nv2u );
	  Double_t v2uee = GetError(v2u, rppsires[4], v2ue, rppsires[5]);
	  
	  v2u_sum += v2u / rppsires[4] * nv2u;
	  v2u_n   += nv2u;
	  v2u_ste += pow(v2uee*nv2u,2);

	}
      }

      Double_t v2_y  = v_y / yn;
      Double_t v2_ye = sqrt(v_ye2)/yn;

      gv_v2->SetPoint( iyv, ymean, v2_y);
      gv_v2->SetPointError( iyv, ystdv, v2_ye);

      Double_t v2u_ave = v2u_sum / v2u_n;
      Double_t v2u_err = sqrt(v2u_ste)/v2u_n;

      gu_v2->SetPoint( iyv, ymean, v2u_ave);
      gu_v2->SetPointError( iyv, ystdv, v2u_err );

      iyv++;
    }
  
    gPt_v2[iy]->Write();
    gUt_v2[iy]->Write();
  }
  ///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  gmpx->Write();  
  gv_v1->Write();
  gv_v2->Write();
  gu_v1->Write();
  gu_v2->Write();
  hmult->Write();
  hyawpitch->Write();
  hmass    ->Write();
  hyptacp  ->Write();
  hirap1   ->Write();
  hirap2   ->Write();
  huy      ->Write();
  hpsi     ->Write();
  hpsi1    ->Write();
  hpsindx  ->Write();
  hiphi    ->Write();
  hPsi     ->Write();
  hRPPsi    ->Write();
  hRPPsipsi ->Write();
  hphi      ->Write();
  //--------------------------------------------------
  //--- Plotting
  //--------------------------------------------------
  id=1;
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),2000,500);
  cc->Divide(10);
  for( UInt_t iy = 0; iy < ybin1-1; iy++ ) {
    cc->cd(id); id++;
    // gPt_v1[iy] -> SetMarkerStyle(20);
    // gPt_v1[iy] -> SetMarkerColor(4);
    // gPt_v1[iy] -> Draw("ALP");

    gUt_v1[iy] -> SetMarkerStyle(20);
    gUt_v1[iy] -> SetMarkerColor(4);
    gUt_v1[iy] -> Draw("ALP");
  }

  if( bEffCorr && corrPt_v1[0] != NULL ){
    id=1;
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),2000,500);
    cc->Divide(10);
    for( UInt_t iy = 0; iy < ybin1-1; iy++ ) {
      cc->cd(id); id++;
      corrPt_v1[iy] -> SetMarkerStyle(20);
      corrPt_v1[iy] -> SetMarkerColor(4);
      corrPt_v1[iy] -> Draw();
    }
  }

  id=1;
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),2000,500);
  cc->Divide(7);
  for( UInt_t iy = 0; iy < ybin2-1; iy++ ) {
    cc->cd(id); id++;
    // gPt_v2[iy] -> SetMarkerStyle(20);
    // gPt_v2[iy] -> SetMarkerColor(4);
    // gPt_v2[iy] -> Draw("ALP");

    gUt_v2[iy] -> SetMarkerStyle(20);
    gUt_v2[iy] -> SetMarkerColor(4);
    gUt_v2[iy] -> Draw("ALP");
  }



  if( bEffCorr && corrPt_v2[0] != NULL ){
    id=1;
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),2000,500);
    cc->Divide(7);
    for( UInt_t iy = 0; iy < ybin2-1; iy++ ) {
      cc->cd(id); id++;
      corrPt_v2[iy] -> SetMarkerStyle(20);
      corrPt_v2[iy] -> SetMarkerColor(4);
      corrPt_v2[iy] -> Draw();
    }
  }

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  gv_v1->Draw("ALP");

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  gu_v1->Draw("ALP");
  if( isys == 5 ) {
    SimFunction();
    fv1y->Draw("same");
  }

  auto LineV = new TLine(0.,gv_v1->GetYaxis()->GetXmin(), 0., gv_v1->GetYaxis()->GetXmax());
  auto LineH = new TLine(gv_v1->GetXaxis()->GetXmin(),    0., gv_v1->GetXaxis()->GetXmax(), 0.);
  LineV->SetLineStyle(3);
  LineH->SetLineStyle(3);
  LineV->Draw();
  LineH->Draw();

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  gv_v2->Draw("ALP");

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  gu_v2->Draw("ALP");
  if(isys == 5)
    fv2y->Draw("same");


  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  hpsi->Draw("colz");


  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  huy->Draw("colz");

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  hpsindx->Draw("colz");


  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  cc->Divide(2,2);
  cc->cd(1);
  hRPPsi->Draw();
  cc->cd(2);
  hPsi->Draw();
  cc->cd(3);
  hRPPsipsi->Draw("colz");
  cc->cd(4);
  hphi->Draw();

  //--------------------------------------------------
  //--- Debug plotting
  // id=1;
  // ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),2000,500);
  // cc->Divide(ybin2, psibin);
  //   for( UInt_t iy = 0; iy < ybin2; iy++ ) {
     for( UInt_t ips = 0; ips < psibin; ips++) {
       //       cc->cd(id); id++;
       //       hdydutcos2[iy][2][ips]->Draw();
       auto mm = hdydutcos2[3][2][ips]->GetMean();
       //       cout << " ips " << ips << " " << mm << endl;;


     }
     //   }


  // id=1;
  // ic++; cc = new TCanvas(Form("cc%d",ic),"y : cos(2#psi)",2000,1500);
  // cc->Divide(12,ybin1);
  // UInt_t ipt = 4;
  // for( UInt_t iy = 0; iy < ybin1; iy++ ) {
  //   for( UInt_t ips = 0; ips < psibin; ips++ ) {
  //     cc->cd(id); id++;
  //     hdydutcos1[iy][ipt][ips]->Draw();
  //   }
  // }

  // id=1;
  // ic++; cc = new TCanvas(Form("cc%d",ic),"pt : cos(#psi) 2",2000,1500);
  // cc->Divide(12,ptbin1);
  // UInt_t iy = 2;

  // std::cout << " Raidity range " << rangeLabel1[iy] << std::endl;
  // for( UInt_t ipt = 0; ipt < ptbin1; ipt++ ) {
  //   for( UInt_t ips = 0; ips < psibin; ips++ ) {
  //     cc->cd(id); id++;
  //     hdydutcos1[iy][ipt][ips]->Draw();
  //   }
  // }

  // id=1;
  // ic++; cc = new TCanvas(Form("cc%d",ic),"pt : #delta #phi 2",2000,1500);
  // cc->Divide(12,ptbin1);
  // for( UInt_t ipt = 0; ipt < ptbin1; ipt++ ) {
  //   for( UInt_t ips = 0; ips < psibin; ips++ ) {
  //     cc->cd(id); id++;
  //     hdydutdphi1[iy][ipt][ips]->Draw();
  //   }
  // }

  // id=1;
  // ic++; cc = new TCanvas(Form("cc%d",ic),"pt : #delta #phi 7",2000,1500);
  // cc->Divide(12,ptbin1);
  // iy = 7;
  // std::cout << " Raidity range " << rangeLabel1[iy] << std::endl;
  // for( UInt_t ipt = 0; ipt < ptbin1; ipt++ ) {
  //   for( UInt_t ips = 0; ips < psibin; ips++ ) {
  //     cc->cd(id); id++;
  //     hdydutdphi1[iy][ipt][ips]->Draw();
  //   }
  // }

  gSystem->cd("..");
} //



void AzimuthalAngleDependence()            //%% Executable :
{
  TString fHeader = "azm_"+ sysName + ".v"+sVer+".";
  auto fName = SetupOutputFile( fHeader );
  auto GraphSave = new TFile(fName,"recreate");

  TH1I *hmult = new TH1I("hmult","multiplicity",100,0,100);

  const UInt_t nphi = 12;
  TH1I *hphibin[nphi];
  TH1D *hphi0_180[nphi];
  TH1D *hphi90_180[nphi];
  
  for(UInt_t k = 0; k < mbin; k++){
    TString htitle = Form("hphi0_180_%d",k);
    hphi0_180[k]  = new TH1D(htitle, "",100,0.,3.2);
    htitle = Form("hphi90_180_%d",k);
    hphi90_180[k] = new TH1D(htitle,"",100,0.,3.2);
  }

  auto hiphi = new TH1I("hiphi","hiphi",15,0,15);


  //-------------------------------------------------- 
  //--- Event Loop 
  //--------------------------------------------------                                                                                                  
  Long64_t nEntry = SetBranch();
 
  for(Int_t i = 0; i < nEntry; i++){

    rChain->GetEntry(i);
    STFlowInfo *aflow = (STFlowInfo*)aFlowArray->At(0);

    /// for reaction plane resolution  
    Bool_t bFill = kFALSE;
    Bool_t bRes  = kFALSE;

    /// centrality selection   
    if(aflow->mtrack1 > Ucent || aflow->mtrack1 <= Lcent || aflow->mtrack4 < 6) continue;

    hmult->Fill( aflow->mtrack2 );
    bRes = kTRUE; //@@@ 
    TIter next(aArray);
    STParticle *aPart = NULL;

    //--------------------------------------------------
    //----- Main loop
    //--------------------------------------------------       
    while( (aPart = (STParticle*)next()) ) {

      auto pid   = aPart->GetPIDTight();  // = GetPID()   

      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
      // ---- quality selection ---
      if( !aPart->GetNDFFlag() ) continue;
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
      
      bFill = kTRUE;

      auto yaw   = aPart->GetYawAngle();
      auto pitch = aPart->GetPitchAngle();
      auto phi   = TVector2::Phi_0_2pi(aPart->GetRotatedMomentum().Phi());
      auto theta = aPart->GetRotatedMomentum().Theta();

      UInt_t iphi = UInt_t(phi/(TMath::Pi()/6.));

      hiphi->Fill(iphi);

      Double_t subdphi = abs(TVector2::Phi_mpi_pi((aflow->unitP_1fc).Phi() - (aflow->unitP_2fc).Phi()));
      hphi0_180[iphi] ->Fill(subdphi);
      if( subdphi > TMath::Pi()/2. )
	hphi90_180[iphi]->Fill(subdphi);
    }
  }

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  hiphi->Draw();


  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),500,1000);
  cc->Divide(2,nphi/2);

  auto gv_phi1 = new TGraphErrors();
  gv_phi1->SetName("gv_phi1");
  gv_phi1->SetTitle("; #phi; <cos(#Delta #Psi)>");
   
  UInt_t id = 1;
  UInt_t ip = 0;
  for(UInt_t i = 0; i < nphi; i++) {
    cc->cd(id); id++;

    if( hphi0_180[i]->GetEntries() < 5 ) continue;

    hphi0_180[i] ->Draw();
    hphi90_180[i]->Draw("same");

    Double_t *rpres = new Double_t[4];
    rpres = GetRPResolutionwChi(hphi0_180[i], hphi90_180[i]);
     
    if( hphi90_180[i]->GetEntries() > 15 && !std::isnan(rpres[0])) {
      gv_phi1->SetPoint(ip, TVector2::Phi_mpi_pi(i*TMath::Pi()/6.), rpres[0]);
      gv_phi1->SetPointError(ip, 0., rpres[1]);
      ip++;
    }
  }

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));  
  gv_phi1->Draw("AP");

  hmult->Write();
  hiphi->Write();
  gv_phi1->Write();


  std::cout << " FILE output is " << fName << std::endl;
}



void PtDistribution(UInt_t selid)
{
  TH1D *hut[ybin1][2];
  TString rangeLabel1[ybin1];
  Double_t u_p = 0.355151 * 1.06974; //p+p(268.9MeV/u)
  
  for( UInt_t iy = 0; iy < ybin1; iy++ ) {
    rangeLabel1[iy] = Form("%5.2f <= y < %5.2f"      ,yrange1[iy],yrange1[iy+1]);
    hut[iy][0]    = new TH1D(Form("hut0_%d",iy)      ,rangeLabel1[iy]+"_yaw>0;U_{t};",100,0., 2.5);
    hut[iy][1]    = new TH1D(Form("hut1_%d",iy)      ,rangeLabel1[iy]+"_yaw<0;U_{t};",100,0., 2.5);
  }

  Int_t nevt = SetBranch();

  //--------------------------------------------------
  //--- Event Loop
  //--------------------------------------------------
  for(Int_t i = 0; i < nevt; i++){

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
///######################################################################
Double_t *GetRPResolution(UInt_t vn, UInt_t xn)
{

  
  Double_t *rpcor = new Double_t[3];
  rpcor[0] = 0.;
  rpcor[1] = 0.;
  rpcor[2] = 0.;

  if( xn < v1x.size() ) {

    if( vn == 1 ){
      rpcor[0] = v1x.at(xn);   // multiplicity
      rpcor[1] = v1y.at(xn);   // v1 corr 
      rpcor[2] = v1ye.at(xn);  // v1 corr error
    }
    else { //v2
      rpcor[0] = v2x.at(xn);
      rpcor[1] = v2y.at(xn);
      rpcor[2] = v2ye.at(xn);
    }
  }

  return rpcor;

}

Double_t *GetPsiRPResolution(UInt_t ival)
{
  Double_t *rpcor = new Double_t[6];
  for( UInt_t i = 0; i < 7; i++ ) 
    rpcor[i] = 0.;

  rpcor[0] = v1psix.at(ival);   // Psi angle
  rpcor[1] = v1psiy.at(ival);   // v1 corr 
  rpcor[2] = v1psiye.at(ival);  // v1 corr error

  rpcor[3] = v2psix.at(ival);
  rpcor[4] = v2psiy.at(ival);    // v2 corr
  rpcor[5] = v2psiye.at(ival);   // v2 corr error

  // rpcor[0] = v1psix.at(ival);  // Psi angle
  // rpcor[1] = 1.;  // v1 corr 
  // rpcor[2] = 0.;  // v1 corr error
	     
  // rpcor[3] = v2psix.at(ival);
  // rpcor[4] = 1.;   // v2 corr
  // rpcor[5] = 0.;   // v2 corr error

  return rpcor;
}

//**************************************************
TString GetPsiRPLoadFileName()
{
  TString fn = "data/bpsi_"+ sysName + ".v" + dVer;
  
  cout << " dver " << dVer << endl;

  if( dVer == "")
    fn = "data/bpsi_"+ sysName + ".v" + sVer;

  if( RPBase < 3 && RPBase > 0)
    fn += Form("_%d", RPBase);

  fn += Form(".m%02dto%02d.root",Lcent,Ucent);

  return fn;
}

//**************************************************
Bool_t SetPsiRPResolution()
{
  //  itrpvpsix = new ROOT::Math::Interpolator(20, ROOT::Math::Interpolation::kPOLYNOMIAL);

  TString fname = GetPsiRPLoadFileName(); 
  TFile *fOpen = TFile::Open(fname);

  TGraphErrors* hgv_psi1;
  TGraphErrors* hgv_psi2;


  if( fOpen != NULL ){
    LOG(INFO) << fname << " is opened. " << FairLogger::endl;
    hgv_psi1 = (TGraphErrors*)fOpen->Get("gv_psi1");
    hgv_psi2 = (TGraphErrors*)fOpen->Get("gv_psi2");
    bCorr = kTRUE;
  }
  else {
    LOG(INFO) << fname << " is NOT opened. " << FairLogger::endl;

    std::cout << " Do you continue ? (y or n) " << std::endl;
    TString aAns;
    std::cin >> aAns;
    if( aAns == "n" ) return kFALSE;
  }



  std::vector< Double_t > ix;
    
  Double_t x, y, xe, ye;
  UInt_t k = 0;

  if( bCorr ) {
    for( Int_t i = 0; i < (Int_t)hgv_psi1->GetN(); i++ ) {

      hgv_psi1->GetPoint(i, x, y);
      xe = hgv_psi1->GetErrorX(i);
      ye = hgv_psi1->GetErrorY(i);

      ix.push_back(k); 
      v1psix.push_back(x);
      v1psiy.push_back(y);
      v1psixe.push_back(xe);
      v1psiye.push_back(ye);

      k++;
    }
  }

  else{
    for( Int_t i = 0; i < npsi; i++ ) {
      // for no correction
      ix.push_back(k); 
      Double_t dphi = 2.*TMath::Pi()/ (Double_t)npsi ;
      v1psix.push_back( TVector2::Phi_mpi_pi( (Double_t)(i - npsi/2 ) * dphi ) );
      v1psiy.push_back(1);
      v1psixe.push_back(0.);
      v1psiye.push_back(0.);

      k++;
    }
  }
  for(UInt_t j = 0; j < (UInt_t)v1psix.size(); j++) {
    LOG(INFO) << "v1 resolution:"<< j << " th " << v1psix.at(j) << " vs " << v1psiy.at(j) << " +- " << v1psiye.at(j) << FairLogger::endl; 
  }

  //  itrpvpsix->SetData(v1psix,ix);

  if( bCorr ) {
    for( Int_t i = 0; i < (Int_t)hgv_psi2->GetN(); i++ ) {
      hgv_psi2->GetPoint(i, x, y);
      xe = hgv_psi2->GetErrorX(i);
      ye = hgv_psi2->GetErrorY(i);

      v2psix.push_back(x);
      v2psiy.push_back(y);
      v2psixe.push_back(xe);
      v2psiye.push_back(ye);
    }
  }
  else {
    for( Int_t i = 0; i < npsi; i++ ) {
      // for no correction
      Double_t dphi = 2.*TMath::Pi()/ (Double_t)npsi ;
      v2psix.push_back( TVector2::Phi_mpi_pi( (Double_t)i * dphi ) );
      v2psiy.push_back(1);
      v2psixe.push_back(0.);
      v2psiye.push_back(0.);
    }
  }

  for(UInt_t j = 0; j < (UInt_t)v2psix.size(); j++) {
    LOG(INFO) << "v2 resolution: "<< j << " th " << v2psix.at(j) << " vs " << v2psiy.at(j) << " +- " << v2psiye.at(j) << FairLogger::endl; 
  }

  if( fOpen != NULL )  fOpen->Close();

  return kTRUE;
}
//**************************************************
Bool_t LoadRPResolution()
{
  itrpvx = new ROOT::Math::Interpolator(20, ROOT::Math::Interpolation::kPOLYNOMIAL);

  TString fname = "data/mlt_"+ sysName + ".v" + dVer + ".root";
  if( dVer == "" )
    fname = "data/mlt_"+ sysName + ".v" + sVer + ".root";

  TFile *fOpen = TFile::Open(fname);
  if( fOpen == NULL ) {
    LOG(ERROR) << "Please do it doflow -2 " << FairLogger::endl;
    return kFALSE;
  }
  
  std::cout << fname << " is opened. " << std::endl;

  auto hgv_mcos1 = (TGraphErrors*)fOpen->Get("gv_mcos1");
  auto hgv_mcos2 = (TGraphErrors*)fOpen->Get("gv_mcos2");

  std::vector< Double_t > ix;
    
  Double_t x, y, xe, ye;
  UInt_t k = 0;
  ix.push_back(k); 
  v1x.push_back(0.); v1xe.push_back(0.);
  v1y.push_back(0.); v1ye.push_back(0.);
  k++;
  for( Int_t i = (Int_t)hgv_mcos1->GetN()-1; i > -1; i-- ) {
    hgv_mcos1->GetPoint(i, x, y);
    xe = hgv_mcos1->GetErrorX(i);
    ye = hgv_mcos1->GetErrorY(i);


    ix.push_back(k); 
    v1x.push_back(x);
    v1y.push_back(y);
    v1xe.push_back(xe);
    v1ye.push_back(ye);

    k++;
  }
  if( v1x.size() > 1 )
    itrpvx->SetData(v1x,ix);


  v2x.push_back(0.); v2xe.push_back(0.);
  v2y.push_back(0.); v2ye.push_back(0.);

  k = 0;
  for( Int_t i = (Int_t)hgv_mcos2->GetN()-1; i > -1; i-- ) {
    hgv_mcos2->GetPoint(i, x, y);
    xe = hgv_mcos2->GetErrorX(i);
    ye = hgv_mcos2->GetErrorY(i);


    v2x.push_back(x);
    v2y.push_back(y);
    v2xe.push_back(xe);
    v2ye.push_back(ye);

    LOG(DEBUG) << "v2 resolution "<< k << " th " << v2x.at(k) << " vs " << v2y.at(k) << " +- " << v2ye.at(k) << FairLogger::endl; 
    k++;
  }

  fOpen->Close();

  return kTRUE;
}

UInt_t GetV1RapidityIndex(Double_t y)
{
  UInt_t irapid1 = ybin1-1;
  for( UInt_t i = 1; i < ybin1; i++){
    if( y < yrange1[i]){
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
    if( y < yrange1[i]/y_cm[isys] ){
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
    if( y < yrange2[i]){
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
    if( y < yrange2[i]/y_cm[isys]){
      irapid2 = i-1;
      break;
    }
  }
  return irapid2;
}

UInt_t GetV1PtIndex(Double_t val)
{
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
  TString oVer  = gSystem->Getenv("OUTVER");
  TString fName = fopt + oVer;

  if( oVer == "" ) {
    UInt_t kVer = 0; 
    TString checkName;

    while( kTRUE ) {
      checkName = Form(fName + "%04d.root", kVer);
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
      std::cout << checkName << " is existing. Do you recreate? (y/n)" << std::endl;
      TString sAns;
      std::cin >> sAns;
      if(sAns == "N" || sAns == "n") {
	std::cout << " Retry" << std::endl;
	exit(0);
      }
    }
    fName += ".root";
  }

  std::cout << "File " << fName << " will be created. " << std::endl;  

  
  return fName;
}
//--------------------------------------------------

Bool_t SetupEffCorrection()
{
  fefffile = TFile::Open("correctedPt.phi60.nomultcut.root");

  if( fefffile != NULL ) return kTRUE;
  else
    return kFALSE;
}


//---------------------------------------------------

void DetectorBias()
{
  UInt_t m = 0;

  SetBranch();
  Int_t nEntry = rChain->GetEntries();


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
  
  std::cout << " <cos Phi> " << cosphi/count 
	    << " <sin phi> " << sinphi/count << std::endl;

  std::cout << " <cos 2Phi> " << cos2phi/count 
	    << " <sin 2phi> " << sin2phi/count << std::endl;

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


///######################################################################


