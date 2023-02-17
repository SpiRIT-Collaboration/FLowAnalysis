#include "DoRPRes.C" //#include "openRunAna.C" #include "SimFunction.C"
#include "STLorentzBoostVector.hh"


//drawing
UInt_t selReactionPlanef = 10000;
UInt_t  phicutID = 5;

TFile* fefffile;

Bool_t bccPsi;


// functions
void     makePhysicsPlot(UInt_t selid = 0);
TString  SetupOutputFile(TString fopt);
Bool_t   SetupEffCorrection();

Bool_t   LoadAcceptanceCorrection(UInt_t selid);
Double_t *GetPsiRPResolution (Double_t *rpcor, UInt_t ival);
Double_t *Get2PsiRPResolution(Double_t *rpcor, UInt_t ival);
Double_t *GetMultRPResolution(Double_t *rpcor, UInt_t vn, UInt_t mult);

UInt_t   GetV1RapidityIndex(Double_t y);
UInt_t   GetV2RapidityIndex(Double_t y);
UInt_t   GetV1cmRapidityIndex(Double_t y);
UInt_t   GetV2cmRapidityIndex(Double_t y);
UInt_t   GetV1PtIndex(Double_t val);
UInt_t   GetV2PtIndex(Double_t val);

void     DrawUt(TH2D *h2, UInt_t sel, TString opt="", Color_t scol=1);
TString  GetPsiRPLoadFileName();
void     GetdNdydPt(Double_t ycm, Double_t pt);
void     GetdNdydUt( Double_t ycm, Double_t ut, TH2D* h2c);
void     GetdNdydUt( Int_t ih, Int_t iy, Int_t ipt, TH2D *h2c);
TH2D*    GetAcpCorrection(TH2D &h2, TH2D *hcorr);
TH1D*    GetdNdy(TH2D *h2c, Double_t ntotal);
TH1D*    GetdNdXt(TH2D* h2c, TString htitle, Double_t ylow, Double_t yup, Double_t ntotal);

TFile*        GraphSave;
TH2D*         hAcpYPtCorr;
TH2D*         hAcpYUtCorr;
TH2D*         hAcpyytCorr;
TH2D*         hAcpyyxCorr;
TH2D*         hAcpyxhCorr;

Double_t      fmass;
Double_t      *acpcorr = new Double_t[2];
Double_t      phiAcp = 1;
Double_t      AtmNumber;

//-------------------//
void DoAna_phys(Int_t isel = -2) 
{
  gROOT->Reset();
  
  openRunAna();

  if(rChain != NULL) 
    LOG(INFO) << " DoAna_phys: System " << isys << "  -> " << sysName << FairLogger::endl; 

  else
    exit(0);


  gROOT->ProcessLine(".! grep -i void DoAna_phys.C ");


  //@@@@@ Configure ++++++++++++++++++++++++++++++++++++++++++++++++++
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

  if( isys == 0 || isys == 2 ) {
    Lcent += 1;
    Ucent += 1;
  }
  LOG(INFO) << "Multiplicity :: " << Lcent << " <= M <=" << Ucent  <<  FairLogger::endl;

  //------  Reaction plane definition --------------------------------------------
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

  //@@@@@ phi angle cut ++++++++++++++++++++++++++++++++++++++++++++++++++
  TString PhiCut = gSystem -> Getenv("PHICUT"); 
  phicutID = PhiCut != "" ? atoi(PhiCut) : 4;
  
  TString PhiCutDef[] = {"|phi|<45 or |phi|>135",
			 "|phi|<45",
			 "|phi|>140",
			 "45<|phi|<135",
			 "all",
			 "-30<phi<20",
			 "-30<phi<20 or |phi|>140"  
};

  LOG(INFO) << "Phi angle cut " << phicutID << " : "  << PhiCutDef[phicutID]  << FairLogger::endl;
  
  acpcorr[0] = 1.;
  acpcorr[1] = 0.;

  // dpt1 = ut_max[isel-2]/(Double_t)ptbin1;
  // dpt2 = ut_max[isel-2]/(Double_t)ptbin2;

  dpt1 = 2.0/(Double_t)ptbin1;
  dpt2 = 2.0/(Double_t)ptbin2;

  if( isel >= 0 ) 
    makePhysicsPlot((UInt_t)isel);  
  else if( isel == -1 )
    makePhysicsPlot(0);  

}


//@@@@@ main ++++++++++++++++++++++++++++++++++++++++++++++++++
void makePhysicsPlot(UInt_t selid = 0)       //%% Executable :
{
  LOG(INFO) << "makePhysicsPlot(" << selid << ")" << FairLogger::endl;

  TDatime beginTime;
  gStyle->SetOptStat(0);

  std::map< TString, UInt_t > pToIDX = {{"1H",0},{"2H",1},{"3H",2},{"3He",3},{"4He",4},{"H",5},{"HHe",6},{"n",7}};
  std::map< TString, Long64_t > ptoPDG = {{"1H",2212},{"1H",1000010020},{"3H",1000010030},{"3He",1000020030},
					  {"4He",1000020040},{"H",2212},{"HHe",2212},{"n",2112}};
  std::map< Long64_t, Int_t > PDGtoIDX{{2212, 0},{1000010020,1},{1000010030,2},{1000020030,3},{1000020040,4},{2112,7}};
  std::vector<Long64_t> nclPDG = {2212,1000010020,1000010030,1000020030,1000020040,2122};

  // beamAToSys[132] = 0;
  // beamAToSys[108] = 1;
  // beamAToSys[124] = 2;
  // beamAToSys[112] = 3;

  ///------>>>> Acceptance Correction Option: --------
  ///$$$$$/
  //Bool_t bAcc_corr = kFALSE;
  Bool_t bAcc_corr = kTRUE;

  if( bAcc_corr && LoadAcceptanceCorrection(selid) ) {
    LOG(INFO) << " Acceptance correction is found. " << FairLogger::endl;
  }
  else {
    LOG(INFO) << " Acceptance correction is OFF." << FairLogger::endl;
    bAcc_corr = kFALSE;
  }


  LOG(INFO) << sysName << " Particle is " << ncls[selid].name << " " << sVer << FairLogger::endl;
  //-- Define Output file name 
  TString fHeader = "phys_"+ sysName + "_" + ncls[selid].name+".v"+sVer+".";
  auto fName = SetupOutputFile( fHeader );
  LOG(INFO) << " output file is " << fName << FairLogger::endl;
  
  GraphSave  = new TFile(fName,"recreate");

  LOG(INFO) << " Rapidity binning " << ybin1 << FairLogger::endl;


  //@@@@@ Booking ++++++++++++++++++++++++++++++++++++++++++++++++++
  LOG(INFO) << " v1 BIN y " << ybin1 << " pt " << ptbin1 << " mult " << mbin << FairLogger::endl;
  LOG(INFO) << " v2 BIN y " << ybin2 << " pt " << ptbin2 << " mult " << mbin << FairLogger::endl;

  TH1I *hmult        = new TH1I("hmult",";Multiplicity",80,0,80);
  TH1I *hmultfill    = new TH1I("hmultfill",";Multiplicity with filling",80,0,80);
  TH1D *hdphi0_180   = new TH1D("hdphi0_180"  , " 0to180",100,0.,3.2);
  TH1D *hdphi90_180  = new TH1D("hdphi90_180" , "90to180",100,0.,3.2);

  auto hpid       = new TH2D("hpid",   "All particles    ; P/Z[MeV/c]; dEdx[ADC/mm]",400,-800.,3000.,300,0.,1500);
  auto hpidsel    = new TH2D("hpidsel","selected particle; P/Z[MeV/c]; dEdx[ADC/mm]",400,-800.,3000.,300,0.,1500);
  auto hmass      = new TH2D("hmass",   ";P/Q; Mass [MeV/c^2]"     ,200,  0.,2500., 200, 0.,4000);
  auto hmassHe    = new TH2D("hmassHe", ";P/Q; MassHe[MeV/c^2]"    ,200,  0.,3000., 200, 0.,7000);

  TString hlabel  = (TString)Form("mtrack1 %2d to %2d",Lcent,Ucent);
  auto hypt       = new TH2D("hypt"   ,hlabel+";y/y_{nn}-1; Pt[GeV/c]",40, -2., 2., 50, 0., 2000.);
  auto hyut       = new TH2D("hyut"   ,hlabel+";y/y_{nn}-1; Ut"       ,40, -2., 2., 50, 0., 2.);   
  auto hyyx       = new TH2D("hyyx"   ,hlabel+";Rapidity;Rapidity_x"  ,40, -2., 2., 40,-2.,2.);    
  auto hyyxrnd    = new TH2D("hyyxrnd",hlabel+";Rapidity;Rapidity_x(rand)",40, -2., 2., 40,-2.,2.);    
  auto hypx       = new TH2D("hypx"   ,hlabel+";y/y_{nn}-1;px/A"      ,40,-2.,2.,60,-1200.,1200.); 
  auto hypxrnd    = new TH2D("hypxrnd",hlabel+";y/y_{nn}-1;px(rnd)/A" ,40,-2.,2.,60,-1200.,1200.); 
  auto hyEt       = new TH2D("hyEt",";;Transverse Ek"                 ,40, -2., 2., 50,    0., 500.);
  auto hyEtrnd    = new TH2D("hyEtrnd",";;Transverse Ek on randon"    ,40, -2., 2.,100, -500., 500.);
  auto hybtgm     = new TH2D("hybtgm",";;#beta*#gamma"                ,40, -2., 2., 50,    0., 1.2);
  auto hybtgmrnd    = new TH2D("hybtgmrnd",";;Transverse beta*gamma"      ,40, -2., 2.,100,  -1.2, 1.2);

  hypt->Sumw2();
  hyut->Sumw2();
  hyyx->Sumw2();
  hyyxrnd->Sumw2();
  hypx->Sumw2();
  hypxrnd->Sumw2();
  hyEt->Sumw2();
  hyEtrnd->Sumw2();
  hybtgm->Sumw2();
  hybtgmrnd->Sumw2();

  TH2D* hyutcrr;

  auto hdndydut   = new TH2D("hdndydut",";y;ut", 40, -2., 2., 50, 0., 2.);
  auto gdndycrr   = new TGraphErrors();   gdndycrr -> SetName("gdndycrr");

  auto huy        = new TH2D("huy","; y_{cm}/y_{proj}; u_{t0} ", 200,-1., 1.8, 200,0.,2.5);
  auto hpsi       = new TH2D("hpsi",";#Psi;",15, 0., 15., 200,0.,1.);
  auto hpsi1      = new TH1D("hpsi1",";#Psi",500,0.3,0.8);
  auto hpsindx    = new TH2D("hpsindx",";#Psi ;index ",100,-TMath::Pi(),TMath::Pi(),20,0.,20.);
  auto hiphi      = new TH1I("hiphi","hiphi",13,0,13);
  auto hPsi       = new TH1D("hPsi"   ,";RP #Psi"  ,100,-TMath::Pi(),TMath::Pi());
  auto hPsinc     = new TH1D("hPsinc" ,"w/o Cut;RP #Psi"  ,100,-TMath::Pi(),TMath::Pi());
  auto hRPPsi     = new TH1D("hRPPsi",";Indiv. #Psi"  ,100,-TMath::Pi(),TMath::Pi());
  auto hRPPsipsi  = new TH2D("hRPPsipsi",";RP #Psi; Indiv #Psi"  ,100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi());
  auto hphi       = new TH1D("hphi",     ";#phi"  ,100,-180.,180.);
  auto hcostm     = new TH1D("hcostm", ";cosphi", 100,-TMath::Pi(),TMath::Pi());
  auto h_vartl    = new TH1D("h_vartl",";vartl;",100,0.,5.);
  auto hycos1     = new TH2D("hycos1", ";y_{cm}/y_{proj}; cos(dphi)" ,100,-1.,1.8, 100.,-1., 1.);
  auto hycos2     = new TH2D("hycos2", ";y_{cm}/y_{proj}; cos(2dphi)",100,-1.,1.8, 100.,-1., 1.);
  auto hycos2am   = new TH2D("hycos2am",";y_{cm}/y_{proj}; cos(2dphi)AM",100,-1.,1.8, 100.,-2.2, 2.2);
  auto hyanum     = new TH2D("hyanum", ";y_{cm}/y_{proj}; Atmic Number" ,100,-1.,1.8, 100.,-2.2, 2.2);

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

  TH1D *hdut1[ptbin1][2];
  TH1D *hdut2[ptbin2];
  TH1D *hdutcos1[ptbin1][2];
  TH1D *hdutcos2[ptbin2];

  //-----  Booking
  TString rangeLabel1[ybin1];
  TString rangeLabel2[ybin2];

  // replace y_cm -> y_bm

  // booking for v1
  for( UInt_t iy = 0; iy < ybin1; iy++ ) {
    rangeLabel1[iy] = Form("%5.2f <= y < %5.2f"      ,yrange1nrm[iy],yrange1nrm[iy+1]);

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
    rangeLabel2[iy] = Form("%5.2f <= y < %5.2f",yrange2nrm[iy],yrange2nrm[iy+1]);
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

  for( UInt_t ipt = 0; ipt < ptbin1; ipt++ ) {
    hdut1[ipt][0]    = new TH1D(Form("hdut1_ut%d_0",ipt), "0.2<y_{0}<0.6", 500,0., 2.5);
    hdutcos1[ipt][0] = new TH1D(Form("hdutcos1_%d_0",ipt),"<cos(#phi)>", 100, -1, 1.);
    hdut1[ipt][1]    = new TH1D(Form("hdut1_ut%d_1",ipt), "0.4<y_{0}<0.8", 500,0., 2.5);
    hdutcos1[ipt][1] = new TH1D(Form("hdutcos1_%d_1",ipt),"<cos(#phi)>", 100, -1, 1.);
  }

  for( UInt_t ipt = 0; ipt < ptbin2; ipt++ ) {
    hdut2[ipt]    = new TH1D(Form("hdut2_ut%d",ipt),"-0.4<y_{0}<0.4", 500,0., 2.5);
    hdutcos2[ipt] = new TH1D(Form("hdutcos2_%d",ipt), "<cos(2#phi)>",100, -1, 1.);
  }


  //------------------------------
  Double_t u_p = 0.355151 * 1.06974; //p+p(268.9MeV/u)

  // Double_t *rppsires  = new Double_t[3];  // Psi dependent correction
  // Double_t *rpres     = new Double_t[4];  // 

  TDatime dtime;
  TRandom3 grnd(dtime.GetSecond());
  gRandom->SetSeed(dtime.GetSecond());

  Long64_t nEntry = SetBranch();
  
  //--------------------------------------------------
  //--- Event Loop
  //--------------------------------------------------

  for(Long64_t i = 0; i < nEntry; i++){

    ShowProcess(i);

    rChain->GetEntry(i);
    //    STFlowInfo *aflow = (STFlowInfo*)aFlowArray->At(0);


    /// for reaction plane resolution
    Bool_t bFill = kFALSE;
    Bool_t bRes  = kFALSE;


    //@@@@@ centrality selection ++++++++++++++++++++++++++++++++++++++++++++++++++
    //    Int_t trackselection = aFlowInfo->mtrack2;
    Int_t trackselection = aFlowInfo->mtrack1;
    
    if(trackselection > Ucent || trackselection < Lcent ) continue; //v56.0.0

    hmult->Fill( trackselection );

    ///    auto RPangle = GetRPBaseAngle(aFlowInfo);
    auto RPangle = aFlowInfo->unitP_fc.Phi();
    hPsinc->Fill(RPangle);

    bRes = kTRUE; //@1

    Double_t subevt_phi = abs(TVector2::Phi_mpi_pi((aFlowInfo->unitP_1fc).Phi()-
						   (aFlowInfo->unitP_2fc).Phi()));
   
 
    Double_t rnd_Phi = 2.*TMath::Pi()*(grnd.Rndm() - 0.5);

    Double_t sumPz = 0;
    Double_t sumPt = 0;

    TIter next(aArray);
    STKParticle *aPart = NULL;

    //--------------------------------------------------
    //----- Main loop 
    //--------------------------------------------------
    UInt_t mtk = 0;

    while( (aPart = (STKParticle*)next()) ) {

      //&&&&&
      if( !( isys == 5 || aPart->GetGoodTrackFlag() >= 111110 ) ) continue;

      mtk++;

      // phi cut
      auto phi    = aPart->GetRotatedMomentum().Phi()*TMath::RadToDeg();
      // w phicut                                                                                                                                              
      switch(phicutID) {
      case 0:
	if( abs(phi) > 45 && abs(phi) < 135) continue;
	phiAcp = 90./180.;
	break;
      case 1:
	if( abs(phi) > 45) continue;
	phiAcp = 45./180.;
	break;
      case 2:
	if( abs(phi) < 140 ) continue;
	phiAcp = 40./180.;
	break;
      case 3:
	if( abs(phi) < 45 || abs(phi) > 135) continue;
	phiAcp = 90./180.;
	break;
      case 5:
	if( phi > 20 || phi < -30 ) continue;
	phiAcp = (30.+20.)/360.;
	break;
      case 6:
	if(( phi > 20 && phi < 140) || (phi > -140 && phi < -30 ) ) continue;
	phiAcp = (30.+20.+80)/360.;
	break;
      }

      hphi     ->Fill(phi);

      //------------------------------

      auto bmass = aPart->GetBBMass();
      auto theta = aPart->GetRotatedMomentum().Theta();
      auto charge= aPart->GetCharge();

      auto ipid   = aPart->GetPID();  

      //----- Particle Selection -----
      // selid:0>1H, 1:2H, 2:3H, 3:3He, 4:4He, 5:H, 6:HHe, 7:n,
      Bool_t bpart = kFALSE;

      if( ipid == 0 ) continue;

      if( selid < 5 && selid == PDGtoIDX[ipid] ) 
	bpart = kTRUE;

      else if( selid == 5 && PDGtoIDX[ipid] < 3 ) {
	AtmNumber  = ncls[PDGtoIDX[ipid]].Z;
	bpart = kTRUE;
      }
      else if( selid == 6 && PDGtoIDX[ipid] < 5 ) {
	AtmNumber  = ncls[PDGtoIDX[ipid]].Z;
	bpart = kTRUE;
      }

      hpid->Fill(aPart->GetRotatedMomentum().Mag()*aPart->GetCharge(), aPart->GetdEdx());

      //-------------------
      if( !bpart && isys != 5 ) continue; //default
      //-----------------

      hpidsel->Fill(aPart->GetRotatedMomentum().Mag()*aPart->GetCharge(), aPart->GetdEdx());

      bFill = kTRUE;

      y_norm = yNN[isys];

      TLorentzVector lrnzVec = aPart->GetLorentzVector();
      TVector3 boostVec = STLorentzBoostVector::GetBoostVector(4+isys);

      if( selid == 5 || selid == 6 ) {
	auto bg = lrnzVec.Beta() * lrnzVec.Gamma();
	auto mom = aPart->GetRotatedMomentum().Unit();
	mom.SetMag( ncls[selid].mass * bg );
	lrnzVec.SetVectM( mom, ncls[selid].mass );
      }

      auto rapidl = aPart->GetRapidity();

      //      cout << pid << " : rapidity1 " << rapidl << " vec " << lrnzVec.Rapidity() << endl;
      //      auto rapid  = aPart->GetRapiditycm();;	
      auto rapidn = lrnzVec.Rapidity() / y_norm - 1.;
      auto pt     = lrnzVec.Pt();
      auto rpphi  = aPart->GetIndividualRPAngle();
      auto dphi   = aPart->GetAzmAngle_wrt_RP();
      auto dphi2  = aPart->GetAzmAngle2_wrt_RP();
      fmass  = aPart->GetMass();
      Double_t u_t0  = aPart->GetRotatedMomentum().Pt()/fmass/u_p;
      Double_t ou_t0 = aPart->GetMomentumAtTarget().Pt()/fmass/u_p;


      lrnzVec.Boost(-boostVec);
      TLorentzVector transvec (lrnzVec.Z(),0,lrnzVec.Pt(),lrnzVec.E());
      TLorentzVector transvecx(lrnzVec.Z(),lrnzVec.Pt()*sin(dphi), lrnzVec.Pt()*cos(dphi), lrnzVec.E());
      TLorentzVector transvecxrnd(lrnzVec.Z(),lrnzVec.Pt()*sin(rnd_Phi), lrnzVec.Pt()*cos(rnd_Phi), lrnzVec.E());
      if( abs(lrnzVec.Rapidity()/y_norm) < 1 && abs(transvecx.Rapidity()/y_norm) < 1. ) {
	sumPz += abs(lrnzVec.Rapidity()/y_norm);
	sumPt += abs(transvecx.Rapidity()/y_norm);
      }

      hmass    ->Fill(aPart->GetRotatedMomentum().Mag(), bmass);
      if( aPart->GetPID() > 100002000 )
	hmassHe  ->Fill(aPart->GetRotatedMomentum().Mag(), aPart->GetBBMass());

      hPsi     ->Fill(RPangle);
      hRPPsi   ->Fill(rpphi);
      hRPPsipsi->Fill(RPPsi, rpphi);
      hypt     ->Fill(rapidn, pt);
      hyut     ->Fill(rapidn, u_t0);

      hyyx     ->Fill(rapidn, transvecx.Rapidity()/y_norm);
      hyyxrnd  ->Fill(rapidn, transvecxrnd.Rapidity()/y_norm);
      hypx     ->Fill(rapidn, pt*cos(dphi)/ncls[selid].A);
      hypxrnd  ->Fill(rapidn, pt*cos(rnd_Phi)/ncls[selid].A);

      hycos1   ->Fill(rapidn, cos(dphi));
      hycos2   ->Fill(rapidn, cos(2.*dphi)); 

      hycos2am ->Fill(rapidn, cos(2.*dphi)*(Double_t)ncls[PDGtoIDX[ipid]].Z);
      hyanum   ->Fill(rapidn, (Double_t)ncls[PDGtoIDX[ipid]].Z);

      //      cout << "ncls[PDGtoIDX[ipid]].Z " << ncls[PDGtoIDX[ipid]].Z << " PID " << ipid << " " << cos(2.*dphi)*(Double_t)ncls[PDGtoIDX[ipid]].Z <<endl;

      Double_t Ek   = lrnzVec.E()-fmass;
      Double_t btgm = lrnzVec.Pt()/fmass; 

      hyEt     ->Fill(rapidn, Ek * sin(lrnzVec.Theta()));
      hyEtrnd  ->Fill(rapidn, Ek * sin(lrnzVec.Theta()*cos(rnd_Phi)));
      hybtgm   ->Fill(rapidn, btgm );
      hybtgmrnd  ->Fill(rapidn, btgm * cos(rnd_Phi));

      if( selid == 6 && AtmNumber == 2. ) {
	hypt     ->Fill(rapidn, pt);
	hyyx     ->Fill(rapidn, transvecx.Rapidity()/y_norm);
	hyyxrnd  ->Fill(rapidn, transvecxrnd.Rapidity()/y_norm);
	hypx     ->Fill(rapidn, pt*cos(dphi)/ncls[selid].A);
	hypxrnd  ->Fill(rapidn, pt*cos(rnd_Phi)/ncls[selid].A);

	hyEt     ->Fill(rapidn, Ek * sin(lrnzVec.Theta()));
	hyEtrnd  ->Fill(rapidn, Ek * sin(lrnzVec.Theta()*cos(rnd_Phi)));
	hybtgm   ->Fill(rapidn, btgm );
	hybtgmrnd  ->Fill(rapidn, btgm * cos(rnd_Phi));

      }


      UInt_t irapid1 = GetV1RapidityIndex(rapidn);
      UInt_t ipt1    = GetV1PtIndex(u_t0);

      UInt_t irapid2 = GetV2RapidityIndex(rapidn);
      UInt_t ipt2    = GetV2PtIndex(u_t0);

      // cout << "irapid1 "<< irapid1 << " irapid2  "<< irapid2 << endl;

      // Comparison with Tommy
      if( rapidn >0.7 && rapidn < 0.8 )
	hcostm -> Fill( cos(dphi) );

      hirap1 -> Fill( rapidn, (Double_t)irapid1 );
      hirap2 -> Fill( rapidn, (Double_t)irapid2 );

      if( rapidn >= 0.2 && rapidn <= 0.6 ) {
	hdut1[ipt1][0]     -> Fill( u_t0 );
	hdutcos1[ipt1][0]  -> Fill( cos(dphi) );
      }
      if( rapidn >= 0.4 && rapidn <= 0.8 ) {
	hdut1[ipt1][1]     -> Fill( u_t0 );
	hdutcos1[ipt1][1]  -> Fill( cos(dphi) );
      }
      
      if( abs(rapidn) <= 0.4 ) {
	hdut2[ipt2]     -> Fill( u_t0 );
	hdutcos2[ipt2]  -> Fill( cos(2.*dphi) );
      }
      

      //@@@@@
      //v1 -----------
      hutphi10[irapid1] -> Fill( subevt_phi );

      if( kTRUE ){
      
	hdy1[irapid1]    -> Fill( rapidn );
	hdycos1[irapid1] -> Fill( cos(dphi) );
	
	hdydut1[irapid1][ipt1]    -> Fill( u_t0 );
	hdydutcos1[irapid1][ipt1] -> Fill( cos(dphi) );

	//v2 ------------
	
	hutphi20[irapid2] -> Fill( subevt_phi );
	
	hdy2[irapid2]     -> Fill( rapidn );
	hdycos2[irapid2]  -> Fill( cos(2.*dphi) );
	
	hdydut2[irapid2][ipt2] -> Fill( u_t0 );
	hdydutcos2[irapid2][ipt2] -> Fill( cos(2.*dphi) );
      }

      ///$$$$$////

      if( u_t0 > 0.) {

	huy->Fill( rapidn, u_t0 );

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

    Double_t vartl = 2.*sumPt/(TMath::Pi()*sumPz);
    h_vartl -> Fill(vartl);


    if( mtk > 0 )
      hmultfill->Fill( trackselection );
  }

  //--------------------------------------------------
  //--- Enf of event Loop
  //--------------------------------------------------
  LOG(INFO) << " End of Event Loop " << FairLogger::endl;
  LOG(INFO) << " Number of Events " << hmult->GetEntries() 
	    << " with good tracks -> " << hmultfill->GetEntries() <<  FairLogger::endl; 

  //@@@@@ acceptance correction ++++++++++++++++++++++++++++++++++++++++++++++++++
  Double_t scalefactor = 1./(hyut->GetXaxis()->GetBinWidth(0)*hyut->GetYaxis()->GetBinWidth(0)*hmultfill->GetEntries());
  hyut -> Scale(scalefactor/phiAcp);

  scalefactor = 1./(hypt->GetXaxis()->GetBinWidth(0)*hypt->GetYaxis()->GetBinWidth(0)*hmultfill->GetEntries());
  hypt -> Scale(scalefactor/phiAcp);

  scalefactor = 1./(hyyx->GetXaxis()->GetBinWidth(0)*hyyx->GetYaxis()->GetBinWidth(0)*hmultfill->GetEntries());
  hyyx -> Scale(scalefactor/phiAcp);

  scalefactor = 1./(hyyxrnd->GetXaxis()->GetBinWidth(0)*hyyxrnd->GetYaxis()->GetBinWidth(0)*hmultfill->GetEntries());
  hyyxrnd -> Scale(scalefactor/phiAcp);

  scalefactor = 1./(hypx->GetXaxis()->GetBinWidth(0)*hypx->GetYaxis()->GetBinWidth(0)*hmultfill->GetEntries());
  hypx -> Scale(scalefactor/phiAcp);

  scalefactor = 1./(hypxrnd->GetXaxis()->GetBinWidth(0)*hypxrnd->GetYaxis()->GetBinWidth(0)*hmultfill->GetEntries());
  hypxrnd -> Scale(scalefactor/phiAcp);

  scalefactor = 1./(hyEt->GetXaxis()->GetBinWidth(0)*hyEt->GetYaxis()->GetBinWidth(0)*hmultfill->GetEntries());
  hyEt -> Scale(scalefactor/phiAcp);

  scalefactor = 1./(hyEtrnd->GetXaxis()->GetBinWidth(0)*hyEtrnd->GetYaxis()->GetBinWidth(0)*hmultfill->GetEntries());
  hyEtrnd -> Scale(scalefactor/phiAcp);

  scalefactor = 1./(hybtgm->GetXaxis()->GetBinWidth(0)*hybtgm->GetYaxis()->GetBinWidth(0)*hmultfill->GetEntries());
  hybtgm -> Scale(scalefactor/phiAcp); 

  scalefactor = 1./(hybtgmrnd->GetXaxis()->GetBinWidth(0)*hybtgmrnd->GetYaxis()->GetBinWidth(0)*hmultfill->GetEntries());
  hybtgmrnd -> Scale(scalefactor/phiAcp); 

  if( bAcc_corr ) {
    hyutcrr = GetAcpCorrection(*hyut, hAcpYUtCorr);
    hyutcrr -> SetName("hyutcrr");
  }

  //-- For Tommy's comparison
  //  hycos2am->Divide(hyanum);

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

  TGraphErrors *g_utv1[2];
  for( UInt_t in = 0; in < 2; in++ ) {
    g_utv1[in] = new TGraphErrors();
    g_utv1[in] -> SetName(Form("g_utv1_%d",in));
  }

  TGraphErrors *g_utv2 = new TGraphErrors();
  g_utv2 -> SetName("g_utv2");

  TGraphErrors *gUt_v1[ybin1];
  TGraphErrors *gUt_v2[ybin2];

  for(UInt_t kn = 0; kn < ybin1 ; kn++){
    gUt_v1[kn] = new TGraphErrors();
    gUt_v1[kn]->SetName((TString)Form("gUt_v1%d",kn));
    TString sname = ncls[selid].sName+" "+rangeLabel1[kn]+"; u_{t0}; v1";
    gUt_v1[kn]->SetTitle(sname); 
  }

  for(UInt_t kn = 0; kn < ybin2 ; kn++){
    gUt_v2[kn] = new TGraphErrors();
    gUt_v2[kn]->SetName((TString)Form("gUt_v2%d",kn));
    TString sname = ncls[selid].sName+" "+rangeLabel2[kn]+"; u_{t0}; v2";
    gUt_v2[kn]->SetTitle(sname);
  }


  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++@@
  // This function returen resolution from subevent correlation with overall events.
  Double_t *rpresall = new Double_t[2];
  GetRPResolutionwChi(rpresall, hdphi0_180, hdphi90_180, 1.);
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //--- v1 ------------------------------
  for( UInt_t in = 0; in < 2; in++ ) {
    UInt_t iutv = 0;

    for( UInt_t ipt = 0; ipt < ptbin1; ipt++ ) { 
      Double_t utn    = (Double_t)hdut1[ipt][in]->GetEntries();
      if( utn == 0 ) continue;

      Double_t utmean = hdut1[ipt][in]->GetMean();
      Double_t ute    = hdut1[ipt][in]->GetStdDev()/sqrt(utn);

      Double_t v1u   = hdutcos1[ipt][in]->GetMean();
      Double_t v1ue  = hdutcos1[ipt][in]->GetStdDev()/sqrt(utn);
      
      Double_t v1ut  = v1u / rpresall[0];
      Double_t v1ute = GetError(v1u, rpresall[0], v1ue, rpresall[1]);
      
      g_utv1[in] -> SetPoint( iutv, utmean, v1ut );
      g_utv1[in] -> SetPointError( iutv, ute, v1ute );
      iutv++;
    }

    g_utv1[in] -> SetTitle( (TString)hdut1[0][in]->GetTitle()+"; u_{t0}; v1" );
  }


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

    if( ymean != 0 ) {
      gy_v1->SetPoint(iyy, ymean, cosc);
      gy_v1->SetPointError(iyy, ystdv, mcosee );
      iyy++;
    }
    
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

      if( !std::isnan(utmean) && !std::isnan(v1ut) && v1ute < 1. ) {
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

	  GetdNdydUt(1, iy, ipt, hyutcrr);
	  ///$$$$$////
	  hdndydut->Fill(yrange1nrm[iy], dpt1*ipt, acpcorr[0]);
	}
	else {
	  acpcorr[0] = utn;
	  acpcorr[1] = 0.;
	}
	//---------------------

	v1u   = hdydutcos1cut[iy][ipt]->GetMean();
	v1ue  = hdydutcos1cut[iy][ipt]->GetStdDev()/sqrt(utn);
	v1ut  = v1u / rpresall[0];
	v1ute = GetError(v1u, rpresall[0], v1ue, rpresall[1]);

	v1ut *= acpcorr[0];
	v1ute = pow(acpcorr[0]*v1ute,2) ;

	v1u_sum += v1ut;	
	v1u_ste += v1ute ;
	v1u_n   += acpcorr[0];
      
      }
    }

    if( yn > 0 && v1u_n > 0 && ystdv != 0) {
      Double_t v1u_ave = v1u_sum / v1u_n;
      Double_t v1u_err = sqrt(v1u_ste) / v1u_n;
      
      gu_v1->SetPoint( iyv, ymean, v1u_ave);
      gu_v1->SetPointError( iyv, ystdv, v1u_err);
      
      iyv++;
    }

    gUt_v1[iy]->Write();
  }
  

  LOG(INFO) <<" Comparison with Tommy " << hcostm->GetMean() << " => " <<hcostm->GetMean()/ rpresall[0]  << " : " << hmult->GetMean() << FairLogger::endl;

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //--- <px>
  
  auto gy_px = new TGraphErrors();
  gy_px->SetName("gy_px");


  //@@@@@ meanpx
  if( hypx->GetEntries() > 0) {       
    UInt_t i = 0;
    for(auto iy : ROOT::TSeqI(hypx->GetXaxis()->GetNbins()) ) {
      //      for( auto ipt : ROOT::TSeqI(hypx->GetYaxis()->GetNbins()) ) {
	
      auto rap = hypx->GetXaxis()->GetBinCenter(iy);
      auto hpx = hypx->ProjectionY("hpx",iy,iy);
      auto mean = hpx->GetMean();;
      mean /= rpresall[0];
      auto meanerr = hpx->GetMeanError();
      meanerr = GetError(mean, rpresall[0], meanerr, rpresall[1]);
	
      if( meanerr > 0 ) {
	gy_px -> SetPoint(i, rap, mean);
	gy_px -> SetPointError(i, 0, meanerr );
	i++;
      }
      //      }
    }
    
    ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
    cc->SetGrid();
    cc->Divide(2,2);
    cc->cd(1);
    hyEt->Draw("colz");
    cc->cd(2);
    hyEtrnd->Draw("colz");
    cc->cd(3);
    hybtgm->Draw("colz");
    cc->cd(4);
    hybtgmrnd->Draw("colz");
  }

  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //--- v2

  GetRPResolutionwChi(rpresall, hdphi0_180, hdphi90_180, 2.);


  
  UInt_t iutv = 0;
  for( UInt_t ipt = 0; ipt < ptbin2; ipt++ ) { 
    Double_t utn    = (Double_t)hdut2[ipt]->GetEntries();
    
    if( utn == 0 ) continue;

    Double_t utmean = hdut2[ipt]->GetMean();
    Double_t ute    = hdut2[ipt]->GetStdDev()/sqrt(utn);

    Double_t v2u   = hdutcos2[ipt]->GetMean();
    Double_t v2ue  = hdutcos2[ipt]->GetStdDev()/sqrt(utn);

    Double_t v2ut  = v2u / rpresall[0];
    Double_t v2ute = GetError(v2u, rpresall[0], v2ue, rpresall[1]);


    g_utv2 -> SetPoint( iutv,utmean, v2ut );
    g_utv2 -> SetPointError( iutv, ute, v2ute );
    iutv++;
  }
  g_utv2 -> SetTitle((TString)hdut2[0]->GetTitle()+"; u_{t0}; v2" );


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

      if( !std::isnan(utmean) && !std::isnan(ute) && v2ute < 1) {
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

	  if( !hyutcrr -> IsZombie() ) {
	    GetdNdydUt(2, iy, ipt, hyutcrr);
	  }
	}
	else {
	  acpcorr[0] = utn;
	  acpcorr[1] = 0.;
	}

	v2u   = hdydutcos2cut[iy][ipt] -> GetMean();
	v2ue  = hdydutcos2cut[iy][ipt] -> GetStdDev()/ sqrt(utn);          
	v2ut  = v2u / rpresall[0];
	v2ute = GetError(v2u, rpresall[0], v2ue, rpresall[1]) ;	


	v2ut *= acpcorr[0];
	v2ute = pow(acpcorr[0]*v2ute,2) ;
	
	v2u_sum += v2ut;
	v2u_ste += v2ute;
	v2u_n   += acpcorr[0];

      }

      hdydutcos2[iy][ipt]->Write();
      hdydutcos2cut[iy][ipt]->Write();
    }

    if( v2u_n > 0 ) {
               
      Double_t v2u_ave = v2u_sum / v2u_n;
      Double_t v2u_err = sqrt(v2u_ste)/v2u_n;
      
      if( ymean > -1. && ystdv != 0) {
	gu_v2->SetPoint( iyv, ymean, v2u_ave);
	gu_v2->SetPointError( iyv, ystdv, v2u_err );
	
	iyv++;
      }
    }
    gUt_v2[iy]->Write();
  }



  //@@@@@ Write ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  gu_v1->Write();
  gu_v2->Write();
  gy_v1->Write();
  gy_v2->Write();
  gy_px->Write();
  g_utv1[0]    ->Write();
  g_utv1[1]    ->Write();
  g_utv2    ->Write();

  GraphSave -> Write();

  return;
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
    hyutcrr->Draw("colz");

    cc->cd(4);
    hAcpyytCorr->Draw("colz");
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
  
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  g_utv1[0]->Draw("ALP");

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  g_utv2->Draw("ALP");
  g_utv2->Print();


  gSystem->cd("..");
}



//**************************************************
Bool_t   LoadAcceptanceCorrection(UInt_t selid)
{
  //  if(selid == 10 ) selid = 2;

  TFile *fopen = GetAcceptanceCorrectionFile(0,phicutID,"womc");
  if( !fopen ) 
    return kFALSE;

  //  cout << selid << " " << ncls[selid].pdg << " " << ncls[selid].Name  << endl;

  TString hname = "h2UtYEff_108Sn_" + ncls[selid].Name+"_mbin0_iter1_GaussBlur";
  hAcpYUtCorr = (TH2D*)fopen->Get(hname);  
  if( !hAcpYUtCorr ) {
    LOG(ERROR) << hname << " is not found. " << FairLogger::endl;
    fopen->Close();
    return kFALSE;
  }
  LOG(INFO) << hAcpYUtCorr->GetName() << " is loaded. " << FairLogger::endl;
  hAcpYUtCorr -> SetDirectory( 0 );

  hname = "h2yytEff_108Sn_" + ncls[selid].Name +"_mbin0_iter1_GaussBlur";
  hAcpyytCorr = (TH2D*)fopen->Get(hname);

  if( !hAcpyytCorr ) {
    LOG(ERROR) << hname << " is not found " << FairLogger::endl;
    fopen->Close();
    return kFALSE;
  }
  LOG(INFO) << hAcpyytCorr->GetName() << " is loaded. " << FairLogger::endl;
  hAcpyytCorr -> SetDirectory(0);

  hname = "h2yyxEff_108Sn_" + ncls[selid].Name+"_mbin0_iter1_GaussBlur";
  hAcpyyxCorr = (TH2D*)fopen->Get(hname);

  if( !hAcpyyxCorr ) {
    LOG(ERROR) << hname << " is not found " << FairLogger::endl;
    fopen->Close();
    return kFALSE;
  }
  LOG(INFO) << hAcpyyxCorr->GetName() << " is loaded. " << FairLogger::endl;
  hAcpyyxCorr -> SetDirectory(0);

  //-- yyxhalf
  hname = "h2yxhEff_108Sn_" + ncls[selid].Name+"_mbin0_iter1_GaussBlur";
  hAcpyxhCorr = (TH2D*)fopen->Get(hname);

  if( !hAcpyxhCorr ) {
    LOG(ERROR) << hname << " is not found " << FairLogger::endl;
    fopen->Close();
    return kFALSE;
  }
  LOG(INFO) << hAcpyxhCorr->GetName() << " is loaded. " << FairLogger::endl;
  hAcpyxhCorr -> SetDirectory(0);

  fopen->Close();
  return kTRUE;
}

void DrawUt(TH2D *h2, UInt_t sel, TString opt, Color_t scol)
{
  TString hname = h2->GetName();
  cout << " DrawUt -> " << hname << endl;

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
    cc->cd(i+1);

    auto firstbin = h2->GetXaxis()->FindBin(*(yrange+i));
    auto lastbin  = h2->GetXaxis()->FindBin(*(yrange+i+1));

    //    cout << " ybin " << i << " " << *(yrange+i) << " y " << firstbin << " ~ " << *(yrange+i+1) << endl;

    hacp = (TH1D*) h2->ProjectionY(hname+Form("_%d_%d",sel,i), firstbin, lastbin);
    hacp -> SetTitle(Form("%3.2f < y < %3.2f" , *(yrange+i), *(yrange+i+1))); 
    if( hacp->GetEntries() > 0) {
      hacp -> SetNormFactor(50);
      //      hacp -> Scale(1./hacp->GetEntries() );

      hacp->SetLineColor(scol);
      
      hacp -> Draw(opt);
    }
    else 
      cout << " NO entry in ybin " << i << endl;

  }
}

void GetdNdydPt( Double_t ycm, Double_t pt)
{
  auto ybin = hAcpYPtCorr->GetXaxis()->FindBin(ycm);
  auto xbin = hAcpYPtCorr->GetYaxis()->FindBin(pt);

  acpcorr[0] = hAcpYPtCorr->GetBinContent(ybin, xbin);
  acpcorr[1] = hAcpYPtCorr->GetBinError(ybin, xbin);  
}


TH1D* GetdNdy(TH2D* h2c, Double_t ntotal)
{
  TH1D* h_dndy = new TH1D();
  h_dndy -> Sumw2();

  h_dndy  = h2c->ProjectionX("h_dndy");
  if( h_dndy == NULL ) return NULL;

  Double_t rnorm  = phiAcp * (h_dndy->GetXaxis()->GetBinWidth(0));

  cout << h2c->GetName()
       << " phiacp = " << phiAcp 
       << " binwidth = " << h_dndy->GetXaxis()->GetBinWidth(0)
       << " total = " << ntotal
       << " y_cm = " << y_cm[isys+6]
       << endl;

  
  h_dndy -> Scale(1./(rnorm*ntotal));
  return h_dndy;
}


TH1D* GetdNdXt(TH2D* h2c, TString htitle, Double_t ylow, Double_t yup, Double_t ntotal)
{
  TH1D* h_dndt = new TH1D();
  h_dndt -> Sumw2();

  Double_t lbin = h2c->GetXaxis()->FindBin(ylow);
  Double_t ubin = h2c->GetXaxis()->FindBin(yup);
  h_dndt = h2c->ProjectionY(htitle,lbin,ubin);

  Double_t rnorm  = phiAcp* (yup-ylow)*h_dndt->GetXaxis()->GetBinWidth(0);
  h_dndt -> Scale(1./(rnorm*ntotal));

  return h_dndt;
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

  if( val > (ptbin1-1)  * dpt1 )
    return ptbin1-1;

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
  if( val > (ptbin2-1)  * dpt2 )
    return ptbin2-1;

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
  //  fefffile = TFile::Open("/home/mizuki/EffCorrection/rootfiles/LCPSpectra/UnfoldedLCPSpectra.slim.root");
  //  fefffile = TFile::Open("/home/mizuki/EffCorrection/rootfiles/LCPSpectra47to52/UnfoldedLCPSpectra.slim.root");
  fefffile = TFile::Open("/home/mizuki/EffCorrection/rootfiles/LCPSpectra47to52/UnfoldedLCPSpectra.slim.root");

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
    
    if(aFlowInfo->mtrack1 > Ucent || aFlowInfo->mtrack1 < Lcent ) continue;
    count++;

    cosphi += cos( (aFlowInfo->unitP_fc).Phi());      
    sinphi += sin( (aFlowInfo->unitP_fc).Phi());

    cos2phi += cos(2.*(aFlowInfo->unitP_fc).Phi());      
    sin2phi += sin(2.*(aFlowInfo->unitP_fc).Phi());

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

  TCut mcrot = Form("mtrack4>%f&&mtrack4<%f",mlt[0]*2.,mlt[1]*2.);
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


