#include "openRunAna.C"
#include "DoFlow.h"

auto *fcos1 = new TF1("fcos1","[0]+2.*[1]*cos(x)"   ,-1.*TMath::Pi(),TMath::Pi());

//drawing

UInt_t selReactionPlanef = 10000;

// Retrivew tree
TClonesArray *aArray; 
TClonesArray *aFlowArray;
TClonesArray *aNLClusterArray;



// functions
void     GetFittingParameters(TH1D &h1, Double_t pp[6]);
void     GetFittingParameters(TH1D &h1, Double_t pp[6], Double_t corr[2]);
Double_t GetRapidity_cm(TVector3 p, Double_t mass, TVector3 bvec);
UInt_t   SetBranch();
void     PlotPtDependence(UInt_t selid);     
Double_t *GetRPResolutionwChi(TH1D *hphi0_180, TH1D *hphi90_180);
void     PlotCosPtDependence(UInt_t selid);

Double_t  GetError(Double_t x, Double_t y, Double_t xe, Double_t ye)
{
  if( y == 0 && ye == 0 ) return 0.;
  Double_t xc = abs(x/y);
  return  xc * sqrt(pow(xe/x,2) + pow(ye/y,2));
}

//Double_t rapoffset[] = {0., 0.018, -0.008, 0.028};
Double_t rapoffset[] = {0.01113, -0.00645, 0.01898, -0.01726};
UInt_t   npb = 30;

UInt_t  icent = 0;
UInt_t  jcent = 3;

//-------------------//
void DoFlow(UInt_t isel = 2) 
{
  gROOT->Reset();

  openRunAna();

  if(rChain != NULL)     
    std::cout << " System " << isys << "  -> " << sysName << std::endl; 

  else
    exit(0);


  gROOT->ProcessLine(".! grep -i void DoFlow.C ");


  // Multiplicity cut ================================
  TString su = gSystem -> Getenv("UC");
  if( su != "" ) {
    UInt_t ucent = (UInt_t)atoi(su);
    if( ucent < sizeof(cent)/sizeof(UInt_t)-1)
      icent = ucent;
  }


  su = gSystem -> Getenv("LC");
  if( su != "" ) {
    UInt_t lcent = (UInt_t)atoi(su);
    if(  lcent < sizeof(cent)/sizeof(UInt_t))
      jcent = lcent;
  }

  std::cout << " Multiplicity :: " << cent[icent] << " to " << cent[jcent] << std::endl;

  //==================================================

  PlotCosPtDependence(isel);  
}


//pdt
void PlotCosPtDependence(UInt_t selid = 2)       //%% Executable :
{
  gStyle->SetOptStat(0);

  std::cout << "PlotCosPtDependence(" << selid << ")" << std::endl;
  // PT binning
  Double_t pt_max = 800.;
  UInt_t   nbin1  = 10; //16
  UInt_t   nbin2  = 8; //10
  Double_t dpt1   = pt_max/(Double_t)nbin1;
  Double_t dpt2   = pt_max/(Double_t)nbin2;

  auto cutfile = new TFile("db/RegionCut.root");
  TCutG *goodThetaPhi = (TCutG*)cutfile->Get("goodThetaPhi");
  cutfile->Close();    


  //--------------------------------------------------
  //----- SETUP OUTPUT FILE --------------------------
  //--------------------------------------------------
  gSystem->cd("data");
  TString oVer = gSystem->Getenv("OUTVER");
  TString fName = "cosYPt_"+ sysName + "_" + partname[selid]+".v"+sVer+"."+oVer;

  if( oVer == "" ) {
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

  std::cout << "File " << fName << " is created. " << std::endl;  

  //--------------------------------------------------


  auto GraphSave = new TFile(fName,"recreate");
  auto hphi0_180  = new TH1D("hphi0_180" ,"#Phi from  0 to 180; #Phi",100,0.,3.2);
  auto hphi90_180 = new TH1D("hphi90_180","#Phi from 90 to 180; #Phi",100,0.,3.2);

  Int_t pcharge = 1;
  if(selid == 0)  pcharge = -1;

  std::cout << " Rapidity binning " << ybin1 << std::endl;
  TString rangeLabel1[ybin1];
  for(UInt_t i = 0; i < ybin1; i++ ){
    if( i == 0 )
      rangeLabel1[0] = Form(" y < %5.2f ",yrange1[0]);
    else if ( i == ybin1 -1 )
      rangeLabel1[i] = Form("%5.2f <= y "    ,yrange1[ybin1-1]);
    else 
      rangeLabel1[i] = Form("%5.2f <= y < %5.2f",yrange1[i-1],yrange1[i]);
  }
  TString rangeLabel2[ybin2];
  for(UInt_t i = 0; i < ybin2; i++ ){
    if( i == 0 )
      rangeLabel2[0] = Form(" y < %5.2f ",yrange2[0]);
    else if ( i == ybin1 -1 )
      rangeLabel2[i] = Form("%5.2f <= y "    ,yrange2[ybin2-1]);
    else 
      rangeLabel2[i] = Form("%5.2f <= y < %5.2f",yrange2[i-1],yrange2[i]);
  }
  //------------------------------
  //-----  booking
  auto hmass   = new TH2D("hmass",   ";P/Q; Mass [MeV/c^2]"     ,200,  0.,2500., 200, 0.,7000);
  auto hmassy  = new TH2D("hmassy",  "; Y_{cm}; Mass [MeV/C^2]" ,200, -0.4, 0.5, 200, 0.,7000);
  auto hpy     = new TH2D("hpy",     "; Y_{cm}; P/Q [MeV/c]"    ,200, -0.4, 0.5, 200, 0.,2500.);
  auto hptmass = new TH2D("hptmass", "; Pt [MeV/c]; Mass [MeV/c^2]",200, 0., 800,200, 0.,7000);
  auto hazm    = new TH1D("hazm",    "; #phi"                   ,100, -3.2, 3.2);
  auto hmidv2  = new TH1D("hmidv2",  "; #Delta#phi"             ,100, -3.2, 3.2);

  TString hlabel  = (TString)Form("mtrack4 %2d to %d ; Y_{cm}; Pt [MeV/c]",
				  (UInt_t)cent[icent],(UInt_t)cent[jcent]);
  auto hyptacp    = new TH2D("hyptacp", hlabel ,200, -0.4, 0.5, 200, 0., 800);
  auto hntrack    = new TH1I("hntrack",  "; Number of good track", 60, 0, 60);
  auto hyawpitch  = new TH2D("hyawpitch","; Yaw angle; Pitch angle",200,-1.5,1.5,200,-1.5,1.5);
  auto hyaw       = new TH1D("hyaw",     "; Yaw angle"  , 100, -90., 90);
  auto hpitch     = new TH1D("hpitch",   "; Pitch angle", 100, -90., 90);
  TH1D *hyphi1[ybin1];
  TH1D *hyphi2[ybin2];
  TH1D *hyptphi1[ybin1][nbin1];
  TH1D *hyptphi2[ybin2][nbin2];
  TH1D *hypt1[ybin1];
  TH1D *hypt2[ybin1];
  //------------------------------

  std::vector< std::vector< Double_t > > cosv1x(ybin1);
  std::vector< std::vector< Double_t > > cosv2x(ybin2);
  Double_t  cosv1[ybin1];
  Double_t  cosv2[ybin2];
  Double_t  sinv1[ybin1];
  Double_t  sinv2[ybin2];


  //---> v1 vs rapidity
  std::vector< std::vector< std::vector< Double_t> > > cosv1ptx(ybin1);
  Double_t  cosv1pt[ybin1][nbin1];
  Double_t  sinv1pt[ybin1][nbin1];
  for( UInt_t i = 0; i < ybin1; i++) {
    cosv1[i]   = 0.;
    sinv1[i]   = 0.;
    cosv1x[i].clear();
    cosv1ptx[i].resize(nbin1);

    hyphi1[i] = new TH1D( Form("hyphi1_%d",i),rangeLabel1[i]+"#De #Phi"   , npb, -3.15, 3.15);
    hypt1[i] = new TH1D( Form("hypt1_%d",i), rangeLabel1[i]+"; Pt", 100, 0., 800.);

    for( UInt_t j = 0; j < nbin1; j++ ){
      cosv1pt[i][j] = 0.;
      sinv1pt[i][j] = 0.;

      hyptphi1[i][j] = new TH1D( Form("hyptphi1_%d%d",i,j),rangeLabel1[i]+"#Delta #phi"   , 
				 npb, -3.15, 3.15);
    }
  }

  //---> v2 vs rapidity
  std::vector< std::vector< std::vector< Double_t> > > cosv2ptx(ybin2);
  Double_t  cosv2pt[ybin2][nbin2];
  Double_t  sinv2pt[ybin2][nbin2];
  for( UInt_t i = 0; i < ybin2; i++) {
    cosv2[i]   = 0.;
    sinv2[i]   = 0.;
    cosv2x[i].clear();
    cosv2ptx[i].resize(nbin2);
    hyphi2[i] = new TH1D( Form("hyphi2_%d",i),rangeLabel2[i]+"2x #Delta #phi", npb, -3.15, 3.15);
    hypt2[i]  = new TH1D( Form("hypt2_%d",i), rangeLabel2[i]+"; Pt", 100, 0., 800.);

    for( UInt_t j = 0; j < nbin2; j++ ){
      cosv2pt[i][j] = 0.;
      sinv2pt[i][j] = 0.;
      hyptphi2[i][j] = new TH1D( Form("hyptphi2_%d%d",i,j),rangeLabel2[i]+"2x #Delta #phi"   , 
				 npb, -3.15, 3.15);
    }
  }
    
  


  Int_t nevt = SetBranch();
  cout << " Number of events " << nevt << endl;
  
  //--------------------------------------------------
  //--- Event Loop
  //--------------------------------------------------
  for(Int_t i = 0; i < nevt; i++){

    rChain->GetEntry(i);
    STFlowInfo *aflow = (STFlowInfo*)aFlowArray->At(0);

    /// for reaction plane resolution
    Bool_t bFill = kFALSE;


    /// centrality selection
    if(aflow->mtrack4 > cent[icent] || aflow->mtrack4 <= cent[jcent] ) continue;


    TIter next(aArray);
    STParticle *aPart = NULL;

    //--------------------------------------------------
    //----- Main loop 
    //--------------------------------------------------
    while( (aPart = (STParticle*)next()) ) {
	
      auto rpf   = aPart->GetReactionPlaneFlag();
      auto pid   = aPart->GetPID();
      auto yaw   = aPart->GetYawAngle();
      auto pitch = aPart->GetPitchAngle();
      auto ndf   = aPart->GetNDF();


      //------------------------------
      //----- Particle Selection -----
      if( pid == partid[selid] && rpf > 10000 ){ //&& yaw < 0 && pitch < 0) {
      
	bFill = kTRUE;
      
	auto pt    = aPart->GetRotatedMomentum().Pt();
	auto dphi  = aPart->GetAzmAngle_wrt_RP();
	auto rapid = aPart->GetRapiditycm();;	
	auto bmass = aPart->GetBBMass();
	auto phi   = aPart->GetRotatedMomentum().Phi();
	auto theta = aPart->GetRotatedMomentum().Theta();

	hazm     ->Fill( phi );	
	hmass    ->Fill( aPart->GetRotatedMomentum().Mag(), bmass);
	
	hyptacp  ->Fill( rapid, pt );
	hmassy   ->Fill( rapid, bmass );
	hpy      ->Fill( rapid, aPart->GetRotatedMomentum().Mag() );
	hptmass  ->Fill( pt,    bmass );
	hyawpitch->Fill( yaw,   pitch );
	hyaw     ->Fill( yaw/TMath::DegToRad() );
	hpitch   ->Fill( pitch/TMath::DegToRad() );

	UInt_t irapid = ybin1 - 1;
	for( UInt_t i = 0; i < ybin1; i++){
	  if(rapid < yrange1[i]){
	    irapid = i;
	    break;
	  }
	}
	
	hyphi1[irapid]->Fill(dphi);
	hypt1[irapid] ->Fill(aPart->GetRotatedMomentum().Pt());
	cosv1x[irapid].push_back( rapid );
	cosv1[irapid] += cos(dphi);
	sinv1[irapid] += sin(dphi);

	UInt_t ipt = nbin1 - 1;
	for(UInt_t i = 0; i < nbin1; i++){
	  if( pt < dpt1*(i+1)) {
	    ipt = i;
	    break;
	  }
	}

	hyptphi1[irapid][ipt]->Fill(dphi);
	cosv1ptx[irapid][ipt].push_back( pt );
	cosv1pt[irapid][ipt] += cos(dphi);
	sinv1pt[irapid][ipt] += sin(dphi);

	irapid = ybin2 - 1;

	for( UInt_t i = 0; i < ybin2; i++){
	  if(rapid < yrange2[i]){
	    irapid = i;
	    break;
	  }
	}

	hyphi2[irapid]->Fill(TVector2::Phi_mpi_pi(2.*dphi));
	hypt2[irapid] ->Fill(aPart->GetRotatedMomentum().Pt());
	cosv2x[irapid].push_back( rapid );
	cosv2[irapid] += cos(2.*dphi);
	sinv2[irapid] += sin(2.*dphi);

	ipt = nbin2 - 1;
	for(UInt_t i = 0; i < nbin2; i++){
	  if( pt < dpt2*(i+1)) {
	    ipt = i;
	    break;
	  }
	}

	hyptphi2[irapid][ipt]->Fill(TVector2::Phi_mpi_pi(2.*dphi));
	cosv2ptx[irapid][ipt].push_back( pt );
	cosv2pt[irapid][ipt] += cos(2.*dphi);
	sinv2pt[irapid][ipt] += sin(2.*dphi);
      }
    }

    if( bFill ) {
      hntrack->Fill( aflow->mtrack4 );

      Double_t subevt_phi = abs(TVector2::Phi_mpi_pi((aflow->unitP_1fc).Phi()-
						     (aflow->unitP_2fc).Phi())); 

      hphi0_180->Fill( subevt_phi );
      if( subevt_phi > TMath::Pi()/2. )
	hphi90_180->Fill( subevt_phi );
    }
  }
  //--------------------------------------------------
  //--- Enf of event Loop
  //--------------------------------------------------
    

  TGraphErrors *gv_v1 = new TGraphErrors();
  gv_v1->SetName("gv_v1");
  TGraphErrors *gv_v2 = new TGraphErrors();
  gv_v2->SetName("gv_v2");
  
  TGraphErrors *gPt_v1[ybin1];
  TGraphErrors *gPt_v2[ybin1];
    
  for(UInt_t kn = 0; kn < ybin1 ; kn++){      
    gPt_v1[kn] = new TGraphErrors();
    gPt_v1[kn]->SetName((TString)Form("gPt_v1%d",kn));
    TString sname = partname[selid]+"; Pt [MeV/c]; v1";
    gPt_v1[kn]->SetTitle(sname);
  }
  
  for(UInt_t kn = 0; kn < ybin2 ; kn++){      
    gPt_v2[kn] = new TGraphErrors();
    gPt_v2[kn]->SetName((TString)Form("gPt_v2%d",kn));
    TString sname = partname[selid]+"; Pt [MeV/c]; v2";
    gPt_v2[kn]->SetTitle(sname);
  }


  Double_t *rpres = new Double_t[4];
  rpres = GetRPResolutionwChi(hphi0_180, hphi90_180);
  std::cout << " <cos(Phi)> = " << rpres[0] << " +- " << rpres[1] 
	    << " <cos(2Phi)> = "<< rpres[2] << " +- " << rpres[3] 
	    << std::endl;
  std::cout << " --------------------------------- " << std::endl;

  std::cout << " ---- Resutls ---------------------" << std::endl;
  
  UInt_t kl = 0;
  UInt_t id1 = 0;
  UInt_t id2 = 0;
  Double_t para[6];

    
  for(UInt_t kn = 0; kn < ybin1; kn++){
      
    if( cosv1x[kn].size() == 0 ) continue;

    //v1 rapidity dependence
    Double_t rapm  = TMath::Mean(  cosv1x[kn].begin(), cosv1x[kn].end());
    Double_t rape  = TMath::StdDev(cosv1x[kn].begin(), cosv1x[kn].end());
      
    Double_t yv1  = cosv1[kn] / (Double_t)cosv1x[kn].size();
    Double_t yv1e = sinv1[kn] / (Double_t)cosv1x[kn].size();

    Double_t yv1c  = yv1/rpres[0];
    Double_t yv1ce = GetError(yv1, rpres[0], yv1e, rpres[1]);

    if( !std::isnan(rapm) ) {
      gv_v1->SetPoint( kl,      rapm, yv1c);
      gv_v1->SetPointError( kl, rape, yv1ce);
      kl++;
    

      // pt dependence 
      UInt_t il = 0; 
      for(UInt_t jn = 0; jn < nbin1; jn++){
	
	if( cosv1ptx[kn][jn].size() > 0 ) {	
	  
	  Double_t ptc  = TMath::Mean(cosv1ptx[kn][jn].begin(), cosv1ptx[kn][jn].end());
	  Double_t ptce = TMath::StdDev(cosv1ptx[kn][jn].begin(), cosv1ptx[kn][jn].end());
	  
	  Double_t ptv1  = cosv1pt[kn][jn] / (Double_t)cosv1ptx[kn][jn].size();
	  Double_t ptv1e = sinv1pt[kn][jn] / (Double_t)cosv1ptx[kn][jn].size();
	  
	  yv1c  = ptv1/rpres[0];
	  yv1ce = GetError(ptv1, rpres[0], ptv1e, rpres[1]);
	  
	  if( kn == 0 ){
	    cout << jn << " erro " << sinv1pt[kn][jn] << " yvc1  " << yv1ce <<  endl;
	    cout << sinv1pt[kn][jn] << " --  " << (Double_t)cosv1ptx[kn][jn].size() << endl;
	  }

	  if( ptce < 100 ) {
	    gPt_v1[kn]->SetPoint(il, ptc, yv1c);
	    gPt_v1[kn]->SetPointError(il, ptce, yv1ce);
	    
	    il++;
	  }
	}
      }
    }
    gPt_v1[kn]->Write();
    hypt1[kn]->Write();
  }

  kl = 0;
  for(UInt_t kn = 0; kn < ybin2; kn++){

    if( cosv2x[kn].size() == 0 ) continue;

    Double_t rapm  = TMath::Mean(cosv2x[kn].begin(), cosv2x[kn].end());
    Double_t rape  = TMath::StdDev(cosv2x[kn].begin(), cosv2x[kn].end());

    Double_t yv2  = cosv2[kn] / (Double_t)cosv2x[kn].size();
    Double_t yv2e = sinv2[kn] / (Double_t)cosv2x[kn].size();

    Double_t yv2c  = yv2 / rpres[2];
    Double_t yv2ce = GetError(yv2, rpres[2], yv2e, rpres[3]);

    if( !std::isnan(rapm) ) {
      gv_v2->SetPoint( kl,    rapm, yv2c);
      gv_v2->SetPointError( kl, rape, yv2ce);      
      kl++;
    

      UInt_t il = 0; 
      for(UInt_t jn = 0; jn < nbin2; jn++){
	
	if( cosv2ptx[kn][jn].size() > 0 ){

	  Double_t ptc  = TMath::Mean(cosv2ptx[kn][jn].begin(), cosv2ptx[kn][jn].end());
	  Double_t ptce = TMath::StdDev(cosv2ptx[kn][jn].begin(), cosv2ptx[kn][jn].end());

	  Double_t ptv2  = cosv2pt[kn][jn] / (Double_t)cosv2ptx[kn][jn].size();
	  Double_t ptv2e = sinv2pt[kn][jn] / (Double_t)cosv2ptx[kn][jn].size();
	  
	  yv2c  = ptv2/ rpres[2];
	  yv2ce = GetError(ptv2, rpres[2], ptv2e, rpres[3]);

	  if( ptce < 100 ) {
	    gPt_v2[kn]->SetPoint(il, ptc, yv2c);
	    gPt_v2[kn]->SetPointError(il, ptce, yv2ce);
	    
	    il++;
	  }
	}
      }
    }
    gPt_v2[kn]->Write();
    hypt2[kn]->Write();
  }

  gv_v1->Write();
  gv_v2->Write();


  //--------------------------------------------------
  //--- Ploting
  //--------------------------------------------------
  id=1;
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  cc->Divide(2,2);

  cc->cd(id); id++;
  hmass->Draw("colz");

  cc->cd(id); id++;
  hntrack->Draw();

  cc->cd(id); id++;
  hpy->Draw("colz");

  // cc->cd(id); id++;
  // hptmass->Draw("colz");

  cc->cd(id); id++;
  hazm->Draw();
  
  id = 1;
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  cc->Divide(2,2);
  cc->cd(id); id++; id++;
  hyawpitch->Draw("colz");
  cc->cd(id); id++;
  hyaw->Draw();
  cc->cd(id); id++;
  hpitch->Draw();

  
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));
  hyptacp->Draw("colz");
  


  id = 1;
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),2000,500);
  cc->Divide(ybin1,2);
  for(UInt_t kn = 0; kn < ybin1; kn++){
    cc->cd(id); id++;
    gPt_v1[kn]->Draw("ALP");

    cc->cd(id+ybin1-1);
    hypt1[kn]->Draw();
  }

  id = 1;
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),1600,500);
  cc->Divide(ybin2,2);
  for(UInt_t kn = 0; kn < ybin2; kn++){
    cc->cd(id); id++;
    gPt_v2[kn]->Draw("ALP");

    cc->cd(id+ybin2-1);
    hypt2[kn]->Draw();
    
  }


  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  gv_v1->Draw("ALP");

  auto LineV = new TLine(0.,gv_v1->GetYaxis()->GetXmin(), 0., gv_v1->GetYaxis()->GetXmax());
  auto LineH = new TLine(gv_v1->GetXaxis()->GetXmin(),    0., gv_v1->GetXaxis()->GetXmax(), 0.);
  LineV->SetLineStyle(3);
  LineH->SetLineStyle(3);
  LineV->Draw();
  LineH->Draw();

  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),600,400);
  gv_v2->Draw("ALP");

  ic++; cc   = new TCanvas("dyphi","dphi1 y bin",500,1200);
  cc->Divide(2, ybin1);

  //--- y-bin and pt-bin
  auto cc1   = new TCanvas("dphi1","dphi1 y bin and pt bin",2000,1200);
  cc1->Divide(nbin1, ybin1);

  auto cc2 = new TCanvas("dphi2","dphi2 ybin and pt bin",1400,1200);
  cc2->Divide(nbin2, ybin2);

  id = 1;
  for(UInt_t i = 0; i < ybin1; i++) {
    cc->cd(2*(i+1)-1);
    hyphi1[i]->SetNormFactor(npb);
    hyphi1[i]->Draw("e");

    for(UInt_t j = 0; j < nbin1; j++){
      cc1->cd(id);    id++;
      
      if(hyptphi1[i][j]->GetEntries() == 0) continue;
      
      hyptphi1[i][j]->SetNormFactor(npb);
      hyptphi1[i][j]->Draw("e");
    }
  } 

  id = 1;
  for(UInt_t i = 0; i < ybin2; i++) {
    cc->cd(2*(i+1));
    hyphi2[i]->SetNormFactor(60);
    hyphi2[i]->Draw("e");

    for(UInt_t j = 0; j < nbin2; j++){
      cc2->cd(id); id++;
      
      if(hyptphi2[i][j]->GetEntries() == 0) continue;
      
      hyptphi2[i][j]->SetNormFactor(npb);
      hyptphi2[i][j]->Draw("e");

    }
  } 

  hphi0_180->Write();
  hyptacp->Write();
  hntrack->Write();

  gSystem->cd("..");

  SaveCanvas(fName);



  //  ic++; 
  //  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic));

  //  gROOT->cd();
}


void PlotAzimuthalDistribution(UInt_t selid = 2)       //%% Executable :
{

  std::cout << "----->  PlotPtDependence(" << selid << ")" << std::cout;

  gStyle->SetOptStat(0);

  gSystem->cd("data");

  TString fName = Form("YPT_n%2dto%d_",cent[icent],cent[jcent]);
  fName  +=  sysName + "_xs" + partname[selid]+"_1.root";
  auto GraphSave = new TFile(fName,"recreate");
  std::cout << "File " << fName << " is created. " << std::endl;

  Int_t pcharge = 1;
  if(selid == 0)  pcharge = -1;

  Int_t RPflag = 0;
  if(selid >= 2) RPflag = 10;

  cout << " Particle " << partname[selid] << endl;

  Double_t yrange[] = {-4,-0.06,0.06,4.};

  std::cout << " Rapidity binning " << ybin1 << std::endl;
  TH1D *hphi[3];

  for(UInt_t i = 0; i < 3; i++ ){
    TString rangeLabel = Form("%5.2f <= y < %5.2f",yrange[i],yrange[i+1]);
    hphi[i]  = new TH1D(Form("hphi%d",i) ,rangeLabel+"; #Phi",80,-3.14,3.14);
  }

  
  Int_t nevt = SetBranch();
  cout << " Number of events " << nevt << endl;
  
  for(Int_t i = 0; i < nevt; i++){
    rChain->GetEntry(i);
      
    TIter next(aArray);
    STParticle *aPart = NULL;

    while( (aPart = (STParticle*)next()) ) {
	
      auto rpf   = aPart->GetReactionPlaneFlag();
      auto pid   = aPart->GetPID();
      auto bmass = aPart->GetBBMass();
      
      if( pid == partid[selid] ) { //&& 

	auto pt    = aPart->GetRotatedMomentum().Pt();
	auto dphi  = aPart->GetAzmAngle_wrt_RP();
	auto rapid = aPart->GetRapiditycm();;


	if( rapid >= yrange[0] && rapid < yrange[1] )
	  hphi[0]->Fill( dphi );
	else if( rapid >= yrange[1] && rapid < yrange[2] )
	  hphi[1]->Fill( dphi );
	else if( rapid >= yrange[2] && rapid < yrange[3] )
	  hphi[2]->Fill( dphi );
      }
    }
  }

  for(UInt_t i = 0; i < 3; i++) 
    hphi[i]->SetNormFactor(80);


  

    

  ic++; 
  ic++; cc   = new TCanvas("cc0","dphi1", 1000, 500);
  cc->Divide(3,1);
  cc->cd(1);

  auto *fv1 = new TF1("fv1","[0]+2.*[1]*cos(x)"   ,-1.*TMath::Pi(),TMath::Pi());
  hphi[0]->Fit("fv1","","",-3.1,3.1);
  auto p0 = fv1->GetParameter(0);
  auto p1 = fv1->GetParameter(1);
  fv1->SetParameter(0, 1.);
  fv1->SetParameter(1, p1/p0);
  hphi[0]->Draw("e");
  fv1->Draw("same");

  //  ic++; 
  //  cc   = new TCanvas("cc1","dphi1");

  cc->cd(2);
  auto *fv2 = new TF1("fv2","[0]+2.*[1]*cos(x) + 2.*[2]*cos(2.* x)"   ,-1.*TMath::Pi(),TMath::Pi());
  hphi[1]->Fit("fv2","","",-3.1,3.1);
  p0 = fv2->GetParameter(0);
  p1 = fv2->GetParameter(1);
  auto p2 = fv2->GetParameter(2);
  fv2->SetParameter(0, 1.);
  fv2->SetParameter(1, p1/p0);
  fv2->SetParameter(2, p2/p0);
  hphi[1]->Draw("e");
  fv2->Draw("same");

  //  ic++; 
  //  cc   = new TCanvas("cc2","dphi1");
  cc->cd(3);
  auto *fv1p = new TF1("fv1p","[0]+2.*[1]*cos(x)"   ,-1.*TMath::Pi(),TMath::Pi());
  hphi[2]->Fit("fv1p","","",-3.1,3.1);
  p0 = fv1p->GetParameter(0);
  p1 = fv1p->GetParameter(1);
  fv1p->SetParameter(0, 1.);
  fv1p->SetParameter(1, p1/p0);
  hphi[2]->Draw("e");
  fv1p->Draw("same");

}

//**************************************************
//**************************************************

Double_t *GetRPResolutionwChi(TH1D *hphi0_180, TH1D *hphi90_180)            //%% Executable : 
{
  hphi90_180->SetLineColor(2);

  Double_t m0 = hphi0_180->GetEntries();
  Double_t m1 = hphi90_180->GetEntries();

  Double_t chi = sqrt(-2.* log(2.* m1/m0));
  
  Double_t m01e = m1/m0*sqrt(1./m1+1./m0);
  Double_t chie = chi - sqrt(-2.* log(2. * (m1/m0+m01e)));

  Double_t *rpres = new Double_t[4];

  rpres[0] = sqrt(TMath::Pi())/(2.*TMath::Gamma(1))*chi;
  rpres[1] = sqrt(TMath::Pi())/(2.*TMath::Gamma(1))*(chi+chie) - rpres[0];
  rpres[2] = sqrt(TMath::Pi())/(2.*2.*TMath::Gamma(1.5))*pow(chi,2);
  rpres[3] = sqrt(TMath::Pi())/(2.*2.*TMath::Gamma(1.5))*pow(chi+chie,2) - rpres[2];

  // ic++;
  // cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  // hphi0_180->Draw();
  // hphi90_180->Draw("same");

  //  Double_t para[6];
  //  GetFittingParameters( *hphi0_180, para);

  std::cout << " Resolution : v1 " 
	    << std::setw(14) << rpres[0] << " +- " << std::setw(10) << rpres[1]
	    << " v2 "
	    << std::setw(14) << rpres[2] << " +- " << std::setw(10) << rpres[3]
	    << std::endl;

  return rpres;
}


void FlatteningCheck()            
{
  //----- Parametres
  Int_t ntrack[7];

  //----- Booking
  TH2D *hphitheta;
  TH2D *hphimtrck;

  hphitheta = new TH2D("hphitheta", sysName+"; #Theta ; #Phi",100,0,0.8, 100,-3.2, 3.2);
  hphimtrck = new TH2D("hphimtrck", sysName+"; Number of Track ; #Phi",60,0,60, 100,-3.2, 3.2);

  //----- Filling
  Int_t nEntry = SetBranch();

  for(Int_t i = 0; i < nEntry; i++){

    rChain->GetEntry(i);
    STFlowInfo *aflow = (STFlowInfo*)aFlowArray->At(0);
    if( aflow == NULL ) continue;
    

    TIter next(aArray);
    STParticle *aPart = NULL;
    
    while( (aPart = (STParticle*)next()) ) {

      auto phiv  = aPart->GetIndividualRPVector();
      auto theta = aPart->GetRotatedMomentum().Theta();
      auto flag  = aPart->GetReactionPlaneFlag();

      if(flag >= selReactionPlanef ){
	hphitheta->Fill( theta, phiv.Phi() );
	hphimtrck->Fill( aflow->mtrack4, phiv.Phi() ); 
      }
    }
  }


  //----- Drawing 
  ic++;
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  hphitheta->Draw("colz");

  ic++;
  ic++; cc = new TCanvas(Form("cc%d",ic),Form("cc%d",ic),700,500);
  hphimtrck->Draw("colz");

}



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


void GetFittingParameters(TH1D &h1, Double_t pp[6])
{
  for(UInt_t i = 0; i < 6; i++)
    pp[i] = 0.;


  Double_t nbin = h1.GetNbinsX();
  Double_t scf  = h1.GetEntries();
  h1.SetNormFactor(nbin);
  // h1.SetMaximum(1.2);
  // h1.SetMinimum(0.8);
  fcos1->SetParameter(0, 1.);
  fcos1->SetParameter(1, 0.1);

  h1.Fit("fcos1","","",-3.1,3.1);

  pp[0] = fcos1->GetParameter(0);
  pp[1] = fcos1->GetParameter(1);
  pp[2] = fcos1->GetParError(0);
  pp[3] = fcos1->GetParError(1);
  pp[3] = pp[1]*sqrt(pp[2]/pp[0]*pp[2]/pp[0] + pp[3]/pp[1]*pp[3]/pp[1]);
  pp[1] = pp[1]/pp[0];
  
  

  h1.SetNormFactor(-1);
}
void GetFittingParameters(TH1D &h1, Double_t pp[6], Double_t corr[2])
{
  GetFittingParameters(h1,pp);
  
  if( corr[0] != 0 && corr[1] != 0 ) {
    pp[4] = pp[1]/corr[0];
    pp[5] = sqrt( pow(pp[4],2)*( pow(pp[3]/pp[1],2) + pow(corr[1]/corr[0], 2) ));
  }
  else {
    pp[4] = 0.;
    pp[5] = 0.;
  }
}

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

    
    if(aflow->mtrack4 > cent[icent] || aflow->mtrack4 <= cent[jcent] ) continue;      
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

void CentralityDependence()            //%% Executable :
{
  //  UInt_t mrange[] = {70, 50, 45, 40, 35, 25, 20, 0};
  UInt_t mrange[] = {70, 35, 28, 0};

  const UInt_t mbin = sizeof(mrange)/sizeof(UInt_t) - 1;
  TH1D *hphi0_180;
  TH1D *hphi90_180;
  TH1I *hmult;

  auto gv_mcos1 = new TGraphErrors();
  gv_mcos1->SetTitle(";Multiplicity; <cos(#Psi)>");
  auto gv_mcos2 = new TGraphErrors();
  gv_mcos2->SetTitle(";Multiplicity; <cos(2#Psi)>");


  auto cc80 = new TCanvas("cc80","cc80",500,1000);
  cc80->Divide(2,7);
  
  UInt_t id = 1;
  for(UInt_t i = 0; i < mbin; i++) {
    cc80->cd(id); id++;
    TCut hcut = Form("mtrack4>%u&&mtrack4<=%u",mrange[i+1],mrange[i]);

    TString htitle = Form("hphi0_%d",i);
    hphi0_180  = new TH1D(htitle, "",100,0.,3.2);
    rChain->Project(htitle,"abs(TVector2::Phi_mpi_pi(unitP_1->Phi()-unitP_2->Phi()))",hcut);
    hphi0_180->Draw();

    htitle = Form("hphi90_%d",i);
    hphi90_180 = new TH1D(htitle,"",100,0.,3.2);
    TCut phicut = "abs(TVector2::Phi_mpi_pi(unitP_1->Phi()-unitP_2->Phi()))>1.5707963";
    rChain->Project(htitle,"abs(TVector2::Phi_mpi_pi(unitP_1->Phi()-unitP_2->Phi()))",hcut&&phicut);
    hphi90_180->Draw("same");

    hmult = new TH1I(Form("hmult_%u",i),"",70,0.,70.);
    rChain->Project(Form("hmult_%u",i),"mtrack4",hcut);
    cc80->cd(id); id++;
    hmult->Draw();
    

     Double_t *rpres = new Double_t[4];
     rpres = GetRPResolutionwChi(hphi0_180, hphi90_180);
     std::cout << mrange[i+1] << " ~ " << mrange[i] 
     	      << " <cos(Phi)> = " << rpres[0] << " +- " << rpres[1] 
     	      << " <cos(2Phi)> = "<< rpres[2] << " +- " << rpres[3] 
     	      << std::endl;
     
     hmult->GetMean();

     gv_mcos1->SetPoint(i, hmult->GetMean(), rpres[0]);
     gv_mcos1->SetPointError(i, hmult->GetStdDev(), rpres[1]);

     gv_mcos2->SetPoint(i, hmult->GetMean(), rpres[2]);
     gv_mcos2->SetPointError(i, hmult->GetStdDev(), rpres[3]);
  }

  hmult = new TH1I("hmult",";Multiplicity",70,0,70);
  rChain->Project("hmult","mtrack4");

  auto cc79 = new TCanvas("cc79","cc79");
  hmult->Draw("");

  auto cc81 = new TCanvas("cc81","cc81");
  gv_mcos1->Draw("ALP");

  auto cc82 = new TCanvas("cc82","cc82");
  gv_mcos2->Draw("ALP");
}


UInt_t SetBranch()
{
  if(rChain == NULL) {
    std::cout << " no file is loaded " << std::endl;
    return 0;
  }

  std::cout << " Nentry ->  " << rChain->GetEntries() << std::endl;

  if(aArray != NULL)
    aArray->Clear();

  rChain->SetBranchAddress("STParticle",&aArray);
  rChain->SetBranchAddress("STFlow"    ,&aFlowArray);

  return rChain->GetEntries();
}

